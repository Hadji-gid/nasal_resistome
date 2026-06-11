# =============================================================================
# Nasal Resistome Mali – Snakemake Workflow
# Author  : Alhadji A. Dicko | ACE-B / INSP Mali
# Version : 1.5  (2026)
#
# Full pipeline from raw reads to Clinical Risk Score:
#   1. fastp        – QC and adapter trimming
#   2. host_depletion  – Human DNA removal (Bowtie2 vs GRCh38)
#   3. spades_assembly – Metagenomic assembly (metaSPAdes)
#   4. quast           – Assembly quality assessment
#   5. kraken2         – Taxonomic classification
#   6. bracken         – Species abundance re-estimation
#   7. amrfinder       – ARG annotation (CARD)
#   8. abricate_vfdb   – Virulence factor annotation (VFDB)
#   9. abricate_plasmidfinder – Plasmid detection
#  10. hgt_and_crs_scoring    – HGT potential + Clinical Risk Score (R)
#  11. visualize              – Bubble plot + ARG heatmap (R)
#
# Usage:
#   snakemake --cores 8 --jobs 2 --rerun-incomplete \
#             --rerun-triggers mtime --latency-wait 30
# =============================================================================

from pathlib import Path

configfile: "config/config.yaml"

with open(config["samples_list"]) as fh:
    SAMPLES = [s.strip() for s in fh if s.strip() and not s.startswith("#")]

D           = config["outdir"]
RAW         = config["raw_dir"]
BRACKEN_BIN = config["bracken_bin"]

# ── Target rule ───────────────────────────────────────────────────────────────
rule all:
    input:
        # QC
        expand("{d}/multiqc_report.html",                  d=D["qc"]),
        # Host depletion stats
        expand("{d}/{sample}/depletion_stats.txt",         d=D["host_depl"], sample=SAMPLES),
        # Assembly
        expand("{d}/{sample}/contigs.fasta",               d=D["assembly"],   sample=SAMPLES),
        # QUAST
        expand("{d}/{sample}/report.tsv",                  d=D["asm_qc"],     sample=SAMPLES),
        # Taxonomy
        expand("{d}/{sample}/{sample}.bracken.txt",        d=D["taxonomy"],   sample=SAMPLES),
        # Annotation
        expand("{d}/{sample}/amrfinder.tsv",               d=D["annotation"], sample=SAMPLES),
        expand("{d}/{sample}/abricate_vfdb.tsv",           d=D["annotation"], sample=SAMPLES),
        expand("{d}/{sample}/abricate_plasmidfinder.tsv",  d=D["annotation"], sample=SAMPLES),
        # Scoring
        f"{D['scoring']}/final_risk_report.csv",
        f"{D['scoring']}/hgt_potential.csv",
        # Visualizations
        f"{D['reports']}/crs_bubble_plot.pdf",
        f"{D['reports']}/arg_heatmap.pdf",


# =============================================================================
# STEP 1 – QUALITY CONTROL  (fastp + MultiQC)
# =============================================================================
rule fastp_trim:
    """
    Adapter trimming, quality filtering, per-sample QC report.
    Clean reads are temp() – deleted once host_depletion finishes.
    """
    input:
        r1 = f"{RAW}/{{sample}}_R1_001.fastq.gz",
        r2 = f"{RAW}/{{sample}}_R2_001.fastq.gz",
    output:
        r1   = temp(f"{D['qc']}/{{sample}}/{{sample}}_R1_clean.fastq.gz"),
        r2   = temp(f"{D['qc']}/{{sample}}/{{sample}}_R2_clean.fastq.gz"),
        json = f"{D['qc']}/{{sample}}/{{sample}}_fastp.json",
        html = f"{D['qc']}/{{sample}}/{{sample}}_fastp.html",
    params:
        min_len = config["fastp"]["min_length"],
        qual    = config["fastp"]["qualified_quality"],
    threads: config["fastp"]["threads"]
    log:       f"{D['logs']}/fastp/{{sample}}.log"
    benchmark: f"{D['benchmarks']}/fastp/{{sample}}.txt"
    shell:
        """
        mkdir -p $(dirname {output.r1}) $(dirname {log})
        fastp \
            --in1 {input.r1} --in2 {input.r2} \
            --out1 {output.r1} --out2 {output.r2} \
            --json {output.json} --html {output.html} \
            --length_required {params.min_len} \
            --qualified_quality_phred {params.qual} \
            --detect_adapter_for_pe \
            --correction \
            --thread {threads} \
            2> {log}
        """

rule multiqc:
    """
    Aggregate all fastp JSON reports into one MultiQC summary.
    Depends only on JSON files (not temp) so always has its inputs.
    """
    input:
        expand("{d}/{sample}/{sample}_fastp.json", d=D["qc"], sample=SAMPLES)
    output:
        f"{D['qc']}/multiqc_report.html"
    params:
        indir  = D["qc"],
        outdir = D["qc"],
    log: f"{D['logs']}/multiqc.log"
    shell:
        """
        mkdir -p $(dirname {log})
        multiqc {params.indir} \
            --outdir {params.outdir} \
            --force \
            2> {log}
        """


# =============================================================================
# STEP 2 – HOST (HUMAN) DNA DEPLETION  (Bowtie2 + Samtools)
# =============================================================================
rule host_depletion:
    """
    Align reads to GRCh38; keep only unmapped (microbial) read pairs.
    Microbial FASTQs are temp() – deleted once spades + kraken2 both finish.
    Clean input FASTQs (temp) are deleted after this rule completes.
    """
    input:
        r1 = rules.fastp_trim.output.r1,
        r2 = rules.fastp_trim.output.r2,
    output:
        r1    = temp(f"{D['host_depl']}/{{sample}}/{{sample}}_R1_microbial.fastq.gz"),
        r2    = temp(f"{D['host_depl']}/{{sample}}/{{sample}}_R2_microbial.fastq.gz"),
        stats = f"{D['host_depl']}/{{sample}}/depletion_stats.txt",
    params:
        index  = config["host_index"],
        sens   = config["bowtie2"]["sensitivity"],
        tmpbam = f"{D['host_depl']}/{{sample}}/tmp_coord.bam",
    threads: config["bowtie2"]["threads"]
    log:       f"{D['logs']}/host_depletion/{{sample}}.log"
    benchmark: f"{D['benchmarks']}/host_depletion/{{sample}}.txt"
    shell:
        """
        mkdir -p $(dirname {output.r1}) $(dirname {log})

        # 1. Align → coordinate-sorted BAM
        bowtie2 \
            {params.sens} \
            -x {params.index} \
            -1 {input.r1} -2 {input.r2} \
            --threads {threads} \
            2>> {log} \
        | samtools view -bS -@ {threads} - \
        | samtools sort -@ {threads} \
            -T {params.tmpbam}_sort \
            -o {params.tmpbam} \
            2>> {log} || true

        # 2. Alignment stats
        samtools index {params.tmpbam}
        samtools flagstat {params.tmpbam} > {output.stats}

        # 3. Extract unmapped pairs only
        samtools view -b -f 12 -F 256 {params.tmpbam} \
        | samtools sort -n -@ {threads} \
            -T {params.tmpbam}_nsort \
        | samtools fastq \
            -1 {output.r1} \
            -2 {output.r2} \
            -0 /dev/null \
            -s /dev/null \
            -n \
            2>> {log} || true

        # 4. Verify outputs were created
        [ -s {output.r1} ] || {{ echo "ERROR: R1 empty"; exit 1; }}
        [ -s {output.r2} ] || {{ echo "ERROR: R2 empty"; exit 1; }}

        # 5. Cleanup BAM
        rm -f {params.tmpbam} {params.tmpbam}.bai
        """


# =============================================================================
# STEP 3 – METAGENOMIC ASSEMBLY  (metaSPAdes)
# =============================================================================
rule spades_assembly:
    """
    De-novo metagenomic assembly.
    Microbial FASTQ inputs (temp) deleted after both spades AND kraken2 finish.
    """
    input:
        r1 = rules.host_depletion.output.r1,
        r2 = rules.host_depletion.output.r2,
    output:
        contigs = f"{D['assembly']}/{{sample}}/contigs.fasta",
    params:
        outdir = f"{D['assembly']}/{{sample}}",
        mem    = config["spades"]["memory"],
    threads: config["spades"]["threads"]
    log:       f"{D['logs']}/spades/{{sample}}.log"
    benchmark: f"{D['benchmarks']}/spades/{{sample}}.txt"
    shell:
        """
        mkdir -p $(dirname {log})
        rm -rf {params.outdir}

        # Check minimum reads before assembly
        READ_COUNT=$(zcat {input.r1} | wc -l | awk '{{print $1/4}}')
        if [ "$READ_COUNT" -lt 1000 ]; then
            echo "SKIP: only $READ_COUNT reads — insufficient for assembly" > {log}
            mkdir -p {params.outdir}
            echo ">insufficient_reads_placeholder" > {output.contigs}
            echo "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" >> {output.contigs}
        else
            spades.py \
            --meta \
            -1 {input.r1} \
            -2 {input.r2} \
            -o {params.outdir} \
            --threads {threads} \
            --memory {params.mem} \
            -k 21,33,55,77 \
            --cov-cutoff auto \
            2> {log}
        fi
        """


# =============================================================================
# STEP 4 – ASSEMBLY QUALITY ASSESSMENT  (QUAST)
# =============================================================================
rule quast:
    """
    Assembly QC: N50, # contigs, largest contig. No reference (metagenome mode).
    Non-fatal: creates placeholder report for low-biomass samples with no contigs.
    """
    input:
        contigs = rules.spades_assembly.output.contigs,
    output:
        report = f"{D['asm_qc']}/{{sample}}/report.tsv",
    params:
        outdir     = f"{D['asm_qc']}/{{sample}}",
        min_contig = config["quast"]["min_contig"],
    threads: config["quast"]["threads"]
    log:       f"{D['logs']}/quast/{{sample}}.log"
    benchmark: f"{D['benchmarks']}/quast/{{sample}}.txt"
    shell:
        """
        mkdir -p $(dirname {log})
        quast.py \
            {input.contigs} \
            --output-dir {params.outdir} \
            --min-contig {params.min_contig} \
            --threads {threads} \
            --no-html \
            2> {log} || true

        [ -f {output.report} ] || \
            echo -e "Assembly\t{wildcards.sample}\nNote\tInsufficient microbial reads after host depletion" \
            > {output.report}
        """


# =============================================================================
# STEP 5 – TAXONOMIC PROFILING  (Kraken2 + Bracken)
# =============================================================================
rule kraken2:
    """
    Classify microbial reads against Kraken2 standard database.
    Raw output is temp() – deleted after Bracken finishes.
    Microbial FASTQs (temp) deleted after both spades AND kraken2 finish.
    """
    input:
        r1 = rules.host_depletion.output.r1,
        r2 = rules.host_depletion.output.r2,
    output:
        report = f"{D['taxonomy']}/{{sample}}/{{sample}}_kraken2_report.txt",
        out    = f"{D['taxonomy']}/{{sample}}/{{sample}}_kraken2_output.txt",
    params:
        db         = config["kraken2_db"],
        confidence = config["kraken2"]["confidence"],
        min_hits   = config["kraken2"]["min_hit_groups"],
    threads: config["kraken2"]["threads"]
    log:       f"{D['logs']}/kraken2/{{sample}}.log"
    benchmark: f"{D['benchmarks']}/kraken2/{{sample}}.txt"
    shell:
        """
        mkdir -p $(dirname {output.report}) $(dirname {log})
        kraken2 \
            --db {params.db} \
            --paired {input.r1} {input.r2} \
            --report {output.report} \
            --output {output.out} \
            --confidence {params.confidence} \
            --minimum-hit-groups {params.min_hits} \
            --gzip-compressed \
            --threads {threads} \
            2> {log} || true

        # Create empty outputs for zero-read samples
        [ -f {output.report} ] || touch {output.report}
        [ -f {output.out} ]    || touch {output.out}
        """

rule bracken:
    """
    Re-estimate species-level abundances from Kraken2 report.
    Kraken2 raw output (temp) deleted after this rule finishes.
    Uses explicit bracken binary path from config (bracken_bin).
    """
    input:
        report = rules.kraken2.output.report,
    output:
        bracken = f"{D['taxonomy']}/{{sample}}/{{sample}}.bracken.txt",
        report  = f"{D['taxonomy']}/{{sample}}/{{sample}}.bracken_report.txt",
    params:
        db        = config["bracken_db"],
        read_len  = config["bracken"]["read_len"],
        level     = config["bracken"]["level"],
        threshold = config["bracken"]["threshold"],
        bin       = BRACKEN_BIN,
    log:       f"{D['logs']}/bracken/{{sample}}.log"
    benchmark: f"{D['benchmarks']}/bracken/{{sample}}.txt"
    shell:
        """
        mkdir -p $(dirname {log})
        {params.bin} \
            -d {params.db} \
            -i {input.report} \
            -o {output.bracken} \
            -w {output.report} \
            -r {params.read_len} \
            -l {params.level} \
            -t {params.threshold} \
            2> {log} || true

        # Create empty outputs for zero-read samples
        [ -f {output.bracken} ] || touch {output.bracken}
        [ -f {output.report} ]  || touch {output.report}
        """


# =============================================================================
# STEP 6 – FUNCTIONAL ANNOTATION
# =============================================================================
rule amrfinder:
    """
    Detect ARGs in assembled contigs using AMRFinderPlus (nucleotide mode).
    """
    input:
        contigs = rules.spades_assembly.output.contigs,
    output:
        tsv = f"{D['annotation']}/{{sample}}/amrfinder.tsv",
    params:
        db    = config["amrfinder_db"],
        ident = config["amrfinder"]["ident_min"],
        cov   = config["amrfinder"]["coverage_min"],
    threads: config["amrfinder"]["threads"]
    log:       f"{D['logs']}/amrfinder/{{sample}}.log"
    benchmark: f"{D['benchmarks']}/amrfinder/{{sample}}.txt"
    shell:
        """
        mkdir -p $(dirname {output.tsv}) $(dirname {log})
        amrfinder \
            --nucleotide {input.contigs} \
            --database {params.db} \
            --ident_min {params.ident} \
            --coverage_min {params.cov} \
            --threads {threads} \
            --output {output.tsv} \
            2> {log}
        """

rule abricate_vfdb:
    """
    Screen contigs against VFDB (virulence factors).
    """
    input:
        contigs = rules.spades_assembly.output.contigs,
    output:
        tsv = f"{D['annotation']}/{{sample}}/abricate_vfdb.tsv",
    params:
        min_id  = config["abricate"]["min_id"],
        min_cov = config["abricate"]["min_cov"],
    log:       f"{D['logs']}/abricate_vfdb/{{sample}}.log"
    benchmark: f"{D['benchmarks']}/abricate_vfdb/{{sample}}.txt"
    shell:
        """
        mkdir -p $(dirname {log})
        abricate \
            --db vfdb \
            --minid {params.min_id} \
            --mincov {params.min_cov} \
            {input.contigs} \
            > {output.tsv} \
            2> {log}
        """

rule abricate_plasmidfinder:
    """
    Detect plasmid replicons in contigs using PlasmidFinder.
    """
    input:
        contigs = rules.spades_assembly.output.contigs,
    output:
        tsv = f"{D['annotation']}/{{sample}}/abricate_plasmidfinder.tsv",
    params:
        min_id  = config["abricate"]["min_id"],
        min_cov = config["abricate"]["min_cov"],
    log:       f"{D['logs']}/abricate_plasmidfinder/{{sample}}.log"
    benchmark: f"{D['benchmarks']}/abricate_plasmidfinder/{{sample}}.txt"
    shell:
        """
        mkdir -p $(dirname {log})
        abricate \
            --db plasmidfinder \
            --minid {params.min_id} \
            --mincov {params.min_cov} \
            {input.contigs} \
            > {output.tsv} \
            2> {log}
        """


# =============================================================================
# STEP 7 – HGT POTENTIAL + CLINICAL RISK SCORING  (R)
# =============================================================================
rule hgt_and_crs_scoring:
    """
    Merge ARG / virulence / plasmid annotations, compute proximity-based
    HGT Potential, and calculate the weighted Clinical Risk Score (CRS).
    """
    input:
        amrfinder     = expand("{d}/{sample}/amrfinder.tsv",
                               d=D["annotation"], sample=SAMPLES),
        vfdb          = expand("{d}/{sample}/abricate_vfdb.tsv",
                               d=D["annotation"], sample=SAMPLES),
        plasmidfinder = expand("{d}/{sample}/abricate_plasmidfinder.tsv",
                               d=D["annotation"], sample=SAMPLES),
        bracken       = expand("{d}/{sample}/{sample}.bracken.txt",
                               d=D["taxonomy"], sample=SAMPLES),
    output:
        hgt  = f"{D['scoring']}/hgt_potential.csv",
        risk = f"{D['scoring']}/final_risk_report.csv",
    params:
        ann_dir    = D["annotation"],
        tax_dir    = D["taxonomy"],
        out_dir    = D["scoring"],
        proximity  = config["crs"]["hgt_proximity_bp"],
        w_plasmid  = config["crs"]["plasmid_weight"],
        w_chrom    = config["crs"]["chrom_weight"],
        w_mdr      = config["crs"]["w_mdr"],
        w_hgt      = config["crs"]["w_hgt"],
        w_pathogen = config["crs"]["w_pathogen"],
        high_risk  = config["crs"]["high_risk"],
        low_risk   = config["crs"]["low_risk"],
        samples    = ",".join(SAMPLES),
    log:    f"{D['logs']}/scoring.log"
    script: "scripts/hgt_crs_scoring.R"


# =============================================================================
# STEP 8 – VISUALIZATIONS  (R)
# =============================================================================
rule visualize:
    """
    1. Bubble plot  – pathogen abundance x MDR index, bubble = HGT potential
    2. ARG heatmap  – sample x drug-class matrix, annotated with CRS tier
    """
    input:
        risk    = rules.hgt_and_crs_scoring.output.risk,
        hgt     = rules.hgt_and_crs_scoring.output.hgt,
        bracken = expand("{d}/{sample}/{sample}.bracken.txt",
                         d=D["taxonomy"], sample=SAMPLES),
        amr     = expand("{d}/{sample}/amrfinder.tsv",
                         d=D["annotation"], sample=SAMPLES),
    output:
        bubble  = f"{D['reports']}/crs_bubble_plot.pdf",
        heatmap = f"{D['reports']}/arg_heatmap.pdf",
    params:
        tax_dir = D["taxonomy"],
        ann_dir = D["annotation"],
        out_dir = D["reports"],
        samples = ",".join(SAMPLES),
    log:    f"{D['logs']}/visualize.log"
    script: "scripts/visualize.R"
