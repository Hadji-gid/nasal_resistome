# =============================================================================
# Nasal Resistome Mali – Snakemake Workflow
# Author  : Alhadji A. Dicko | ACE-B / INSP Mali
# Version : 1.6  (2026)
#
# Full pipeline from raw reads to Clinical Risk Score + MultiQC report:
#   1.  fastp              – QC and adapter trimming
#   2.  multiqc            – Aggregate QC reports
#   3.  host_depletion     – Human DNA removal (Bowtie2 vs GRCh38)
#   4.  spades_assembly    – Metagenomic assembly (metaSPAdes)
#   5.  quast              – Assembly quality assessment (non-fatal)
#   6.  kraken2            – Taxonomic classification
#   7.  bracken            – Species abundance re-estimation
#   8.  amrfinder          – ARG annotation (CARD)
#   9.  abricate_vfdb      – Virulence factor annotation (VFDB)
#   10. abricate_plasmidfinder – Plasmid detection (PlasmidFinder)
#   11. hgt_and_crs_scoring    – HGT potential + Clinical Risk Score (R)
#   12. visualize              – Bubble plot + ARG heatmap (R)
#   13. multiqc_custom         – Custom MultiQC content (taxonomy + resistome + virulome)
#   14. multiqc_full           – Full aggregated MultiQC report
#
# Key fixes in v1.6 vs v1.5:
#   - kraken2 output no longer temp() — prevents MissingOutputException
#   - kraken2 and bracken non-fatal for zero-read samples (|| true + touch)
#   - quast non-fatal for low-biomass samples (|| true + placeholder)
#   - SPAdes: k=21,33,55,77; memory=28GB; checks min reads before assembly
#   - host_depletion: || true on samtools pipes; output size validation
#   - All rules: mkdir -p for log and output dirs
#   - bracken uses explicit binary path from config (bracken_bin)
#   - Added multiqc_custom and multiqc_full rules
#
# Usage (with resistome2 conda env active):
#   snakemake --cores 8 --jobs 2 --rerun-incomplete \
#             --rerun-triggers mtime --latency-wait 60
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
        expand("{d}/{sample}/contigs.fasta",               d=D["assembly"],  sample=SAMPLES),
        # QUAST
        expand("{d}/{sample}/report.tsv",                  d=D["asm_qc"],    sample=SAMPLES),
        # Taxonomy
        expand("{d}/{sample}/{sample}.bracken.txt",        d=D["taxonomy"],  sample=SAMPLES),
        # Annotation
        expand("{d}/{sample}/amrfinder.tsv",               d=D["annotation"],sample=SAMPLES),
        expand("{d}/{sample}/abricate_vfdb.tsv",           d=D["annotation"],sample=SAMPLES),
        expand("{d}/{sample}/abricate_plasmidfinder.tsv",  d=D["annotation"],sample=SAMPLES),
        # Scoring
        f"{D['scoring']}/final_risk_report.csv",
        f"{D['scoring']}/hgt_potential.csv",
        # Visualizations
        f"{D['reports']}/crs_bubble_plot.pdf",
        f"{D['reports']}/arg_heatmap.pdf",
        # MultiQC full report
        f"{D['reports']}/multiqc_full_report.html",


# =============================================================================
# STEP 1 – QUALITY CONTROL  (fastp + MultiQC)
# =============================================================================
rule fastp_trim:
    """
    Adapter trimming and quality filtering.
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
    Aggregate all fastp JSON reports into one MultiQC QC summary.
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
    Uses || true on samtools pipes to handle SIGPIPE gracefully.
    Output size validation prevents silent empty file production.
    Microbial FASTQs are temp() — deleted once spades + kraken2 both finish.
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

        # 1. Align to human genome → coordinate-sorted BAM
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

        # 2. Alignment stats (% human reads removed)
        samtools index {params.tmpbam}
        samtools flagstat {params.tmpbam} > {output.stats}

        # 3. Extract unmapped pairs only (-f 12 = both reads unmapped)
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

        # 4. Validate outputs exist (allow empty for zero-microbe samples)
        [ -f {output.r1} ] || touch {output.r1}
        [ -f {output.r2} ] || touch {output.r2}

        # 5. Cleanup BAM
        rm -f {params.tmpbam} {params.tmpbam}.bai
        """


# =============================================================================
# STEP 3 – METAGENOMIC ASSEMBLY  (metaSPAdes)
# =============================================================================
rule spades_assembly:
    """
    De-novo metagenomic assembly. Skips assembly for samples with <1000 reads
    (low-biomass / viral-only samples) and creates a placeholder contig.
    k-mers limited to 21,33,55,77 — appropriate for 150 bp reads and
    significantly faster than default (avoids k=99,127 which caused timeouts).
    Memory set to 28 GB to handle large samples (up to 20M read pairs).
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

        # Check minimum reads before attempting assembly
        READ_COUNT=$(zcat {input.r1} 2>/dev/null | wc -l | awk '{{print int($1/4)}}')

        if [ "$READ_COUNT" -lt 1000 ]; then
            echo "SKIP: only $READ_COUNT reads — insufficient for assembly" \
                > {log}
            mkdir -p {params.outdir}
            printf ">low_biomass_placeholder\\nNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\\n" \
                > {output.contigs}
        else
            spades.py \
                --meta \
                -1 {input.r1} \
                -2 {input.r2} \
                -o {params.outdir} \
                --threads {threads} \
                --memory {params.mem} \
                -k 21,33,55,77 \
                2> {log}
        fi
        """


# =============================================================================
# STEP 4 – ASSEMBLY QUALITY ASSESSMENT  (QUAST)
# =============================================================================
rule quast:
    """
    Assembly QC: N50, # contigs, largest contig.
    Non-fatal: creates placeholder report for low-biomass samples.
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

        # Create placeholder if QUAST failed (low-biomass sample)
        [ -f {output.report} ] || \
            printf "Assembly\\t{wildcards.sample}\\nNote\\tInsufficient microbial reads\\n" \
            > {output.report}
        """


# =============================================================================
# STEP 5 – TAXONOMIC PROFILING  (Kraken2 + Bracken)
# =============================================================================
rule kraken2:
    """
    Classify microbial reads against Kraken2 standard database.
    Non-fatal for zero-read samples (creates empty output files).
    Note: output NOT marked as temp() to prevent MissingOutputException
    on macOS filesystem with --latency-wait.
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
    Non-fatal for empty Kraken2 reports (zero-read samples).
    Uses explicit binary path from config (bracken_bin) — required because
    bracken may not be in the active conda environment PATH.
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
    Detect ARGs in assembled contigs (nucleotide mode).
    Uses CARD database via AMRFinderPlus.
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
    Detects adherence, immune modulation, motility, toxin, and iron
    acquisition factors. PRODUCT column [13] contains VF category info.
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
    Detect plasmid replicons using PlasmidFinder.
    Plasmid-positive contigs receive higher HGT vehicle weight (1.0)
    in the CRS calculation.
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
    Merge ARG / virulence / plasmid annotations per sample.
    Compute proximity-based HGT Potential (ARG-MGE co-localization, <5 kb,
    vehicle-weighted, VF-boosted).
    Compute weighted CRS = 0.4*MDR + 0.3*HGT + 0.3*Pathogen_Abundance.
    Stratify into High/Moderate/Low risk tiers.
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
    Figure 1: Bubble plot — pathogen abundance x MDR index, bubble = HGT potential
    Figure 2: ARG heatmap — sample x drug-class presence/absence, ordered by CRS
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


# =============================================================================
# STEP 9 – MULTIQC CUSTOM CONTENT (taxonomy + resistome + virulome)
# =============================================================================
rule multiqc_custom:
    """
    Generate 5 MultiQC custom content modules:
      1. Resistome, Virulome & Taxonomy Summary table
      2. Top 15 Taxa Abundance heatmap (Bracken)
      3. ARG Drug Class Presence/Absence heatmap (AMRFinderPlus)
      4. Virulome Summary table (VFDB)
      5. VF Functional Category heatmap (VFDB)
    Also writes cohort_stats.json with key aggregate statistics.
    """
    input:
        amrfinder = expand("{d}/{sample}/amrfinder.tsv",
                           d=D["annotation"], sample=SAMPLES),
        vfdb      = expand("{d}/{sample}/abricate_vfdb.tsv",
                           d=D["annotation"], sample=SAMPLES),
        bracken   = expand("{d}/{sample}/{sample}.bracken.txt",
                           d=D["taxonomy"], sample=SAMPLES),
        risk      = rules.hgt_and_crs_scoring.output.risk,
    output:
        summary   = f"{D['reports']}/multiqc_custom_data/sample_summary_mqc.tsv",
        taxonomy  = f"{D['reports']}/multiqc_custom_data/taxonomy_abundance_mqc.tsv",
        drugclass = f"{D['reports']}/multiqc_custom_data/drug_class_matrix_mqc.tsv",
        virulome  = f"{D['reports']}/multiqc_custom_data/virulome_summary_mqc.tsv",
        vf_cats   = f"{D['reports']}/multiqc_custom_data/vf_category_matrix_mqc.tsv",
        stats     = f"{D['reports']}/multiqc_custom_data/cohort_stats.json",
    params:
        tax_dir    = D["taxonomy"],
        ann_dir    = D["annotation"],
        risk_path  = rules.hgt_and_crs_scoring.output.risk,
        outdir     = f"{D['reports']}/multiqc_custom_data",
        samples    = config["samples_list"],
    log: f"{D['logs']}/multiqc_custom.log"
    shell:
        """
        mkdir -p {params.outdir} $(dirname {log})
        python3 scripts/build_custom_report.py \
            --samples {params.samples} \
            --taxonomy-dir {params.tax_dir} \
            --annotation-dir {params.ann_dir} \
            --risk-report {params.risk_path} \
            --outdir {params.outdir} \
            2> {log}
        """


# =============================================================================
# STEP 10 – FULL MULTIQC REPORT
# =============================================================================
rule multiqc_full:
    """
    Aggregate all tool outputs + custom content into one MultiQC HTML report.
    Includes: fastp QC, QUAST assembly stats, Kraken2 taxonomy,
    plus custom resistome/virulome/taxonomy sections.
    """
    input:
        rules.multiqc_custom.output.summary,
        rules.multiqc_custom.output.taxonomy,
        rules.multiqc_custom.output.drugclass,
        rules.multiqc_custom.output.virulome,
        rules.multiqc_custom.output.vf_cats,
        rules.visualize.output.bubble,
    output:
        f"{D['reports']}/multiqc_full_report.html"
    params:
        qc_dir      = D["qc"],
        asm_qc_dir  = D["asm_qc"],
        tax_dir     = D["taxonomy"],
        custom_dir  = f"{D['reports']}/multiqc_custom_data",
        outdir      = D["reports"],
        config      = "multiqc_config.yaml",
    log: f"{D['logs']}/multiqc_full.log"
    shell:
        """
        mkdir -p $(dirname {log})
        multiqc \
            {params.qc_dir} \
            {params.asm_qc_dir} \
            {params.tax_dir} \
            {params.custom_dir} \
            --config {params.config} \
            --outdir {params.outdir} \
            --filename multiqc_full_report \
            --force \
            2> {log}
        """
