# =============================================================================
# Nasal Resistome Mali – Snakemake Workflow
# Author  : Alhadji A. Dicko | ICER-Mali / INRSP
# Version : 1.4  (2026)
#
# Changes in v1.4:
#   - Removed fastp and host_depletion rules (already completed externally)
#   - Pipeline starts directly from metaSPAdes assembly
#   - Input: results/02_host_depletion/{sample}_cleaned_R1/R2.fastq.gz
#   - Assembler: metaSPAdes (replaces MEGAHIT)
#   - Bracken uses explicit binary path from config
#
# Usage:
#   snakemake --cores 8 --jobs 2 --rerun-incomplete
# =============================================================================

from pathlib import Path

configfile: "config/config.yaml"

with open(config["samples_list"]) as fh:
    SAMPLES = [s.strip() for s in fh if s.strip() and not s.startswith("#")]

D           = config["outdir"]
DEPL        = D["host_depl"]
BRACKEN_BIN = config["bracken_bin"]

# ── Target rule ───────────────────────────────────────────────────────────────
rule all:
    input:
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
# STEP 1 – METAGENOMIC ASSEMBLY  (metaSPAdes)
# Input: pre-processed host-depleted reads (flat directory)
# =============================================================================
rule spades_assembly:
    input:
        r1 = f"{DEPL}/{{sample}}_cleaned_R1.fastq.gz",
        r2 = f"{DEPL}/{{sample}}_cleaned_R2.fastq.gz",
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

        spades.py \
            --meta \
            -1 {input.r1} \
            -2 {input.r2} \
            -o {params.outdir} \
            --threads {threads} \
            --memory {params.mem} \
            2> {log}
        """


# =============================================================================
# STEP 2 – ASSEMBLY QUALITY ASSESSMENT  (QUAST)
# =============================================================================
rule quast:
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

        # Create empty report if quast failed (low biomass sample)
        [ -f {output.report} ] || echo -e "Assembly\t{wildcards.sample}\nNote\tNo contigs >= {params.min_contig}bp" > {output.report}
        """


# =============================================================================
# STEP 3 – TAXONOMIC PROFILING  (Kraken2 + Bracken)
# =============================================================================
rule kraken2:
    input:
        r1 = f"{DEPL}/{{sample}}_cleaned_R1.fastq.gz",
        r2 = f"{DEPL}/{{sample}}_cleaned_R2.fastq.gz",
    output:
        report = f"{D['taxonomy']}/{{sample}}/{{sample}}_kraken2_report.txt",
        out    = temp(f"{D['taxonomy']}/{{sample}}/{{sample}}_kraken2_output.txt"),
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
            2> {log}
        """

rule bracken:
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
            2> {log}
        """


# =============================================================================
# STEP 4 – FUNCTIONAL ANNOTATION
# =============================================================================
rule amrfinder:
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
# STEP 5 – HGT POTENTIAL + CLINICAL RISK SCORING  (R)
# =============================================================================
rule hgt_and_crs_scoring:
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
# STEP 6 – VISUALIZATIONS  (R)
# =============================================================================
rule visualize:
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
