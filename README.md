# Nasal-Resistome-Mali

**Characterization of the Nasal Resistome and Virulome in Suspected Grippal Cases (Bamako, Mali)**
*MSc Bioinformatics Thesis – Alhadji A. Dicko | ACE-B / INSP Mali | 2026*

---

## Project Overview

This project provides an automated, reproducible Snakemake pipeline for analyzing
shotgun metagenomic data from human nasal swabs. The study focuses on a cohort of
**79 samples from Bamako, Mali**, specifically investigating:

- **Taxonomic Composition:** Species-level profiling of the nasal microbiota (Kraken2 + Bracken)
- **The Resistome:** Identification of Antimicrobial Resistance Genes (ARGs) via AMRFinderPlus (CARD)
- **The Virulome:** Identification of virulence factors via Abricate (VFDB)
- **Mobilization Risk:** Proximity-based HGT potential (ARG–MGE co-localization on assembled contigs)
- **Clinical Risk Score (CRS):** Weighted composite score to stratify patients by AMR transmission risk

---

## Repository Structure

```
nasal-resistome-mali/
├── Snakefile                      ← Main workflow (11 steps, raw reads → CRS)
├── config/
│   ├── config.yaml                ← All parameters, paths, and CRS weights
│   └── samples.txt                ← One sample ID per line (79 samples)
├── envs/
│   ├── resistome.yaml             ← Conda env: fastp, Bowtie2, SPAdes, Kraken2 ...
│   └── r_scoring.yaml             ← Conda env: R + tidyverse, pheatmap, ggrepel
├── scripts/
│   ├── setup_databases.sh         ← One-time database download helper
│   ├── hgt_crs_scoring.R          ← HGT potential + Clinical Risk Score algorithm
│   └── visualize.R                ← Bubble plot + ARG heatmap (ggplot2, pheatmap)
└── resources/                     ← Databases (not tracked by git – see .gitignore)
    ├── human_index/               ← GRCh38 Bowtie2 index (host depletion)
    ├── kraken2_db/                ← Kraken2 standard database
    └── amrfinderplus_db/          ← AMRFinderPlus / CARD database
```

---

## Pipeline Steps

| Step | Rule | Tool(s) | Output |
|------|------|---------|--------|
| 1 | `fastp_trim` | fastp | Trimmed reads + per-sample QC report |
| 2 | `multiqc` | MultiQC | Aggregated QC report (all samples) |
| 3 | `host_depletion` | Bowtie2, Samtools (vs GRCh38) | Microbial reads (human DNA removed) |
| 4 | `spades_assembly` | metaSPAdes `--meta` | `contigs.fasta` per sample |
| 5 | `quast` | QUAST | Assembly quality report (N50, # contigs) |
| 6 | `kraken2` | Kraken2 | Read-level taxonomic classification |
| 7 | `bracken` | Bracken | Species-level abundance table |
| 8 | `amrfinder` | AMRFinderPlus (CARD) | ARG annotation table |
| 9 | `abricate_vfdb` | Abricate (VFDB) | Virulence factor table |
| 10 | `abricate_plasmidfinder` | Abricate (PlasmidFinder) | Plasmid replicon table |
| 11 | `hgt_and_crs_scoring` | R | `hgt_potential.csv` + `final_risk_report.csv` |
| 12 | `visualize` | R (ggplot2, pheatmap) | Bubble plot + ARG heatmap (PDF) |

### Disk space optimization
Large intermediate files are automatically deleted via Snakemake `temp()`:

| File | Deleted after |
|------|--------------|
| Clean reads (`_R1/R2_clean.fastq.gz`) | `host_depletion` finishes |
| Microbial reads (`_R1/R2_microbial.fastq.gz`) | `spades_assembly` AND `kraken2` both finish |
| Kraken2 raw output (`_kraken2_output.txt`) | `bracken` finishes |
| SPAdes intermediate files | Deleted inside assembly rule |

---

## Quick Start

### Prerequisites
- macOS or Linux
- Conda / Mamba
- Snakemake ≥ 7.0

### Setup
```bash
git clone https://github.com/Hadji-gid/nasal_resistome.git
cd nasal_resistome

# Create main environment
conda env create -f envs/resistome.yaml
conda activate resistome2

# Build databases (one-time, ~4-6 hrs)
bash scripts/setup_databases.sh
```

### Configure
Edit `config/config.yaml` to set your database paths:
```yaml
host_index:   "resources/human_index/GRCh38"
kraken2_db:   "resources/kraken2_db"
amrfinder_db: "/path/to/amrfinderplus/data/latest"
bracken_bin:  "/path/to/conda/envs/resistome_env/bin/bracken"
```

Add your sample IDs to `config/samples.txt` (one per line):
```bash
ls data/raw/*_R1_001.fastq.gz \
  | xargs -n1 basename \
  | sed 's/_R1_001.fastq.gz//' \
  > config/samples.txt
```

### Run
```bash
# Dry run first
snakemake --cores 8 --jobs 2 --dry-run

# Full run (79 samples, ~20-25 hrs on MacBook Pro i7)
snakemake --cores 8 --jobs 2 \
          --rerun-incomplete \
          --rerun-triggers mtime \
          --latency-wait 30
```

---

## Clinical Risk Score (CRS)

```
CRS = (0.4 × MDR_Index) + (0.3 × HGT_Potential) + (0.3 × Pathogen_Abundance)
```

All components normalized to [0, 1] across the cohort before weighting.

| Tier | CRS Range | Interpretation | Action |
|------|-----------|----------------|--------|
| **High** | > 0.7 | MDR pathogens with mobilizable resistance | Alert clinician; stewardship intervention |
| **Moderate** | 0.4 – 0.7 | ARGs present with limited mobility | Standard monitoring |
| **Low** | < 0.4 | Predominantly commensal flora | No immediate action |

### HGT Potential Scoring

Per-ARG score based on four criteria applied at contig level:

| Criterion | Description | Weight |
|-----------|-------------|--------|
| Co-localization | ARG and MGE on the same contig | Required |
| Distance | ARG within 5 kb of an MGE | Higher risk |
| Vehicle | Plasmid-borne vs chromosomal | 1.0 vs 0.6 |
| Pathogenicity | Co-localized with high-priority VF | +20% boost |

---

## Key Outputs

| File | Description |
|------|-------------|
| `results/07_scoring/final_risk_report.csv` | Per-sample CRS, risk tier, drug classes |
| `results/07_scoring/hgt_potential.csv` | Per-sample HGT scores and co-localization counts |
| `results/reports/crs_bubble_plot.pdf` | Bubble plot: pathogen abundance × MDR × HGT |
| `results/reports/arg_heatmap.pdf` | Sample × drug-class presence/absence heatmap |
| `results/01_qc/multiqc_report.html` | Aggregated QC report for all samples |

---

## Study Cohort

- **Site:** Bamako, Republic of Mali
- **Samples:** 79 nasal swabs from suspected grippal (ILI) cases
- **Sequencing:** Illumina NextSeq (shotgun metagenomics, paired-end)
- **Institution:** African Center of Excellence in Bioinformatics (ACE-B) / INSP, Mali

---

## Citation

If you use this pipeline, please cite:

> Dicko, A. A. (2026). *Characterization of the Nasal Resistome and Virulome
> in Suspected Grippal Cases (Bamako, Mali)*. Master's Thesis,
> African Center of Excellence in Bioinformatics (ACE-B),
> Institut National de Santé Publique (INSP), Mali.

---

## Contact

**Alhadji A. Dicko** | alhadji-a.dicko@icermali.org
African Center of Excellence in Bioinformatics (ACE-B) / INSP, Mali

*Project Status: In Progress – 2026*
