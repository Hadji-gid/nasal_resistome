# Nasal-Resistome-Mali

**Characterization of the Nasal Resistome and Virulome in Suspected Grippal Cases (Bamako, Mali)**
*MSc Bioinformatics Thesis – Alhadji A. Dicko | ACE-B / INSP | 2026*

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
├── Snakefile                   ← Main workflow (6 steps, starts from assembly)
├── config/
│   ├── config.yaml             ← All parameters, paths, and weights
│   └── samples.txt             ← One sample ID per line (79 samples)
├── envs/
│   ├── resistome.yaml          ← Conda env: SPAdes, Kraken2, AMRFinderPlus ...
│   └── r_scoring.yaml          ← Conda env: R + tidyverse, pheatmap, ggrepel
├── scripts/
│   ├── setup_databases.sh      ← One-time database download helper
│   ├── hgt_crs_scoring.R       ← HGT potential + Clinical Risk Score algorithm
│   └── visualize.R             ← Bubble plot + ARG heatmap (ggplot2, pheatmap)
└── resources/                  ← Databases (not tracked by git)
    ├── human_index/            ← GRCh38 Bowtie2 index (host depletion)
    ├── kraken2_db/             ← Kraken2 standard database
    └── amrfinderplus_db/       ← AMRFinderPlus / CARD database
```

---

## Pipeline Steps

| Step | Rule | Tool(s) | Output |
|------|------|---------|--------|
| 1 | `spades_assembly` | metaSPAdes `--meta` | `contigs.fasta` per sample |
| 2 | `quast` | QUAST | Assembly quality report (N50, # contigs) |
| 3 | `kraken2` + `bracken` | Kraken2, Bracken | Species abundance table |
| 4a | `amrfinder` | AMRFinderPlus (CARD) | ARG annotation table |
| 4b | `abricate_vfdb` | Abricate (VFDB) | Virulence factor table |
| 4c | `abricate_plasmidfinder` | Abricate (PlasmidFinder) | Plasmid replicon table |
| 5 | `hgt_and_crs_scoring` | R | `hgt_potential.csv` + `final_risk_report.csv` |
| 6 | `visualize` | R (ggplot2, pheatmap) | Bubble plot + ARG heatmap (PDF) |

> **Note:** Quality control (fastp) and host depletion (Bowtie2 vs GRCh38)
> were performed prior to this pipeline. Input files are host-depleted
> paired-end FASTQs named `{sample}_cleaned_R1/R2.fastq.gz`.

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

# Build databases (one-time, see script for instructions)
bash scripts/setup_databases.sh
```

### Run
```bash
# Dry run
snakemake --cores 8 --jobs 2 --dry-run

# Full run (79 samples, ~20-25 hrs on MacBook Pro)
snakemake --cores 8 --jobs 2 --rerun-incomplete \
          --rerun-triggers mtime --latency-wait 30
```

---

## Clinical Risk Score (CRS)

```
CRS = (0.4 × MDR_Index) + (0.3 × HGT_Potential) + (0.3 × Pathogen_Abundance)
```

All components normalized to [0, 1] across the cohort before weighting.

| Tier | CRS Range | Interpretation |
|------|-----------|----------------|
| **High** | > 0.7 | MDR pathogens with mobilizable resistance — alert clinician |
| **Moderate** | 0.4 – 0.7 | ARGs present with limited mobility |
| **Low** | < 0.4 | Predominantly commensal flora |

### HGT Potential Scoring

Per-ARG score based on four criteria:
- **Co-localization:** ARG and MGE on the same contig
- **Distance weighting:** ARG within 5 kb of an MGE (higher risk)
- **Vehicle weighting:** Plasmid-borne (1.0) vs chromosomal (0.6)
- **Pathogenicity flag:** Co-localized with high-priority virulence factor (+20%)

---

## Key Outputs

| File | Description |
|------|-------------|
| `results/07_scoring/final_risk_report.csv` | Per-sample CRS, risk tier, drug classes |
| `results/07_scoring/hgt_potential.csv` | Per-sample HGT scores and co-localization counts |
| `results/reports/crs_bubble_plot.pdf` | Bubble plot: pathogen abundance × MDR × HGT potential |
| `results/reports/arg_heatmap.pdf` | Sample × drug-class presence/absence heatmap |

---

## Study Cohort

- **Site:** Bamako, Republic of Mali
- **Samples:** 79 nasal swabs from suspected grippal (ILI) cases
- **Sequencing:** Illumina NextSeq (shotgun metagenomics)
- **Institution:** African Center of Excellence in Bioinformatics (ACE-B) / INSP

---

## Citation

If you use this pipeline, please cite:

> Dicko, A. A. (2026). *Characterization of the Nasal Resistome and Virulome
> in Suspected Grippal Cases (Bamako, Mali)*. Master's Thesis,
> Institut National de Santé Publique (INSP), Mali.

---

## Contact

**Alhadji A. Dicko** | alhadji-a.dicko@icermali.org
African Center of Excellence in Bioinformatics (ACE-B) / INSP, Mali

*Project Status: In Progress – 2026*
