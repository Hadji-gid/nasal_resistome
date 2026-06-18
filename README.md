# Nasal-Resistome-Mali

**Characterization of the Nasal Resistome and Virulome in Suspected Grippal Cases (Bamako, Mali)**
*MSc Bioinformatics Thesis – Alhadji A. Dicko | ACE-B / INSP Mali | 2026*

---

## Project Overview

This project provides an automated, reproducible Snakemake pipeline for analyzing shotgun metagenomic data from human nasal swabs collected from patients presenting with WHO-defined influenza-like illness (ILI) symptoms in Bamako, Mali. The study investigates:

- **Taxonomic Composition:** Species-level profiling of the nasal microbiota (Kraken2 + Bracken)
- **The Resistome:** Identification of ARGs via AMRFinderPlus (CARD database)
- **The Virulome:** Identification of virulence factors via Abricate (VFDB)
- **Mobilization Risk:** Proximity-based HGT potential (ARG–MGE co-localization on contigs)
- **Clinical Risk Score (CRS):** Weighted composite score stratifying patients by AMR transmission risk

---

## Cohort Summary (n=78)

| Metric | Value |
|---|---|
| Total samples analyzed | 78 |
| Samples with ARGs detected | 33 (42.3%) |
| Samples MDR (≥3 drug classes) | 24 (30.8%) |
| Unique drug classes detected | 17 |
| Samples with virulence factors | 20 (25.6%) |
| Total VF hits (VFDB) | 708 |
| High Risk (CRS >0.70) | 0 (0%) |
| Moderate Risk (CRS 0.40–0.70) | 2 (2.6%) |
| Low Risk (CRS <0.40) | 76 (97.4%) |
| Mean CRS | 0.089 |

---

## Repository Structure

```
nasal-resistome-mali/
├── Snakefile                        ← Main workflow (10 steps, v1.6)
├── config/
│   ├── config.yaml                  ← All parameters, paths, CRS weights
│   └── samples.txt                  ← 78 sample IDs
├── envs/
│   ├── resistome.yaml               ← Conda: fastp, Bowtie2, SPAdes, Kraken2...
│   └── r_scoring.yaml               ← Conda: R + tidyverse, pheatmap, ggrepel
├── scripts/
│   ├── setup_databases.sh           ← One-time database download helper
│   ├── hgt_crs_scoring.R            ← HGT potential + CRS algorithm (R)
│   ├── visualize.R                  ← Bubble plot + ARG heatmap (R)
│   └── build_custom_report.py       ← MultiQC custom content generator
├── multiqc_config.yaml              ← MultiQC configuration
├── outputs/
│   ├── final_risk_report.csv        ← Per-sample CRS results
│   ├── hgt_potential.csv            ← Per-sample HGT scores
│   └── cohort_stats.json            ← Aggregate cohort statistics
└── resources/                       ← Databases (not tracked — see .gitignore)
    ├── human_index/                 ← GRCh38 Bowtie2 index
    ├── kraken2_db/                  ← Kraken2 standard database
    └── amrfinderplus_db/            ← AMRFinderPlus / CARD database
```

---

## Pipeline Steps (v1.6)

| Step | Rule | Tool(s) | Output |
|------|------|---------|--------|
| 1 | `fastp_trim` | fastp | Clean paired reads + QC JSON |
| 2 | `multiqc` | MultiQC | Aggregated QC report |
| 3 | `host_depletion` | Bowtie2 + Samtools (GRCh38) | Microbial reads |
| 4 | `spades_assembly` | metaSPAdes (k=21,33,55,77) | `contigs.fasta` |
| 5 | `quast` | QUAST | Assembly quality report |
| 6 | `kraken2` | Kraken2 | Taxonomic classification |
| 7 | `bracken` | Bracken | Species abundance table |
| 8 | `amrfinder` | AMRFinderPlus (CARD) | ARG annotation |
| 9 | `abricate_vfdb` | Abricate (VFDB) | Virulence factor annotation |
| 10 | `abricate_plasmidfinder` | Abricate (PlasmidFinder) | Plasmid replicons |
| 11 | `hgt_and_crs_scoring` | R | HGT potential + CRS report |
| 12 | `visualize` | R (ggplot2, pheatmap) | Bubble plot + ARG heatmap |
| 13 | `multiqc_custom` | Python | 5 custom MultiQC modules |
| 14 | `multiqc_full` | MultiQC | Full aggregated HTML report |

### Disk Space Optimization

`temp()` files deleted automatically when no longer needed:

| File | Deleted after |
|------|--------------|
| Clean reads (`_R1/R2_clean.fastq.gz`) | `host_depletion` finishes |
| Microbial reads (`_R1/R2_microbial.fastq.gz`) | `spades_assembly` AND `kraken2` both finish |

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

conda env create -f envs/resistome.yaml
conda activate resistome2

# Build databases (one-time, ~4-6 hrs)
bash scripts/setup_databases.sh
```

### Configure
Edit `config/config.yaml`:
```yaml
host_index:   "resources/human_index/GRCh38"
kraken2_db:   "resources/kraken2_db"
amrfinder_db: "/path/to/amrfinderplus/data/latest"
bracken_bin:  "/path/to/conda/envs/resistome_env/bin/bracken"
```

Generate `config/samples.txt`:
```bash
ls data/raw/*_R1_001.fastq.gz \
  | xargs -n1 basename \
  | sed 's/_R1_001.fastq.gz//' \
  > config/samples.txt
```

### Run
```bash
# Dry run
snakemake --cores 8 --jobs 2 --dry-run

# Full run
snakemake --cores 8 --jobs 2 \
          --rerun-incomplete \
          --rerun-triggers mtime \
          --latency-wait 60
```

---

## Clinical Risk Score (CRS)

```
CRS = (0.4 × MDR_Index) + (0.3 × HGT_Potential) + (0.3 × Pathogen_Abundance)
```

Weights assigned exploratorily based on AMR surveillance literature
(Magiorakos et al., 2012; Tacconelli et al., 2018) as a proof-of-concept
framework subject to future validation.

| Tier | CRS | Interpretation |
|------|-----|----------------|
| **High** | >0.70 | MDR pathogens with mobilizable resistance |
| **Moderate** | 0.40–0.70 | ARGs present with limited mobility |
| **Low** | <0.40 | Predominantly commensal flora |

### HGT Potential Scoring

| Criterion | Description | Weight |
|-----------|-------------|--------|
| Co-localization | ARG and MGE on same contig | Required |
| Distance | Within 5 kb of MGE | Higher risk |
| Vehicle | Plasmid vs chromosomal | 1.0 vs 0.6 |
| Pathogenicity | Co-localized with VF | +20% boost |

---

## Key Outputs

| File | Description |
|------|-------------|
| `outputs/final_risk_report.csv` | Per-sample CRS, risk tier, drug classes |
| `outputs/hgt_potential.csv` | Per-sample HGT scores, co-localization counts |
| `outputs/cohort_stats.json` | Aggregate cohort statistics |
| `results/reports/multiqc_full_report.html` | Full MultiQC report |
| `results/reports/crs_bubble_plot.pdf` | Bubble plot: abundance × MDR × HGT |
| `results/reports/arg_heatmap.pdf` | Sample × drug-class heatmap |

---

## Known Limitations

- Virological confirmation of ILI was not systematically available; cases defined by WHO clinical criteria
- CRS weights are exploratory and require prospective clinical validation
- One sample (22_S61, 20M microbial reads) exceeded assembly capacity and is excluded from contig-level analyses
- Assembly quality varies with microbial biomass; 12 samples had <1,000 reads after host depletion

---

## Citation

> Dicko, A. A. (2026). *Characterization of the Nasal Resistome and Virulome in Suspected Grippal Cases (Bamako, Mali)*. Master's Thesis, African Center of Excellence in Bioinformatics (ACE-B), Institut National de Santé Publique (INSP), Mali.

## References

- Magiorakos et al. (2012). CMI. MDR definition.
- Tacconelli et al. (2018). Lancet ID. WHO priority pathogens.
- Partridge et al. (2018). CMR. Mobile genetic elements and AMR.
- Saurith-Coronell et al. (2026). ESPR. HGT potential in the airborne resistome.
- Feldgarden et al. (2021). Sci Rep. AMRFinderPlus.

---

## Contact

**Alhadji A. Dicko** | alhadji-a.dicko@icermali.org
African Center of Excellence in Bioinformatics (ACE-B) / INSP, Mali
*Project Status: MSc Thesis — In Progress 2026*
