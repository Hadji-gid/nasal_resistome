# Nasal-Resistome-Mali Pipeline

**Characterization of the Nasal Resistome and Virulome in Suspected Grippal Cases (Bamako, Mali)**
*MSc Bioinformatics Thesis – Alhadji A. Dicko | ACE-Bamako / INSP | 2026*

---

## Repository Structure

```
nasal-resistome-mali/
├── Snakefile                   ← Main workflow (all 8 steps)
├── config/
│   ├── config.yaml             ← All parameters and paths
│   └── samples.txt             ← One sample ID per line
├── envs/
│   ├── resistome.yaml          ← Conda env: fastp, bowtie2, megahit, kraken2 …
│   └── r_scoring.yaml          ← Conda env: R + tidyverse, pheatmap …
├── scripts/
│   ├── setup_databases.sh      ← One-time DB download helper
│   ├── hgt_crs_scoring.R       ← HGT potential + Clinical Risk Score
│   └── visualize.R             ← Bubble plot + ARG heatmap
├── data/
│   └── raw/                    ← Place raw FASTQs here: {sample}_R1/R2.fastq.gz
└── resources/                  ← Databases built by setup_databases.sh
    ├── human_index/
    ├── kraken2_db/
    └── amrfinderplus_db/
```

---

## Quick Start

### Step 1 – Clone and set up environment

```bash
git clone https://github.com/hadji/nasal-resistome-mali.git
cd nasal-resistome-mali

# Install Mamba (faster than conda) if not already installed
conda install -n base -c conda-forge mamba

# Create environments (Snakemake will do this automatically with --use-conda,
# but you can pre-build them manually)
mamba env create -f envs/resistome.yaml
mamba env create -f envs/r_scoring.yaml
```

### Step 2 – Build databases (one-time, ~4–6 hrs)

```bash
# Activate the main environment first
conda activate resistome_env

# Follow the instructions printed by this script:
bash scripts/setup_databases.sh
```

Update all paths in `config/config.yaml` to match where databases were built.

### Step 3 – Prepare your samples

1. Place all raw paired FASTQ files in `data/raw/`:
   - `data/raw/SAMPLE_001_R1.fastq.gz`
   - `data/raw/SAMPLE_001_R2.fastq.gz`
2. Edit `config/samples.txt` — one sample ID per line (no header, no extensions).

### Step 4 – Run the pipeline

```bash
# Dry run first (checks everything without running)
snakemake --use-conda --cores 8 --dry-run

# Full run on 83 samples
snakemake --use-conda --cores 8

# Run on a cluster (SLURM example)
snakemake --use-conda --cores 64 \
  --cluster "sbatch -c {threads} --mem=32G -t 12:00:00" \
  --jobs 20
```

---

## Pipeline Steps

| Step | Rule | Tool(s) | Output |
|------|------|---------|--------|
| 1 | `fastp_trim` | fastp, MultiQC | Clean reads + QC report |
| 2 | `host_depletion` | Bowtie2, Samtools | Microbial reads (human DNA removed) |
| 3 | `megahit_assembly` | MEGAHIT | Assembled contigs (N50 ≥ 5,000 bp target) |
| 4 | `quast` | QUAST | Assembly quality report |
| 5 | `kraken2` + `bracken` | Kraken2, Bracken | Species abundance table |
| 6a | `amrfinder` | AMRFinderPlus (CARD) | ARG annotation table |
| 6b | `abricate_vfdb` | Abricate (VFDB) | Virulence factor table |
| 6c | `abricate_plasmidfinder` | Abricate (PlasmidFinder) | Plasmid replicon table |
| 7 | `hgt_and_crs_scoring` | R | `hgt_potential.csv` + `final_risk_report.csv` |
| 8 | `visualize` | R (ggplot2, pheatmap) | Bubble plot + ARG heatmap (PDF) |

---

## Clinical Risk Score Formula

```
CRS = (0.4 × MDR_Index) + (0.3 × HGT_Potential) + (0.3 × Pathogen_Abundance)
```

| Tier | CRS Range | Interpretation |
|------|-----------|----------------|
| **High** | > 0.7 | MDR pathogens with mobilizable resistance — alert clinician |
| **Moderate** | 0.4 – 0.7 | ARGs present with limited mobility |
| **Low** | < 0.4 | Predominantly commensal flora |

---

## Key Outputs

| File | Description |
|------|-------------|
| `results/07_scoring/final_risk_report.csv` | Per-sample CRS, risk tier, drug classes |
| `results/07_scoring/hgt_potential.csv` | Per-sample HGT scores and co-localization counts |
| `results/reports/crs_bubble_plot.pdf` | Bubble plot: abundance × MDR × HGT |
| `results/reports/arg_heatmap.pdf` | Sample × drug-class presence/absence heatmap |
| `results/01_qc/multiqc_report.html` | Aggregated QC report for all 83 samples |

---

## Citation

If you use this pipeline, please cite:

> Dicko, A. A. (2026). *Characterization of the Nasal Resistome and Virulome
> in Suspected Grippal Cases (Bamako, Mali)*. Master's Thesis,
> Institut National de Santé Publique (INSP), Mali.

---

## Contact

**Alhadji A. Dicko** | alhadji-a.dicko@icermali.org
International Center for Excellence in Research (ICER-Mali) / INRSP, Mali
