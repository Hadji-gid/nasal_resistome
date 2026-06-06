#!/usr/bin/env bash
# =============================================================================
# setup_databases.sh
# One-time database download and index building script.
# Run ONCE before executing the Snakemake pipeline.
#
# Usage:  bash scripts/setup_databases.sh
# =============================================================================
set -euo pipefail

THREADS=8
RESOURCES="resources"

echo "========================================"
echo " Nasal Resistome Mali – Database Setup"
echo "========================================"

mkdir -p ${RESOURCES}/{human_index,kraken2_db,amrfinderplus_db}

# ── 1. Human genome index (GRCh38) ──────────────────────────────────────────
echo ""
echo "[1/3] Building Bowtie2 index for GRCh38 ..."
echo "      Download GRCh38 FASTA first if not present:"
echo "      wget -P ${RESOURCES} https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/"
echo ""
echo "      Then run:"
echo "      bowtie2-build --threads ${THREADS} \\"
echo "        ${RESOURCES}/GRCh38.fa \\"
echo "        ${RESOURCES}/human_index/GRCh38"
echo ""

# ── 2. Kraken2 + Bracken standard database ──────────────────────────────────
echo "[2/3] Building Kraken2 standard database (~60 GB, ~2–4 hrs) ..."
echo "      kraken2-build --standard \\"
echo "        --db ${RESOURCES}/kraken2_db \\"
echo "        --threads ${THREADS}"
echo ""
echo "      Then build Bracken database (read length must match your data):"
echo "      bracken-build \\"
echo "        -d ${RESOURCES}/kraken2_db \\"
echo "        -t ${THREADS} \\"
echo "        -l 150"
echo ""

# ── 3. AMRFinderPlus database ────────────────────────────────────────────────
echo "[3/3] Downloading AMRFinderPlus database ..."
amrfinder --update --database ${RESOURCES}/amrfinderplus_db \
  && echo "      AMRFinderPlus database updated successfully." \
  || echo "      ERROR: amrfinder update failed. Install amrfinderplus first."

echo ""
echo "========================================"
echo " Database setup instructions complete."
echo " Update paths in config/config.yaml"
echo " before running the pipeline."
echo "========================================"
