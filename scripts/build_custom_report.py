#!/usr/bin/env python3
"""
build_custom_report.py
Generates a MultiQC-compatible custom content module summarizing:
  1. Taxonomic profiling (top pathogens per sample, from Bracken)
  2. Resistome annotation (ARG counts, drug classes, from AMRFinderPlus)
  3. Virulome annotation (VF counts, from Abricate)
  4. Clinical Risk Score summary

Output: results/reports/multiqc_custom_data/ (consumed automatically by MultiQC
        when placed alongside other search paths, OR run standalone as an
        HTML/TSV summary if MultiQC custom content parsing isn't desired).

Usage:
    python3 build_custom_report.py \
        --samples config/samples.txt \
        --taxonomy-dir results/05_taxonomy \
        --annotation-dir results/06_annotation \
        --risk-report results/07_scoring/final_risk_report.csv \
        --outdir results/reports/multiqc_custom_data
"""

import argparse
import csv
import json
import os
from pathlib import Path
from collections import defaultdict

PRIORITY_PATHOGENS = {
    "Staphylococcus aureus", "Streptococcus pneumoniae",
    "Haemophilus influenzae", "Moraxella catarrhalis",
    "Klebsiella pneumoniae", "Acinetobacter baumannii",
    "Pseudomonas aeruginosa", "Escherichia coli",
    "Enterococcus faecalis", "Enterococcus faecium",
    "Staphylococcus epidermidis"
}


def read_bracken(path):
    """Return list of (taxon, fraction) sorted descending, plus top taxon."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return [], None
    rows = []
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            try:
                frac = float(row.get("fraction_total_reads", 0))
            except (ValueError, TypeError):
                continue
            rows.append((row.get("name", "Unknown"), frac))
    rows.sort(key=lambda x: x[1], reverse=True)
    top = rows[0][0] if rows else None
    return rows, top


def read_amrfinder(path):
    """Return (n_args, set_of_drug_classes)."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return 0, set()
    classes = set()
    n = 0
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        # Find the class column flexibly
        if reader.fieldnames is None:
            return 0, set()
        class_col = next((c for c in reader.fieldnames if c.strip().lower() == "class"), None)
        if class_col is None:
            class_col = next((c for c in reader.fieldnames if "class" in c.lower()), None)
        for row in reader:
            n += 1
            if class_col and row.get(class_col):
                classes.add(row[class_col].strip())
    return n, classes


def read_abricate(path):
    """Return number of hits in an abricate TSV."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return 0
    n = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            n += 1
    return n


def read_risk_report(path):
    """Return dict sample -> row dict."""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            out[row["sample"]] = row
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True)
    ap.add_argument("--taxonomy-dir", required=True)
    ap.add_argument("--annotation-dir", required=True)
    ap.add_argument("--risk-report", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    with open(args.samples) as fh:
        samples = [s.strip() for s in fh if s.strip() and not s.startswith("#")]

    risk = read_risk_report(args.risk_report)

    tax_dir = Path(args.taxonomy_dir)
    ann_dir = Path(args.annotation_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    summary_rows = []
    taxon_abundance = defaultdict(dict)   # taxon -> {sample: fraction}
    drug_class_presence = defaultdict(dict)  # class -> {sample: 0/1}

    for s in samples:
        bracken_path = tax_dir / s / f"{s}.bracken.txt"
        amr_path = ann_dir / s / "amrfinder.tsv"
        vfdb_path = ann_dir / s / "abricate_vfdb.tsv"
        plasmid_path = ann_dir / s / "abricate_plasmidfinder.tsv"

        tax_rows, top_taxon = read_bracken(bracken_path)
        n_args, classes = read_amrfinder(amr_path)
        n_vf = read_abricate(vfdb_path)
        n_plasmid = read_abricate(plasmid_path)

        priority_hits = [(t, f) for t, f in tax_rows if t in PRIORITY_PATHOGENS]
        priority_str = "; ".join(f"{t} ({f*100:.1f}%)" for t, f in priority_hits[:3])

        for t, f in tax_rows[:5]:
            taxon_abundance[t][s] = round(f * 100, 3)

        for c in classes:
            drug_class_presence[c][s] = 1

        r = risk.get(s, {})
        summary_rows.append({
            "Sample": s,
            "Top_Taxon": top_taxon or "N/A",
            "Priority_Pathogens_Detected": priority_str or "None",
            "N_ARGs": n_args,
            "N_Drug_Classes": len(classes),
            "Drug_Classes": "; ".join(sorted(classes)) or "None",
            "N_Virulence_Factors": n_vf,
            "N_Plasmid_Hits": n_plasmid,
            "CRS": r.get("crs", "NA"),
            "Risk_Tier": r.get("risk_tier", "NA"),
        })

    # ── Write master summary TSV ──────────────────────────────────────────────
    summary_tsv = outdir / "sample_summary_mqc.tsv"
    with open(summary_tsv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=summary_rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"Wrote {summary_tsv}")

    # ── MultiQC custom content header for general stats table ────────────────
    header_lines = [
        "# id: 'nasal_resistome_summary'",
        "# section_name: 'Resistome & Taxonomy Summary'",
        "# description: 'Per-sample ARG counts, drug classes, top pathogens, and Clinical Risk Score'",
        "# plot_type: 'table'",
        "# pconfig:",
        "#     id: 'resistome_summary_table'",
        "#     namespace: 'Cohort Summary'",
    ]
    with open(summary_tsv, "r") as fh:
        body = fh.read()
    with open(outdir / "sample_summary_mqc.tsv", "w") as fh:
        fh.write("\n".join(header_lines) + "\n" + body)

    # ── Top taxa abundance matrix (for heatmap-style custom content) ──────────
    all_top_taxa = sorted(taxon_abundance.keys(),
                           key=lambda t: -sum(taxon_abundance[t].values()))[:15]
    taxa_matrix_path = outdir / "taxonomy_abundance_mqc.tsv"
    with open(taxa_matrix_path, "w") as fh:
        fh.write("# id: 'nasal_resistome_taxonomy'\n")
        fh.write("# section_name: 'Top Taxa Relative Abundance (%)'\n")
        fh.write("# description: 'Top 15 most abundant taxa across cohort (Bracken fraction_total_reads x100)'\n")
        fh.write("# plot_type: 'heatmap'\n")
        fh.write("Taxon\t" + "\t".join(samples) + "\n")
        for t in all_top_taxa:
            row = [str(taxon_abundance[t].get(s, 0)) for s in samples]
            fh.write(f"{t}\t" + "\t".join(row) + "\n")
    print(f"Wrote {taxa_matrix_path}")

    # ── Drug class presence/absence matrix ─────────────────────────────────────
    drugclass_matrix_path = outdir / "drug_class_matrix_mqc.tsv"
    all_classes = sorted(drug_class_presence.keys())
    with open(drugclass_matrix_path, "w") as fh:
        fh.write("# id: 'nasal_resistome_drugclass'\n")
        fh.write("# section_name: 'ARG Drug Class Presence/Absence'\n")
        fh.write("# description: 'Antibiotic drug classes with detected resistance genes (AMRFinderPlus/CARD)'\n")
        fh.write("# plot_type: 'heatmap'\n")
        fh.write("Drug_Class\t" + "\t".join(samples) + "\n")
        for c in all_classes:
            row = [str(drug_class_presence[c].get(s, 0)) for s in samples]
            fh.write(f"{c}\t" + "\t".join(row) + "\n")
    print(f"Wrote {drugclass_matrix_path}")

    # ── Cohort-level stats JSON ────────────────────────────────────────────────
    cohort_stats = {
        "total_samples": len(samples),
        "samples_with_args": sum(1 for r in summary_rows if r["N_ARGs"] > 0),
        "samples_mdr": sum(1 for r in summary_rows if r["N_Drug_Classes"] >= 3),
        "total_unique_drug_classes": len(all_classes),
        "high_risk": sum(1 for r in summary_rows if r["Risk_Tier"] == "High"),
        "moderate_risk": sum(1 for r in summary_rows if r["Risk_Tier"] == "Moderate"),
        "low_risk": sum(1 for r in summary_rows if r["Risk_Tier"] == "Low"),
    }
    with open(outdir / "cohort_stats.json", "w") as fh:
        json.dump(cohort_stats, fh, indent=2)
    print(f"Wrote cohort_stats.json: {cohort_stats}")


if __name__ == "__main__":
    main()
