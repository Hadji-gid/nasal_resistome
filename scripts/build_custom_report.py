#!/usr/bin/env python3
"""
build_custom_report.py  v1.1
Generates MultiQC-compatible custom content modules:
  1. Resistome & Taxonomy Summary table (per sample)
  2. Top Taxa Abundance heatmap (Bracken)
  3. Drug Class Presence/Absence heatmap (AMRFinderPlus)
  4. Virulome Summary table (Abricate VFDB)
  5. VF Functional Category heatmap (per sample x category)

Fix v1.1: PRODUCT column is col[13], parse VF category from bracket notation.

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
import re
import os
import json
import glob
from pathlib import Path
from collections import Counter, defaultdict

PRIORITY_PATHOGENS = {
    "Staphylococcus aureus", "Streptococcus pneumoniae",
    "Haemophilus influenzae", "Moraxella catarrhalis",
    "Klebsiella pneumoniae", "Acinetobacter baumannii",
    "Pseudomonas aeruginosa", "Escherichia coli",
    "Enterococcus faecalis", "Enterococcus faecium",
    "Staphylococcus epidermidis"
}

# VFDB functional category codes -> readable names
VFC_NAMES = {
    "VFC0001": "Adherence",
    "VFC0204": "Motility",
    "VFC0258": "Immune modulation",
    "VFC0272": "Nutritional/Metabolic",
    "VFC0282": "Stress survival",
    "VFC0325": "Antimicrobial activity",
    "VFC0001": "Adherence",
    "VFC0200": "Invasion",
    "VFC0203": "Exotoxin",
    "VFC0205": "Exoenzyme",
    "VFC0206": "Iron uptake",
}


def parse_vfdb_product(product_str):
    """
    Parse VFDB PRODUCT field like:
    '(flhA) flagellar biosynthesis protein FlhA [Flagella (VF1400) - Motility (VFC0204)] [Enterobacter...]'
    Returns: (vf_name, vf_category, vfc_code, source_org)
    Note: Abricate truncates at ~80 chars so source_org may be cut off.
    """
    if not product_str or not product_str.strip():
        return "Unknown", "Unknown", "Unknown", "Unknown"

    # Extract VF name from first bracket e.g. [Flagella (VF1400) - Motility (VFC0204)]
    vf_name = "Unknown"
    vf_category = "Unknown"
    vfc_code = "Unknown"
    source_org = "Unknown"

    # Find VF name: text before " - " in first bracket group
    m = re.search(r'\[([^\[\]]+)\s*-\s*([^\[\]]+)\]', product_str)
    if m:
        vf_name = m.group(1).strip()
        # Remove VF ID code e.g. "(VF1400)"
        vf_name = re.sub(r'\s*\(VF\d+\)', '', vf_name).strip()
        cat_part = m.group(2).strip()
        # Extract VFC code
        vfc_m = re.search(r'\(VFC(\d+)\)', cat_part)
        if vfc_m:
            vfc_code = f"VFC{vfc_m.group(1)}"
        # Category name = text before the VFC code
        vf_category = re.sub(r'\s*\(VFC\d+\)', '', cat_part).strip()

    # Source organism: last bracket content (often truncated)
    brackets = re.findall(r'\[([^\[\]]+)\]', product_str)
    if len(brackets) >= 2:
        source_org = brackets[-1].strip()
    elif len(brackets) == 1:
        source_org = brackets[0].strip()

    return vf_name, vf_category, vfc_code, source_org


def read_vfdb(path):
    """Return list of dicts with VF info per hit."""
    hits = []
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return hits
    with open(path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 14:
                continue
            gene = cols[5]
            product = cols[13] if len(cols) > 13 else ""
            pct_cov = cols[9] if len(cols) > 9 else "0"
            pct_id = cols[10] if len(cols) > 10 else "0"
            vf_name, vf_cat, vfc_code, source_org = parse_vfdb_product(product)
            hits.append({
                "gene": gene,
                "vf_name": vf_name,
                "vf_category": vf_cat,
                "vfc_code": vfc_code,
                "source_org": source_org,
                "pct_cov": pct_cov,
                "pct_id": pct_id,
            })
    return hits


def read_bracken(path):
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
    return rows, (rows[0][0] if rows else None)


def read_amrfinder(path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return 0, set()
    classes = set()
    n = 0
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            return 0, set()
        class_col = next((c for c in reader.fieldnames
                          if c.strip().lower() == "class"), None)
        if class_col is None:
            class_col = next((c for c in reader.fieldnames
                              if "class" in c.lower()), None)
        for row in reader:
            n += 1
            if class_col and row.get(class_col):
                classes.add(row[class_col].strip())
    return n, classes


def read_risk_report(path):
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
        samples = [s.strip() for s in fh
                   if s.strip() and not s.startswith("#")]

    risk     = read_risk_report(args.risk_report)
    tax_dir  = Path(args.taxonomy_dir)
    ann_dir  = Path(args.annotation_dir)
    outdir   = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    summary_rows = []
    taxon_abundance   = defaultdict(dict)
    drug_class_pres   = defaultdict(dict)
    vf_category_pres  = defaultdict(dict)  # category -> {sample: count}
    vf_org_counter    = Counter()
    vf_gene_counter   = Counter()
    vf_cat_counter    = Counter()

    for s in samples:
        bracken_path = tax_dir / s / f"{s}.bracken.txt"
        amr_path     = ann_dir / s / "amrfinder.tsv"
        vfdb_path    = ann_dir / s / "abricate_vfdb.tsv"

        tax_rows, top_taxon = read_bracken(bracken_path)
        n_args, classes     = read_amrfinder(amr_path)
        vf_hits             = read_vfdb(vfdb_path)

        priority_hits = [(t, f) for t, f in tax_rows if t in PRIORITY_PATHOGENS]
        priority_str  = "; ".join(
            f"{t} ({f*100:.1f}%)" for t, f in priority_hits[:3]) or "None"

        # Taxonomy heatmap (top 5 per sample)
        for t, f in tax_rows[:5]:
            taxon_abundance[t][s] = round(f * 100, 3)

        # Drug class heatmap
        for c in classes:
            drug_class_pres[c][s] = 1

        # VF stats
        n_vf = len(vf_hits)
        vf_cats_this = Counter()
        vf_orgs_this = Counter()
        for h in vf_hits:
            vf_gene_counter[h["gene"]] += 1
            vf_cat_counter[h["vf_category"]] += 1
            vf_org_counter[h["source_org"]] += 1
            vf_cats_this[h["vf_category"]] += 1
            vf_orgs_this[h["source_org"]] += 1
            vf_category_pres[h["vf_category"]][s] = \
                vf_category_pres[h["vf_category"]].get(s, 0) + 1

        top_vf_cat = vf_cats_this.most_common(1)[0][0] \
            if vf_cats_this else "None"
        top_vf_org = vf_orgs_this.most_common(1)[0][0] \
            if vf_orgs_this else "None"

        r = risk.get(s, {})
        summary_rows.append({
            "Sample":                    s,
            "CRS":                       r.get("crs", "NA"),
            "Risk_Tier":                 r.get("risk_tier", "NA"),
            "Top_Taxon":                 top_taxon or "N/A",
            "Priority_Pathogens":        priority_str,
            "N_ARGs":                    n_args,
            "N_Drug_Classes":            len(classes),
            "Drug_Classes":              "; ".join(sorted(classes)) or "None",
            "N_Virulence_Factors":       n_vf,
            "Top_VF_Category":           top_vf_cat,
            "Top_VF_Source_Organism":    top_vf_org,
        })

    # ── 1. Master summary table ───────────────────────────────────────────────
    summary_path = outdir / "sample_summary_mqc.tsv"
    header_lines = [
        "# id: 'nasal_resistome_summary'",
        "# section_name: 'Resistome, Virulome & Taxonomy Summary'",
        "# description: 'Per-sample ARG counts, drug classes, virulence factors, top pathogens and Clinical Risk Score'",
        "# plot_type: 'table'",
        "# pconfig:",
        "#     id: 'resistome_summary_table'",
        "#     namespace: 'Cohort Summary'",
    ]
    with open(summary_path, "w", newline="") as fh:
        fh.write("\n".join(header_lines) + "\n")
        writer = csv.DictWriter(fh, fieldnames=summary_rows[0].keys(),
                                delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"Wrote {summary_path}")

    # ── 2. Taxonomy heatmap ───────────────────────────────────────────────────
    all_top_taxa = sorted(taxon_abundance.keys(),
                          key=lambda t: -sum(taxon_abundance[t].values()))[:15]
    taxa_path = outdir / "taxonomy_abundance_mqc.tsv"
    with open(taxa_path, "w") as fh:
        fh.write("# id: 'nasal_resistome_taxonomy'\n")
        fh.write("# section_name: 'Top Taxa Relative Abundance (%)'\n")
        fh.write("# description: 'Top 15 most abundant taxa (Bracken fraction_total_reads x100)'\n")
        fh.write("# plot_type: 'heatmap'\n")
        fh.write("Taxon\t" + "\t".join(samples) + "\n")
        for t in all_top_taxa:
            row = [str(taxon_abundance[t].get(s, 0)) for s in samples]
            fh.write(f"{t}\t" + "\t".join(row) + "\n")
    print(f"Wrote {taxa_path}")

    # ── 3. Drug class heatmap ─────────────────────────────────────────────────
    dc_path = outdir / "drug_class_matrix_mqc.tsv"
    all_classes = sorted(drug_class_pres.keys())
    with open(dc_path, "w") as fh:
        fh.write("# id: 'nasal_resistome_drugclass'\n")
        fh.write("# section_name: 'ARG Drug Class Presence/Absence'\n")
        fh.write("# description: 'Antibiotic drug classes with detected resistance genes (AMRFinderPlus/CARD)'\n")
        fh.write("# plot_type: 'heatmap'\n")
        fh.write("Drug_Class\t" + "\t".join(samples) + "\n")
        for c in all_classes:
            row = [str(drug_class_pres[c].get(s, 0)) for s in samples]
            fh.write(f"{c}\t" + "\t".join(row) + "\n")
    print(f"Wrote {dc_path}")

    # ── 4. Virulome summary table ─────────────────────────────────────────────
    vf_summary_path = outdir / "virulome_summary_mqc.tsv"
    vf_header = [
        "# id: 'nasal_resistome_virulome'",
        "# section_name: 'Virulome Summary (VFDB)'",
        "# description: 'Per-sample virulence factor hits, top functional category and source organism'",
        "# plot_type: 'table'",
        "# pconfig:",
        "#     id: 'virulome_table'",
        "#     namespace: 'Virulome'",
    ]
    vf_rows = [r for r in summary_rows if r["N_Virulence_Factors"] > 0]
    vf_cols = ["Sample", "N_Virulence_Factors", "Top_VF_Category",
               "Top_VF_Source_Organism", "CRS", "Risk_Tier"]
    with open(vf_summary_path, "w", newline="") as fh:
        fh.write("\n".join(vf_header) + "\n")
        writer = csv.DictWriter(fh, fieldnames=vf_cols,
                                delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(vf_rows)
    print(f"Wrote {vf_summary_path}")

    # ── 5. VF functional category heatmap ────────────────────────────────────
    vf_cat_path = outdir / "vf_category_matrix_mqc.tsv"
    all_vf_cats = sorted(vf_category_pres.keys(),
                         key=lambda c: -sum(vf_category_pres[c].values()))
    with open(vf_cat_path, "w") as fh:
        fh.write("# id: 'nasal_resistome_vf_categories'\n")
        fh.write("# section_name: 'Virulence Factor Functional Categories'\n")
        fh.write("# description: 'Number of VF hits per functional category per sample (VFDB)'\n")
        fh.write("# plot_type: 'heatmap'\n")
        fh.write("VF_Category\t" + "\t".join(samples) + "\n")
        for cat in all_vf_cats:
            row = [str(vf_category_pres[cat].get(s, 0)) for s in samples]
            fh.write(f"{cat}\t" + "\t".join(row) + "\n")
    print(f"Wrote {vf_cat_path}")

    # ── 6. Cohort stats ───────────────────────────────────────────────────────
    cohort_stats = {
        "total_samples":           len(samples),
        "samples_with_args":       sum(1 for r in summary_rows if r["N_ARGs"] > 0),
        "samples_mdr":             sum(1 for r in summary_rows if r["N_Drug_Classes"] >= 3),
        "total_unique_drug_classes": len(all_classes),
        "samples_with_vf":         sum(1 for r in summary_rows if r["N_Virulence_Factors"] > 0),
        "total_vf_hits":           sum(r["N_Virulence_Factors"] for r in summary_rows),
        "top_vf_categories":       [c for c, _ in vf_cat_counter.most_common(5)],
        "top_vf_genes":            [g for g, _ in vf_gene_counter.most_common(10)],
        "high_risk":               sum(1 for r in summary_rows if r["Risk_Tier"] == "High"),
        "moderate_risk":           sum(1 for r in summary_rows if r["Risk_Tier"] == "Moderate"),
        "low_risk":                sum(1 for r in summary_rows if r["Risk_Tier"] == "Low"),
    }
    stats_path = outdir / "cohort_stats.json"
    with open(stats_path, "w") as fh:
        json.dump(cohort_stats, fh, indent=2)
    print(f"\nCohort stats:")
    for k, v in cohort_stats.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
