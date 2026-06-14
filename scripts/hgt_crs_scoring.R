#!/usr/bin/env Rscript
# =============================================================================
# hgt_crs_scoring.R  v1.1
# HGT Potential scoring + Clinical Risk Score (CRS) calculation
# Author : Alhadji A. Dicko | ACE-B / INSP Mali
# Fix    : include all hgt_results columns in final join
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

# ── Helper: safe read of potentially empty TSVs ──────────────────────────────
safe_read <- function(path, ...) {
  if (!file.exists(path) || file.size(path) == 0) return(NULL)
  tryCatch(read_tsv(path, show_col_types = FALSE, ...), error = function(e) NULL)
}

# ── Parameters from Snakemake ────────────────────────────────────────────────
ann_dir    <- snakemake@params$ann_dir
tax_dir    <- snakemake@params$tax_dir
out_dir    <- snakemake@params$out_dir
proximity  <- as.integer(snakemake@params$proximity)
w_plasmid  <- as.numeric(snakemake@params$w_plasmid)
w_chrom    <- as.numeric(snakemake@params$w_chrom)
w_mdr      <- as.numeric(snakemake@params$w_mdr)
w_hgt      <- as.numeric(snakemake@params$w_hgt)
w_pathogen <- as.numeric(snakemake@params$w_pathogen)
high_thresh<- as.numeric(snakemake@params$high_risk)
low_thresh <- as.numeric(snakemake@params$low_risk)
samples    <- strsplit(snakemake@params$samples, ",")[[1]]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# SECTION 1 – HGT POTENTIAL SCORING
# =============================================================================
message("[HGT Scoring] Processing ", length(samples), " samples ...")

hgt_results <- map_dfr(samples, function(s) {
  message("  -> ", s)

  amr  <- safe_read(file.path(ann_dir, s, "amrfinder.tsv"))
  vfdb <- safe_read(file.path(ann_dir, s, "abricate_vfdb.tsv"))
  plas <- safe_read(file.path(ann_dir, s, "abricate_plasmidfinder.tsv"))

  if (is.null(amr) || nrow(amr) == 0) {
    return(tibble(
      sample = s, n_args = 0L, n_mges = 0L,
      n_colocalized = 0L, n_proximal = 0L,
      n_plasmid_borne = 0L, n_vf_colocalized = 0L,
      hgt_raw = 0, hgt_potential = 0
    ))
  }

  # Standardize AMRFinder column names
  amr <- amr %>%
    rename_with(~ gsub(" ", "_", .x)) %>%
    rename_with(tolower)

  # Find contig, start, stop, class columns flexibly
  contig_col <- grep("contig", names(amr), value = TRUE)[1]
  start_col  <- grep("^start$", names(amr), value = TRUE)[1]
  stop_col   <- grep("^stop$", names(amr), value = TRUE)[1]
  class_col  <- grep("^class$", names(amr), value = TRUE)[1]
  if (is.na(class_col)) class_col <- grep("class", names(amr), value = TRUE)[1]

  if (is.na(contig_col)) {
    return(tibble(
      sample = s, n_args = nrow(amr), n_mges = 0L,
      n_colocalized = 0L, n_proximal = 0L,
      n_plasmid_borne = 0L, n_vf_colocalized = 0L,
      hgt_raw = 0, hgt_potential = 0
    ))
  }

  amr <- amr %>%
    rename(contig_id = all_of(contig_col)) %>%
    mutate(
      start   = as.integer(.data[[start_col]]),
      stop    = as.integer(.data[[stop_col]]),
      arg_mid = (start + stop) / 2
    )

  # Plasmid-positive contigs
  plasmid_contigs <- character(0)
  mge_positions   <- tibble(contig_id = character(), mge_mid = numeric())
  if (!is.null(plas) && nrow(plas) > 0) {
    plas <- plas %>% rename_with(tolower)
    seq_col   <- grep("sequence|contig", names(plas), value = TRUE)[1]
    start_p   <- grep("^start$", names(plas), value = TRUE)[1]
    end_p     <- grep("^end$|^stop$", names(plas), value = TRUE)[1]
    if (!is.na(seq_col)) {
      plasmid_contigs <- unique(plas[[seq_col]])
      if (!is.na(start_p) && !is.na(end_p)) {
        mge_positions <- plas %>%
          mutate(mge_mid = (as.integer(.data[[start_p]]) + as.integer(.data[[end_p]])) / 2) %>%
          select(contig_id = all_of(seq_col), mge_mid)
      }
    }
  }

  # VF contigs
  vf_contigs <- character(0)
  if (!is.null(vfdb) && nrow(vfdb) > 0) {
    vfdb <- vfdb %>% rename_with(tolower)
    seq_col_v <- grep("sequence|contig", names(vfdb), value = TRUE)[1]
    if (!is.na(seq_col_v)) vf_contigs <- unique(vfdb[[seq_col_v]])
  }

  n_args        <- nrow(amr)
  n_mges        <- nrow(mge_positions)
  n_colocalized <- 0L
  n_proximal    <- 0L
  n_plasmid_borne <- 0L
  n_vf_coloc    <- 0L
  hgt_score_sum <- 0

  for (i in seq_len(nrow(amr))) {
    arg_contig <- amr$contig_id[i]
    arg_pos    <- amr$arg_mid[i]

    vehicle_w <- if (arg_contig %in% plasmid_contigs) w_plasmid else w_chrom
    if (arg_contig %in% plasmid_contigs) n_plasmid_borne <- n_plasmid_borne + 1L

    mge_same <- mge_positions %>% filter(contig_id == arg_contig)
    colocalized <- nrow(mge_same) > 0
    if (colocalized) n_colocalized <- n_colocalized + 1L

    dist_w <- 0
    if (colocalized) {
      min_dist <- min(abs(mge_same$mge_mid - arg_pos))
      dist_w   <- if (min_dist <= proximity) 1.0 else max(0, 1 - min_dist / (proximity * 2))
      if (min_dist <= proximity) n_proximal <- n_proximal + 1L
    }

    vf_flag <- if (arg_contig %in% vf_contigs) 1.2 else 1.0
    if (arg_contig %in% vf_contigs) n_vf_coloc <- n_vf_coloc + 1L

    hgt_score_sum <- hgt_score_sum + vehicle_w * dist_w * vf_flag
  }

  max_possible  <- n_args * w_plasmid * 1.0 * 1.2
  hgt_potential <- if (max_possible > 0) min(1, hgt_score_sum / max_possible) else 0

  tibble(
    sample          = s,
    n_args          = n_args,
    n_mges          = n_mges,
    n_colocalized   = n_colocalized,
    n_proximal      = n_proximal,
    n_plasmid_borne = n_plasmid_borne,
    n_vf_colocalized= n_vf_coloc,
    hgt_raw         = round(hgt_score_sum, 4),
    hgt_potential   = round(hgt_potential, 4)
  )
})

write_csv(hgt_results, snakemake@output$hgt)
message("[HGT Scoring] Written to ", snakemake@output$hgt)

# =============================================================================
# SECTION 2 – MDR INDEX
# =============================================================================
message("[MDR Index] Calculating ...")

mdr_results <- map_dfr(samples, function(s) {
  amr <- safe_read(file.path(ann_dir, s, "amrfinder.tsv"))

  if (is.null(amr) || nrow(amr) == 0) {
    return(tibble(sample = s, n_drug_classes = 0L, mdr_flag = FALSE,
                  drug_classes = "", mdr_index_raw = 0))
  }

  amr <- amr %>%
    rename_with(~ gsub(" ", "_", .x)) %>%
    rename_with(tolower)

  class_col <- grep("^class$", names(amr), value = TRUE)[1]
  if (is.na(class_col)) class_col <- grep("class", names(amr), value = TRUE)[1]

  classes   <- unique(na.omit(amr[[class_col]]))
  n_classes <- length(classes)

  tibble(
    sample         = s,
    n_drug_classes = n_classes,
    mdr_flag       = n_classes >= 3,
    drug_classes   = paste(sort(classes), collapse = "; "),
    mdr_index_raw  = n_classes
  )
})

# =============================================================================
# SECTION 3 – PATHOGEN ABUNDANCE
# =============================================================================
message("[Pathogen Abundance] Calculating ...")

PRIORITY_PATHOGENS <- c(
  "Staphylococcus aureus", "Streptococcus pneumoniae",
  "Haemophilus influenzae", "Moraxella catarrhalis",
  "Klebsiella pneumoniae", "Acinetobacter baumannii",
  "Pseudomonas aeruginosa", "Escherichia coli",
  "Enterococcus faecalis", "Enterococcus faecium",
  "Staphylococcus epidermidis"
)

pathogen_results <- map_dfr(samples, function(s) {
  bracken_path <- file.path(tax_dir, s, paste0(s, ".bracken.txt"))
  brac <- safe_read(bracken_path)

  if (is.null(brac) || nrow(brac) == 0) {
    return(tibble(sample = s, pathogen_abundance_raw = 0,
                  top_pathogen = NA_character_, n_priority_taxa = 0L))
  }

  brac <- brac %>% rename_with(tolower)
  frac_col <- grep("fraction", names(brac), value = TRUE)[1]

  path_rows <- brac %>%
    filter(name %in% PRIORITY_PATHOGENS) %>%
    arrange(desc(.data[[frac_col]]))

  tibble(
    sample                = s,
    pathogen_abundance_raw= round(sum(path_rows[[frac_col]], na.rm = TRUE), 4),
    top_pathogen          = if (nrow(path_rows) > 0) path_rows$name[1] else NA_character_,
    n_priority_taxa       = nrow(path_rows)
  )
})

# =============================================================================
# SECTION 4 – NORMALISE & COMPUTE CRS
# =============================================================================
message("[CRS] Computing final scores ...")

normalize_01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(0, length(x)))
  (x - rng[1]) / diff(rng)
}

# Merge ALL hgt columns (not just hgt_potential)
combined <- mdr_results %>%
  left_join(hgt_results, by = "sample") %>%
  left_join(pathogen_results, by = "sample") %>%
  mutate(
    mdr_norm      = normalize_01(mdr_index_raw),
    hgt_norm      = hgt_potential,
    pathogen_norm = normalize_01(pathogen_abundance_raw),
    crs = round(
      w_mdr      * mdr_norm +
      w_hgt      * hgt_norm +
      w_pathogen * pathogen_norm,
      4
    ),
    risk_tier = case_when(
      crs > high_thresh ~ "High",
      crs >= low_thresh ~ "Moderate",
      TRUE              ~ "Low"
    ),
    risk_tier = factor(risk_tier, levels = c("High", "Moderate", "Low"))
  ) %>%
  arrange(desc(crs))

# Write final report — all columns now available
final_report <- combined %>%
  select(
    sample, crs, risk_tier,
    mdr_norm, n_drug_classes, mdr_flag, drug_classes,
    hgt_norm = hgt_potential,
    n_args, n_mges, n_colocalized, n_proximal,
    n_plasmid_borne, n_vf_colocalized, hgt_raw,
    pathogen_norm, pathogen_abundance_raw,
    top_pathogen, n_priority_taxa
  )

write_csv(final_report, snakemake@output$risk)
message("[CRS] Written to ", snakemake@output$risk)

message("\n=== Cohort Summary ===")
message("Total samples : ", nrow(final_report))
message("High Risk     : ", sum(final_report$risk_tier == "High"))
message("Moderate Risk : ", sum(final_report$risk_tier == "Moderate"))
message("Low Risk      : ", sum(final_report$risk_tier == "Low"))
message("Mean CRS      : ", round(mean(final_report$crs), 3))
message("======================\n")
