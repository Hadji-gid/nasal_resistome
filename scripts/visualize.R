#!/usr/bin/env Rscript
# =============================================================================
# visualize.R  v1.1
# Generates two publication-quality figures:
#   1. Bubble plot  – taxon abundance vs MDR index, bubble = HGT potential
#   2. ARG heatmap  – sample x drug-class matrix annotated with CRS tier
#
# Author : Alhadji A. Dicko | ACE-B / INSP Mali
# Fix    : flexible drug_class column detection for AMRFinder output
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(scales)
})

# ── Parameters ────────────────────────────────────────────────────────────────
risk_path <- snakemake@input$risk
tax_dir   <- snakemake@params$tax_dir
ann_dir   <- snakemake@params$ann_dir
out_dir   <- snakemake@params$out_dir
samples   <- strsplit(snakemake@params$samples, ",")[[1]]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

safe_read <- function(path, ...) {
  if (!file.exists(path) || file.size(path) == 0) return(NULL)
  tryCatch(read_tsv(path, show_col_types = FALSE, ...), error = function(e) NULL)
}

# ── Load risk report ──────────────────────────────────────────────────────────
risk <- read_csv(risk_path, show_col_types = FALSE) %>%
  mutate(risk_tier = factor(risk_tier, levels = c("High", "Moderate", "Low")))

TIER_COLORS <- c(High = "#C0392B", Moderate = "#E67E22", Low = "#27AE60")

# =============================================================================
# FIGURE 1 – BUBBLE PLOT
# =============================================================================
message("[Viz] Creating bubble plot ...")

bubble_data <- risk %>%
  mutate(
    bubble_size = rescale(hgt_norm, to = c(2, 14)),
    label = if_else(risk_tier == "High", sample, NA_character_)
  )

p_bubble <- ggplot(bubble_data,
    aes(x = pathogen_abundance_raw, y = n_drug_classes,
        size = bubble_size, color = risk_tier, label = label)) +
  geom_point(alpha = 0.75, stroke = 0.4) +
  geom_label_repel(
    size = 2.8, fontface = "bold", max.overlaps = 20,
    box.padding = 0.4, point.padding = 0.3,
    show.legend = FALSE, na.rm = TRUE
  ) +
  scale_size_identity() +
  scale_color_manual(values = TIER_COLORS, name = "CRS Risk Tier") +
  scale_x_continuous(labels = percent_format(accuracy = 0.1),
                     name = "Priority Pathogen Relative Abundance") +
  scale_y_continuous(breaks = pretty_breaks(),
                     name = "MDR Index (# Antibiotic Drug Classes)") +
  labs(
    title    = "Clinical Risk Landscape - Nasal Resistome (Bamako, Mali, n=78)",
    subtitle = "Bubble size = HGT Potential  |  Color = CRS Risk Tier  |  Labels on High-Risk samples",
    caption  = "Dicko AA (2026). MSc Bioinformatics Thesis, ACE-B / INSP Mali."
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9.5, color = "grey40"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave(snakemake@output$bubble, p_bubble,
       width = 11, height = 7.5, device = "pdf")
message("[Viz] Bubble plot saved.")

# =============================================================================
# FIGURE 2 – ARG HEATMAP
# =============================================================================
message("[Viz] Creating ARG heatmap ...")

# Build sample x drug-class matrix from AMRFinder outputs
amr_list <- map(samples, function(s) {
  f  <- file.path(ann_dir, s, "amrfinder.tsv")
  df <- safe_read(f)
  if (is.null(df) || nrow(df) == 0) return(NULL)

  # Standardize column names
  df <- df %>%
    rename_with(~ gsub(" ", "_", .x)) %>%
    rename_with(tolower)

  # Find the drug class column flexibly (exact match preferred)
  class_col <- grep("^class$", names(df), value = TRUE)[1]
  if (is.na(class_col)) class_col <- grep("^drug_class$", names(df), value = TRUE)[1]
  if (is.na(class_col)) class_col <- grep("class", names(df), value = TRUE)[1]

  if (is.na(class_col)) return(NULL)

  df %>%
    select(drug_class = all_of(class_col)) %>%
    filter(!is.na(drug_class), drug_class != "") %>%
    distinct() %>%
    mutate(sample = s, present = 1L)
}) %>% bind_rows()

if (!is.null(amr_list) && nrow(amr_list) > 0) {
  # Wide matrix: rows = drug_class, cols = sample
  mat <- amr_list %>%
    pivot_wider(names_from = sample, values_from = present,
                values_fill = 0L) %>%
    column_to_rownames("drug_class") %>%
    as.matrix()

  # Order columns by descending CRS
  sample_order <- risk %>% arrange(desc(crs)) %>% pull(sample)
  sample_order <- intersect(sample_order, colnames(mat))
  mat <- mat[, sample_order, drop = FALSE]

  # Annotation data frame
  ann_col <- risk %>%
    filter(sample %in% sample_order) %>%
    arrange(match(sample, sample_order)) %>%
    select(sample, `Risk Tier` = risk_tier) %>%
    column_to_rownames("sample")

  ann_colors <- list(`Risk Tier` = TIER_COLORS)

  pdf(snakemake@output$heatmap,
      width = max(12, ncol(mat) * 0.18 + 4), height = 8)
  pheatmap(
    mat,
    color             = c("#F7F7F7", "#2166AC"),
    breaks            = c(-0.5, 0.5, 1.5),
    annotation_col    = ann_col,
    annotation_colors = ann_colors,
    cluster_rows      = TRUE,
    cluster_cols      = FALSE,
    fontsize_row      = 9,
    fontsize_col      = 6,
    border_color      = "grey85",
    main              = "ARG Drug-Class Matrix (Bamako Cohort, n=78)\nOrdered by CRS (High to Low)",
    legend_breaks     = c(0, 1),
    legend_labels     = c("Absent", "Present"),
    silent            = TRUE
  )
  dev.off()
  message("[Viz] ARG heatmap saved.")

} else {
  # Placeholder if no ARGs detected
  pdf(snakemake@output$heatmap)
  plot.new()
  text(0.5, 0.5, "No ARGs detected across cohort.", cex = 1.5)
  dev.off()
  message("[Viz] No ARG data - placeholder heatmap written.")
}
