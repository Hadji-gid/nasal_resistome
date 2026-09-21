#!/usr/bin/env Rscript
# =============================================================================
# publication_figures.R  --  v2.1
#
# Publication-quality figures and tables for:
#   "Characterization of the Nasal Resistome and Virulome in Suspected
#    Grippal Cases (Bamako, Mali)"
#
# Author : Alhadji A. Dicko | ACE-B / INSP Mali | 2026
#
# CHANGES IN v2.1
#   - Captions are wrapped with wrap_cap() so they never run off the page.
#   - cap_base shortened (~30 characters saved on every figure).
#   - Bar-label clipping fixed: wider y-axis expansion on all coord_flip()
#     panels that carry outside labels.
#   - plot.caption gains lineheight so wrapped captions are legible.
#
# CHANGES IN v2.0
#   - Removed ComplexHeatmap (Bioconductor, not CRAN). Heatmaps are built
#     with ggplot2: fully vector and journal-safe.
#   - CRAN mirror set explicitly at the top of the script.
#   - Every figure exported as PDF (vector), TIFF (600 dpi LZW) and PNG.
#
# USAGE (run from the project root, ~/nasal_resistome):
#   conda activate resistome_r
#   Rscript scripts/publication_figures.R
#
# OUTPUT: results/figures_publication/
# =============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

# ---- Dependencies (CRAN only) ----------------------------------------------
required <- c("tidyverse", "patchwork", "viridis", "scales",
              "cowplot", "ggrepel", "forcats")

for (p in required) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message("Installing ", p, " ...")
    install.packages(p, quiet = TRUE)
  }
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(viridis)
  library(scales)
  library(cowplot)
  library(ggrepel)
  library(forcats)
})

# =============================================================================
# CONFIGURATION
# =============================================================================
PROJECT  <- getwd()
RISK_CSV <- file.path(PROJECT, "results/07_scoring/final_risk_report.csv")
HGT_CSV  <- file.path(PROJECT, "results/07_scoring/hgt_potential.csv")
ANN_DIR  <- file.path(PROJECT, "results/06_annotation")
TAX_DIR  <- file.path(PROJECT, "results/05_taxonomy")
SAMPLES  <- file.path(PROJECT, "config/samples.txt")
OUT      <- file.path(PROJECT, "results/figures_publication")

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- Palette (colour-blind safe, print safe) -------------------------------
PAL <- list(
  navy = "#1F4E79", blue = "#2E75B6", teal = "#0F6E56", lteal = "#1D9E75",
  green = "#3B6D11", amber = "#BA7517", orange = "#E67E22", red = "#C0392B",
  purple = "#6A4C93", grey = "#7F7F7F", lgrey = "#E8ECEF"
)

TIER_COL <- c(High = PAL$red, Moderate = PAL$orange, Low = PAL$lteal)

DC_COL <- c(
  "BETA-LACTAM"                                = "#C0392B",
  "TETRACYCLINE"                               = "#E67E22",
  "LINCOSAMIDE/MACROLIDE/STREPTOGRAMIN"        = "#6A4C93",
  "AMINOGLYCOSIDE"                             = "#2E75B6",
  "QUATERNARY AMMONIUM"                        = "#0F6E56",
  "MACROLIDE"                                  = "#D35400",
  "MACROLIDE/STREPTOGRAMIN"                    = "#E74C3C",
  "FOSFOMYCIN"                                 = "#1D9E75",
  "LINCOSAMIDE"                                = "#9B59B6",
  "TRIMETHOPRIM"                               = "#16A085",
  "PHENICOL"                                   = "#F39C12",
  "SULFONAMIDE"                                = "#3498DB",
  "STREPTOTHRICIN"                             = "#2ECC71",
  "NITROFURAN/PHENICOL/QUINOLONE/TETRACYCLINE" = "#AD1457",
  "LINCOSAMIDE/STREPTOGRAMIN"                  = "#795548",
  "AMINOGLYCOSIDE/QUINOLONE"                   = "#546E7A",
  "PHENICOL/QUINOLONE"                         = "#FF5722"
)

VF_COL <- c(
  "Immune modulation"                            = "#1F4E79",
  "Adherence"                                    = "#0F6E56",
  "Motility"                                     = "#BA7517",
  "Nutritional/Metabolic factor"                 = "#2E75B6",
  "Stress survival"                              = "#3B6D11",
  "Effector delivery system"                     = "#6A4C93",
  "Exotoxin"                                     = "#C0392B",
  "Exoenzyme"                                    = "#8D6E63",
  "Invasion"                                     = "#16A085",
  "Regulation"                                   = "#9E9E9E",
  "Antimicrobial activity/Competitive advantage" = "#E74C3C",
  "Biofilm"                                      = "#00838F",
  "Other"                                        = "#BDBDBD"
)

# =============================================================================
# JOURNAL THEME
# =============================================================================
theme_pub <- function(base_size = 9, base_family = "Helvetica") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      plot.title    = element_text(size = base_size + 2, face = "bold",
                                   colour = PAL$navy, hjust = 0,
                                   margin = margin(b = 3)),
      plot.subtitle = element_text(size = base_size - 0.5, colour = PAL$grey,
                                   hjust = 0, margin = margin(b = 7)),
      plot.caption  = element_text(size = base_size - 2.5, colour = PAL$grey,
                                   hjust = 0, lineheight = 1.25,
                                   margin = margin(t = 7)),
      axis.title    = element_text(size = base_size, colour = "black"),
      axis.text     = element_text(size = base_size - 1, colour = "black"),
      panel.grid.major.y = element_line(colour = PAL$lgrey, linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.border  = element_blank(),
      axis.line     = element_line(colour = "black", linewidth = 0.35),
      axis.ticks    = element_line(colour = "black", linewidth = 0.35),
      axis.ticks.length = unit(2, "pt"),
      legend.title  = element_text(size = base_size - 1, face = "bold"),
      legend.text   = element_text(size = base_size - 1),
      legend.key.size   = unit(8, "pt"),
      legend.background = element_blank(),
      legend.margin = margin(0, 0, 0, 0),
      strip.text    = element_text(size = base_size, face = "bold",
                                   colour = PAL$navy),
      strip.background = element_blank(),
      plot.margin   = margin(6, 8, 6, 6)
    )
}

theme_pub_flip <- function(...) {
  theme_pub(...) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(colour = PAL$lgrey,
                                            linewidth = 0.3))
}

theme_pub_heat <- function(...) {
  theme_pub(...) +
    theme(panel.grid = element_blank(),
          axis.line  = element_blank(),
          axis.ticks = element_blank())
}

theme_set(theme_pub())

# ---- Caption wrapper (prevents text running off the page) ------------------
# width is in characters; tune per figure if a caption is unusually long.
wrap_cap <- function(txt, width = 130) {
  paste(strwrap(txt, width = width), collapse = "\n")
}

# ---- Shared caption theme for plot_annotation() ----------------------------
cap_theme <- function(size = 6.4) {
  theme(plot.caption = element_text(size = size, colour = PAL$grey,
                                    hjust = 0, lineheight = 1.25))
}

# ---- Export helper ---------------------------------------------------------
save_fig <- function(plot, name, width, height) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), plot,
         width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(file.path(OUT, paste0(name, ".tiff")), plot,
         width = width, height = height, units = "in",
         dpi = 600, compression = "lzw")
  ggsave(file.path(OUT, paste0(name, ".png")), plot,
         width = width, height = height, units = "in", dpi = 300, bg = "white")
  message("  [ok] ", name, "  (", width, " x ", height, " in)")
}

safe_read_tsv <- function(path, ...) {
  if (!file.exists(path) || file.size(path) == 0) return(NULL)
  tryCatch(suppressWarnings(read_tsv(path, show_col_types = FALSE,
                                     progress = FALSE, ...)),
           error = function(e) NULL)
}

# =============================================================================
# LOAD CORE DATA
# =============================================================================
message("\n=== Loading data ===")

stopifnot(file.exists(RISK_CSV))
risk <- read_csv(RISK_CSV, show_col_types = FALSE) %>%
  mutate(risk_tier = factor(risk_tier, levels = c("High", "Moderate", "Low")))

hgt <- read_csv(HGT_CSV, show_col_types = FALSE)

samples <- readLines(SAMPLES)
samples <- samples[nzchar(samples) & !startsWith(samples, "#")]

N <- nrow(risk)
message("  Samples: ", N)

tier_lookup <- setNames(as.character(risk$risk_tier), risk$sample)

hgt_cols <- hgt %>%
  select(sample, any_of(c("n_args", "n_mges", "n_colocalized",
                          "n_proximal", "n_plasmid_borne", "hgt_potential")))
dat <- risk %>%
  select(-any_of(setdiff(names(hgt_cols), "sample"))) %>%
  left_join(hgt_cols, by = "sample")

# =============================================================================
# PARSE PER-SAMPLE ANNOTATIONS
# =============================================================================
message("=== Parsing annotations ===")

# ---- AMRFinderPlus ----------------------------------------------------------
amr <- map_dfr(samples, function(s) {
  df <- safe_read_tsv(file.path(ANN_DIR, s, "amrfinder.tsv"))
  if (is.null(df) || nrow(df) == 0) return(NULL)
  nm      <- names(df)
  gene_c  <- nm[str_detect(nm, regex("^(Element symbol|Gene symbol)$", TRUE))][1]
  class_c <- nm[str_detect(nm, regex("^Class$", TRUE))][1]
  if (is.na(gene_c) || is.na(class_c)) return(NULL)
  tibble(sample     = s,
         gene       = as.character(df[[gene_c]]),
         drug_class = as.character(df[[class_c]]))
}) %>%
  filter(!is.na(gene), gene != "") %>%
  mutate(risk_tier = factor(tier_lookup[sample],
                            levels = c("High", "Moderate", "Low")))

message("  ARG hits: ", nrow(amr), " across ", n_distinct(amr$sample), " samples")

# ---- Abricate reader (VFDB / PlasmidFinder) --------------------------------
ABR_COLS <- c("FILE","SEQUENCE","START","END","STRAND","GENE","COVERAGE",
              "COVERAGE_MAP","GAPS","PCT_COVERAGE","PCT_IDENTITY","DATABASE",
              "ACCESSION","PRODUCT","RESISTANCE")

read_abricate <- function(path) {
  if (!file.exists(path) || file.size(path) == 0) return(NULL)
  df <- tryCatch(
    suppressWarnings(read_tsv(path, comment = "#", col_names = ABR_COLS,
                              skip = 1, show_col_types = FALSE,
                              progress = FALSE)),
    error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df %>% filter(!is.na(GENE), GENE != "GENE")
}

# ---- VFDB -------------------------------------------------------------------
vf <- map_dfr(samples, function(s) {
  df <- read_abricate(file.path(ANN_DIR, s, "abricate_vfdb.tsv"))
  if (is.null(df)) return(NULL)
  df %>%
    mutate(
      sample = s,
      vf_category = str_trim(str_match(PRODUCT,
                      "-\\s*([^()\\[\\]]+?)\\s*\\(VFC\\d+\\)")[, 2]),
      vf_system   = str_trim(str_match(PRODUCT,
                      "\\[([^\\[\\]]+?)\\s*\\(VF\\d+\\)")[, 2]),
      source_org  = str_trim(str_match(PRODUCT,
                      "\\[([^\\[\\]]+)\\]\\s*$")[, 2])
    ) %>%
    select(sample, gene = GENE, pct_id = PCT_IDENTITY,
           pct_cov = PCT_COVERAGE, vf_category, vf_system, source_org)
}) %>%
  mutate(vf_category = ifelse(is.na(vf_category) | vf_category == "",
                              "Other", vf_category),
         risk_tier   = factor(tier_lookup[sample],
                              levels = c("High", "Moderate", "Low")))

message("  VF hits: ", nrow(vf), " across ", n_distinct(vf$sample), " samples")

# ---- PlasmidFinder ----------------------------------------------------------
pf <- map_dfr(samples, function(s) {
  df <- read_abricate(file.path(ANN_DIR, s, "abricate_plasmidfinder.tsv"))
  if (is.null(df)) return(NULL)
  df %>% mutate(sample = s) %>%
    select(sample, gene = GENE, pct_id = PCT_IDENTITY, pct_cov = PCT_COVERAGE)
})
if (nrow(pf) > 0) {
  pf <- pf %>% mutate(risk_tier = factor(tier_lookup[sample],
                                         levels = c("High","Moderate","Low")))
}
message("  Plasmid hits: ", nrow(pf), " across ", n_distinct(pf$sample), " samples")

# ---- Bracken taxonomy -------------------------------------------------------
PRIORITY <- c("Staphylococcus aureus", "Streptococcus pneumoniae",
              "Haemophilus influenzae", "Moraxella catarrhalis",
              "Klebsiella pneumoniae", "Acinetobacter baumannii",
              "Pseudomonas aeruginosa", "Escherichia coli",
              "Enterococcus faecalis", "Enterococcus faecium",
              "Staphylococcus epidermidis")

tax <- map_dfr(samples, function(s) {
  df <- safe_read_tsv(file.path(TAX_DIR, s, paste0(s, ".bracken.txt")))
  if (is.null(df) || nrow(df) == 0) return(NULL)
  if (!all(c("name", "fraction_total_reads") %in% names(df))) return(NULL)
  df %>%
    transmute(sample = s, taxon = name,
              frac = as.numeric(fraction_total_reads)) %>%
    filter(!is.na(frac))
})
message("  Taxonomy rows: ", nrow(tax), " across ", n_distinct(tax$sample), " samples")

# =============================================================================
# DERIVED SUMMARY STATISTICS
# =============================================================================
arg_per_sample <- amr %>% count(sample, name = "n_arg_hits")
vf_per_sample  <- vf  %>% count(sample, name = "n_vf_hits")

dat <- dat %>%
  left_join(arg_per_sample, by = "sample") %>%
  left_join(vf_per_sample,  by = "sample") %>%
  mutate(across(c(n_arg_hits, n_vf_hits), ~replace_na(., 0)))

n_arg_pos <- sum(dat$n_arg_hits > 0)
n_mdr     <- sum(dat$n_drug_classes >= 3, na.rm = TRUE)
n_vf_pos  <- sum(dat$n_vf_hits > 0)
n_class   <- n_distinct(amr$drug_class[!is.na(amr$drug_class)])

# Shortened base caption (v2.1): ~30 characters lighter than v2.0
cap_base <- sprintf(
  "Bamako nasal resistome cohort (n = %d). CRS = 0.4 x MDR + 0.3 x HGT + 0.3 x pathogen abundance (weights exploratory).", N)

samp_order <- dat %>% arrange(desc(crs)) %>% pull(sample)

# =============================================================================
# FIGURE 1 -- CRS distribution and risk-tier composition
# =============================================================================
message("\n=== Building figures ===")

p1a <- ggplot(dat, aes(x = crs)) +
  annotate("rect", xmin = -Inf, xmax = 0.40, ymin = -Inf, ymax = Inf,
           fill = PAL$lteal, alpha = 0.05) +
  annotate("rect", xmin = 0.40, xmax = 0.70, ymin = -Inf, ymax = Inf,
           fill = PAL$orange, alpha = 0.07) +
  annotate("rect", xmin = 0.70, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = PAL$red, alpha = 0.07) +
  geom_histogram(binwidth = 0.025, fill = PAL$navy,
                 colour = "white", linewidth = 0.2, boundary = 0) +
  geom_rug(aes(colour = risk_tier), sides = "b",
           linewidth = 0.4, length = unit(3, "pt"), show.legend = FALSE) +
  geom_vline(xintercept = c(0.40, 0.70), linetype = "dashed",
             colour = c(PAL$orange, PAL$red), linewidth = 0.4) +
  scale_colour_manual(values = TIER_COL, drop = FALSE) +
  scale_x_continuous(limits = c(-0.02, 1), breaks = seq(0, 1, 0.2),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Clinical Risk Score (CRS)", y = "Samples",
       title = "a", subtitle = sprintf("Mean CRS = %.3f", mean(dat$crs)))

tier_df <- dat %>%
  count(risk_tier, .drop = FALSE) %>%
  mutate(pct = n / sum(n) * 100,
         lab = sprintf("%d (%.1f%%)", n, pct))

p1b <- ggplot(tier_df, aes(x = fct_rev(risk_tier), y = n, fill = risk_tier)) +
  geom_col(width = 0.62, colour = "white", linewidth = 0.25) +
  geom_text(aes(label = lab), hjust = -0.12, size = 2.9, colour = PAL$navy) +
  scale_fill_manual(values = TIER_COL, guide = "none", drop = FALSE) +
  # v2.1: 0.30 -> 0.42 so the widest label ("77 (97.5%)") is not clipped
  scale_y_continuous(expand = expansion(mult = c(0, 0.42))) +
  coord_flip() +
  labs(x = NULL, y = "Samples", title = "b",
       subtitle = "CRS risk-tier composition") +
  theme_pub_flip()

fig1 <- (p1a | p1b) +
  plot_layout(widths = c(1.7, 1)) +
  plot_annotation(
    caption = wrap_cap(paste0(cap_base,
      " Thresholds: High > 0.70, Moderate 0.40-0.70, Low < 0.40.")),
    theme = cap_theme(6.5))
save_fig(fig1, "Figure1_CRS_distribution", 7.2, 3.2)

# =============================================================================
# FIGURE 2 -- Ranked CRS with component decomposition
# =============================================================================
comp <- dat %>%
  transmute(sample, risk_tier, crs,
            `MDR index (x0.4)`          = 0.4 * replace_na(mdr_norm, 0),
            `HGT potential (x0.3)`      = 0.3 * replace_na(hgt_norm, 0),
            `Pathogen abundance (x0.3)` = 0.3 * replace_na(pathogen_norm, 0)) %>%
  pivot_longer(cols = c("MDR index (x0.4)", "HGT potential (x0.3)",
                        "Pathogen abundance (x0.3)"),
               names_to = "component", values_to = "value") %>%
  mutate(sample    = factor(sample, levels = samp_order),
         component = factor(component,
                            levels = c("MDR index (x0.4)",
                                       "HGT potential (x0.3)",
                                       "Pathogen abundance (x0.3)")))

lab_df <- dat %>%
  filter(risk_tier %in% c("High", "Moderate")) %>%
  mutate(sample = factor(sample, levels = samp_order))

fig2 <- ggplot(comp, aes(x = sample, y = value, fill = component)) +
  geom_col(width = 0.85) +
  geom_hline(yintercept = c(0.40, 0.70), linetype = "dashed",
             colour = c(PAL$orange, PAL$red), linewidth = 0.35) +
  geom_text(data = lab_df, aes(x = sample, y = crs + 0.03, label = sample),
            inherit.aes = FALSE, angle = 90, hjust = 0, size = 2.3,
            colour = PAL$orange, fontface = "bold") +
  scale_fill_manual(values = c(PAL$navy, PAL$teal, PAL$amber), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14)),
                     breaks = seq(0, 1, 0.2)) +
  labs(x = sprintf("Samples ordered by CRS (n = %d)", N),
       y = "CRS contribution",
       title = "Clinical Risk Score decomposition across the cohort",
       caption = wrap_cap(paste0(cap_base,
         " Dashed lines mark the Moderate (0.40) and High (0.70) thresholds."))) +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "top",
        panel.grid.major.y = element_line(colour = PAL$lgrey, linewidth = 0.3))
save_fig(fig2, "Figure2_CRS_components", 7.2, 3.6)

# =============================================================================
# FIGURE 3 -- Clinical risk landscape (bubble plot)
# =============================================================================
bub <- dat %>%
  mutate(hgt_plot = replace_na(hgt_norm, 0),
         lab = ifelse(risk_tier %in% c("High", "Moderate") | crs >= 0.35,
                      sample, NA_character_))

fig3 <- ggplot(bub, aes(x = pathogen_abundance_raw, y = n_drug_classes)) +
  geom_point(aes(size = hgt_plot, fill = risk_tier),
             shape = 21, colour = "white", stroke = 0.35, alpha = 0.9) +
  geom_text_repel(aes(label = lab, colour = risk_tier),
                  size = 2.4, fontface = "bold", na.rm = TRUE,
                  min.segment.length = 0, segment.size = 0.25,
                  box.padding = 0.4, max.overlaps = 30, show.legend = FALSE) +
  scale_fill_manual(values = TIER_COL, name = "CRS risk tier", drop = FALSE) +
  scale_colour_manual(values = TIER_COL, drop = FALSE) +
  scale_size_continuous(range = c(1.4, 9), name = "HGT potential") +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0.03, 0.07))) +
  scale_y_continuous(breaks = pretty_breaks()) +
  guides(fill = guide_legend(order = 1, override.aes = list(size = 3.2)),
         size = guide_legend(order = 2)) +
  labs(x = "Priority-pathogen relative abundance",
       y = "MDR index (no. of antibiotic drug classes)",
       title = "Clinical risk landscape of the nasal resistome",
       subtitle = "Point area scales with contig-level HGT potential",
       caption = wrap_cap(cap_base, width = 108))
save_fig(fig3, "Figure3_risk_landscape", 6.8, 5.2)

# =============================================================================
# FIGURE 4 -- Resistome overview (four panels)
# =============================================================================
dc_prev <- amr %>%
  filter(!is.na(drug_class), drug_class != "") %>%
  distinct(sample, drug_class) %>%
  count(drug_class, name = "n_samples") %>%
  mutate(pct = n_samples / N * 100) %>%
  arrange(n_samples)

p4a <- ggplot(dc_prev, aes(x = fct_inorder(drug_class), y = n_samples,
                           fill = drug_class)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.2) +
  geom_text(aes(label = n_samples), hjust = -0.25, size = 2.4,
            colour = PAL$navy) +
  scale_fill_manual(values = DC_COL, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  coord_flip() +
  labs(x = NULL, y = "Samples", title = "a",
       subtitle = sprintf("Drug-class prevalence (%d classes)", n_class)) +
  theme_pub_flip() +
  theme(axis.text.y = element_text(size = 5.6))

top_arg <- amr %>%
  distinct(sample, gene, drug_class) %>%
  count(gene, drug_class, name = "n_samples") %>%
  slice_max(n_samples, n = 20, with_ties = FALSE) %>%
  arrange(n_samples)

p4b <- ggplot(top_arg, aes(x = fct_inorder(gene), y = n_samples,
                           fill = drug_class)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.2) +
  geom_text(aes(label = n_samples), hjust = -0.3, size = 2.4,
            colour = PAL$navy) +
  scale_fill_manual(values = DC_COL, name = "Drug class",
                    labels = label_wrap(28)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  coord_flip() +
  labs(x = NULL, y = "Samples", title = "b",
       subtitle = "Twenty most frequent ARGs") +
  theme_pub_flip() +
  theme(axis.text.y  = element_text(face = "italic", size = 6.2),
        legend.position = "right",
        legend.text  = element_text(size = 5.2),
        legend.title = element_text(size = 6),
        legend.key.size = unit(6, "pt"))

arg_stack <- amr %>%
  filter(!is.na(drug_class)) %>%
  count(sample, drug_class, name = "hits") %>%
  mutate(sample = factor(sample, levels = samp_order)) %>%
  filter(!is.na(sample))

p4c <- ggplot(arg_stack, aes(x = sample, y = hits, fill = drug_class)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = DC_COL, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
  labs(x = "ARG-positive samples, ordered by CRS", y = "ARG hits",
       title = "c", subtitle = "Per-sample ARG burden by drug class") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                   size = 4.2))

mdr_df <- dat %>%
  mutate(status = case_when(
    n_drug_classes >= 3 ~ "MDR (>= 3 classes)",
    n_arg_hits > 0      ~ "ARG-positive, non-MDR",
    TRUE                ~ "No ARG detected")) %>%
  count(status) %>%
  mutate(status = factor(status, levels = c("No ARG detected",
                                            "ARG-positive, non-MDR",
                                            "MDR (>= 3 classes)")),
         pct = n / sum(n) * 100,
         lab = sprintf("%d\n(%.1f%%)", n, pct))

p4d <- ggplot(mdr_df, aes(x = status, y = n, fill = status)) +
  geom_col(width = 0.6, colour = "white", linewidth = 0.25) +
  geom_text(aes(label = lab), vjust = -0.25, size = 2.5, colour = PAL$navy,
            lineheight = 0.9) +
  scale_fill_manual(values = c(PAL$lgrey, PAL$blue, PAL$red), guide = "none") +
  scale_x_discrete(labels = label_wrap(14)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.24))) +
  labs(x = NULL, y = "Samples", title = "d",
       subtitle = "Multidrug-resistance status")

fig4 <- (p4a | p4b) / (p4c | p4d) +
  plot_layout(heights = c(1, 0.85)) +
  plot_annotation(
    caption = wrap_cap(sprintf(
      "%s ARGs detected in %d/%d samples (%.1f%%); %d/%d (%.1f%%) met the MDR definition (>= 3 drug classes) across %d classes. AMRFinderPlus/CARD, >= 90%% identity, >= 50%% coverage.",
      cap_base, n_arg_pos, N, n_arg_pos / N * 100,
      n_mdr, N, n_mdr / N * 100, n_class), width = 155),
    theme = cap_theme(6.2))
save_fig(fig4, "Figure4_resistome_overview", 9.0, 8.2)

# =============================================================================
# FIGURE 5 -- ARG gene x sample heatmap (ggplot2, replaces ComplexHeatmap)
# =============================================================================
heat_genes <- amr %>%
  distinct(sample, gene) %>% count(gene, name = "n") %>%
  slice_max(n, n = 25, with_ties = FALSE) %>% pull(gene)

heat_samples <- samp_order[samp_order %in% amr$sample]

heat_df <- expand_grid(gene = heat_genes, sample = heat_samples) %>%
  left_join(amr %>% distinct(sample, gene) %>% mutate(present = 1L),
            by = c("gene", "sample")) %>%
  mutate(present = replace_na(present, 0L),
         sample  = factor(sample, levels = heat_samples),
         gene    = factor(gene, levels = rev(heat_genes)))

ann_df <- tibble(sample = heat_samples) %>%
  mutate(risk_tier = factor(tier_lookup[sample],
                            levels = c("High", "Moderate", "Low")),
         sample = factor(sample, levels = heat_samples))

p5_ann <- ggplot(ann_df, aes(x = sample, y = 1, fill = risk_tier)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  scale_fill_manual(values = TIER_COL, name = "CRS tier", drop = FALSE) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = NULL,
       title = "ARG gene detection across ARG-positive samples",
       subtitle = "Twenty-five most frequent genes; samples ordered by CRS") +
  theme_pub_heat() +
  theme(axis.text = element_blank(),
        legend.position = "top",
        plot.margin = margin(6, 8, 0, 6))

p5_main <- ggplot(heat_df, aes(x = sample, y = gene, fill = factor(present))) +
  geom_tile(colour = "white", linewidth = 0.25) +
  scale_fill_manual(values = c(`0` = PAL$lgrey, `1` = PAL$navy),
                    labels = c("Absent", "Present"), name = NULL) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Samples (ordered by CRS)", y = NULL,
       caption = wrap_cap(cap_base, width = 125)) +
  theme_pub_heat() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                   size = 4.6),
        axis.text.y = element_text(face = "italic", size = 6.2),
        legend.position = "bottom",
        plot.margin = margin(0, 8, 6, 6))

fig5 <- p5_ann / p5_main + plot_layout(heights = c(0.055, 1))
save_fig(fig5, "Figure5_ARG_heatmap", 8.0, 5.8)

# =============================================================================
# FIGURE 6 -- Virulome overview
# =============================================================================
if (nrow(vf) > 0) {

  vf_cat <- vf %>% count(vf_category, name = "hits") %>%
    slice_max(hits, n = 10, with_ties = FALSE) %>% arrange(hits)

  p6a <- ggplot(vf_cat, aes(x = fct_inorder(vf_category), y = hits,
                            fill = vf_category)) +
    geom_col(width = 0.72, colour = "white", linewidth = 0.2) +
    geom_text(aes(label = hits), hjust = -0.22, size = 2.4, colour = PAL$navy) +
    scale_fill_manual(values = VF_COL, guide = "none") +
    scale_x_discrete(labels = label_wrap(26)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
    coord_flip() +
    labs(x = NULL, y = "VF hits", title = "a",
         subtitle = "Virulence-factor functional categories") +
    theme_pub_flip() +
    theme(axis.text.y = element_text(size = 6))

  top_vf <- vf %>%
    distinct(sample, gene, vf_category) %>%
    count(gene, vf_category, name = "n_samples") %>%
    slice_max(n_samples, n = 20, with_ties = FALSE) %>%
    arrange(n_samples)

  p6b <- ggplot(top_vf, aes(x = fct_inorder(gene), y = n_samples,
                            fill = vf_category)) +
    geom_col(width = 0.72, colour = "white", linewidth = 0.2) +
    geom_text(aes(label = n_samples), hjust = -0.3, size = 2.4,
              colour = PAL$navy) +
    scale_fill_manual(values = VF_COL, name = "VF category",
                      labels = label_wrap(24)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
    coord_flip() +
    labs(x = NULL, y = "Samples", title = "b",
         subtitle = "Twenty most frequent VF genes") +
    theme_pub_flip() +
    theme(axis.text.y  = element_text(face = "italic", size = 6.2),
          legend.text  = element_text(size = 5.2),
          legend.title = element_text(size = 6),
          legend.key.size = unit(6, "pt"))

  vf_burden <- vf %>% count(sample, name = "hits") %>%
    left_join(dat %>% select(sample, risk_tier), by = "sample") %>%
    slice_max(hits, n = 20, with_ties = FALSE) %>% arrange(hits)

  p6c <- ggplot(vf_burden, aes(x = fct_inorder(sample), y = hits,
                               fill = risk_tier)) +
    geom_col(width = 0.72, colour = "white", linewidth = 0.2) +
    geom_text(aes(label = hits), hjust = -0.25, size = 2.4, colour = PAL$navy) +
    scale_fill_manual(values = TIER_COL, name = "CRS tier", drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
    coord_flip() +
    labs(x = NULL, y = "VF hits", title = "c",
         subtitle = "Per-sample virulence burden") +
    theme_pub_flip() +
    theme(axis.text.y = element_text(size = 5.8))

  p6d <- ggplot(dat, aes(x = n_arg_hits, y = n_vf_hits)) +
    geom_point(aes(fill = risk_tier), shape = 21, size = 2.1,
               colour = "white", stroke = 0.3, alpha = 0.9) +
    geom_text_repel(data = dat %>%
                      filter(risk_tier %in% c("High", "Moderate") |
                             n_vf_hits > quantile(n_vf_hits, 0.95)),
                    aes(label = sample), size = 2.2, colour = PAL$navy,
                    min.segment.length = 0, segment.size = 0.2,
                    max.overlaps = 20) +
    scale_fill_manual(values = TIER_COL, guide = "none", drop = FALSE) +
    labs(x = "ARG hits", y = "VF hits", title = "d",
         subtitle = "Resistome-virulome co-occurrence")

  fig6 <- (p6a | p6b) / (p6c | p6d) +
    plot_annotation(caption = wrap_cap(sprintf(
      "%s Virulence factors detected in %d/%d samples (%.1f%%); %d total VFDB hits (>= 80%% identity, >= 60%% coverage).",
      cap_base, n_vf_pos, N, n_vf_pos / N * 100, nrow(vf)), width = 155),
      theme = cap_theme(6.2))
  save_fig(fig6, "Figure6_virulome_overview", 9.0, 7.8)
} else {
  message("  [skip] Figure 6 -- no VF hits")
}

# =============================================================================
# FIGURE 7 -- HGT potential and mobilisation evidence
# =============================================================================
p7a <- ggplot(dat, aes(x = replace_na(hgt_norm, 0))) +
  geom_histogram(binwidth = 0.01, fill = PAL$teal,
                 colour = "white", linewidth = 0.2, boundary = 0) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.05))) +
  labs(x = "HGT potential", y = "Samples", title = "a",
       subtitle = sprintf("%d/%d samples with HGT potential > 0",
                          sum(replace_na(dat$hgt_norm, 0) > 0), N))

p7b <- ggplot(dat, aes(x = n_arg_hits, y = replace_na(hgt_norm, 0))) +
  geom_point(aes(fill = risk_tier), shape = 21, size = 2.2,
             colour = "white", stroke = 0.3, alpha = 0.9) +
  geom_text_repel(data = dat %>% filter(replace_na(hgt_norm, 0) > 0.02),
                  aes(label = sample), size = 2.2, colour = PAL$navy,
                  min.segment.length = 0, segment.size = 0.2,
                  max.overlaps = 20) +
  scale_fill_manual(values = TIER_COL, name = "CRS tier", drop = FALSE) +
  labs(x = "ARG hits", y = "HGT potential", title = "b",
       subtitle = "ARG burden is decoupled from mobilisation potential")

mob <- tibble(
  evidence = c("ARG-MGE co-localised", "Proximal (< 5 kb)",
               "Plasmid-borne ARG", "Plasmid replicon"),
  n = c(sum(replace_na(dat$n_colocalized, 0)   > 0),
        sum(replace_na(dat$n_proximal, 0)      > 0),
        sum(replace_na(dat$n_plasmid_borne, 0) > 0),
        n_distinct(pf$sample))
) %>%
  mutate(evidence = fct_inorder(evidence),
         lab = sprintf("%d (%.1f%%)", n, n / N * 100))

p7c <- ggplot(mob, aes(x = fct_rev(evidence), y = n)) +
  geom_col(width = 0.6, fill = PAL$teal, colour = "white", linewidth = 0.25) +
  geom_text(aes(label = lab), hjust = -0.12, size = 2.5, colour = PAL$navy) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.42))) +
  scale_x_discrete(labels = label_wrap(18)) +
  coord_flip() +
  labs(x = NULL, y = "Samples", title = "c",
       subtitle = "Contig-level mobilisation evidence") +
  theme_pub_flip()

fig7 <- (p7a | p7b) / p7c +
  plot_layout(heights = c(1, 0.72)) +
  plot_annotation(caption = wrap_cap(paste0(cap_base,
    " HGT potential integrates ARG-MGE co-localisation on shared contigs, distance weighting (< 5 kb) and vehicle classification (plasmid 1.0 vs chromosomal 0.6)."),
    width = 135),
    theme = cap_theme(6.2))
save_fig(fig7, "Figure7_HGT_potential", 7.6, 6.2)

# =============================================================================
# FIGURE 8 -- Taxonomy
# =============================================================================
if (nrow(tax) > 0) {

  top_taxa <- tax %>% group_by(taxon) %>%
    summarise(total = sum(frac), .groups = "drop") %>%
    slice_max(total, n = 12, with_ties = FALSE) %>% pull(taxon)

  tax_plot <- tax %>%
    mutate(grp = ifelse(taxon %in% top_taxa, taxon, "Other taxa")) %>%
    group_by(sample, grp) %>%
    summarise(frac = sum(frac), .groups = "drop") %>%
    mutate(sample = factor(sample, levels = samp_order),
           grp    = factor(grp, levels = c(top_taxa, "Other taxa"))) %>%
    filter(!is.na(sample))

  p8a <- ggplot(tax_plot, aes(x = sample, y = frac, fill = grp)) +
    geom_col(width = 0.9) +
    scale_fill_manual(values = c(viridis(length(top_taxa), option = "D",
                                         end = 0.92), PAL$lgrey),
                      name = "Taxon", labels = label_wrap(28)) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Samples (ordered by CRS)", y = "Relative abundance",
         title = "a", subtitle = "Twelve most abundant species") +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          legend.text  = element_text(size = 5.4, face = "italic"),
          legend.title = element_text(size = 6.2),
          legend.key.size = unit(6, "pt"))

  prio <- tax %>%
    filter(taxon %in% PRIORITY, frac > 0.0001) %>%
    distinct(sample, taxon) %>%
    count(taxon, name = "n_samples") %>%
    mutate(pct = n_samples / N * 100) %>%
    arrange(n_samples)

  p8b <- ggplot(prio, aes(x = fct_inorder(taxon), y = n_samples)) +
    geom_col(width = 0.66, fill = PAL$red, colour = "white", linewidth = 0.25) +
    geom_text(aes(label = sprintf("%d (%.1f%%)", n_samples, pct)),
              hjust = -0.12, size = 2.4, colour = PAL$navy) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.42))) +
    coord_flip() +
    labs(x = NULL, y = "Samples", title = "b",
         subtitle = "WHO priority pathogen detection") +
    theme_pub_flip() +
    theme(axis.text.y = element_text(face = "italic", size = 6.2))

  fig8 <- p8a / p8b +
    plot_layout(heights = c(1, 0.85)) +
    plot_annotation(caption = wrap_cap(paste0(cap_base,
      " Taxonomy: Kraken2 and Bracken, species level, threshold 10 reads."),
      width = 135),
      theme = cap_theme(6.2))
  save_fig(fig8, "Figure8_taxonomy", 7.6, 6.6)
} else {
  message("  [skip] Figure 8 -- no Bracken data")
}

# =============================================================================
# PUBLICATION TABLES
# =============================================================================
message("\n=== Writing tables ===")

tab1 <- tribble(
  ~Metric, ~Value,
  "Samples analysed",                       as.character(N),
  "Samples with ARGs detected",             sprintf("%d (%.1f%%)", n_arg_pos, n_arg_pos / N * 100),
  "MDR samples (>= 3 drug classes)",        sprintf("%d (%.1f%%)", n_mdr, n_mdr / N * 100),
  "Unique antibiotic drug classes",         as.character(n_class),
  "Total ARG hits",                         as.character(nrow(amr)),
  "Samples with virulence factors",         sprintf("%d (%.1f%%)", n_vf_pos, n_vf_pos / N * 100),
  "Total VFDB hits",                        as.character(nrow(vf)),
  "Samples with plasmid replicons",         sprintf("%d (%.1f%%)", n_distinct(pf$sample), n_distinct(pf$sample) / N * 100),
  "Samples with ARG-MGE co-localisation",   as.character(sum(replace_na(dat$n_colocalized, 0) > 0)),
  "High risk (CRS > 0.70)",                 sprintf("%d (%.1f%%)", sum(dat$risk_tier == "High"), sum(dat$risk_tier == "High") / N * 100),
  "Moderate risk (CRS 0.40-0.70)",          sprintf("%d (%.1f%%)", sum(dat$risk_tier == "Moderate"), sum(dat$risk_tier == "Moderate") / N * 100),
  "Low risk (CRS < 0.40)",                  sprintf("%d (%.1f%%)", sum(dat$risk_tier == "Low"), sum(dat$risk_tier == "Low") / N * 100),
  "Mean CRS",                               sprintf("%.3f", mean(dat$crs)),
  "Median CRS (IQR)",                       sprintf("%.3f (%.3f-%.3f)", median(dat$crs), quantile(dat$crs, .25), quantile(dat$crs, .75))
)
write_csv(tab1, file.path(OUT, "Table1_cohort_summary.csv"))

tab2 <- dat %>%
  filter(risk_tier %in% c("High", "Moderate") | crs >= 0.35) %>%
  arrange(desc(crs)) %>%
  transmute(Sample = sample,
            CRS = round(crs, 3),
            `Risk tier` = as.character(risk_tier),
            `ARG hits` = n_arg_hits,
            `Drug classes` = n_drug_classes,
            `VF hits` = n_vf_hits,
            `HGT potential` = round(replace_na(hgt_norm, 0), 4),
            `ARG-MGE co-localised` = replace_na(n_colocalized, 0),
            `Plasmid-borne ARG` = replace_na(n_plasmid_borne, 0),
            `Dominant priority pathogen` = top_pathogen,
            `Pathogen abundance` = sprintf("%.2f%%",
                                    replace_na(pathogen_abundance_raw, 0) * 100))
write_csv(tab2, file.path(OUT, "Table2_elevated_risk_samples.csv"))

tab3 <- dc_prev %>% arrange(desc(n_samples)) %>%
  transmute(`Drug class` = drug_class,
            `Samples (n)` = n_samples,
            `Prevalence (%)` = round(pct, 1))
write_csv(tab3, file.path(OUT, "Table3_drug_class_prevalence.csv"))

tab4 <- amr %>%
  group_by(gene, drug_class) %>%
  summarise(n_samples = n_distinct(sample), total_hits = n(), .groups = "drop") %>%
  arrange(desc(n_samples), gene) %>%
  transmute(Gene = gene, `Drug class` = drug_class,
            `Samples (n)` = n_samples,
            `Prevalence (%)` = round(n_samples / N * 100, 1),
            `Total hits` = total_hits)
write_csv(tab4, file.path(OUT, "Table4_ARG_inventory.csv"))

if (nrow(vf) > 0) {
  tab5 <- vf %>%
    group_by(gene, vf_category) %>%
    summarise(n_samples = n_distinct(sample), total_hits = n(), .groups = "drop") %>%
    arrange(desc(n_samples), gene) %>%
    transmute(Gene = gene, `VF category` = vf_category,
              `Samples (n)` = n_samples,
              `Prevalence (%)` = round(n_samples / N * 100, 1),
              `Total hits` = total_hits)
  write_csv(tab5, file.path(OUT, "Table5_VF_inventory.csv"))
}

supp <- dat %>%
  arrange(desc(crs)) %>%
  transmute(Sample = sample, CRS = round(crs, 4),
            `Risk tier` = as.character(risk_tier),
            `MDR index (normalised)` = round(replace_na(mdr_norm, 0), 4),
            `Drug classes` = n_drug_classes,
            `ARG hits` = n_arg_hits,
            `VF hits` = n_vf_hits,
            `HGT potential` = round(replace_na(hgt_norm, 0), 4),
            MGEs = replace_na(n_mges, 0),
            `ARG-MGE co-localised` = replace_na(n_colocalized, 0),
            `Proximal (<5kb)` = replace_na(n_proximal, 0),
            `Plasmid-borne` = replace_na(n_plasmid_borne, 0),
            `Pathogen abundance (normalised)` = round(replace_na(pathogen_norm, 0), 4),
            `Dominant priority pathogen` = top_pathogen,
            `Priority taxa detected` = n_priority_taxa,
            `Drug class detail` = drug_classes)
write_csv(supp, file.path(OUT, "TableS1_per_sample_matrix.csv"))

write_csv(amr, file.path(OUT, "TableS2_ARG_hits_long.csv"))
if (nrow(vf) > 0) write_csv(vf, file.path(OUT, "TableS3_VF_hits_long.csv"))
if (nrow(pf) > 0) write_csv(pf, file.path(OUT, "TableS4_plasmid_hits_long.csv"))

message("  [ok] 5 main tables + 4 supplementary tables")

# =============================================================================
# SESSION INFO (reproducibility)
# =============================================================================
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))

message("\n=============================================")
message(" Publication figures and tables complete")
message(" Output: ", OUT)
message(" Formats: PDF (vector) | TIFF (600 dpi) | PNG (300 dpi)")
message("=============================================\n")
