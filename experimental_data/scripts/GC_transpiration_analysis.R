# ==============================================================================
# GC Transpiration Analysis — Growth Chamber VPD vs Transpiration
# Experiment: GC Experiment 3 (Part 3), Second time GC
# Target journal: Agricultural Water Management
# ==============================================================================

setwd("/home/user/LT-trait")

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(RColorBrewer)
  library(scales)
})

dir.create("experimental_data/GC_analysis", recursive = TRUE, showWarnings = FALSE)

# ── Colour / theme constants ───────────────────────────────────────────────────
COL_PLASTIC    <- "#d95f02"
COL_NONPLASTIC <- "#31a354"
GENO_COLS      <- brewer.pal(4, "Set1")          # 4 genotype colours
VPD_THRESHOLD  <- 2.0                            # Sinclair 2008 reference

base_theme <- theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey92", colour = "grey70"),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"
  )

# ==============================================================================
# 1. LOAD & CLEAN
# ==============================================================================
cat("Loading data...\n")

DATA_PATH <- "experimental_data/data-analysis/data-analysis/Second time GC/Part 3/GC experiment 3.xlsx"

raw <- read_excel(DATA_PATH, sheet = "Final")

cat("Raw dimensions:", nrow(raw), "x", ncol(raw), "\n")

df <- raw %>%
  mutate(
    # VPD target labels from Treatment code
    VPD_target = factor(
      case_when(
        Treatment == 1 ~ "1.0",
        Treatment == 2 ~ "1.5",
        Treatment == 3 ~ "2.0",
        Treatment == 4 ~ "2.7",
        Treatment == 5 ~ "3.3"
      ),
      levels = c("1.0", "1.5", "2.0", "2.7", "3.3")
    ),
    # Parse Expected_WUE into two columns
    Plasticity = case_when(
      grepl("^Plastic,",    Expected_WUE) ~ "Plastic",
      grepl("^Nonplastic,", Expected_WUE) ~ "Nonplastic"
    ),
    WUE_class = case_when(
      grepl("High WUE", Expected_WUE) ~ "High WUE",
      grepl("Low WUE",  Expected_WUE) ~ "Low WUE"
    ),
    # Short genotype labels (keep PI numbers)
    Genotype_short = Genotype,
    # Facet label: Genotype + Expected_WUE
    Facet_label = paste0(Genotype, "\n", Expected_WUE),
    Day_label = paste("Day", Day),
    # Ensure Day and Treatment are factors where needed
    Day_f      = factor(Day, labels = c("Day 1", "Day 2")),
    Plasticity = factor(Plasticity, levels = c("Plastic", "Nonplastic")),
    WUE_class  = factor(WUE_class,  levels = c("High WUE", "Low WUE"))
  )

# Fix Genotype factor order to match Quadrant
geno_order <- c("PI398617", "PI070229", "PI423802", "PI603465C")
df <- df %>%
  mutate(
    Genotype = factor(Genotype, levels = geno_order),
    Facet_label = factor(Facet_label,
      levels = df %>%
        distinct(Genotype, Facet_label) %>%
        arrange(match(Genotype, geno_order)) %>%
        pull(Facet_label)
    )
  )

cat("Clean data dimensions:", nrow(df), "x", ncol(df), "\n")
cat("Genotypes:", levels(df$Genotype), "\n")
cat("VPD_target levels:", levels(df$VPD_target), "\n")
cat("Plasticity levels:", levels(df$Plasticity), "\n\n")

# Save clean CSV
write.csv(df, "experimental_data/GC_analysis/GC3_clean.csv", row.names = FALSE)
cat("Saved GC3_clean.csv\n\n")

# ==============================================================================
# 2. FIGURES
# ==============================================================================

save_fig <- function(p, name, width_cm, height_cm) {
  base <- file.path("experimental_data/GC_analysis", name)
  ggsave(paste0(base, ".png"),  plot = p, width = width_cm, height = height_cm,
         units = "cm", dpi = 150)
  ggsave(paste0(base, ".tiff"), plot = p, width = width_cm, height = height_cm,
         units = "cm", dpi = 300, compression = "lzw")
  cat("  Saved", name, ".png/.tiff\n")
}

# ── Fig A — VPD vs Relative Transpiration ─────────────────────────────────────
cat("Plotting Fig A...\n")

figA <- ggplot(df, aes(x = VPD, y = RelTranspiration, colour = Plasticity)) +
  geom_vline(xintercept = VPD_THRESHOLD, linetype = "dashed",
             colour = "grey40", linewidth = 0.7) +
  geom_point(alpha = 0.4,
             position = position_jitter(width = 0, height = 0.01, seed = 42),
             size = 1.8) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  scale_colour_manual(
    values = c(Plastic = COL_PLASTIC, Nonplastic = COL_NONPLASTIC),
    name   = "Plasticity"
  ) +
  facet_wrap(~ Facet_label, ncol = 2) +
  scale_x_continuous(breaks = c(1.0, 1.5, 2.0, 2.7, 3.3)) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(
    x     = "VPD (kPa)",
    y     = "Relative Transpiration",
    title = "VPD vs Relative Transpiration by Genotype",
    caption = "Dashed line: VPD = 2.0 kPa (Sinclair 2008 threshold)"
  ) +
  base_theme

save_fig(figA, "GC3_Fig_A_VPD_RelTransp", width_cm = 20, height_cm = 18)

# ── Fig B — Mean +/- SE summary by Genotype x Day ────────────────────────────
cat("Plotting Fig B...\n")

sum_B <- df %>%
  group_by(Genotype, Day_f, VPD_target) %>%
  summarise(
    mean_RT = mean(RelTranspiration, na.rm = TRUE),
    se_RT   = sd(RelTranspiration, na.rm = TRUE) / sqrt(n()),
    n       = n(),
    .groups = "drop"
  )

figB <- ggplot(sum_B, aes(x = VPD_target, y = mean_RT,
                           colour = Genotype, group = Genotype)) +
  geom_hline(yintercept = 1.0, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_RT - se_RT, ymax = mean_RT + se_RT),
                width = 0.15, linewidth = 0.7) +
  scale_colour_manual(values = setNames(GENO_COLS, geno_order)) +
  facet_wrap(~ Day_f) +
  labs(
    x     = "Target VPD (kPa)",
    y     = "Mean Relative Transpiration +/- SE",
    title = "Mean Relative Transpiration by Genotype and Day",
    colour = "Genotype"
  ) +
  base_theme

save_fig(figB, "GC3_Fig_B_summary", width_cm = 20, height_cm = 12)

# ── Fig C — Transpiration per Leaf Area vs VPD ────────────────────────────────
cat("Plotting Fig C...\n")

df_C <- df %>% filter(!is.na(TranspLeafArea))
cat("  Rows after filtering NA TranspLeafArea:", nrow(df_C), "\n")

figC <- ggplot(df_C, aes(x = VPD, y = TranspLeafArea, colour = Plasticity)) +
  geom_vline(xintercept = VPD_THRESHOLD, linetype = "dashed",
             colour = "grey40", linewidth = 0.7) +
  geom_point(alpha = 0.4,
             position = position_jitter(width = 0, height = 0, seed = 42),
             size = 1.8) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  scale_colour_manual(
    values = c(Plastic = COL_PLASTIC, Nonplastic = COL_NONPLASTIC),
    name   = "Plasticity"
  ) +
  facet_wrap(~ Facet_label, ncol = 2) +
  scale_x_continuous(breaks = c(1.0, 1.5, 2.0, 2.7, 3.3)) +
  labs(
    x     = "VPD (kPa)",
    y     = "Transpiration per Leaf Area (ml/cm2)",
    title = "Transpiration per Leaf Area vs VPD by Genotype",
    caption = "Dashed line: VPD = 2.0 kPa (Sinclair 2008 threshold)"
  ) +
  base_theme

save_fig(figC, "GC3_Fig_C_TranspLeafArea", width_cm = 20, height_cm = 18)

# ── Fig D — WUE Quadrant plot ─────────────────────────────────────────────────
cat("Plotting Fig D...\n")

sum_D <- df %>%
  group_by(Quadrant, Plasticity, WUE_class, VPD_target) %>%
  summarise(
    mean_RT = mean(RelTranspiration, na.rm = TRUE),
    .groups = "drop"
  )

# Colour each quadrant distinctly
quad_cols <- c(Q1 = "#1b9e77", Q2 = "#d95f02", Q3 = "#7570b3", Q4 = "#e7298a")

figD <- ggplot(sum_D, aes(x = VPD_target, y = mean_RT,
                           colour = Quadrant, group = Quadrant)) +
  geom_hline(yintercept = 1.0, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = quad_cols,
                      labels = c(Q1 = "Q1: PI398617 (Plastic, High WUE)",
                                 Q2 = "Q2: PI070229 (Plastic, Low WUE)",
                                 Q3 = "Q3: PI423802 (Nonplastic, Low WUE)",
                                 Q4 = "Q4: PI603465C (Nonplastic, High WUE)")) +
  facet_grid(Plasticity ~ WUE_class) +
  labs(
    x      = "Target VPD (kPa)",
    y      = "Mean Relative Transpiration",
    title  = "VPD Response by Plasticity x WUE Class (Chamarthi et al. 2022 quadrants)",
    colour = "Quadrant"
  ) +
  base_theme +
  theme(legend.position = "bottom",
        legend.direction = "vertical")

save_fig(figD, "GC3_Fig_D_quadrant", width_cm = 20, height_cm = 14)

# ── Fig E — Dry matter by genotype x treatment ────────────────────────────────
cat("Plotting Fig E...\n")

# Unique pot-level rows (DryMatter_total is pot-level, not per-measurement)
dm_df <- df %>%
  distinct(Pot, Genotype, Expected_WUE, VPD_target, .keep_all = TRUE) %>%
  filter(!is.na(DryMatter_total))

cat("  Unique pot rows for dry matter:", nrow(dm_df), "\n")

figE <- ggplot(dm_df, aes(x = Genotype, y = DryMatter_total, fill = Expected_WUE)) +
  geom_boxplot(outlier.size = 1.5, alpha = 0.8) +
  scale_fill_brewer(palette = "Set2", name = "Expected WUE") +
  facet_wrap(~ VPD_target, labeller = labeller(VPD_target = function(x) paste("VPD:", x, "kPa"))) +
  labs(
    x     = "Genotype",
    y     = "Total Dry Matter (g)",
    title = "Dry Matter by Genotype and VPD Treatment"
  ) +
  base_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

save_fig(figE, "GC3_Fig_E_drymatter", width_cm = 22, height_cm = 14)

# ==============================================================================
# 3. STATISTICS
# ==============================================================================
cat("\n=== STATISTICS ===\n\n")

# ── 3a. Segmented (or piecewise) regression by Genotype x Day ─────────────────
has_segmented <- requireNamespace("segmented", quietly = TRUE)
if (has_segmented) {
  library(segmented)
  cat("Using 'segmented' package for breakpoint estimation\n\n")
} else {
  cat("'segmented' package not available — fitting piecewise lm at VPD = 2.0\n\n")
}

stats_list <- list()

for (geno in levels(df$Genotype)) {
  for (d in c(1, 2)) {
    sub <- df %>% filter(Genotype == geno, Day == d, !is.na(RelTranspiration), !is.na(VPD))

    if (nrow(sub) < 6) {
      cat("Skipping", geno, "Day", d, "-- insufficient data\n")
      next
    }

    row_base <- data.frame(Genotype = geno, Day = d,
                           n = nrow(sub),
                           stringsAsFactors = FALSE)

    seg_done <- FALSE
    if (has_segmented) {
      tryCatch({
        lm0 <- lm(RelTranspiration ~ VPD, data = sub)
        seg  <- segmented(lm0, seg.Z = ~VPD, psi = list(VPD = 2.0))
        bp   <- seg$psi["VPD", "Est."]
        sl   <- slope(seg)$VPD
        sl1  <- sl[1, "Est."]
        sl2  <- sl[2, "Est."]
        cat(sprintf("%-12s Day %d | Breakpoint: %.2f kPa | Slope below: %.3f | Slope above: %.3f\n",
                    geno, d, bp, sl1, sl2))
        row_base$method      <- "segmented"
        row_base$breakpoint  <- round(bp, 3)
        row_base$slope_below <- round(sl1, 4)
        row_base$slope_above <- round(sl2, 4)
        seg_done <- TRUE
      }, error = function(e) {
        cat(sprintf("%-12s Day %d | segmented failed: %s -- falling back to piecewise\n",
                    geno, d, conditionMessage(e)))
      })
    }

    if (!seg_done) {
      # Piecewise at fixed VPD = 2.0
      lo  <- sub %>% filter(VPD <= VPD_THRESHOLD)
      hi  <- sub %>% filter(VPD >  VPD_THRESHOLD)
      sl_lo <- if (nrow(lo) >= 2) coef(lm(RelTranspiration ~ VPD, data = lo))["VPD"] else NA
      sl_hi <- if (nrow(hi) >= 2) coef(lm(RelTranspiration ~ VPD, data = hi))["VPD"] else NA
      cat(sprintf("%-12s Day %d | Breakpoint: 2.00 kPa (fixed) | Slope below: %.3f | Slope above: %.3f\n",
                  geno, d,
                  ifelse(is.na(sl_lo), NA_real_, sl_lo),
                  ifelse(is.na(sl_hi), NA_real_, sl_hi)))
      row_base$method      <- "piecewise_fixed"
      row_base$breakpoint  <- 2.0
      row_base$slope_below <- round(sl_lo, 4)
      row_base$slope_above <- round(sl_hi, 4)
    }

    stats_list[[paste(geno, d)]] <- row_base
  }
}

stats_df <- bind_rows(stats_list)
write.csv(stats_df, "experimental_data/GC_analysis/GC3_stats_summary.csv", row.names = FALSE)
cat("\nSaved GC3_stats_summary.csv\n\n")

# ── 3b. ANOVA: RelTranspiration ~ Genotype x VPD + Day ───────────────────────
cat("=== ANOVA: RelTranspiration ~ Genotype * VPD + Day ===\n")

aov_mod <- aov(RelTranspiration ~ Genotype * VPD + Day, data = df)
aov_tab <- summary(aov_mod)
print(aov_tab)

# Also fit as lm for cleaner F-table
lm_mod  <- lm(RelTranspiration ~ Genotype * VPD + factor(Day), data = df)
cat("\n=== LM summary (RelTranspiration ~ Genotype * VPD + Day) ===\n")
print(summary(lm_mod))

# ==============================================================================
# DONE
# ==============================================================================
cat("\n=== All figures and statistics complete ===\n")
cat("Output directory: experimental_data/GC_analysis/\n")
cat("Files produced:\n")
cat("  GC3_clean.csv\n")
cat("  GC3_Fig_A_VPD_RelTransp.png/.tiff\n")
cat("  GC3_Fig_B_summary.png/.tiff\n")
cat("  GC3_Fig_C_TranspLeafArea.png/.tiff\n")
cat("  GC3_Fig_D_quadrant.png/.tiff\n")
cat("  GC3_Fig_E_drymatter.png/.tiff\n")
cat("  GC3_stats_summary.csv\n")
