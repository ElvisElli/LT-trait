## GC_transpiration_analysis.R
## Growth-chamber VPD vs. transpiration analysis — GC Experiment 3
## Target journal: Agricultural Water Management
## -------------------------------------------------------------

setwd("/home/user/LT-trait")

# ── 0. Packages ──────────────────────────────────────────────
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(RColorBrewer)
  library(segmented)
})

dir.create("experimental_data/GC_analysis", recursive = TRUE, showWarnings = FALSE)

# ── 1. Load & clean ──────────────────────────────────────────
DATA_PATH <- "experimental_data/data-analysis/data-analysis/Second time GC/Part 3/GC experiment 3.xlsx"

raw <- read_excel(DATA_PATH, sheet = "Final")

# VPD target labels
vpd_labels <- c("1" = "1.0", "2" = "1.5", "3" = "2.0", "4" = "2.7", "5" = "3.3")

gc <- raw %>%
  mutate(
    VPD_target = factor(vpd_labels[as.character(Treatment)],
                        levels = c("1.0", "1.5", "2.0", "2.7", "3.3")),
    # Parse Expected_WUE into Plasticity and WUE_class
    Plasticity = case_when(
      grepl("^Plastic,", Expected_WUE)    ~ "Plastic",
      grepl("^Nonplastic,", Expected_WUE) ~ "Nonplastic",
      TRUE ~ NA_character_
    ),
    WUE_class = case_when(
      grepl("High WUE", Expected_WUE) ~ "High WUE",
      grepl("Low WUE",  Expected_WUE) ~ "Low WUE",
      TRUE ~ NA_character_
    ),
    # Short genotype labels (used as secondary axis info)
    Genotype_label = Genotype,
    # Facet label: Genotype + Expected_WUE
    Facet_label = paste0(Genotype, "\n(", Expected_WUE, ")"),
    Day_label = paste("Day", Day),
    # Ensure numeric
    VPD             = as.numeric(VPD),
    RelTranspiration = as.numeric(RelTranspiration),
    TranspLeafArea  = as.numeric(TranspLeafArea),
    DryMatter_total = as.numeric(DryMatter_total),
    Replication     = as.integer(Replication),
    Day             = as.integer(Day)
  )

# Fix factor ordering for Facet_label (match Quadrant Q1-Q4 order)
quad_order <- c("Q1", "Q2", "Q3", "Q4")
geno_order <- gc %>%
  distinct(Genotype, Quadrant) %>%
  arrange(match(Quadrant, quad_order)) %>%
  pull(Genotype)

gc <- gc %>%
  mutate(
    Genotype = factor(Genotype, levels = geno_order),
    Facet_label = factor(Facet_label,
                         levels = gc %>%
                           distinct(Genotype, Facet_label) %>%
                           arrange(match(Genotype, geno_order)) %>%
                           pull(Facet_label))
  )

# Save clean CSV
write.csv(gc, "experimental_data/GC_analysis/GC3_clean.csv", row.names = FALSE)
cat("Clean data saved: experimental_data/GC_analysis/GC3_clean.csv\n")
cat("Rows:", nrow(gc), " | Cols:", ncol(gc), "\n\n")

# ── Common theme ─────────────────────────────────────────────
my_theme <- theme_bw(base_size = 12) +
  theme(
    strip.background   = element_rect(fill = "grey92", colour = "grey70"),
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom",
    legend.title       = element_text(face = "bold"),
    strip.text         = element_text(size = 10)
  )

plasticity_cols <- c("Plastic" = "#d95f02", "Nonplastic" = "#31a354")
geno_cols       <- brewer.pal(4, "Set1")
names(geno_cols) <- geno_order

save_fig <- function(plot, name, width_cm, height_cm) {
  base <- file.path("experimental_data/GC_analysis", name)
  ggsave(paste0(base, ".png"),  plot, width = width_cm, height = height_cm,
         units = "cm", dpi = 150)
  ggsave(paste0(base, ".tiff"), plot, width = width_cm, height = height_cm,
         units = "cm", dpi = 300, compression = "lzw")
  cat("Saved:", base, ".png/.tiff\n")
}

# ── Fig A — VPD vs Relative Transpiration ────────────────────
figA <- ggplot(gc, aes(x = VPD, y = RelTranspiration, colour = Plasticity)) +
  geom_jitter(alpha = 0.4, height = 0.01, width = 0, size = 1.8) +
  geom_smooth(method = "loess", se = TRUE, span = 0.85, linewidth = 0.9) +
  geom_vline(xintercept = 2.0, linetype = "dashed", colour = "grey40", linewidth = 0.6) +
  scale_colour_manual(values = plasticity_cols) +
  scale_x_continuous(breaks = c(1.0, 1.5, 2.0, 2.7, 3.3)) +
  scale_y_continuous(limits = c(0, NA)) +
  facet_wrap(~ Facet_label, ncol = 2) +
  labs(
    x      = "Measured VPD (kPa)",
    y      = "Relative Transpiration",
    colour = "Plasticity",
    title  = "Fig A — VPD vs Relative Transpiration by Genotype",
    caption = "Dashed line: VPD = 2.0 kPa reference threshold (Sinclair 2008)"
  ) +
  my_theme

save_fig(figA, "GC3_Fig_A_VPD_RelTransp", width_cm = 20, height_cm = 18)

# ── Fig B — Mean ± SE by VPD treatment × Genotype × Day ─────
sumB <- gc %>%
  group_by(Genotype, Day_label, VPD_target) %>%
  summarise(
    mean_RT = mean(RelTranspiration, na.rm = TRUE),
    se_RT   = sd(RelTranspiration,   na.rm = TRUE) / sqrt(sum(!is.na(RelTranspiration))),
    .groups = "drop"
  )

figB <- ggplot(sumB, aes(x = VPD_target, y = mean_RT,
                          colour = Genotype, group = Genotype)) +
  geom_hline(yintercept = 1.0, linetype = "dotted", colour = "grey50") +
  geom_errorbar(aes(ymin = mean_RT - se_RT, ymax = mean_RT + se_RT),
                width = 0.15, linewidth = 0.5) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = geno_cols) +
  scale_y_continuous(limits = c(0, NA)) +
  facet_wrap(~ Day_label, ncol = 2) +
  labs(
    x      = "Target VPD (kPa)",
    y      = "Mean Relative Transpiration ± SE",
    colour = "Genotype",
    title  = "Fig B — Summary: Mean ± SE by VPD Treatment × Genotype × Day"
  ) +
  my_theme

save_fig(figB, "GC3_Fig_B_summary", width_cm = 20, height_cm = 12)

# ── Fig C — Transpiration per leaf area vs VPD ───────────────
gcC <- gc %>% filter(!is.na(TranspLeafArea))

figC <- ggplot(gcC, aes(x = VPD, y = TranspLeafArea, colour = Plasticity)) +
  geom_jitter(alpha = 0.4, height = 0, width = 0, size = 1.8) +
  geom_smooth(method = "loess", se = TRUE, span = 0.85, linewidth = 0.9) +
  geom_vline(xintercept = 2.0, linetype = "dashed", colour = "grey40", linewidth = 0.6) +
  scale_colour_manual(values = plasticity_cols) +
  scale_x_continuous(breaks = c(1.0, 1.5, 2.0, 2.7, 3.3)) +
  facet_wrap(~ Facet_label, ncol = 2) +
  labs(
    x      = "Measured VPD (kPa)",
    y      = expression("Transpiration / Leaf Area (ml cm"^{-2}*")"),
    colour = "Plasticity",
    title  = "Fig C — VPD vs Transpiration per Leaf Area by Genotype",
    caption = "Dashed line: VPD = 2.0 kPa reference threshold"
  ) +
  my_theme

save_fig(figC, "GC3_Fig_C_TranspLeafArea", width_cm = 20, height_cm = 18)

# ── Fig D — WUE Quadrant plot ─────────────────────────────────
sumD <- gc %>%
  group_by(Quadrant, Plasticity, WUE_class, VPD_target) %>%
  summarise(mean_RT = mean(RelTranspiration, na.rm = TRUE), .groups = "drop")

quad_cols <- c("Q1" = "#e41a1c", "Q2" = "#377eb8", "Q3" = "#4daf4a", "Q4" = "#984ea3")

figD <- ggplot(sumD, aes(x = VPD_target, y = mean_RT,
                          colour = Quadrant, group = Quadrant)) +
  geom_hline(yintercept = 1.0, linetype = "dotted", colour = "grey50") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = quad_cols) +
  scale_y_continuous(limits = c(0, NA)) +
  facet_grid(Plasticity ~ WUE_class) +
  labs(
    x      = "Target VPD (kPa)",
    y      = "Mean Relative Transpiration",
    colour = "Quadrant",
    title  = "Fig D — WUE Quadrant: Plasticity × WUE class (Chamarthi et al. 2022 framework)"
  ) +
  my_theme +
  theme(legend.position = "right")

save_fig(figD, "GC3_Fig_D_quadrant", width_cm = 20, height_cm = 14)

# ── Fig E — Dry matter by genotype × treatment ───────────────
gcE <- gc %>%
  distinct(Pot, Genotype, Expected_WUE, VPD_target, DryMatter_total) %>%
  mutate(Expected_WUE = factor(Expected_WUE,
                                levels = c("Plastic, High WUE", "Plastic, Low WUE",
                                           "Nonplastic, High WUE", "Nonplastic, Low WUE")))

figE <- ggplot(gcE, aes(x = Genotype, y = DryMatter_total, fill = Expected_WUE)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1.5) +
  scale_fill_brewer(palette = "Set2", name = "Expected WUE") +
  facet_wrap(~ VPD_target, ncol = 5,
             labeller = labeller(VPD_target = function(x) paste("VPD", x, "kPa"))) +
  labs(
    x     = "Genotype",
    y     = "Total Dry Matter (g)",
    title = "Fig E — Total Dry Matter by Genotype × VPD Treatment"
  ) +
  my_theme +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

save_fig(figE, "GC3_Fig_E_drymatter", width_cm = 24, height_cm = 12)

# ── 3. Statistics ─────────────────────────────────────────────
cat("\n═══════════════════════════════════════════\n")
cat("STATISTICAL RESULTS\n")
cat("═══════════════════════════════════════════\n\n")

## A) Segmented regression: RelTranspiration ~ VPD per Genotype × Day
seg_results <- list()

for (geno in levels(gc$Genotype)) {
  for (day in c(1, 2)) {
    sub <- gc %>%
      filter(Genotype == geno, Day == day) %>%
      filter(!is.na(VPD), !is.na(RelTranspiration))

    if (nrow(sub) < 6) next

    cat(sprintf("─── %s  Day %d (n=%d) ───\n", geno, day, nrow(sub)))

    seg_fit <- tryCatch({
      lm0 <- lm(RelTranspiration ~ VPD, data = sub)
      seg  <- segmented(lm0, seg.Z = ~VPD,
                        psi = list(VPD = 2.0),
                        control = seg.control(display = FALSE, it.max = 100))
      bp   <- seg$psi[, "Est."]
      sl   <- slope(seg)$VPD
      cat(sprintf("  Breakpoint: %.3f kPa\n", bp))
      cat(sprintf("  Slope below BP: %.4f  |  above BP: %.4f\n",
                  sl[1, "Est."], sl[2, "Est."]))
      list(method = "segmented", genotype = geno, day = day,
           breakpoint = bp,
           slope_below = sl[1, "Est."], slope_above = sl[2, "Est."],
           se_below = sl[1, "Std. Err."], se_above = sl[2, "Std. Err."])
    }, error = function(e) {
      # Fallback: two separate linear models split at VPD = 2.0
      cat(sprintf("  [segmented failed: %s] — using two-lm fallback\n", conditionMessage(e)))
      sub_lo <- filter(sub, VPD <= 2.0)
      sub_hi <- filter(sub, VPD >  2.0)
      sl_lo  <- if (nrow(sub_lo) > 1) coef(lm(RelTranspiration ~ VPD, data = sub_lo))["VPD"] else NA
      sl_hi  <- if (nrow(sub_hi) > 1) coef(lm(RelTranspiration ~ VPD, data = sub_hi))["VPD"] else NA
      cat(sprintf("  Slope ≤ 2.0 kPa: %.4f  |  > 2.0 kPa: %.4f\n",
                  ifelse(is.na(sl_lo), NA, sl_lo),
                  ifelse(is.na(sl_hi), NA, sl_hi)))
      list(method = "two_lm", genotype = geno, day = day,
           breakpoint = 2.0,
           slope_below = sl_lo, slope_above = sl_hi,
           se_below = NA, se_above = NA)
    })

    seg_results[[paste(geno, day)]] <- seg_fit
    cat("\n")
  }
}

seg_df <- do.call(rbind, lapply(seg_results, function(x) {
  data.frame(Genotype    = x$genotype,
             Day         = x$day,
             Method      = x$method,
             Breakpoint  = round(x$breakpoint, 3),
             Slope_below = round(x$slope_below, 4),
             Slope_above = round(x$slope_above, 4),
             SE_below    = round(x$se_below, 4),
             SE_above    = round(x$se_above, 4),
             stringsAsFactors = FALSE)
}))
rownames(seg_df) <- NULL

cat("Segmented regression summary:\n")
print(seg_df)

## B) ANOVA: RelTranspiration ~ Genotype × VPD + Day
cat("\n═══════════════════════════════════════════\n")
cat("ANOVA: RelTranspiration ~ Genotype * VPD + Day\n")
cat("═══════════════════════════════════════════\n")
aov_mod  <- aov(RelTranspiration ~ Genotype * VPD + as.factor(Day), data = gc)
aov_tbl  <- summary(aov_mod)
print(aov_tbl)

## Save stats CSV
write.csv(seg_df, "experimental_data/GC_analysis/GC3_stats_summary.csv", row.names = FALSE)
cat("\nStats saved: experimental_data/GC_analysis/GC3_stats_summary.csv\n")

# ── 4. Descriptive summary ────────────────────────────────────
cat("\n═══════════════════════════════════════════\n")
cat("DESCRIPTIVE SUMMARY\n")
cat("═══════════════════════════════════════════\n\n")
desc <- gc %>%
  group_by(Genotype, Expected_WUE) %>%
  summarise(
    n          = n(),
    mean_RT    = round(mean(RelTranspiration, na.rm = TRUE), 3),
    sd_RT      = round(sd(RelTranspiration,   na.rm = TRUE), 3),
    mean_Tla   = round(mean(TranspLeafArea,   na.rm = TRUE), 2),
    mean_DM    = round(mean(DryMatter_total,  na.rm = TRUE), 2),
    .groups = "drop"
  )
print(desc)

cat("\n── Analysis complete ──\n")
