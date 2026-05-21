## GC_critical_analysis.R
## Critical analysis of Growth Chamber transpiration experiments (GC1 + GC3)
## and field stomatal conductance data.
## Outputs saved to: experimental_data/GC_analysis/

## ── Path configuration ────────────────────────────────────────────────────────
## Windows (local): edit and uncomment the line below
# project_root <- "C:/Users/YourName/Projects/LT-trait"
if (.Platform$OS.type == "windows") setwd(project_root) else setwd("/home/user/LT-trait")
## ─────────────────────────────────────────────────────────────────────────────

# ─── packages ─────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(purrr)
})

outdir <- "experimental_data/GC_analysis"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ─── style constants ──────────────────────────────────────────────────────────
col_exp  <- c(GC1 = "#1f78b4", GC3 = "#e31a1c")
col_day  <- c(`1` = "grey50",  `2` = "#2171b5")
col_quad <- c(Q1 = "#d95f02", Q2 = "#1b9e77", Q3 = "#7570b3", Q4 = "#e7298a")
base_theme <- theme_bw(base_size = 12) +
  theme(
    strip.background  = element_rect(fill = "grey92", colour = NA),
    panel.grid.minor  = element_blank(),
    legend.position   = "bottom"
  )

save_fig <- function(plot_obj, stem, w_cm, h_cm) {
  ggsave(file.path(outdir, paste0(stem, ".png")),
         plot_obj, width = w_cm, height = h_cm, units = "cm", dpi = 300)
  ggsave(file.path(outdir, paste0(stem, ".tiff")),
         plot_obj, width = w_cm, height = h_cm, units = "cm", dpi = 300,
         compression = "lzw")
  invisible(NULL)
}

# ─── helper: clean Expected_WUE ───────────────────────────────────────────────
clean_wue <- function(x) {
  # forward-fill NA within group (some rows have NA after first row)
  x[is.na(x)] <- NA_character_
  zoo_fill <- function(v) {
    last_val <- NA_character_
    for (i in seq_along(v)) {
      if (!is.na(v[i])) last_val <- v[i] else v[i] <- last_val
    }
    v
  }
  zoo_fill(x)
}

# ══════════════════════════════════════════════════════════════════════════════
# STEP 1 — Load and combine GC1 + GC3
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── Step 1: Loading data ──────────────────────────────────────────────────\n")

gc1_raw <- read_excel(
  "experimental_data/data-analysis/data-analysis/transpiration/Data analysis.xlsx",
  sheet = "Final-ee"
)
gc3_raw <- read_excel(
  "experimental_data/data-analysis/data-analysis/Second time GC/Part 3/GC experiment 3.xlsx",
  sheet = "Final"
)

cat("GC1:", nrow(gc1_raw), "rows x", ncol(gc1_raw), "cols\n")
cat("GC3:", nrow(gc3_raw), "rows x", ncol(gc3_raw), "cols\n")

# Mean VPD per treatment (for labelling)
vpd_label_gc1 <- gc1_raw %>%
  group_by(Treatment) %>%
  summarise(VPD_target = round(mean(VPD, na.rm = TRUE), 2), .groups = "drop")
vpd_label_gc3 <- gc3_raw %>%
  group_by(Treatment) %>%
  summarise(VPD_target = round(mean(VPD, na.rm = TRUE), 2), .groups = "drop")

cat("\nGC1 mean VPD by treatment:\n"); print(vpd_label_gc1)
cat("\nGC3 mean VPD by treatment:\n"); print(vpd_label_gc3)

# Plasticity and WUE class lookup
plasticity_wue <- function(quadrant) {
  case_when(
    quadrant == "Q1" ~ list(plasticity = "Plastic",    wue_class = "High WUE"),
    quadrant == "Q2" ~ list(plasticity = "Plastic",    wue_class = "Low WUE"),
    quadrant == "Q3" ~ list(plasticity = "Nonplastic", wue_class = "Low WUE"),
    quadrant == "Q4" ~ list(plasticity = "Nonplastic", wue_class = "High WUE"),
    TRUE             ~ list(plasticity = NA_character_, wue_class = NA_character_)
  )
}

prep_gc <- function(df, exp_label, vpd_lut) {
  df %>%
    mutate(
      Experiment           = exp_label,
      Day                  = as.character(Day),
      StomatalConductance  = if ("StomatalConductance" %in% names(.)) {
                               sc_v <- suppressWarnings(as.numeric(StomatalConductance))
                               # Plausible range for stomatal conductance: 0–2 mol m-2 s-1
                               # Data entry errors detected at rows with values >10 or <0
                               sc_v[sc_v > 10 | sc_v < 0] <- NA_real_
                               sc_v
                             } else NA_real_,
      E_apparent           = if ("E_apparent" %in% names(.))
                               suppressWarnings(as.numeric(E_apparent))
                             else NA_real_,
      T_leaf_degreecelcius = if ("T_leaf_degreecelcius" %in% names(.))
                               suppressWarnings(as.numeric(T_leaf_degreecelcius))
                             else NA_real_,
      DryMatter_rep        = if ("DryMatter_rep" %in% names(.))
                               DryMatter_rep else NA_real_
    ) %>%
    left_join(vpd_lut, by = "Treatment") %>%
    # forward-fill Expected_WUE where NA
    group_by(Genotype) %>%
    mutate(Expected_WUE = {
      v <- Expected_WUE
      last_v <- NA_character_
      for (i in seq_along(v)) {
        if (!is.na(v[i])) last_v <- v[i] else v[i] <- last_v
      }
      v
    }) %>%
    ungroup() %>%
    mutate(
      Plasticity = case_when(
        Quadrant == "Q1" | Quadrant == "Q2" ~ "Plastic",
        Quadrant == "Q3" | Quadrant == "Q4" ~ "Nonplastic",
        TRUE ~ NA_character_
      ),
      WUE_class = case_when(
        Quadrant == "Q1" | Quadrant == "Q4" ~ "High WUE",
        Quadrant == "Q2" | Quadrant == "Q3" ~ "Low WUE",
        TRUE ~ NA_character_
      )
    )
}

gc1 <- prep_gc(gc1_raw, "GC1", vpd_label_gc1)
gc3 <- prep_gc(gc3_raw, "GC3", vpd_label_gc3)

# Harmonise columns for binding
all_cols <- union(names(gc1), names(gc3))
add_missing <- function(df, cols) {
  missing <- setdiff(cols, names(df))
  for (col in missing) df[[col]] <- NA_real_
  df
}
gc1 <- add_missing(gc1, all_cols)
gc3 <- add_missing(gc3, all_cols)

combined <- bind_rows(gc1, gc3) %>%
  mutate(
    Experiment  = factor(Experiment, levels = c("GC1", "GC3")),
    Day         = factor(Day, levels = c("1", "2")),
    Quadrant    = factor(Quadrant, levels = c("Q1", "Q2", "Q3", "Q4")),
    Treatment   = as.integer(Treatment),
    VPD_target  = as.numeric(VPD_target)
  )

cat("\nCombined data: ", nrow(combined), "rows\n")
cat("Experiments:   ", levels(combined$Experiment), "\n")
cat("Genotypes GC1: ", paste(unique(gc1$Genotype), collapse=", "), "\n")
cat("Genotypes GC3: ", paste(unique(gc3$Genotype), collapse=", "), "\n")

write.csv(combined, file.path(outdir, "GC_combined_clean.csv"), row.names = FALSE)
cat("Saved: GC_combined_clean.csv\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2 — Day 1 investigation and outlier detection
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── Step 2: Outlier detection ─────────────────────────────────────────────\n")

combined <- combined %>%
  group_by(Experiment, Genotype, Day, Treatment) %>%
  mutate(
    grp_mean  = mean(RelTranspiration, na.rm = TRUE),
    grp_sd    = sd(RelTranspiration, na.rm = TRUE),
    is_outlier = !is.na(RelTranspiration) &
                   !is.na(grp_sd) &
                   abs(RelTranspiration - grp_mean) > 2 * grp_sd
  ) %>%
  ungroup()

outlier_summary <- combined %>%
  group_by(Experiment, Day) %>%
  summarise(
    n_total   = n(),
    n_outlier = sum(is_outlier, na.rm = TRUE),
    pct       = round(100 * n_outlier / n_total, 1),
    .groups   = "drop"
  )

cat("\nOutliers per Experiment × Day:\n")
print(outlier_summary)

# Diagnostic plot: Day 1 vs Day 2, outliers in red
p_diag <- ggplot(combined,
                 aes(x = VPD, y = RelTranspiration,
                     colour = is_outlier, shape = is_outlier)) +
  geom_point(alpha = 0.7, size = 1.8) +
  scale_colour_manual(values = c(`FALSE` = "grey40", `TRUE` = "red"),
                      labels = c("Normal", "Outlier (>2 SD)"),
                      name = NULL) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17),
                     labels = c("Normal", "Outlier (>2 SD)"),
                     name = NULL) +
  facet_grid(Day ~ Experiment + Genotype,
             labeller = labeller(Day = c(`1` = "Day 1", `2` = "Day 2"))) +
  geom_vline(xintercept = 2.0, linetype = "dashed", colour = "grey30", linewidth = 0.5) +
  labs(title = "Day 1 vs Day 2: RelTranspiration with outlier flags",
       subtitle = paste("Red = |obs - group mean| > 2 SD |",
                        "Groups: Experiment × Genotype × Day × Treatment"),
       x = "VPD (kPa)", y = "Relative transpiration") +
  base_theme
save_fig(p_diag, "GC_Fig_diag_Day1vsDay2", w_cm = 32, h_cm = 16)
cat("Saved: GC_Fig_diag_Day1vsDay2.png/.tiff\n")

# Day 1 vs Day 2 variability summary
day_var <- combined %>%
  group_by(Experiment, Day) %>%
  summarise(
    mean_T  = mean(RelTranspiration, na.rm = TRUE),
    sd_T    = sd(RelTranspiration, na.rm = TRUE),
    cv_T    = sd_T / mean_T,
    .groups = "drop"
  )
cat("\nDay 1 vs Day 2 variability:\n")
print(day_var)

# ── filter non-outlier data for main analyses
clean <- combined %>% filter(!is_outlier)
cat("\nRows after removing outliers:", nrow(clean), "of", nrow(combined), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3 — Main VPD vs RelTranspiration figure (all genotypes, both experiments)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── Step 3: Main VPD vs RelTranspiration figure ───────────────────────────\n")

# Quadrant label annotation
quad_labels <- clean %>%
  distinct(Experiment, Genotype, Quadrant, WUE_class, Plasticity) %>%
  mutate(label = paste0(Quadrant, "\n", Plasticity, "\n", WUE_class),
         VPD = 1.05, RelTranspiration = 0.98)

p_main <- ggplot(clean,
                 aes(x = VPD, y = RelTranspiration, colour = Day)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(aes(group = Day), method = "lm", se = TRUE,
              linewidth = 0.8, alpha = 0.2) +
  geom_vline(xintercept = 2.0, linetype = "dashed",
             colour = "grey30", linewidth = 0.5) +
  geom_text(data = quad_labels,
            mapping = aes(x = VPD, y = RelTranspiration, label = label),
            colour = "grey25",
            size = 2.8, hjust = 0, vjust = 1, lineheight = 0.9,
            inherit.aes = FALSE) +
  scale_colour_manual(values = col_day, name = "Day",
                      labels = c(`1` = "Day 1", `2` = "Day 2")) +
  facet_grid(Experiment ~ Genotype, scales = "free_x") +
  labs(title = "VPD vs Relative Transpiration — all genotypes, both experiments",
       x = "VPD (kPa)", y = "Relative transpiration") +
  base_theme
save_fig(p_main, "GC_Fig_main_allgeno", w_cm = 26, h_cm = 14)
cat("Saved: GC_Fig_main_allgeno.png/.tiff\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 4 — VPD vs RelTranspiration by Quadrant (pooled, Day 2 only)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── Step 4: Quadrant pooled figure (Day 2 only) ──────────────────────────\n")

day2_clean <- clean %>% filter(Day == "2")

p_quad <- ggplot(day2_clean,
                 aes(x = VPD, y = RelTranspiration, colour = Experiment)) +
  geom_point(alpha = 0.4, size = 1.8) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1, alpha = 0.2,
              aes(group = Experiment, fill = Experiment)) +
  geom_vline(xintercept = 2.0, linetype = "dashed",
             colour = "grey30", linewidth = 0.5) +
  scale_colour_manual(values = col_exp, name = "Experiment") +
  scale_fill_manual(values   = col_exp, name = "Experiment") +
  facet_wrap(~ Quadrant, nrow = 2, ncol = 2,
             labeller = labeller(Quadrant = c(
               Q1 = "Q1: Plastic / High WUE",
               Q2 = "Q2: Plastic / Low WUE",
               Q3 = "Q3: Nonplastic / Low WUE",
               Q4 = "Q4: Nonplastic / High WUE"
             ))) +
  labs(title = "VPD vs Relative Transpiration by Quadrant — Day 2, pooled GC1 & GC3",
       subtitle = "Cross-experiment validation: different genotypes, same quadrant classification",
       x = "VPD (kPa)", y = "Relative transpiration") +
  base_theme
save_fig(p_quad, "GC_Fig_quadrant_pooled", w_cm = 22, h_cm = 18)
cat("Saved: GC_Fig_quadrant_pooled.png/.tiff\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 5 — Stomatal conductance analysis (GC1 only)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── Step 5: Stomatal conductance analysis (GC1) ──────────────────────────\n")

gc1_sc <- clean %>%
  filter(Experiment == "GC1", !is.na(StomatalConductance))

cat("GC1 rows with SC data:", nrow(gc1_sc), "\n")
cat("SC range:", round(range(gc1_sc$StomatalConductance), 4), "\n")

# Plot 1: VPD vs SC by Genotype
p_sc_vpd <- ggplot(gc1_sc,
                   aes(x = VPD, y = StomatalConductance, colour = Day)) +
  geom_point(alpha = 0.6, size = 1.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, alpha = 0.2,
              aes(group = Day)) +
  scale_colour_manual(values = col_day, name = "Day",
                      labels = c(`1` = "Day 1", `2` = "Day 2")) +
  geom_vline(xintercept = 2.0, linetype = "dashed",
             colour = "grey30", linewidth = 0.5) +
  facet_wrap(~ Genotype, nrow = 1) +
  labs(title = "VPD vs Stomatal Conductance — GC1",
       x = "VPD (kPa)", y = "Stomatal conductance (mol m⁻² s⁻¹)") +
  base_theme
save_fig(p_sc_vpd, "GC_Fig_SC_vs_VPD", w_cm = 26, h_cm = 10)
cat("Saved: GC_Fig_SC_vs_VPD.png/.tiff\n")

# Plot 2: SC vs RelTranspiration, colour by VPD treatment
p_sc_transp <- ggplot(gc1_sc,
                      aes(x = StomatalConductance, y = RelTranspiration,
                          colour = VPD_target)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              colour = "grey40", linewidth = 0.5) +
  geom_point(alpha = 0.6, size = 1.8) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey20",
              linewidth = 0.8, alpha = 0.2, aes(group = 1)) +
  scale_colour_viridis_c(name = "VPD target\n(kPa)", option = "plasma") +
  facet_wrap(~ Genotype, nrow = 1) +
  labs(title = "Stomatal Conductance vs Relative Transpiration — GC1",
       subtitle = "Dashed line = 1:1 reference; colour = VPD level",
       x = "Stomatal conductance (mol m⁻² s⁻¹)", y = "Relative transpiration") +
  base_theme
save_fig(p_sc_transp, "GC_Fig_SC_vs_Transp", w_cm = 26, h_cm = 10)
cat("Saved: GC_Fig_SC_vs_Transp.png/.tiff\n")

# Pearson correlation SC vs RelTranspiration per Genotype
cor_table <- gc1_sc %>%
  group_by(Genotype, Quadrant, WUE_class) %>%
  summarise(
    n          = n(),
    r_pearson  = cor(StomatalConductance, RelTranspiration, use = "complete.obs"),
    p_value    = tryCatch(
      cor.test(StomatalConductance, RelTranspiration)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(across(c(r_pearson, p_value), ~ round(.x, 4)))

cat("\nPearson correlation SC vs RelTranspiration per Genotype (GC1):\n")
print(cor_table)

# Also overall correlation
r_all <- cor(gc1_sc$StomatalConductance, gc1_sc$RelTranspiration, use = "complete.obs")
cat("\nOverall Pearson r (GC1 all genotypes):", round(r_all, 4), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 6 — Piecewise regression (Day 2 only, non-outliers)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── Step 6: Piecewise regression (Day 2, non-outliers) ───────────────────\n")

vpd_break <- 2.0

piece_results <- day2_clean %>%
  group_by(Experiment, Genotype, Quadrant) %>%
  group_map(function(df, keys) {
    if (nrow(df) < 3) return(NULL)

    # Overall linear model
    m_all <- lm(RelTranspiration ~ VPD, data = df)
    r2_all   <- summary(m_all)$r.squared
    slope_all <- coef(m_all)[["VPD"]]

    # Below breakpoint
    df_lo <- df %>% filter(VPD < vpd_break)
    df_hi <- df %>% filter(VPD >= vpd_break)

    fit_seg <- function(d, tag) {
      if (nrow(d) < 3) return(data.frame(slope=NA, r2=NA, n=nrow(d)))
      m <- lm(RelTranspiration ~ VPD, data = d)
      data.frame(slope = coef(m)[["VPD"]], r2 = summary(m)$r.squared, n = nrow(d))
    }
    lo <- fit_seg(df_lo, "lo")
    hi <- fit_seg(df_hi, "hi")

    data.frame(
      Experiment  = keys$Experiment,
      Genotype    = keys$Genotype,
      Quadrant    = keys$Quadrant,
      n_total     = nrow(df),
      slope_all   = round(slope_all, 4),
      r2_all      = round(r2_all, 4),
      n_below     = lo$n,
      slope_below = round(lo$slope, 4),
      r2_below    = round(lo$r2, 4),
      n_above     = hi$n,
      slope_above = round(hi$slope, 4),
      r2_above    = round(hi$r2, 4),
      ratio_above_below = round(hi$slope / lo$slope, 4)
    )
  }) %>%
  bind_rows()

cat("\nPiecewise regression results (VPD breakpoint = 2.0 kPa):\n")
print(piece_results)

write.csv(piece_results, file.path(outdir, "GC_stats_piecewise_day2.csv"), row.names = FALSE)
cat("Saved: GC_stats_piecewise_day2.csv\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 7 — Field SC summary
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── Step 7: Field stomatal conductance analysis ───────────────────────────\n")

sc_field <- read_excel(
  "experimental_data/SC_vs_T/SC_vs_T/soybean_dataset_SC.xlsx",
  sheet = "dates_ee"
) %>%
  mutate(
    Date_label = format(Date, "%b %d"),
    Irrig      = factor(Irrig, levels = c("Rainfed", "Irrigated")),
    Genotype   = as.character(Genotype)
  )

cat("Field SC rows:", nrow(sc_field), "\n")
cat("Genotypes:", paste(sort(unique(sc_field$Genotype)), collapse=", "), "\n")
cat("Irrig levels:", levels(sc_field$Irrig), "\n")
cat("Dates:", paste(sort(unique(sc_field$Date_label)), collapse=", "), "\n")
cat("SC range:", round(range(sc_field$StomatalConductance, na.rm=TRUE), 3), "\n")

# Plot 1: SC by Genotype × Irrig × Date (boxplot)
p_field1 <- ggplot(sc_field %>% filter(!is.na(StomatalConductance)),
                   aes(x = Genotype, y = StomatalConductance,
                       fill = Irrig, colour = Irrig)) +
  geom_boxplot(alpha = 0.4, outlier.size = 1.2, outlier.alpha = 0.6,
               position = position_dodge(0.75), width = 0.6) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
              alpha = 0.4, size = 1.2) +
  scale_fill_manual(values   = c(Rainfed = "#fc8d59", Irrigated = "#91bfdb"),
                    name = "Water regime") +
  scale_colour_manual(values = c(Rainfed = "#d73027", Irrigated = "#4575b4"),
                      name = "Water regime") +
  facet_wrap(~ Date_label, nrow = 1) +
  labs(title = "Field Stomatal Conductance by Genotype, Irrigation and Date",
       x = "Genotype", y = "Stomatal conductance (mol m⁻² s⁻¹)") +
  base_theme +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 9))
save_fig(p_field1, "SC_field_Fig_overview", w_cm = 28, h_cm = 14)

# Plot 2: SC vs Yield
p_field2 <- ggplot(sc_field %>% filter(!is.na(StomatalConductance), !is.na(Yield)),
                   aes(x = StomatalConductance, y = Yield,
                       colour = Irrig)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, alpha = 0.2,
              aes(group = Irrig)) +
  scale_colour_manual(values = c(Rainfed = "#d73027", Irrigated = "#4575b4"),
                      name = "Water regime") +
  facet_wrap(~ Date_label, nrow = 1) +
  labs(title = "Stomatal Conductance vs Yield (Field)",
       x = "Stomatal conductance (mol m⁻² s⁻¹)", y = "Yield (kg ha⁻¹)") +
  base_theme
save_fig(p_field2, "SC_field_Fig_SC_vs_Yield", w_cm = 22, h_cm = 12)

# Plot 3: WiltingScore vs SC
p_field3 <- ggplot(sc_field %>% filter(!is.na(StomatalConductance), !is.na(WltingScore)),
                   aes(x = StomatalConductance, y = WltingScore,
                       colour = Genotype)) +
  geom_point(alpha = 0.65, size = 2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7,
              aes(group = 1), colour = "grey30") +
  scale_colour_brewer(palette = "Paired", name = "Genotype") +
  labs(title = "Wilting Score vs Stomatal Conductance (Field)",
       subtitle = "Lower SC associated with higher wilting under drought",
       x = "Stomatal conductance (mol m⁻² s⁻¹)", y = "Wilting score") +
  base_theme
save_fig(p_field3, "SC_field_Fig_Wilting_vs_SC", w_cm = 18, h_cm = 12)

cat("Saved: SC_field_Fig_overview.png/.tiff\n")
cat("Saved: SC_field_Fig_SC_vs_Yield.png/.tiff\n")
cat("Saved: SC_field_Fig_Wilting_vs_SC.png/.tiff\n")

# Field SC summary table
sc_summary <- sc_field %>%
  group_by(Genotype, Irrig, Date_label) %>%
  summarise(
    n       = n(),
    SC_mean = round(mean(StomatalConductance, na.rm=TRUE), 3),
    SC_sd   = round(sd(StomatalConductance, na.rm=TRUE), 3),
    Yield_mean = round(mean(Yield, na.rm=TRUE), 1),
    .groups = "drop"
  )
cat("\nField SC summary by Genotype × Irrig × Date:\n")
print(sc_summary)

# ══════════════════════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ══════════════════════════════════════════════════════════════════════════════
cat("\n══════════════════════════════════════════════════════════════════════════\n")
cat("ANALYSIS COMPLETE — key findings:\n")
cat("══════════════════════════════════════════════════════════════════════════\n\n")

cat("1. OUTLIERS (|obs - group mean| > 2 SD):\n")
print(outlier_summary)

cat("\n2. Day 1 vs Day 2 variability (CV of RelTranspiration):\n")
print(day_var %>% select(Experiment, Day, mean_T, sd_T, cv_T))

cat("\n3. Piecewise regression slopes (Day 2 non-outliers):\n")
print(piece_results %>% select(Experiment, Genotype, Quadrant,
                                slope_below, slope_above, ratio_above_below, r2_all))

cat("\n4. SC vs RelTranspiration correlations (GC1):\n")
print(cor_table)

cat("\n5. Field SC summary (Irrigated vs Rainfed):\n")
print(sc_summary %>% filter(Date_label == "Aug 09") %>%
      select(Genotype, Irrig, SC_mean, Yield_mean))

cat("\nAll files saved to:", outdir, "\n")
cat("Script completed successfully.\n")
