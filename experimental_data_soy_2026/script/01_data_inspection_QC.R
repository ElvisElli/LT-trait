# =============================================================================
# 01_data_inspection_QC.R
# Growth Chamber 2026 — Soybean Limited-Transpiration Experiment
#
# Purpose: Load, inspect, and quality-check data_all.xlsx
# Output:  figures/QC_*.png  (8 diagnostic plots)
#          data/data_clean.csv  (cleaned dataset ready for modelling)
#
# Run from project root (experimental_data_soy_2026/):
#   source("script/01_data_inspection_QC.R")
# =============================================================================

library(readxl)
library(tidyverse)   # dplyr + ggplot2 + tidyr
library(ggbeeswarm)  # jitter-free dot plots
library(patchwork)   # multi-panel layout


# --------------------------------------------------------------------------- #
# 1. LOAD DATA
# --------------------------------------------------------------------------- #

raw <- read_excel("data/data_all.xlsx", sheet = "dataAll")
cat("Raw data: ", nrow(raw), "rows x", ncol(raw), "cols\n\n")


# --------------------------------------------------------------------------- #
# 2. STRUCTURE CHECK
# --------------------------------------------------------------------------- #

# Variable types
str(raw)

# Basic summary
summary(raw)


# --------------------------------------------------------------------------- #
# 3. RENAME / RECODE FOR CLARITY
# --------------------------------------------------------------------------- #

# Treatment → numeric VPD level (1–6, low-to-high)
# avg_vpd is the measured chamber VPD for each observation row

df <- raw %>%
  rename(
    Rep        = RepNo,
    Day        = day,
    TR         = TR_mg_sqm_seg,   # transpiration rate (mg m-2 s-1)
    RelTR      = RelativeTransp,  # transpiration normalized to lowest VPD
    VPD        = avg_vpd,         # measured VPD (kPa)
    Temp_C     = avg_temperature,
    RH_pct     = avg_rh,
    LeafArea   = LeafArea_squarecm
  ) %>%
  mutate(
    # Factor levels
    Rep      = factor(Rep),
    Day      = factor(Day, levels = c(1, 2), labels = c("Day 1", "Day 2")),
    Genotype = factor(Genotype),
    WUE_class = factor(WUE_class, levels = c(
      "slow_wilting", "high_wilting",
      "Plastic, High WUE", "Nonplastic, Low WUE", "MAGIC parent"
    )),
    Treatment = factor(Treatment),
    # Average stomatal conductance across the 3 LICOR readings
    gsw = rowMeans(cbind(gsw1, gsw2, gsw3), na.rm = FALSE)
  ) %>%
  # Remove the trailing NULL column (col 31 from Excel)
  select(Rep, Day, Module, Scale, Genotype, WUE_class, Treatment,
         VPD, TR, RelTR, gsw, gsw1, gsw2, gsw3,
         Tleaf1, Tleaf2, Tleaf3, LeafArea,
         Temp_C, RH_pct, Initial_w, Final_w,
         Transpiration_g, Transpiration_mg, treatment_time_mins)

cat("Cleaned columns:\n")
print(names(df))


# --------------------------------------------------------------------------- #
# 4. MISSING VALUES
# --------------------------------------------------------------------------- #

cat("\n--- Missing value count per column ---\n")
miss <- df %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
  filter(n_missing > 0) %>%
  arrange(desc(n_missing))
print(miss)

# gsw is only available for a subset of plant × measurement combinations
cat("\ngsw data availability by Day:\n")
df %>%
  group_by(Day) %>%
  summarise(
    n_total   = n(),
    n_gsw     = sum(!is.na(gsw)),
    pct_gsw   = round(100 * n_gsw / n_total, 1)
  ) %>%
  print()


# --------------------------------------------------------------------------- #
# 5. SAMPLE DESIGN SUMMARY
# --------------------------------------------------------------------------- #

cat("\n--- Observations per Genotype × Day ---\n")
df %>%
  count(Genotype, Day) %>%
  pivot_wider(names_from = Day, values_from = n, values_fill = 0) %>%
  print(n = 20)

cat("\n--- Genotype WUE class ---\n")
df %>%
  distinct(Genotype, WUE_class) %>%
  arrange(WUE_class, Genotype) %>%
  print(n = 20)

cat("\n--- VPD levels achieved per Treatment ---\n")
df %>%
  group_by(Treatment) %>%
  summarise(
    n      = n(),
    VPD_mean = round(mean(VPD, na.rm = TRUE), 3),
    VPD_min  = round(min(VPD,  na.rm = TRUE), 3),
    VPD_max  = round(max(VPD,  na.rm = TRUE), 3)
  ) %>%
  print()


# --------------------------------------------------------------------------- #
# 6. OUTLIER DETECTION
# --------------------------------------------------------------------------- #

# Flag observations where RelTR > 1 (transpiration exceeds Day 1 T1 baseline)
cat("\n--- RelTR > 1 (possible calibration issue) ---\n")
df %>%
  filter(RelTR > 1) %>%
  select(Rep, Day, Genotype, Treatment, VPD, TR, RelTR) %>%
  print()

# IQR-based outlier flag for TR within each Genotype × Day × Treatment cell
df <- df %>%
  group_by(Genotype, Day, Treatment) %>%
  mutate(
    TR_Q1    = quantile(TR, 0.25, na.rm = TRUE),
    TR_Q3    = quantile(TR, 0.75, na.rm = TRUE),
    TR_IQR   = TR_Q3 - TR_Q1,
    TR_outlier = TR < (TR_Q1 - 3 * TR_IQR) | TR > (TR_Q3 + 3 * TR_IQR)
  ) %>%
  ungroup()

cat("\n--- Potential outliers (TR, 3×IQR rule) ---\n")
df %>%
  filter(TR_outlier) %>%
  select(Rep, Day, Genotype, Treatment, VPD, TR, RelTR) %>%
  print()

# Leaf area sanity check (should be > 0 and < ~2000 cm2 for soybean seedlings)
cat("\n--- Leaf area summary ---\n")
summary(df$LeafArea)
cat("Rows with LeafArea <= 0:", sum(df$LeafArea <= 0, na.rm = TRUE), "\n")


# --------------------------------------------------------------------------- #
# 7. PLOTS
# --------------------------------------------------------------------------- #

dir.create("figures", showWarnings = FALSE)
FIG <- "figures"


# --- QC 1: TR distribution by Genotype × Day ---------------------------------
p_tr_dist <- ggplot(df, aes(x = Genotype, y = TR, fill = Day)) +
  geom_boxplot(outlier.size = 0.8, position = position_dodge(0.8)) +
  scale_fill_manual(values = c("Day 1" = "#2166ac", "Day 2" = "#d6604d")) +
  labs(
    title = "QC1: Transpiration rate distribution by genotype and day",
    x     = NULL,
    y     = expression(TR~(mg~m^{-2}~s^{-1}))
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(FIG, "QC1_TR_distribution.png"),
       p_tr_dist, width = 12, height = 5, dpi = 300)


# --- QC 2: TR vs measured VPD, coloured by Day --------------------------------
p_tr_vpd <- ggplot(df, aes(x = VPD, y = TR, colour = Day)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.8) +
  scale_colour_manual(values = c("Day 1" = "#2166ac", "Day 2" = "#d6604d")) +
  facet_wrap(~Genotype, ncol = 5) +
  labs(
    title = "QC2: TR vs VPD by genotype (raw, both days)",
    x     = "VPD (kPa)",
    y     = expression(TR~(mg~m^{-2}~s^{-1}))
  ) +
  theme_bw(base_size = 9) +
  theme(strip.text = element_text(size = 7))

ggsave(file.path(FIG, "QC2_TR_vs_VPD_by_genotype.png"),
       p_tr_vpd, width = 14, height = 8, dpi = 300)


# --- QC 3: RelTR vs VPD, Day 1 only -----------------------------------------
day1 <- filter(df, Day == "Day 1")

p_reltr <- ggplot(day1, aes(x = VPD, y = RelTR, colour = WUE_class)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_point(alpha = 0.7, size = 1.5) +
  facet_wrap(~Genotype, ncol = 5) +
  scale_colour_brewer(palette = "Dark2") +
  labs(
    title = "QC3: Relative transpiration vs VPD — Day 1",
    x     = "VPD (kPa)",
    y     = "Relative transpiration",
    colour = "WUE class"
  ) +
  theme_bw(base_size = 9) +
  theme(strip.text = element_text(size = 7),
        legend.position = "bottom")

ggsave(file.path(FIG, "QC3_RelTR_vs_VPD_Day1.png"),
       p_reltr, width = 14, height = 8, dpi = 300)


# --- QC 4: Stomatal conductance distribution ---------------------------------
gsw_data <- df %>% filter(!is.na(gsw))

p_gsw <- ggplot(gsw_data, aes(x = Genotype, y = gsw, fill = Day)) +
  geom_boxplot(outlier.size = 0.8, position = position_dodge(0.8)) +
  scale_fill_manual(values = c("Day 1" = "#2166ac", "Day 2" = "#d6604d")) +
  labs(
    title = "QC4: Stomatal conductance distribution (gsw, mol m-2 s-1)",
    x     = NULL,
    y     = expression(g[sw]~(mol~m^{-2}~s^{-1}))
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(FIG, "QC4_gsw_distribution.png"),
       p_gsw, width = 12, height = 5, dpi = 300)


# --- QC 5: gsw vs VPD, Day 1 -------------------------------------------------
day1_gsw <- filter(day1, !is.na(gsw))

p_gsw_vpd <- ggplot(day1_gsw, aes(x = VPD, y = gsw, colour = WUE_class)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
  scale_colour_brewer(palette = "Dark2") +
  facet_wrap(~Genotype, ncol = 5) +
  labs(
    title = "QC5: Stomatal conductance vs VPD — Day 1 (LICOR data)",
    x     = "VPD (kPa)",
    y     = expression(g[sw]~(mol~m^{-2}~s^{-1})),
    colour = "WUE class"
  ) +
  theme_bw(base_size = 9) +
  theme(strip.text = element_text(size = 7), legend.position = "bottom")

ggsave(file.path(FIG, "QC5_gsw_vs_VPD_Day1.png"),
       p_gsw_vpd, width = 14, height = 8, dpi = 300)


# --- QC 6: Leaf area by Rep and Genotype -------------------------------------
p_la <- ggplot(df, aes(x = Rep, y = LeafArea, fill = Rep)) +
  geom_boxplot(show.legend = FALSE) +
  facet_wrap(~Genotype, ncol = 5, scales = "free_y") +
  labs(
    title = "QC6: Leaf area (cm2) by rep and genotype",
    x     = "Replicate",
    y     = expression(Leaf~area~(cm^2))
  ) +
  theme_bw(base_size = 9) +
  theme(strip.text = element_text(size = 7))

ggsave(file.path(FIG, "QC6_leaf_area_by_rep_genotype.png"),
       p_la, width = 14, height = 8, dpi = 300)


# --- QC 7: Chamber conditions (Temp, RH) across treatments -------------------
p_cond <- df %>%
  select(Treatment, Temp_C, RH_pct, VPD) %>%
  pivot_longer(-Treatment, names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = Treatment, y = value)) +
  geom_boxplot() +
  facet_wrap(~variable, scales = "free_y") +
  labs(
    title = "QC7: Chamber conditions by treatment",
    x     = "Treatment",
    y     = "Value"
  ) +
  theme_bw()

ggsave(file.path(FIG, "QC7_chamber_conditions.png"),
       p_cond, width = 9, height = 4, dpi = 300)


# --- QC 8: TR Day 1 vs Day 2 within same Rep × Genotype × Treatment ---------
day_compare <- df %>%
  select(Rep, Genotype, Treatment, Day, TR) %>%
  pivot_wider(names_from = Day, values_from = TR,
              values_fn = mean)  # average if >1 obs per cell

p_day_comp <- ggplot(day_compare, aes(x = `Day 1`, y = `Day 2`, colour = Genotype)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(alpha = 0.6, size = 1.5) +
  labs(
    title = "QC8: TR Day 1 vs Day 2 (diagonal = no day effect)",
    x     = expression(TR~Day~1~(mg~m^{-2}~s^{-1})),
    y     = expression(TR~Day~2~(mg~m^{-2}~s^{-1}))
  ) +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 3))

ggsave(file.path(FIG, "QC8_TR_day1_vs_day2.png"),
       p_day_comp, width = 7, height = 6, dpi = 300)


# --------------------------------------------------------------------------- #
# 8. SAVE CLEANED DATA
# --------------------------------------------------------------------------- #

# Remove helper IQR columns; keep outlier flag for transparency
data_clean <- df %>%
  select(-TR_Q1, -TR_Q3, -TR_IQR)

write_csv(data_clean, "data/data_clean.csv")
cat("\nCleaned data saved to data/data_clean.csv\n")
cat("Rows:", nrow(data_clean), "| Columns:", ncol(data_clean), "\n")
cat("\nNote: 'TR_outlier' column flags potential outliers (3×IQR rule).\n")
cat("      'gsw' is the mean of gsw1/gsw2/gsw3; NA where LICOR data absent.\n")
cat("      Day 1 = first measurement day; Day 2 = second measurement day.\n")
