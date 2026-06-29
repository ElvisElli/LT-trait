# =============================================================================
# 01_data_inspection_QC.R
# Growth Chamber 2026 — Soybean Limited-Transpiration Experiment
#
# Purpose : Load data_all.xlsx, inspect structure, flag QC issues,
#           produce 7 diagnostic PNG figures, save data_clean.csv.
#
# Run from project root  (experimental_data_soy_2026/):
#   source("script/01_data_inspection_QC.R")
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(forcats)
  library(stringr)
  library(tibble)
  library(ggbeeswarm)
  library(patchwork)
})

dir.create("figures", showWarnings = FALSE)

# --------------------------------------------------------------------------- #
# COLOUR SYSTEM — University of Arkansas
# --------------------------------------------------------------------------- #
UA_RED   <- "#9D2235"
UA_WINE  <- "#6E1424"
UA_BLUE  <- "#4393c3"
UA_GRAY  <- "#7A7A7A"

WUE_COLS <- c(
  "slow_wilting"        = UA_WINE,
  "high_wilting"        = "#C8102E",
  "Plastic, High WUE"   = UA_BLUE,
  "Nonplastic, Low WUE" = "#d4824a",
  "MAGIC parent"        = "#333333"
)
DAY_COLS <- c("Day 1" = UA_RED, "Day 2" = UA_BLUE)

UA_THEME <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(linewidth = 0.3, colour = "grey90"),
    strip.background  = element_rect(fill = UA_RED, colour = NA),
    strip.text        = element_text(colour = "white", face = "bold", size = 7),
    axis.title        = element_text(size = 9),
    plot.title        = element_text(face = "bold", size = 10, colour = UA_DARK <- "#1A1A1A"),
    legend.key.size   = unit(0.45, "cm"),
    legend.text       = element_text(size = 8)
  )


# --------------------------------------------------------------------------- #
# 1. LOAD DATA
# --------------------------------------------------------------------------- #
cat("Loading data ...\n")
raw <- read_excel("data/data_all.xlsx", sheet = "dataAll")
cat("  Raw dimensions:", nrow(raw), "rows x", ncol(raw), "cols\n\n")


# --------------------------------------------------------------------------- #
# 2. RENAME / RECODE
# --------------------------------------------------------------------------- #
WUE_ORDER <- c("slow_wilting", "high_wilting",
                "Plastic, High WUE", "Nonplastic, Low WUE", "MAGIC parent")

df <- raw |>
  # Drop empty trailing column
  select(where(~ !all(is.na(.)))) |>
  rename(
    Rep        = RepNo,
    Day        = day,
    TR         = TR_mg_sqm_seg,
    RelTR      = RelativeTransp,
    VPD        = avg_vpd,
    Temp_C     = avg_temperature,
    RH_pct     = avg_rh,
    LeafArea   = LeafArea_squarecm,
    duration_min = treatment_time_mins
  ) |>
  mutate(
    Rep       = factor(Rep),
    Day_label = factor(ifelse(Day == 1, "Day 1", "Day 2"),
                       levels = c("Day 1", "Day 2")),
    Genotype  = factor(Genotype),
    WUE_class = factor(WUE_class, levels = WUE_ORDER),
    Treatment = factor(Treatment),
    # Average the three LICOR stomatal conductance readings
    gsw = rowMeans(cbind(gsw1, gsw2, gsw3), na.rm = FALSE)
  )

# Where ALL three gsw readings are NA, set gsw to NA explicitly
all_gsw_na <- is.na(df$gsw1) & is.na(df$gsw2) & is.na(df$gsw3)
df$gsw[all_gsw_na] <- NA

cat("Columns after renaming:\n")
print(names(df))


# --------------------------------------------------------------------------- #
# 3. MISSING VALUE REPORT
# --------------------------------------------------------------------------- #
cat("\n--- Missing values per column ---\n")
miss <- colSums(is.na(df))
print(miss[miss > 0])

cat("\n--- gsw availability by Day ---\n")
df |>
  group_by(Day_label) |>
  summarise(
    n_total = n(),
    n_gsw   = sum(!is.na(gsw)),
    pct_gsw = sprintf("%.1f%%", 100 * n_gsw / n_total),
    .groups = "drop"
  ) |>
  print()


# --------------------------------------------------------------------------- #
# 4. SAMPLE DESIGN SUMMARY
# --------------------------------------------------------------------------- #
cat("\n--- Observations per Genotype x Day ---\n")
df |>
  count(Genotype, Day_label) |>
  pivot_wider(names_from = Day_label, values_from = n, values_fill = 0L) |>
  print(n = 20)

cat("\n--- Genotype WUE class ---\n")
df |>
  distinct(Genotype, WUE_class) |>
  arrange(WUE_class, Genotype) |>
  print(n = 20)

cat("\n--- VPD achieved per Treatment ---\n")
df |>
  group_by(Treatment) |>
  summarise(
    n        = n(),
    VPD_mean = round(mean(VPD, na.rm = TRUE), 3),
    VPD_min  = round(min(VPD,  na.rm = TRUE), 3),
    VPD_max  = round(max(VPD,  na.rm = TRUE), 3),
    .groups = "drop"
  ) |>
  print()


# --------------------------------------------------------------------------- #
# 5. OUTLIER FLAGS
# --------------------------------------------------------------------------- #
# IQR outlier flag (3× rule within Genotype x Day x Treatment)
df <- df |>
  group_by(Genotype, Day_label, Treatment) |>
  mutate(
    TR_Q1      = quantile(TR, 0.25, na.rm = TRUE),
    TR_Q3      = quantile(TR, 0.75, na.rm = TRUE),
    TR_IQR     = TR_Q3 - TR_Q1,
    TR_outlier = TR < (TR_Q1 - 3 * TR_IQR) | TR > (TR_Q3 + 3 * TR_IQR)
  ) |>
  ungroup()

cat("\n--- TR outliers (3×IQR rule) ---\n")
out <- filter(df, TR_outlier)
if (nrow(out) == 0) {
  cat("  None detected\n")
} else {
  out |> select(Rep, Day_label, Genotype, Treatment, VPD, TR) |> print()
}

cat("\n--- Leaf area summary ---\n")
summary(df$LeafArea)
cat("Rows with LeafArea <= 0:", sum(df$LeafArea <= 0, na.rm = TRUE), "\n")


# --------------------------------------------------------------------------- #
# 6. SAVE CLEAN DATA
# --------------------------------------------------------------------------- #
data_clean <- df |> select(-TR_Q1, -TR_Q3, -TR_IQR)
write_csv(data_clean, "data/data_clean.csv")
cat("\nSaved data/data_clean.csv (", nrow(data_clean), "rows)\n")


# --------------------------------------------------------------------------- #
# 7. FIGURES
# --------------------------------------------------------------------------- #
day1     <- filter(df, Day_label == "Day 1")
day1_gsw <- filter(day1, !is.na(gsw))


# --- QC1: TR distribution by Genotype x Day ----------------------------------
p_qc1 <- ggplot(df, aes(x = Genotype, y = TR, fill = Day_label)) +
  geom_boxplot(outlier.size = 0.7, position = position_dodge(0.8),
               linewidth = 0.5) +
  scale_fill_manual(values = DAY_COLS, name = NULL) +
  labs(title = "QC1  Transpiration rate distribution by genotype and day",
       x = NULL, y = expression(TR~(mg~m^{-2}~s^{-1}))) +
  UA_THEME +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

ggsave("figures/QC1_TR_distribution.png", p_qc1,
       width = 13, height = 4.5, dpi = 300)
cat("Saved QC1\n")


# --- QC2: TR vs VPD, both days, per genotype ---------------------------------
p_qc2 <- ggplot(df, aes(x = VPD, y = TR, colour = Day_label)) +
  geom_point(size = 1.2, alpha = 0.65) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, formula = y ~ x) +
  scale_colour_manual(values = DAY_COLS, name = NULL) +
  facet_wrap(~Genotype, ncol = 5) +
  labs(title = "QC2  TR vs VPD by genotype (both days)",
       x = "VPD (kPa)", y = expression(TR~(mg~m^{-2}~s^{-1}))) +
  UA_THEME

ggsave("figures/QC2_TR_vs_VPD_by_genotype.png", p_qc2,
       width = 14, height = 7, dpi = 300)
cat("Saved QC2\n")


# --- QC4: gsw distribution ---------------------------------------------------
gsw_all <- filter(df, !is.na(gsw))

p_qc4 <- ggplot(gsw_all, aes(x = Genotype, y = gsw, fill = Day_label)) +
  geom_boxplot(outlier.size = 0.7, position = position_dodge(0.8),
               linewidth = 0.5) +
  scale_fill_manual(values = DAY_COLS, name = NULL) +
  labs(title = "QC4  Stomatal conductance distribution",
       x = NULL, y = expression(g[sw]~(mol~m^{-2}~s^{-1}))) +
  UA_THEME +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

ggsave("figures/QC4_gsw_distribution.png", p_qc4,
       width = 11, height = 4.5, dpi = 300)
cat("Saved QC4\n")


# --- QC5: gsw vs VPD, Day 1 --------------------------------------------------
p_qc5 <- ggplot(day1_gsw, aes(x = VPD, y = gsw, colour = WUE_class)) +
  geom_point(size = 2.2, alpha = 0.80) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, formula = y ~ x) +
  scale_colour_manual(values = WUE_COLS, name = "WUE class") +
  facet_wrap(~Genotype, ncol = 5) +
  labs(title = "QC5  Stomatal conductance vs VPD — Day 1 (LICOR)",
       x = "VPD (kPa)", y = expression(g[sw]~(mol~m^{-2}~s^{-1}))) +
  UA_THEME +
  theme(legend.position = "bottom")

ggsave("figures/QC5_gsw_vs_VPD_Day1.png", p_qc5,
       width = 14, height = 6, dpi = 300)
cat("Saved QC5\n")


# --- QC6: Leaf area by Rep x Genotype ----------------------------------------
p_qc6 <- ggplot(df, aes(x = Rep, y = LeafArea, colour = Rep)) +
  geom_beeswarm(size = 1.4, alpha = 0.75, cex = 2) +
  scale_colour_manual(values = c("1"=UA_RED,"2"=UA_BLUE,"3"="#2ca25f","4"="#f4a582"),
                      guide = "none") +
  facet_wrap(~Genotype, ncol = 5, scales = "free_y") +
  labs(title = "QC6  Leaf area by replicate and genotype",
       x = "Replicate", y = expression(Leaf~area~(cm^2))) +
  UA_THEME

ggsave("figures/QC6_leaf_area_by_rep_genotype.png", p_qc6,
       width = 14, height = 7, dpi = 300)
cat("Saved QC6\n")


# --- QC7: Chamber conditions per treatment -----------------------------------
p_qc7 <- df |>
  select(Treatment, Temp_C, RH_pct, VPD) |>
  pivot_longer(-Treatment, names_to = "variable", values_to = "value") |>
  mutate(variable = recode(variable,
    Temp_C = "Temperature (°C)",
    RH_pct = "Relative humidity (%)",
    VPD    = "VPD (kPa)"
  )) |>
  ggplot(aes(x = Treatment, y = value)) +
  geom_boxplot(fill = UA_RED, alpha = 0.55, outlier.size = 0.7,
               linewidth = 0.5) +
  facet_wrap(~variable, scales = "free_y", ncol = 3) +
  labs(title = "QC7  Chamber conditions by treatment",
       x = "Treatment", y = NULL) +
  UA_THEME

ggsave("figures/QC7_chamber_conditions.png", p_qc7,
       width = 11, height = 4, dpi = 300)
cat("Saved QC7\n")


# --- QC8: Day 1 vs Day 2 TR --------------------------------------------------
day_pivot <- df |>
  group_by(Rep, Genotype, Treatment, Day_label) |>
  summarise(TR = mean(TR, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = Day_label, values_from = TR)

p_qc8 <- ggplot(day_pivot, aes(x = `Day 1`, y = `Day 2`,
                                colour = Genotype)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = UA_GRAY, linewidth = 0.6) +
  geom_point(size = 1.8, alpha = 0.65) +
  labs(title = "QC8  TR Day 1 vs Day 2 (diagonal = no day effect)",
       x = expression(TR~Day~1~(mg~m^{-2}~s^{-1})),
       y = expression(TR~Day~2~(mg~m^{-2}~s^{-1}))) +
  UA_THEME +
  theme(legend.position = "right") +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 2)))

ggsave("figures/QC8_TR_day1_vs_day2.png", p_qc8,
       width = 7, height = 5.5, dpi = 300)
cat("Saved QC8\n")

cat("\nQC complete — all figures in figures/\n")
