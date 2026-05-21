## GC_comprehensive_analysis.R
## Comprehensive growth chamber transpiration analysis — GC1, GC2 (pending), GC3
## Shows: (1) per-experiment results  (2) overall combined analysis
## Primary response: TranspLeafArea (ml cm⁻²) — actual transpiration per leaf area
## Secondary:        RelTranspiration (0–1 relative) — for cross-experiment comparison
## Day 2 recommended for modelling (Day 1 acclimation artifact documented)

## ── Path configuration ────────────────────────────────────────────────────────
# project_root <- "C:/Users/YourName/Projects/LT-trait"
if (.Platform$OS.type == "windows") setwd(project_root) else setwd("/home/user/LT-trait")
## ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(readxl); library(scales); library(ggrepel)
})

figdir  <- "experimental_data/figures"
statsdir <- "experimental_data/stats"
dir.create(figdir,   showWarnings = FALSE, recursive = TRUE)
dir.create(statsdir, showWarnings = FALSE, recursive = TRUE)

## ── Shared style ──────────────────────────────────────────────────────────────
col_quad  <- c(Q1 = "#d95f02", Q2 = "#1b9e77", Q3 = "#7570b3", Q4 = "#e7298a")
col_exp   <- c(GC1 = "#1f78b4", GC3 = "#e31a1c")
col_day   <- c(`1` = "grey60",  `2` = "#2171b5")
VPD_BREAK <- 2.0   # Sinclair 2008 reference breakpoint

base_theme <- theme_bw(base_size = 12) +
  theme(
    strip.background  = element_rect(fill = "grey15", colour = "grey15"),
    strip.text        = element_text(colour = "white", face = "bold"),
    panel.grid.minor  = element_blank(),
    legend.position   = "bottom"
  )

save_fig <- function(p, name, w = 22, h = 16) {
  ggsave(file.path(figdir, paste0(name, ".png")),
         p, width = w, height = h, units = "cm", dpi = 150, bg = "white")
  message("Saved: ", name, ".png")
}

## ── Quadrant label helper ─────────────────────────────────────────────────────
quad_label <- function(q) {
  c(Q1 = "Q1 · Plastic / High WUE",
    Q2 = "Q2 · Plastic / Low WUE",
    Q3 = "Q3 · Nonplastic / Low WUE",
    Q4 = "Q4 · Nonplastic / High WUE")[q]
}

## ══════════════════════════════════════════════════════════════════════════════
## 1. LOAD DATA
## ══════════════════════════════════════════════════════════════════════════════
cat("\n── Loading cleaned data ──────────────────────────────────────────────────\n")
df <- read.csv("experimental_data/cleaned_data/GC_combined_clean.csv",
               stringsAsFactors = FALSE) %>%
  mutate(
    Day       = as.character(Day),
    Quadrant  = factor(Quadrant, levels = c("Q1","Q2","Q3","Q4")),
    Experiment = factor(Experiment, levels = c("GC1","GC3")),
    Geno_label = paste0(Genotype, "\n(", quad_label(Quadrant), ")")
  )

cat("Total rows:", nrow(df),
    "| GC1:", sum(df$Experiment == "GC1"),
    "| GC3:", sum(df$Experiment == "GC3"), "\n")

## ── Outlier detection (>2 SD from group mean per Experiment×Genotype×Day×VPD) ──
df <- df %>%
  group_by(Experiment, Genotype, Day, Treatment) %>%
  mutate(
    mean_TLA  = mean(TranspLeafArea, na.rm = TRUE),
    sd_TLA    = sd(TranspLeafArea,   na.rm = TRUE),
    mean_Rel  = mean(RelTranspiration, na.rm = TRUE),
    sd_Rel    = sd(RelTranspiration,   na.rm = TRUE),
    outlier   = !is.na(TranspLeafArea) &
                abs(TranspLeafArea - mean_TLA) > 2 * sd_TLA
  ) %>%
  ungroup()

cat("Outliers flagged:", sum(df$outlier, na.rm = TRUE), "\n\n")

## Day 2 subset (primary analysis)
d2 <- df %>% filter(Day == "2")

## ══════════════════════════════════════════════════════════════════════════════
## 2. GC1 ANALYSIS
## ══════════════════════════════════════════════════════════════════════════════
cat("── GC1 Analysis ─────────────────────────────────────────────────────────\n")
gc1    <- df %>% filter(Experiment == "GC1")
gc1_d2 <- gc1  %>% filter(Day == "2", !outlier)

## GC1 means (Day 2, no outliers)
gc1_means <- gc1_d2 %>%
  group_by(Genotype, Quadrant, VPD) %>%
  summarise(
    mean_TLA  = mean(TranspLeafArea,   na.rm = TRUE),
    se_TLA    = sd(TranspLeafArea,     na.rm = TRUE) / sqrt(n()),
    mean_Rel  = mean(RelTranspiration, na.rm = TRUE),
    se_Rel    = sd(RelTranspiration,   na.rm = TRUE) / sqrt(n()),
    mean_SC   = mean(StomatalConductance, na.rm = TRUE),
    se_SC     = sd(StomatalConductance,   na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

## — Fig GC1-1: Actual transpiration per leaf area vs VPD —
p_gc1_1 <- ggplot(gc1_d2, aes(x = VPD, colour = Quadrant, fill = Quadrant)) +
  geom_vline(xintercept = VPD_BREAK, linetype = "dashed",
             colour = "grey50", linewidth = 0.7) +
  annotate("text", x = VPD_BREAK + 0.05, y = Inf, vjust = 1.5, hjust = 0,
           label = "VPD threshold\n(Sinclair 2008)", size = 2.8, colour = "grey40") +
  geom_jitter(aes(y = TranspLeafArea), width = 0.04, alpha = 0.35, size = 1.5,
              shape = 16, show.legend = FALSE) +
  geom_line(data = gc1_means, aes(y = mean_TLA), linewidth = 1.0) +
  geom_errorbar(data = gc1_means,
                aes(y = mean_TLA, ymin = mean_TLA - se_TLA, ymax = mean_TLA + se_TLA),
                width = 0.06, linewidth = 0.7) +
  geom_point(data = gc1_means, aes(y = mean_TLA),
             size = 3.5, shape = 21, colour = "white", stroke = 0.8) +
  facet_wrap(~ Quadrant + Genotype, nrow = 1,
             labeller = labeller(Quadrant = c(
               Q1 = "Q1 · Plastic / High WUE",
               Q2 = "Q2 · Plastic / Low WUE",
               Q3 = "Q3 · Nonplastic / Low WUE",
               Q4 = "Q4 · Nonplastic / High WUE"))) +
  scale_colour_manual(values = col_quad) +
  scale_fill_manual(values   = col_quad) +
  labs(
    title    = "GC1 — Actual transpiration per leaf area vs VPD (Day 2)",
    subtitle = "Points = individual pots · Line = mean ± SE · Dashed = Sinclair 2008 VPD threshold",
    x        = "VPD (kPa)",
    y        = expression(Transpiration~(ml~cm^{-2})),
    colour = "Quadrant", fill = "Quadrant"
  ) +
  base_theme + theme(legend.position = "none")

save_fig(p_gc1_1, "GC1_Fig1_actual_transpiration", w = 28, h = 14)

## — Fig GC1-2: Relative transpiration vs VPD —
p_gc1_2 <- ggplot(gc1_d2, aes(x = VPD, colour = Quadrant, fill = Quadrant)) +
  geom_vline(xintercept = VPD_BREAK, linetype = "dashed",
             colour = "grey50", linewidth = 0.7) +
  geom_hline(yintercept = 1.0, linetype = "dotted", colour = "grey70") +
  geom_jitter(aes(y = RelTranspiration), width = 0.04, alpha = 0.35, size = 1.5,
              shape = 16, show.legend = FALSE) +
  geom_line(data = gc1_means, aes(y = mean_Rel), linewidth = 1.0) +
  geom_errorbar(data = gc1_means,
                aes(y = mean_Rel, ymin = mean_Rel - se_Rel, ymax = mean_Rel + se_Rel),
                width = 0.06, linewidth = 0.7) +
  geom_point(data = gc1_means, aes(y = mean_Rel),
             size = 3.5, shape = 21, colour = "white", stroke = 0.8) +
  facet_wrap(~ Quadrant + Genotype, nrow = 1,
             labeller = labeller(Quadrant = c(
               Q1 = "Q1 · Plastic / High WUE",
               Q2 = "Q2 · Plastic / Low WUE",
               Q3 = "Q3 · Nonplastic / Low WUE",
               Q4 = "Q4 · Nonplastic / High WUE"))) +
  scale_colour_manual(values = col_quad) +
  scale_fill_manual(values   = col_quad) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(
    title    = "GC1 — Relative transpiration vs VPD (Day 2)",
    subtitle = "Relative to pot maximum · Dashed = Sinclair 2008 VPD threshold",
    x = "VPD (kPa)", y = "Relative transpiration (0–1)"
  ) +
  base_theme + theme(legend.position = "none")

save_fig(p_gc1_2, "GC1_Fig2_relative_transpiration", w = 28, h = 14)

## — Fig GC1-3: Stomatal conductance vs VPD —
gc1_sc <- gc1_d2 %>% filter(!is.na(StomatalConductance), StomatalConductance > 0)

p_gc1_3 <- ggplot(gc1_sc, aes(x = VPD, colour = Quadrant, fill = Quadrant)) +
  geom_vline(xintercept = VPD_BREAK, linetype = "dashed",
             colour = "grey50", linewidth = 0.7) +
  geom_jitter(aes(y = StomatalConductance), width = 0.04, alpha = 0.4,
              size = 1.8, shape = 16, show.legend = FALSE) +
  geom_smooth(aes(y = StomatalConductance), method = "lm",
              se = TRUE, linewidth = 1.0, alpha = 0.15) +
  facet_wrap(~ Quadrant + Genotype, nrow = 1,
             labeller = labeller(Quadrant = c(
               Q1 = "Q1 · Plastic / High WUE",
               Q2 = "Q2 · Plastic / Low WUE",
               Q3 = "Q3 · Nonplastic / Low WUE",
               Q4 = "Q4 · Nonplastic / High WUE"))) +
  scale_colour_manual(values = col_quad) +
  scale_fill_manual(values   = col_quad) +
  labs(
    title    = "GC1 — Stomatal conductance vs VPD (Day 2, LI-600 data)",
    subtitle = "GC1 only — GC3 has no stomatal conductance measurements",
    x = "VPD (kPa)",
    y = expression(Stomatal~conductance~(mol~m^{-2}~s^{-1}))
  ) +
  base_theme + theme(legend.position = "none")

save_fig(p_gc1_3, "GC1_Fig3_stomatal_conductance", w = 28, h = 14)

## — Fig GC1-4: Day 1 vs Day 2 diagnostic —
gc1_day <- gc1 %>%
  group_by(Genotype, Quadrant, Day, VPD) %>%
  summarise(mean_TLA = mean(TranspLeafArea, na.rm = TRUE),
            se_TLA   = sd(TranspLeafArea,   na.rm = TRUE) / sqrt(n()),
            .groups  = "drop")

p_gc1_4 <- ggplot(gc1, aes(x = VPD, y = TranspLeafArea,
                            colour = Day, shape = Day)) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_line(data = gc1_day, aes(y = mean_TLA, group = Day), linewidth = 0.9) +
  geom_point(data = gc1_day, aes(y = mean_TLA), size = 3) +
  geom_point(data = filter(gc1, outlier),
             aes(y = TranspLeafArea), shape = 4, size = 3,
             colour = "red", stroke = 1.2, show.legend = FALSE) +
  facet_wrap(~ Quadrant + Genotype, nrow = 1) +
  scale_colour_manual(values = col_day, name = "Day",
                      labels = c("1" = "Day 1 (acclimation)",
                                 "2" = "Day 2 (stable)")) +
  scale_shape_manual(values = c("1" = 1, "2" = 16), name = "Day",
                     labels = c("1" = "Day 1 (acclimation)",
                                "2" = "Day 2 (stable)")) +
  labs(
    title    = "GC1 — Day 1 vs Day 2 comparison (outliers flagged in red ×)",
    subtitle = "Day 1 shows acclimation artifact — Day 2 used for all modelling",
    x = "VPD (kPa)", y = expression(Transpiration~(ml~cm^{-2}))
  ) +
  base_theme

save_fig(p_gc1_4, "GC1_Fig4_day_comparison_outliers", w = 28, h = 14)

## ── GC1 outlier table ─────────────────────────────────────────────────────────
gc1_out <- gc1 %>% filter(outlier) %>%
  select(Experiment, Genotype, Quadrant, Day, Treatment, VPD,
         Pot, TranspLeafArea, mean_TLA, sd_TLA) %>%
  mutate(z_score = round((TranspLeafArea - mean_TLA) / sd_TLA, 2))
write.csv(gc1_out,
          file.path(statsdir, "GC1_outliers.csv"), row.names = FALSE)
cat("GC1 outliers:", nrow(gc1_out), "\n")

## ══════════════════════════════════════════════════════════════════════════════
## 3. GC3 ANALYSIS
## ══════════════════════════════════════════════════════════════════════════════
cat("\n── GC3 Analysis ─────────────────────────────────────────────────────────\n")
gc3    <- df %>% filter(Experiment == "GC3")
gc3_d2 <- gc3  %>% filter(Day == "2", !outlier)

gc3_means <- gc3_d2 %>%
  group_by(Genotype, Quadrant, VPD) %>%
  summarise(
    mean_TLA = mean(TranspLeafArea,   na.rm = TRUE),
    se_TLA   = sd(TranspLeafArea,     na.rm = TRUE) / sqrt(n()),
    mean_Rel = mean(RelTranspiration, na.rm = TRUE),
    se_Rel   = sd(RelTranspiration,   na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

## — Fig GC3-1: Actual transpiration per leaf area vs VPD —
p_gc3_1 <- ggplot(gc3_d2, aes(x = VPD, colour = Quadrant, fill = Quadrant)) +
  geom_vline(xintercept = VPD_BREAK, linetype = "dashed",
             colour = "grey50", linewidth = 0.7) +
  geom_jitter(aes(y = TranspLeafArea), width = 0.04, alpha = 0.35, size = 1.5,
              shape = 16, show.legend = FALSE) +
  geom_line(data = gc3_means, aes(y = mean_TLA), linewidth = 1.0) +
  geom_errorbar(data = gc3_means,
                aes(y = mean_TLA, ymin = mean_TLA - se_TLA, ymax = mean_TLA + se_TLA),
                width = 0.06, linewidth = 0.7) +
  geom_point(data = gc3_means, aes(y = mean_TLA),
             size = 3.5, shape = 21, colour = "white", stroke = 0.8) +
  facet_wrap(~ Quadrant + Genotype, nrow = 1,
             labeller = labeller(Quadrant = c(
               Q1 = "Q1 · Plastic / High WUE",
               Q2 = "Q2 · Plastic / Low WUE",
               Q3 = "Q3 · Nonplastic / Low WUE",
               Q4 = "Q4 · Nonplastic / High WUE"))) +
  scale_colour_manual(values = col_quad) +
  scale_fill_manual(values   = col_quad) +
  labs(
    title    = "GC3 — Actual transpiration per leaf area vs VPD (Day 2)",
    subtitle = "Points = individual pots · Line = mean ± SE · Dashed = Sinclair 2008 VPD threshold",
    x = "VPD (kPa)", y = expression(Transpiration~(ml~cm^{-2}))
  ) +
  base_theme + theme(legend.position = "none")

save_fig(p_gc3_1, "GC3_Fig1_actual_transpiration", w = 28, h = 14)

## — Fig GC3-2: Relative transpiration vs VPD —
p_gc3_2 <- ggplot(gc3_d2, aes(x = VPD, colour = Quadrant, fill = Quadrant)) +
  geom_vline(xintercept = VPD_BREAK, linetype = "dashed",
             colour = "grey50", linewidth = 0.7) +
  geom_hline(yintercept = 1.0, linetype = "dotted", colour = "grey70") +
  geom_jitter(aes(y = RelTranspiration), width = 0.04, alpha = 0.35, size = 1.5,
              shape = 16, show.legend = FALSE) +
  geom_line(data = gc3_means, aes(y = mean_Rel), linewidth = 1.0) +
  geom_errorbar(data = gc3_means,
                aes(y = mean_Rel, ymin = mean_Rel - se_Rel, ymax = mean_Rel + se_Rel),
                width = 0.06, linewidth = 0.7) +
  geom_point(data = gc3_means, aes(y = mean_Rel),
             size = 3.5, shape = 21, colour = "white", stroke = 0.8) +
  facet_wrap(~ Quadrant + Genotype, nrow = 1,
             labeller = labeller(Quadrant = c(
               Q1 = "Q1 · Plastic / High WUE",
               Q2 = "Q2 · Plastic / Low WUE",
               Q3 = "Q3 · Nonplastic / Low WUE",
               Q4 = "Q4 · Nonplastic / High WUE"))) +
  scale_colour_manual(values = col_quad) +
  scale_fill_manual(values   = col_quad) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(
    title    = "GC3 — Relative transpiration vs VPD (Day 2)",
    x = "VPD (kPa)", y = "Relative transpiration (0–1)"
  ) +
  base_theme + theme(legend.position = "none")

save_fig(p_gc3_2, "GC3_Fig2_relative_transpiration", w = 28, h = 14)

## — Fig GC3-3: Day 1 vs Day 2 diagnostic —
gc3_day <- gc3 %>%
  group_by(Genotype, Quadrant, Day, VPD) %>%
  summarise(mean_TLA = mean(TranspLeafArea, na.rm = TRUE),
            .groups  = "drop")

p_gc3_3 <- ggplot(gc3, aes(x = VPD, y = TranspLeafArea,
                            colour = Day, shape = Day)) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_line(data = gc3_day, aes(y = mean_TLA, group = Day), linewidth = 0.9) +
  geom_point(data = gc3_day, aes(y = mean_TLA), size = 3) +
  geom_point(data = filter(gc3, outlier),
             aes(y = TranspLeafArea), shape = 4, size = 3,
             colour = "red", stroke = 1.2, show.legend = FALSE) +
  facet_wrap(~ Quadrant + Genotype, nrow = 1) +
  scale_colour_manual(values = col_day, name = "Day",
                      labels = c("1" = "Day 1 (acclimation)",
                                 "2" = "Day 2 (stable)")) +
  scale_shape_manual(values = c("1" = 1, "2" = 16), name = "Day",
                     labels = c("1" = "Day 1 (acclimation)",
                                "2" = "Day 2 (stable)")) +
  labs(
    title    = "GC3 — Day 1 vs Day 2 comparison (outliers flagged in red ×)",
    x = "VPD (kPa)", y = expression(Transpiration~(ml~cm^{-2}))
  ) +
  base_theme

save_fig(p_gc3_3, "GC3_Fig3_day_comparison_outliers", w = 28, h = 14)

gc3_out <- gc3 %>% filter(outlier) %>%
  select(Experiment, Genotype, Quadrant, Day, Treatment, VPD,
         Pot, TranspLeafArea, mean_TLA, sd_TLA) %>%
  mutate(z_score = round((TranspLeafArea - mean_TLA) / sd_TLA, 2))
write.csv(gc3_out,
          file.path(statsdir, "GC3_outliers.csv"), row.names = FALSE)
cat("GC3 outliers:", nrow(gc3_out), "\n")

## ══════════════════════════════════════════════════════════════════════════════
## 4. GC2 STATUS NOTE
## ══════════════════════════════════════════════════════════════════════════════
cat("\n── GC2 Status ───────────────────────────────────────────────────────────\n")
cat("GC2 raw data: experimental_data/raw_data/GC2/GC_experiment_2.xlsx\n")
cat("Sheets: Experiment map, Weights, Irrigation, Evaluation targets,\n")
cat("        Tracking VPD, Tracking weights\n")
cat("Status: No 'Final' sheet yet. Raw pot weights recorded; transpiration\n")
cat("        calculation (weight difference × correction) pending.\n")
cat("Genotypes: PI398245, PI567488B, PI377574, PI424397\n")
cat("Action: Process Weights sheet into Final format (like GC1/GC3)\n")

## ══════════════════════════════════════════════════════════════════════════════
## 5. OVERALL COMBINED ANALYSIS
## ══════════════════════════════════════════════════════════════════════════════
cat("\n── Overall Combined Analysis ────────────────────────────────────────────\n")

## — Fig Overall-1: All genotypes, both experiments, actual transpiration —
all_means <- d2 %>% filter(!outlier) %>%
  group_by(Experiment, Genotype, Quadrant, VPD) %>%
  summarise(
    mean_TLA = mean(TranspLeafArea,   na.rm = TRUE),
    se_TLA   = sd(TranspLeafArea,     na.rm = TRUE) / sqrt(n()),
    mean_Rel = mean(RelTranspiration, na.rm = TRUE),
    se_Rel   = sd(RelTranspiration,   na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    Geno_exp = paste0(Genotype, "\n[", Experiment, "]")
  )

## Re-order genotypes by quadrant
geno_order <- d2 %>%
  distinct(Genotype, Quadrant, Experiment) %>%
  arrange(Quadrant, Experiment) %>%
  mutate(Geno_exp = paste0(Genotype, "\n[", Experiment, "]"))

all_means$Geno_exp <- factor(all_means$Geno_exp,
                              levels = geno_order$Geno_exp)

p_all_1 <- ggplot(d2 %>% filter(!outlier) %>%
                    mutate(Geno_exp = factor(
                      paste0(Genotype, "\n[", Experiment, "]"),
                      levels = geno_order$Geno_exp)),
                  aes(x = VPD, colour = Quadrant, fill = Quadrant)) +
  geom_vline(xintercept = VPD_BREAK, linetype = "dashed",
             colour = "grey50", linewidth = 0.6) +
  geom_jitter(aes(y = TranspLeafArea), width = 0.04, alpha = 0.25, size = 1.2,
              shape = 16, show.legend = FALSE) +
  geom_line(data = all_means, aes(y = mean_TLA), linewidth = 1.0) +
  geom_errorbar(data = all_means,
                aes(y = mean_TLA, ymin = mean_TLA - se_TLA,
                    ymax = mean_TLA + se_TLA),
                width = 0.06, linewidth = 0.6) +
  geom_point(data = all_means, aes(y = mean_TLA),
             size = 3.0, shape = 21, colour = "white", stroke = 0.8) +
  facet_wrap(~ Geno_exp, nrow = 2) +
  scale_colour_manual(values = col_quad, name = "Quadrant") +
  scale_fill_manual(values   = col_quad, name = "Quadrant") +
  labs(
    title    = "Overall — Actual transpiration per leaf area vs VPD (GC1 + GC3, Day 2)",
    subtitle = "Genotypes ordered by quadrant (Q1–Q4) · Each column = same quadrant across experiments",
    x = "VPD (kPa)", y = expression(Transpiration~(ml~cm^{-2}))
  ) +
  base_theme +
  guides(colour = guide_legend(nrow = 1),
         fill   = guide_legend(nrow = 1))

save_fig(p_all_1, "Overall_Fig1_actual_transpiration", w = 36, h = 22)

## — Fig Overall-2: Relative transpiration, all genotypes —
p_all_2 <- ggplot(d2 %>% filter(!outlier) %>%
                    mutate(Geno_exp = factor(
                      paste0(Genotype, "\n[", Experiment, "]"),
                      levels = geno_order$Geno_exp)),
                  aes(x = VPD, colour = Quadrant, fill = Quadrant)) +
  geom_vline(xintercept = VPD_BREAK, linetype = "dashed",
             colour = "grey50", linewidth = 0.6) +
  geom_hline(yintercept = 1.0, linetype = "dotted", colour = "grey70") +
  geom_jitter(aes(y = RelTranspiration), width = 0.04, alpha = 0.25, size = 1.2,
              shape = 16, show.legend = FALSE) +
  geom_line(data = all_means, aes(y = mean_Rel), linewidth = 1.0) +
  geom_errorbar(data = all_means,
                aes(y = mean_Rel, ymin = mean_Rel - se_Rel,
                    ymax = mean_Rel + se_Rel),
                width = 0.06, linewidth = 0.6) +
  geom_point(data = all_means, aes(y = mean_Rel),
             size = 3.0, shape = 21, colour = "white", stroke = 0.8) +
  facet_wrap(~ Geno_exp, nrow = 2) +
  scale_colour_manual(values = col_quad, name = "Quadrant") +
  scale_fill_manual(values   = col_quad, name = "Quadrant") +
  scale_y_continuous(limits = c(0, NA)) +
  labs(
    title    = "Overall — Relative transpiration vs VPD (GC1 + GC3, Day 2)",
    subtitle = "Genotypes ordered by quadrant (Q1–Q4)",
    x = "VPD (kPa)", y = "Relative transpiration (0–1)"
  ) +
  base_theme +
  guides(colour = guide_legend(nrow = 1),
         fill   = guide_legend(nrow = 1))

save_fig(p_all_2, "Overall_Fig2_relative_transpiration", w = 36, h = 22)

## — Fig Overall-3: Quadrant-pooled (GC1+GC3 combined per quadrant) —
quad_means <- d2 %>% filter(!outlier) %>%
  group_by(Quadrant, VPD) %>%
  summarise(
    mean_TLA = mean(TranspLeafArea,   na.rm = TRUE),
    se_TLA   = sd(TranspLeafArea,     na.rm = TRUE) / sqrt(n()),
    mean_Rel = mean(RelTranspiration, na.rm = TRUE),
    se_Rel   = sd(RelTranspiration,   na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

p_all_3a <- ggplot(d2 %>% filter(!outlier),
                   aes(x = VPD, colour = Quadrant, fill = Quadrant)) +
  geom_vline(xintercept = VPD_BREAK, linetype = "dashed",
             colour = "grey50", linewidth = 0.7) +
  geom_jitter(aes(y = TranspLeafArea), width = 0.04, alpha = 0.25, size = 1.5,
              shape = 16, show.legend = FALSE) +
  geom_line(data = quad_means, aes(y = mean_TLA), linewidth = 1.2) +
  geom_errorbar(data = quad_means,
                aes(y = mean_TLA, ymin = mean_TLA - se_TLA,
                    ymax = mean_TLA + se_TLA),
                width = 0.07, linewidth = 0.7) +
  geom_point(data = quad_means, aes(y = mean_TLA),
             size = 4, shape = 21, colour = "white", stroke = 1.0) +
  facet_wrap(~ Quadrant, nrow = 1,
             labeller = labeller(Quadrant = c(
               Q1 = "Q1\nPlastic / High WUE",
               Q2 = "Q2\nPlastic / Low WUE",
               Q3 = "Q3\nNonplastic / Low WUE",
               Q4 = "Q4\nNonplastic / High WUE"))) +
  scale_colour_manual(values = col_quad) +
  scale_fill_manual(values   = col_quad) +
  labs(
    title    = "Overall — Quadrant-pooled actual transpiration (GC1 + GC3, Day 2)",
    subtitle = "Each quadrant pools 2 genotypes across 2 experiments (n = 8 per treatment per quadrant)",
    x = "VPD (kPa)", y = expression(Transpiration~(ml~cm^{-2}))
  ) +
  base_theme + theme(legend.position = "none")

save_fig(p_all_3a, "Overall_Fig3_quadrant_actual_transp", w = 28, h = 14)

## — Fig Overall-4: Experiment comparison within each quadrant —
p_all_4 <- ggplot(all_means,
                  aes(x = VPD, y = mean_TLA,
                      colour = Experiment, linetype = Experiment)) +
  geom_vline(xintercept = VPD_BREAK, linetype = "dashed",
             colour = "grey50", linewidth = 0.6) +
  geom_line(linewidth = 1.0) +
  geom_errorbar(aes(ymin = mean_TLA - se_TLA, ymax = mean_TLA + se_TLA),
                width = 0.07, linewidth = 0.6) +
  geom_point(aes(shape = Experiment), size = 3.5) +
  facet_wrap(~ Quadrant, nrow = 1,
             labeller = labeller(Quadrant = c(
               Q1 = "Q1\nPlastic / High WUE",
               Q2 = "Q2\nPlastic / Low WUE",
               Q3 = "Q3\nNonplastic / Low WUE",
               Q4 = "Q4\nNonplastic / High WUE"))) +
  scale_colour_manual(values = col_exp) +
  scale_shape_manual(values  = c(GC1 = 16, GC3 = 17)) +
  labs(
    title    = "Experiment comparison — GC1 vs GC3 per quadrant (Day 2)",
    subtitle = "Consistency across experiments validates quadrant classifications",
    x = "VPD (kPa)", y = expression(Transpiration~(ml~cm^{-2})),
    colour = "Experiment", linetype = "Experiment", shape = "Experiment"
  ) +
  base_theme

save_fig(p_all_4, "Overall_Fig4_GC1_vs_GC3_comparison", w = 28, h = 14)

## ── Summary statistics table ──────────────────────────────────────────────────
summary_stats <- d2 %>% filter(!outlier) %>%
  group_by(Experiment, Genotype, Quadrant, Expected_WUE) %>%
  summarise(
    n              = n(),
    mean_TLA       = round(mean(TranspLeafArea,   na.rm = TRUE), 3),
    sd_TLA         = round(sd(TranspLeafArea,     na.rm = TRUE), 3),
    mean_Rel       = round(mean(RelTranspiration, na.rm = TRUE), 3),
    sd_Rel         = round(sd(RelTranspiration,   na.rm = TRUE), 3),
    VPD_range      = paste0(round(min(VPD, na.rm=TRUE),1), "–",
                             round(max(VPD, na.rm=TRUE),1)),
    n_outliers_all = sum(df$outlier[df$Genotype == unique(Genotype)], na.rm=TRUE),
    .groups = "drop"
  )

write.csv(summary_stats,
          file.path(statsdir, "GC_summary_by_genotype.csv"), row.names = FALSE)

cat("\nSummary statistics saved to stats/GC_summary_by_genotype.csv\n")
cat("\n── Done ─────────────────────────────────────────────────────────────────\n")
cat("Figures saved to:", figdir, "\n")
cat("  GC1_Fig1_actual_transpiration.png\n")
cat("  GC1_Fig2_relative_transpiration.png\n")
cat("  GC1_Fig3_stomatal_conductance.png\n")
cat("  GC1_Fig4_day_comparison_outliers.png\n")
cat("  GC3_Fig1_actual_transpiration.png\n")
cat("  GC3_Fig2_relative_transpiration.png\n")
cat("  GC3_Fig3_day_comparison_outliers.png\n")
cat("  Overall_Fig1_actual_transpiration.png\n")
cat("  Overall_Fig2_relative_transpiration.png\n")
cat("  Overall_Fig3_quadrant_actual_transp.png\n")
cat("  Overall_Fig4_GC1_vs_GC3_comparison.png\n")
