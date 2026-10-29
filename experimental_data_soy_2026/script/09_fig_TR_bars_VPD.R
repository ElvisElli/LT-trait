## 09_fig_TR_bars_VPD.R
## Bar plot of model-predicted transpiration at two VPD levels (VPD_TARGETS,
## currently 1.5 and 3.5 kPa) per accession,
## to show the variability among genotypes.
##   * same accessions and data filters as the presentation figure
##     (08_fig_TR_presentation.R): 9 accessions, rep 5 of PI424231 dropped
##   * per accession, the AIC-selected mixed model (script/model_functions.R)
##     is evaluated at the two target VPDs
##   * error bars = 95% confidence interval of the predicted mean (fixed
##     effects; delta method on the model's fixed-effect covariance matrix)
##   * accessions are ordered by transpiration at the higher VPD
##   * the figure is drawn in several ColorBrewer palettes (PALETTES below);
##     the first one is saved as the main figure
##
## Input : data/data_clean.csv
## Output: figures/presentation/Fig_TR_bars_<lo>_<hi>kPa.png / .tiff (600 dpi)
##         figures/presentation/Fig_TR_bars_<lo>_<hi>kPa_<palette>.png (palette options)
##         figures/presentation/TR_predicted_<lo>_<hi>kPa.csv
## Run from experimental_data_soy_2026/:  Rscript script/09_fig_TR_bars_VPD.R
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(readr) })
source("script/model_functions.R")

DROP_GENO <- c("PI398272", "PI404191", "PI442012A", "AG43XF2", "PI377574")
DROP_REPS <- list(PI424231 = "5")
VPD_TARGETS <- c(1.5, 3.5)
VPD_LABELS  <- paste(VPD_TARGETS, "kPa")
TAG <- paste0(VPD_TARGETS[1], "_", VPD_TARGETS[2], "kPa")
## ColorBrewer pairs (low VPD, high VPD)
PALETTES <- list(
  BrBG    = c("#dfc27d", "#35978f"),   # diverging brown / blue-green
  PuOr    = c("#fdb863", "#8073ac"),   # diverging orange / purple
  RdBu    = c("#f4a582", "#4393c3"),   # diverging red / blue
  Greens  = c("#a1d99b", "#31a354"),   # sequential (previous version)
  BluesTeal = c("#c7eae5", "#01665e")  # light / dark teal from BrBG
)
GROUP <- c("slow_wilting" = "Slow wilting", "high_wilting" = "High wilting",
           "Plastic, High WUE" = "High WUE plasticity", "Nonplastic, Low WUE" = "Low WUE plasticity",
           "MAGIC parent" = "Low WUE plasticity")

d <- read_csv("data/data_clean.csv", show_col_types = FALSE) |>
  mutate(Rep = factor(Rep)) |>
  filter(Day_label == "Day 1", !is.na(TR), !is.na(VPD), !(Genotype %in% DROP_GENO))
for (g in names(DROP_REPS)) d <- d |> filter(!(Genotype == g & Rep %in% DROP_REPS[[g]]))
genos <- sort(unique(d$Genotype))

## prediction and its standard error (delta method) for one fitted model
predict_se <- function(m, model, x) {
  b <- if (inherits(m, c("lme", "nlme"))) fixef(m) else coef(m)
  V <- vcov(m)
  rows <- lapply(x, function(v) {
    if (model == "Linear") {
      est <- b[1] + b[2] * v;                 g <- c(1, v)
    } else if (model == "Lin. plateau") {
      a <- b["a"]; s <- b["b"]; bp <- b["bp"]
      if (v <= bp) { est <- a + s * v;  g <- c(1, v, 0) }
      else         { est <- a + s * bp; g <- c(1, bp, s) }
    } else {                                   # Bilinear
      a <- b["a"]; b1 <- b["b1"]; b2 <- b["b2"]; bp <- b["bp"]
      if (v <= bp) { est <- a + b1 * v;                   g <- c(1, v, 0, 0) }
      else         { est <- a + b1 * bp + b2 * (v - bp);  g <- c(1, bp, v - bp, b1 - b2) }
    }
    c(est = unname(est), se = sqrt(as.numeric(t(g) %*% V %*% g)))
  })
  as.data.frame(do.call(rbind, rows))
}

pred <- bind_rows(lapply(genos, function(g) {
  f  <- fit_all_models(d |> filter(Genotype == g) |> select(TR, VPD, Rep), response = "TR")
  ps <- predict_se(f[[f$best]], f$best, VPD_TARGETS)
  tibble(Genotype = g, Best_model = f$best, Fit = class(f[[f$best]])[1], VPD = VPD_TARGETS,
         TR = pmax(ps$est, 0), SE = ps$se,
         lo = pmax(ps$est - 1.96 * ps$se, 0), hi = ps$est + 1.96 * ps$se)
})) |>
  left_join(d |> distinct(Genotype, WUE_class), by = "Genotype") |>
  mutate(Group = unname(GROUP[WUE_class]), VPD_label = factor(paste(VPD, "kPa"), levels = VPD_LABELS))

wide <- pred |> select(Genotype, Group, Best_model, Fit, VPD, TR, lo, hi) |>
  mutate(VPD = ifelse(VPD == VPD_TARGETS[1], "low", "high")) |>
  pivot_wider(names_from = VPD, values_from = c(TR, lo, hi)) |>
  mutate(increase_pct = 100 * (TR_high - TR_low) / TR_low) |> arrange(TR_high)
write_csv(wide |> mutate(across(where(is.numeric), ~ round(., 2))), sprintf("figures/presentation/TR_predicted_%s.csv", TAG))
print(wide |> mutate(across(where(is.numeric), ~ round(., 1))), n = Inf)
cat(sprintf("\nRange across accessions: %.1f-%.1f at %s, %.1f-%.1f at %s (mg m-2 s-1)\n",
            min(wide$TR_low), max(wide$TR_low), VPD_LABELS[1], min(wide$TR_high), max(wide$TR_high), VPD_LABELS[2]))

## order by transpiration at the higher VPD; x label = accession name only
pred <- pred |> mutate(Genotype = factor(Genotype, levels = wide$Genotype))

make_plot <- function(cols) {
  ggplot(pred, aes(Genotype, TR, fill = VPD_label)) +
    geom_col(position = position_dodge(width = 0.78), width = 0.72, colour = "black", linewidth = 0.3) +
    geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.78),
                  width = 0.22, linewidth = 0.5, colour = "black") +
    scale_fill_manual(values = setNames(cols, VPD_LABELS), name = "VPD") +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0, 0)) +
    labs(x = NULL, y = expression(Transpiration~(mg~m^{-2}~s^{-1}))) +
    theme_bw(base_size = 14) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(colour = "#e6e6e6"),
          panel.border = element_rect(colour = "black", linewidth = 0.7),
          axis.text = element_text(colour = "black"), axis.text.x = element_text(size = 12),
          legend.position = "top", legend.text = element_text(size = 13),
          plot.background = element_rect(fill = "white", colour = NA))
}

for (nm in names(PALETTES)) {
  ggsave(sprintf("figures/presentation/Fig_TR_bars_%s_%s.png", TAG, nm), make_plot(PALETTES[[nm]]),
         width = 12, height = 5.5, dpi = 300, bg = "white")
}
p <- make_plot(PALETTES[[1]])
ggsave(sprintf("figures/presentation/Fig_TR_bars_%s.png", TAG), p, width = 12, height = 5.5, dpi = 300, bg = "white")
ggsave(sprintf("figures/presentation/Fig_TR_bars_%s.tiff", TAG), p, width = 12, height = 5.5, dpi = 600, bg = "white",
       device = "tiff", compression = "lzw")
cat("Saved figures/presentation/Fig_TR_bars_", TAG, ".png / .tiff (", names(PALETTES)[1], ") and palette options\n", sep = "")
