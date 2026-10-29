## 09_fig_TR_bars_1_3kPa.R
## Bar plot of model-predicted transpiration at VPD = 1 and 3 kPa per accession,
## to show the variability among genotypes.
##   * same accessions and data filters as the presentation figure
##     (08_fig_TR_presentation.R): 9 accessions, rep 5 of PI424231 dropped
##   * per accession, the AIC-selected mixed model (script/model_functions.R)
##     is evaluated at 1 and 3 kPa
##   * accessions are ordered by transpiration at 3 kPa; the label under each
##     accession gives its group and the % increase from 1 to 3 kPa
##
## Input : data/data_clean.csv
## Output: figures/presentation/Fig_TR_bars_1_3kPa.png / .tiff (600 dpi)
##         figures/presentation/TR_predicted_1_3kPa.csv
## Run from experimental_data_soy_2026/:  Rscript script/09_fig_TR_bars_1_3kPa.R
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(readr) })
source("script/model_functions.R")

DROP_GENO <- c("PI398272", "PI404191", "PI442012A", "AG43XF2", "PI377574")
DROP_REPS <- list(PI424231 = "5")
VPD_TARGETS <- c(1, 3)
VPD_COLS <- c("1 kPa" = "#a1d99b", "3 kPa" = "#31a354")
GROUP <- c("slow_wilting" = "Slow wilting", "high_wilting" = "High wilting",
           "Plastic, High WUE" = "High WUE plasticity", "Nonplastic, Low WUE" = "Low WUE plasticity",
           "MAGIC parent" = "Low WUE plasticity")

d <- read_csv("data/data_clean.csv", show_col_types = FALSE) |>
  mutate(Rep = factor(Rep)) |>
  filter(Day_label == "Day 1", !is.na(TR), !is.na(VPD), !(Genotype %in% DROP_GENO))
for (g in names(DROP_REPS)) d <- d |> filter(!(Genotype == g & Rep %in% DROP_REPS[[g]]))
genos <- sort(unique(d$Genotype))

pred <- bind_rows(lapply(genos, function(g) {
  f <- fit_all_models(d |> filter(Genotype == g) |> select(TR, VPD, Rep), response = "TR")
  p <- extract_params(f, f$best)
  v <- switch(f$best,
    "Linear"       = (p[["a"]] %||% p[["(Intercept)"]]) + p[["VPD"]] * VPD_TARGETS,
    "Bilinear"     = bilinear_fn(VPD_TARGETS, p$a, p$b1, p$b2, p$bp),
    "Lin. plateau" = lin_plateau_fn(VPD_TARGETS, p$a, p$b, p$bp))
  tibble(Genotype = g, Best_model = f$best, VPD = VPD_TARGETS, TR = pmax(v, 0))
})) |>
  left_join(d |> distinct(Genotype, WUE_class), by = "Genotype") |>
  mutate(Group = unname(GROUP[WUE_class]), VPD_label = factor(paste(VPD, "kPa"), levels = names(VPD_COLS)))

wide <- pred |> select(Genotype, Group, Best_model, VPD, TR) |>
  pivot_wider(names_from = VPD, values_from = TR, names_prefix = "TR_") |>
  mutate(increase_pct = 100 * (TR_3 - TR_1) / TR_1) |> arrange(TR_3)
write_csv(wide |> mutate(across(where(is.numeric), ~ round(., 2))), "figures/presentation/TR_predicted_1_3kPa.csv")
print(wide |> mutate(across(where(is.numeric), ~ round(., 1))), n = Inf)
cat(sprintf("\nRange across accessions: %.1f-%.1f at 1 kPa, %.1f-%.1f at 3 kPa (mg m-2 s-1)\n",
            min(wide$TR_1), max(wide$TR_1), min(wide$TR_3), max(wide$TR_3)))

## order by transpiration at 3 kPa; x label = accession + group + % increase
lab <- wide |> mutate(x = sprintf("%s\n%s\n+%.0f%%", Genotype, sub(" plasticity", "\nplasticity", sub(" wilting", "\nwilting", Group)), increase_pct))
pred <- pred |> left_join(lab |> select(Genotype, x), by = "Genotype") |>
  mutate(x = factor(x, levels = lab$x))

p <- ggplot(pred, aes(x, TR, fill = VPD_label)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7, colour = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.0f", TR)), position = position_dodge(width = 0.75),
            vjust = -0.45, size = 3.6) +
  scale_fill_manual(values = VPD_COLS, name = "VPD") +
  scale_y_continuous(limits = c(0, 90), breaks = seq(0, 80, 20), expand = c(0, 0)) +
  labs(x = NULL, y = expression(Predicted~transpiration~(mg~m^{-2}~s^{-1}))) +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "#e6e6e6"),
        panel.border = element_rect(colour = "black", linewidth = 0.7),
        axis.text = element_text(colour = "black"), axis.text.x = element_text(size = 11, lineheight = 0.95),
        legend.position = "top", legend.text = element_text(size = 13),
        plot.background = element_rect(fill = "white", colour = NA))

ggsave("figures/presentation/Fig_TR_bars_1_3kPa.png", p, width = 13, height = 6, dpi = 300, bg = "white")
ggsave("figures/presentation/Fig_TR_bars_1_3kPa.tiff", p, width = 13, height = 6, dpi = 600, bg = "white",
       device = "tiff", compression = "lzw")
cat("Saved figures/presentation/Fig_TR_bars_1_3kPa.png / .tiff\n")
