## 06_predicted_gsw_TR_by_group.R
## Fit the per-genotype mixed models (linear / bilinear / linear-plateau, Rep
## as random intercept, AIC selection; Day 1 data, as in the pipeline) for
## stomatal conductance (gsw) and transpiration (TR), predict both at VPD = 1
## and 3 kPa from the best model, and compare plasticity groups:
##   Plastic / high WUE, Nonplastic / low WUE (PI398198 included, as in the
##   paper figures), Wilting checks.
##
## Input : data/data_clean.csv (01_data_inspection_QC.R)
## Output: data/predicted_TR_gsw_at_VPD_1_3_by_genotype.csv
##         data/predicted_TR_gsw_at_VPD_1_3_by_group.csv
##         figures/predicted_gsw_TR_by_group.png
## Run from experimental_data_soy_2026/:  Rscript script/06_predicted_gsw_TR_by_group.R
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(readr); library(patchwork); library(ggrepel)
})
source("script/model_functions.R")

VPD_TARGETS <- c(1, 3)
COLS <- c("Plastic / high WUE" = "#4393c3", "Nonplastic / low WUE" = "#d4824a", "Wilting checks" = "#9d2235")

d <- read_csv("data/data_clean.csv", show_col_types = FALSE) |>
  mutate(Rep = factor(Rep)) |>
  filter(Day_label == "Day 1") |>
  mutate(Group = case_when(
    WUE_class == "Plastic, High WUE"                          ~ "Plastic / high WUE",
    WUE_class %in% c("Nonplastic, Low WUE", "MAGIC parent") ~ "Nonplastic / low WUE",
    TRUE                                                     ~ "Wilting checks"),
    Group = factor(Group, levels = names(COLS)))
groups <- d |> distinct(Genotype, Group)

## ── Fit per genotype and predict at the target VPDs ───────────────────────────
fit_predict <- function(df, col, label) {
  genos <- sort(unique(df$Genotype))
  bind_rows(lapply(genos, function(g) {
    dd <- df |> filter(Genotype == g) |> select(all_of(col), VPD, Rep) |> rename(TR = all_of(col))
    fit <- fit_all_models(dd, response = "TR")
    if (is.null(fit)) return(NULL)
    pars <- extract_params(fit, fit$best)
    vals <- switch(fit$best,
      "Linear"       = (pars[["a"]] %||% pars[["(Intercept)"]]) + pars[["VPD"]] * VPD_TARGETS,
      "Bilinear"     = bilinear_fn(VPD_TARGETS, pars$a, pars$b1, pars$b2, pars$bp),
      "Lin. plateau" = lin_plateau_fn(VPD_TARGETS, pars$a, pars$b, pars$bp))
    tibble(variable = label, Genotype = g, Best_model = fit$best, n = nrow(na.omit(dd)),
           VPD = VPD_TARGETS, predicted = pmax(vals, 0))
  }))
}

cat("Fitting TR models ...\n")
tr_pred  <- fit_predict(d, "TR", "TR")
cat("Fitting gsw models ...\n")
gsw_pred <- fit_predict(d |> filter(!is.na(gsw)), "gsw", "gsw")

pred <- bind_rows(tr_pred, gsw_pred) |> left_join(groups, by = "Genotype") |>
  mutate(VPD_label = factor(paste(VPD, "kPa"), levels = c("1 kPa", "3 kPa")))
write_csv(pred, "data/predicted_TR_gsw_at_VPD_1_3_by_genotype.csv")

grp <- pred |> group_by(variable, Group, VPD, VPD_label) |>
  summarise(n_geno = n(), mean = mean(predicted), se = sd(predicted) / sqrt(n()), .groups = "drop")
write_csv(grp, "data/predicted_TR_gsw_at_VPD_1_3_by_group.csv")

chg <- grp |> select(variable, Group, VPD, mean) |> pivot_wider(names_from = VPD, values_from = mean, names_prefix = "kPa_") |>
  mutate(change_pct = 100 * (kPa_3 - kPa_1) / kPa_1)
cat("\nModel-predicted group means (mean of genotype predictions, SE across genotypes):\n")
print(grp |> mutate(across(c(mean, se), ~ round(., 3))), n = Inf)
cat("\nChange from 1 to 3 kPa:\n"); print(chg |> mutate(across(where(is.numeric), ~ round(., 3))))
cat("\nGenotype predictions:\n")
print(pred |> select(variable, Group, Genotype, Best_model, VPD, predicted) |>
        pivot_wider(names_from = VPD, values_from = predicted, names_prefix = "kPa_") |>
        mutate(across(where(is.numeric), ~ round(., 3))) |> arrange(variable, Group), n = Inf)

## ── Plot: rows = variable, left = group mean +/- SE, right = genotype slopes ─
thm <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "#e6e6e6"),
        panel.border = element_rect(colour = "black", linewidth = 0.7),
        axis.text = element_text(colour = "black"), legend.position = "top",
        plot.title = element_text(size = 12), plot.background = element_rect(fill = "white", colour = NA))
pd <- position_dodge(width = 0.4)

panel_group <- function(var, ylab, title) {
  g <- filter(grp, variable == var); p <- filter(pred, variable == var)
  ggplot() +
    geom_point(data = p, aes(VPD_label, predicted, colour = Group),
               position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.4, seed = 1),
               alpha = 0.35, size = 1.8) +
    geom_line(data = g, aes(VPD_label, mean, colour = Group, group = Group), position = pd, linewidth = 0.9) +
    geom_errorbar(data = g, aes(VPD_label, ymin = mean - se, ymax = mean + se, colour = Group),
                  position = pd, width = 0, linewidth = 0.9) +
    geom_point(data = g, aes(VPD_label, mean, fill = Group), position = pd, shape = 21,
               colour = "white", size = 4, stroke = 0.8) +
    geom_text(data = g, aes(VPD_label, mean, label = if (var == "gsw") sprintf("%.2f", mean) else sprintf("%.0f", mean),
                            group = Group), position = pd, vjust = -1.3, size = 3.1, colour = "#333333") +
    scale_colour_manual(values = COLS, name = NULL) + scale_fill_manual(values = COLS, name = NULL) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.12))) +
    labs(x = "VPD", y = ylab, title = title) + thm
}
panel_geno <- function(var, title) {
  p <- filter(pred, variable == var)
  ggplot(p, aes(VPD_label, predicted, group = Genotype, colour = Group)) +
    geom_line(linewidth = 0.7, alpha = 0.9) + geom_point(size = 2.2) +
    geom_text_repel(data = filter(p, VPD == 3), aes(label = Genotype), direction = "y", hjust = 0,
                    nudge_x = 0.12, size = 2.6, segment.size = 0.25, segment.colour = "grey60",
                    box.padding = 0.1, min.segment.length = 0, seed = 1, show.legend = FALSE) +
    scale_colour_manual(values = COLS, name = NULL) +
    scale_x_discrete(expand = expansion(mult = c(0.1, 0.6))) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.12))) +
    labs(x = "VPD", y = NULL, title = title) + thm + theme(legend.position = "none")
}

p <- (panel_group("gsw", expression(g[sw]~(mol~m^{-2}~s^{-1})), "A  Predicted gsw, group mean +/- SE across genotypes") +
      panel_geno("gsw", "B  Predicted gsw per genotype")) /
     (panel_group("TR", expression(Transpiration~(mg~m^{-2}~s^{-1})), "C  Predicted transpiration, group mean +/- SE") +
      panel_geno("TR", "D  Predicted transpiration per genotype")) +
  plot_layout(guides = "collect") & theme(legend.position = "top")
ggsave("figures/predicted_gsw_TR_by_group.png", p, width = 11, height = 10, dpi = 200, bg = "white")
ggsave("figures/predicted_gsw_TR_by_group.tiff", p, width = 11, height = 10, dpi = 600, bg = "white",
       device = "tiff", compression = "lzw", type = "cairo")
cat("\nSaved figures/predicted_gsw_TR_by_group.png\n")
