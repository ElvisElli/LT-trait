## 10_fig_gsw_presentation.R
## Presentation (backup slide) version of the stomatal conductance figure:
## gsw vs VPD for the same 9 accessions and panel order as the TR presentation
## figure (08_fig_TR_presentation.R), Day 1 data, same mixed models and AIC
## selection (script/model_functions.R). Each panel is labelled with the
## predicted change in gsw from 1.5 to 3.5 kPa.
##   * gsw (porometer) was measured on a subset of plants, mostly reps 3-5, so n
##     per accession is smaller than for TR; all Day 1 gsw values are kept
##     (the rep 5 drop for PI424231 in the TR figure is TR-specific)
##   * accessions with fewer than MIN_N gsw values are plotted without a fit
##     (A5959: 5 values from 2 plants; the rep intercepts absorb the VPD effect)
##
## Input : data/data_clean.csv (01_data_inspection_QC.R)
## Also a raw-data version: every individual porometer reading (gsw1-gsw3 per
## plant and step), Day 1 and Day 2, no fits.
##
## Output: figures/presentation/Fig_gsw_response_VPD_presentation.png / .tiff
##         figures/presentation/Fig_gsw_raw_VPD_presentation.png / .tiff
##         figures/presentation/gsw_best_model_presentation.csv
## Run from experimental_data_soy_2026/:  Rscript script/10_fig_gsw_presentation.R
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(readr) })
source("script/model_functions.R")

BLUE  <- "#2c7fb8"
MIN_N <- 6
ORDER <- c("PI398492", "PI398617", "PI424231", "PI603418D",
           "PI603171", "PI398198",
           "A5959", "PI416937", "PI471938")
dir.create("figures/presentation", showWarnings = FALSE)

d <- read_csv("data/data_clean.csv", show_col_types = FALSE) |>
  mutate(Rep = factor(Rep)) |>
  filter(Day_label == "Day 1", !is.na(gsw), !is.na(VPD), Genotype %in% ORDER)
cat("gsw:", n_distinct(d$Genotype), "accessions,", nrow(d), "observations\n")

## fit per accession, predict at 1.5 and 3.5 kPa
FIT <- ORDER[sapply(ORDER, function(g) sum(d$Genotype == g) >= MIN_N)]
fits <- lapply(setNames(FIT, FIT), function(g)
  fit_all_models(d |> filter(Genotype == g) |> select(gsw, VPD, Rep), response = "gsw"))
at <- function(g, v) predict_curve(fits[[g]], fits[[g]]$best, v)$y_hat
best <- bind_rows(lapply(FIT, function(g) {
  f <- fits[[g]]; p <- extract_params(f, f$best)
  tibble(Genotype = g, n = sum(d$Genotype == g), Best_model = f$best,
         Intercept = p[["a"]] %||% p[["(Intercept)"]],
         Slope1 = p[["VPD"]] %||% p[["b1"]] %||% p[["b"]],
         BP_kPa = p[["bp"]] %||% NA_real_,
         gsw_1.5 = at(g, 1.5), gsw_3.5 = at(g, 3.5))
})) |> mutate(change_pct = 100 * (gsw_3.5 / gsw_1.5 - 1),
              across(where(is.numeric), ~ round(., 3)))
print(best, n = Inf, width = 200)
write_csv(best, "figures/presentation/gsw_best_model_presentation.csv")

## curves and panel labels
vpd_seq <- seq(min(d$VPD) * 0.95, max(d$VPD) * 1.05, length.out = 300)
curves <- bind_rows(lapply(FIT, function(g)
  predict_curve(fits[[g]], fits[[g]]$best, vpd_seq) |> mutate(Genotype = g)))
labs_df <- tibble(Genotype = ORDER, n = as.integer(table(factor(d$Genotype, levels = ORDER)))) |>
  left_join(best |> select(Genotype, change_pct), by = "Genotype") |>
  mutate(lab = ifelse(is.na(change_pct), sprintf("not fitted  (n = %d)", n),
                      sprintf("%+.0f%%  (n = %d)", change_pct, n)))
d$Genotype       <- factor(d$Genotype, levels = ORDER)
curves$Genotype  <- factor(curves$Genotype, levels = ORDER)
labs_df$Genotype <- factor(labs_df$Genotype, levels = ORDER)

p <- ggplot(d, aes(VPD, gsw)) +
  geom_point(size = 2.8, shape = 21, colour = "black", fill = BLUE, stroke = 0.7, alpha = 0.9) +
  geom_line(data = curves, aes(VPD, y_hat), colour = BLUE, linewidth = 1.2) +
  geom_text(data = labs_df, aes(x = 4.05, y = 1.12, label = lab), hjust = 1, size = 4.6) +
  facet_wrap(~Genotype, ncol = 3) +
  coord_cartesian(ylim = c(0, 1.2)) +
  scale_x_continuous(breaks = 1:4) +
  labs(x = "Vapor Pressure Deficit (kPa)",
       y = expression(Stomatal~conductance~(mol~m^{-2}~s^{-1}))) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 0.7, fill = NA),
        strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.7),
        strip.text = element_text(size = 13, colour = "black"),
        axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 15),
        plot.background = element_rect(fill = "white", colour = NA))
ggsave("figures/presentation/Fig_gsw_response_VPD_presentation.png", p, width = 9, height = 8.2, dpi = 300, bg = "white")
ggsave("figures/presentation/Fig_gsw_response_VPD_presentation.tiff", p, width = 9, height = 8.2, dpi = 600, bg = "white",
       device = "tiff", compression = "lzw", type = "cairo")
cat("Saved figures/presentation/Fig_gsw_response_VPD_presentation.png\n")

## ── Raw data: all individual readings, both days ─────────────────────────────
DAY_COL <- c("Day 1" = BLUE, "Day 2" = "#fdae61")
raw <- read_csv("data/data_clean.csv", show_col_types = FALSE) |>
  filter(Genotype %in% ORDER, !is.na(VPD)) |>
  pivot_longer(c(gsw1, gsw2, gsw3), names_to = "reading", values_to = "gsw_raw") |>
  filter(!is.na(gsw_raw)) |>
  mutate(Genotype = factor(Genotype, levels = ORDER), Day_label = factor(Day_label, levels = names(DAY_COL)))
n_raw <- raw |> count(Genotype, Day_label) |> pivot_wider(names_from = Day_label, values_from = n, values_fill = 0) |>
  mutate(lab = sprintf("n = %d + %d", `Day 1`, `Day 2`))
cat("Raw gsw readings:", nrow(raw), "\n"); print(count(raw, Day_label))
p_raw <- ggplot(raw, aes(VPD, gsw_raw, fill = Day_label)) +
  geom_point(size = 2.2, shape = 21, colour = "black", stroke = 0.4, alpha = 0.75) +
  geom_text(data = n_raw, aes(x = 4.05, y = 1.43, label = lab), inherit.aes = FALSE, hjust = 1, size = 4.4) +
  facet_wrap(~Genotype, ncol = 3) +
  scale_fill_manual(values = DAY_COL, name = NULL) +
  coord_cartesian(ylim = c(0, 1.5)) +
  scale_x_continuous(breaks = 1:4) +
  labs(x = "Vapor Pressure Deficit (kPa)",
       y = expression(Stomatal~conductance~(mol~m^{-2}~s^{-1}))) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(), legend.position = "top",
        panel.border = element_rect(colour = "black", linewidth = 0.7, fill = NA),
        strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.7),
        strip.text = element_text(size = 13, colour = "black"),
        axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 15), legend.text = element_text(size = 13),
        plot.background = element_rect(fill = "white", colour = NA))
ggsave("figures/presentation/Fig_gsw_raw_VPD_presentation.png", p_raw, width = 9, height = 8.6, dpi = 300, bg = "white")
ggsave("figures/presentation/Fig_gsw_raw_VPD_presentation.tiff", p_raw, width = 9, height = 8.6, dpi = 600, bg = "white",
       device = "tiff", compression = "lzw", type = "cairo")
cat("Saved figures/presentation/Fig_gsw_raw_VPD_presentation.png\n")
