## 08_fig_TR_presentation.R
## Presentation version of Figure 1 (transpiration vs VPD per accession).
##   * accessions with only two replications are dropped (PI398272, PI404191,
##     PI442012A); the commercial check AG43XF2 and PI377574 (very high
##     scatter between plants) are dropped as well
##   * replication 5 of PI424231 is dropped (clearly offset from reps 1-4)
##   * all points in one colour, curves refitted on the retained data with the
##     same mixed models as the paper (script/model_functions.R)
##
## Input : data/data_clean.csv (01_data_inspection_QC.R)
## Output: figures/presentation/Fig_TR_response_VPD_presentation.png
##         figures/presentation/TR_best_model_presentation.csv
## Run from experimental_data_soy_2026/:  Rscript script/08_fig_TR_presentation.R
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(readr) })
source("script/model_functions.R")

DROP_GENO <- c("PI398272", "PI404191", "PI442012A", "AG43XF2", "PI377574")
DROP_REPS <- list(PI424231 = "5")
GREEN <- "#31a354"
dir.create("figures/presentation", showWarnings = FALSE)

WUE_LABELS <- c("slow_wilting" = "Slow wilting", "high_wilting" = "High wilting",
                "Plastic, High WUE" = "High WUE plasticity", "Nonplastic, Low WUE" = "Low WUE plasticity",
                "MAGIC parent" = "Low WUE plasticity")     # PI398198 shown with the low-plasticity lines
ORDER <- c("PI398492", "PI398617", "PI424231", "PI603418D",
           "PI603171", "PI398198",
           "A5959", "PI416937", "PI471938")

d <- read_csv("data/data_clean.csv", show_col_types = FALSE) |>
  mutate(Rep = factor(Rep)) |>
  filter(Day_label == "Day 1", !is.na(TR), !is.na(VPD), !(Genotype %in% DROP_GENO))
for (g in names(DROP_REPS)) d <- d |> filter(!(Genotype == g & Rep %in% DROP_REPS[[g]]))
cat("Retained:", n_distinct(d$Genotype), "accessions,", nrow(d), "observations\n")

## fit per accession
fits <- lapply(setNames(ORDER, ORDER), function(g) {
  fit_all_models(d |> filter(Genotype == g) |> select(TR, VPD, Rep), response = "TR")
})
best <- bind_rows(lapply(ORDER, function(g) {
  f <- fits[[g]]; p <- extract_params(f, f$best)
  tibble(Genotype = g, n = sum(d$Genotype == g), Best_model = f$best,
         Intercept = p[["a"]] %||% p[["(Intercept)"]],
         Slope1 = p[["VPD"]] %||% p[["b1"]] %||% p[["b"]],
         BP_kPa = p[["bp"]] %||% NA_real_, Slope2 = p[["b2"]] %||% NA_real_,
         AIC_linear = f$AIC["Linear"], AIC_bilinear = f$AIC["Bilinear"], AIC_plateau = f$AIC["Lin. plateau"])
})) |> mutate(across(where(is.numeric), ~ round(., 3)))
print(best, n = Inf, width = 200)
write_csv(best, "figures/presentation/TR_best_model_presentation.csv")

## curves and panel labels
vpd_seq <- seq(min(d$VPD) * 0.95, max(d$VPD) * 1.05, length.out = 300)
curves <- bind_rows(lapply(ORDER, function(g) {
  predict_curve(fits[[g]], fits[[g]]$best, vpd_seq) |> mutate(Genotype = g)
}))
wue <- d |> distinct(Genotype, WUE_class)
labels <- tibble(Genotype = ORDER) |> left_join(wue, by = "Genotype") |>
  mutate(Panel = paste0(Genotype, "\n(", WUE_LABELS[WUE_class], ")"))
d      <- d      |> left_join(labels |> select(Genotype, Panel), by = "Genotype")
curves <- curves |> left_join(labels |> select(Genotype, Panel), by = "Genotype")
d$Panel      <- factor(d$Panel, levels = labels$Panel)
curves$Panel <- factor(curves$Panel, levels = labels$Panel)

p <- ggplot(d, aes(VPD, TR)) +
  geom_point(size = 2.6, shape = 21, colour = "black", fill = GREEN, stroke = 0.7, alpha = 0.9) +
  geom_line(data = curves, aes(VPD, y_hat), colour = GREEN, linewidth = 1.2) +
  facet_wrap(~Panel, ncol = 3) +
  coord_cartesian(ylim = c(0, 120)) +
  scale_x_continuous(breaks = 1:4) +
  labs(x = "Vapor Pressure Deficit (kPa)", y = expression(Transpiration~(mg~m^{-2}~s^{-1}))) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 0.7, fill = NA),
        strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.7),
        strip.text = element_text(size = 13, colour = "black"),
        axis.text = element_text(size = 13, colour = "black"),
        axis.title = element_text(size = 15),
        plot.background = element_rect(fill = "white", colour = NA))
ggsave("figures/presentation/Fig_TR_response_VPD_presentation.png", p, width = 12, height = 10.5, dpi = 300, bg = "white")
ggsave("figures/presentation/Fig_TR_response_VPD_presentation.tiff", p, width = 12, height = 10.5, dpi = 600, bg = "white",
       device = "tiff", compression = "lzw", type = "cairo")
cat("Saved figures/presentation/Fig_TR_response_VPD_presentation.png\n")
