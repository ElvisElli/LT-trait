## 07_supp_table_model_selection.R
## Supplementary tables: candidate models (linear, bilinear, linear-plateau)
## fitted per genotype for transpiration (TR, Table S1) and stomatal
## conductance (gsw, Table S2), with parameters, AIC and delta AIC, and the
## selected model flagged. Built from the parameter tables written by
## 02_model_fitting_TR_gs_VPD.R.
##
## Output: figures/paper_figures/Supplementary_Tables_model_selection.docx
##         figures/paper_figures/Table_S1_TR_model_selection.csv
##         figures/paper_figures/Table_S2_gsw_model_selection.csv
## Run from experimental_data_soy_2026/ in a UTF-8 locale (the file holds
## Unicode symbols):  LC_ALL=C.UTF-8 Rscript script/07_supp_table_model_selection.R
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(flextable); library(officer) })
invisible(Sys.setlocale("LC_CTYPE", "C.UTF-8"))

WUE_ORDER  <- c("Plastic, High WUE", "Nonplastic, Low WUE", "MAGIC parent", "slow_wilting", "high_wilting")
WUE_LABELS <- c("Plastic, High WUE" = "Plastic / high WUE", "Nonplastic, Low WUE" = "Nonplastic / low WUE",
                "MAGIC parent" = "MAGIC parent", "slow_wilting" = "Slow-wilting check", "high_wilting" = "High-wilting check")
MODEL_ORDER <- c("Linear", "Bilinear", "Lin. plateau")
MODEL_LABELS <- c("Linear" = "Linear", "Bilinear" = "Bilinear", "Lin. plateau" = "Linear-plateau")

n_obs <- read_csv("data/data_clean.csv", show_col_types = FALSE) |>
  filter(Day_label == "Day 1") |>
  group_by(Genotype) |>
  summarise(n_TR = sum(!is.na(TR) & !is.na(VPD)), n_gsw = sum(!is.na(gsw) & !is.na(VPD)), .groups = "drop")

build <- function(params, n_col, digits = 2) {
  # complete the genotype x model grid so non-converged fits appear as rows
  grid <- params |> distinct(Genotype, WUE_class) |>
    crossing(Model = MODEL_ORDER)
  tab <- grid |>
    left_join(params |> select(-WUE_class), by = c("Genotype", "Model")) |>
    left_join(n_obs |> select(Genotype, n = all_of(n_col)), by = "Genotype") |>
    mutate(WUE_class = factor(WUE_class, levels = WUE_ORDER),
           Model = factor(Model, levels = MODEL_ORDER)) |>
    arrange(WUE_class, Genotype, Model) |>
    mutate(converged = !is.na(AIC),
           Selected  = ifelse(!is.na(Best_model) & Best_model == "yes", "✓", ""),
           `WUE class` = unname(WUE_LABELS[as.character(WUE_class)]),
           Model = unname(MODEL_LABELS[as.character(Model)]),
           Fit = case_when(!converged ~ "", Fit == "lme" ~ "lme", Fit == "nlme" ~ "nlme",
                           Fit == "nls" ~ "nls*", Fit == "lm" ~ "lm*", TRUE ~ Fit),
           Intercept = ifelse(converged, formatC(a, format = "f", digits = digits), "not converged"),
           `Slope 1` = ifelse(converged, formatC(Slope1, format = "f", digits = digits + 1), ""),
           `Slope 2` = ifelse(!is.na(Slope2), formatC(Slope2, format = "f", digits = digits + 1), "–"),
           `Breakpoint (kPa)` = ifelse(!is.na(BP_kPa), formatC(BP_kPa, format = "f", digits = 2), "–"),
           AIC = ifelse(converged, formatC(AIC, format = "f", digits = 1), ""),
           `ΔAIC` = ifelse(converged, formatC(Delta_AIC, format = "f", digits = 1), "")) |>
    select(Genotype, `WUE class`, n, Model, Fit, Intercept, `Slope 1`, `Slope 2`, `Breakpoint (kPa)`, AIC, `ΔAIC`, Selected)
  tab
}

tr  <- build(read_csv("data/TR_model_parameters.csv",  show_col_types = FALSE), "n_TR",  digits = 1)
gsw <- build(read_csv("data/gsw_model_parameters.csv", show_col_types = FALSE), "n_gsw", digits = 3)
write_csv(tr,  "figures/paper_figures/Table_S1_TR_model_selection.csv")
write_csv(gsw, "figures/paper_figures/Table_S2_gsw_model_selection.csv")

make_ft <- function(tab) {
  sel_rows <- which(tab$Selected == "✓")
  # show genotype, class and n only on the first row of each genotype block
  first <- !duplicated(tab$Genotype)
  tab <- tab |> mutate(Genotype = ifelse(first, Genotype, ""), n = ifelse(first, as.character(n), ""),
                       `WUE class` = ifelse(first, `WUE class`, ""))
  ft <- flextable(tab) |>
    bold(i = sel_rows, j = c("Model", "AIC", "Selected")) |>
    align(j = 3:12, align = "center", part = "all") |>
    align(j = 1:2, align = "left", part = "all") |>
    fontsize(size = 8, part = "all") |>
    font(fontname = "Calibri", part = "all") |>
    padding(padding = 2, part = "all") |>
    border_remove() |>
    hline_top(border = fp_border(width = 1), part = "header") |>
    hline_bottom(border = fp_border(width = 1), part = "header") |>
    hline_bottom(border = fp_border(width = 1), part = "body") |>
    set_table_properties(layout = "autofit", width = 1)
  block_end <- which(first) - 1
  block_end <- block_end[block_end > 0]
  ft <- hline(ft, i = block_end, border = fp_border(width = 0.5, color = "#999999"))
  ft
}

cap1 <- paste0(
  "Table S1. Candidate models fitted to plant transpiration (TR, mg m⁻² s⁻¹) as a function of vapor ",
  "pressure deficit (VPD, kPa) for each soybean accession. Linear: TR = a + b₁·VPD; bilinear: slope b₁ up ",
  "to the breakpoint and b₂ above it; linear-plateau: slope b₁ up to the breakpoint, constant above it. ",
  "Models were fitted by maximum likelihood with replicate as a random intercept (nlme) to Day 1 measurements and ",
  "compared by AIC; the model with the lowest AIC (ΔAIC = 0, ✓) was selected. n = number of plant-level ",
  "observations. Fit: lme and nlme are mixed models with the replicate random intercept; nls* (no random ",
  "effect) is the fallback used when the mixed-model fit did not converge from the nls starting values. ",
  "“Not converged”: the non-linear fit did not converge from any starting value. WUE class ",
  "follows Chamarthi et al. (2022); PI398198 is a MAGIC parent.")
cap2 <- paste0(
  "Table S2. Candidate models fitted to stomatal conductance (gsw, mol m⁻² s⁻¹) as a function ",
  "of VPD for the accessions measured with a porometer. Model definitions, fitting and selection as in Table S1. ",
  "Non-linear models did not converge for any accession, so the linear model was retained throughout.")

doc <- read_docx() |>
  body_add_par(cap1, style = "Normal") |>
  body_add_flextable(make_ft(tr)) |>
  body_add_par("", style = "Normal") |>
  body_add_break() |>
  body_add_par(cap2, style = "Normal") |>
  body_add_flextable(make_ft(gsw))
doc <- body_end_section_landscape(doc)
print(doc, target = "figures/paper_figures/Supplementary_Tables_model_selection.docx")
cat("Saved figures/paper_figures/Supplementary_Tables_model_selection.docx\n")
cat(sprintf("Table S1: %d rows (%d genotypes); Table S2: %d rows (%d genotypes)\n",
            nrow(tr), n_distinct(tr$Genotype), nrow(gsw), n_distinct(gsw$Genotype)))
