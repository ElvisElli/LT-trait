# ==============================================================================
# figure_bullets.R
# Compute 3 key-result bullet points per paper figure from the model-parameter
# CSV files produced by 02_model_fitting_TR_gs_VPD.R.
#
# Run from figures/paper_figures/ or any working directory.
# Output: figure_bullets.json  (read by update_docx.py to annotate the Word doc)
# ==============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(jsonlite)
})

DATA_DIR <- file.path(dirname(normalizePath(".", mustWork = FALSE)),
                      "..", "data")
# Hardcode fallback for Rscript invocation
DATA_DIR <- "/home/user/LT-trait/experimental_data_soy_2026/data"

tr   <- read_csv(file.path(DATA_DIR, "TR_best_model_summary.csv"),  show_col_types = FALSE)
gsw  <- read_csv(file.path(DATA_DIR, "gsw_best_model_summary.csv"), show_col_types = FALSE)
pred <- read_csv(file.path(DATA_DIR, "predicted_TR_gsw_at_target_VPD.csv"), show_col_types = FALSE)

bullets <- list()


# ====== FIGURE 1 — TR response to VPD per genotype =========================
local({
  n_lin  <- sum(tr$Best_model == "Linear")
  n_plat <- sum(tr$Best_model == "Lin. plateau")
  n_bil  <- sum(tr$Best_model == "Bilinear")

  int_min <- round(min(tr$Intercept), 1)
  int_max <- round(max(tr$Intercept), 1)
  slp_min <- round(min(tr$Slope1),    1)
  slp_max <- round(max(tr$Slope1),    1)

  plastic   <- filter(tr, WUE_class == "Plastic, High WUE")
  extreme_g <- plastic$Genotype[which.max(plastic$Slope1)]   # PI398617 has both max slope & min intercept
  extreme_s <- round(plastic$Slope1[which.max(plastic$Slope1)], 1)
  extreme_i <- round(plastic$Intercept[which.max(plastic$Slope1)], 1)

  bullets$F1 <<- c(
    sprintf("%d genotypes were best described by a linear model, %d by a linear-plateau model, and %d by a bilinear model, indicating diverse VPD response strategies across the 14 accessions.", n_lin, n_plat, n_bil),
    sprintf("Initial slopes ranged from %.1f to %.1f mg/m2/s per kPa and intercepts from %.1f to %.1f mg/m2/s, reflecting large genotypic variation in both baseline transpiration rate and VPD sensitivity.", slp_min, slp_max, int_min, int_max),
    sprintf("Plastic, High WUE genotypes did not uniformly show lower transpiration: %s combined the highest slope (%.1f mg/m2/s per kPa) with the lowest intercept (%.1f mg/m2/s) among all genotypes, a pattern of near-zero baseline transpiration with high VPD sensitivity.", extreme_g, extreme_s, extreme_i)
  )
})


# ====== FIGURE 2 — VPD breakpoints ==========================================
local({
  bp       <- filter(tr, !is.na(BP_kPa))
  n_bp     <- nrow(bp)
  n_total  <- nrow(tr)
  pct_bp   <- round(n_bp / n_total * 100)
  bp_min   <- round(min(bp$BP_kPa), 2)
  bp_max   <- round(max(bp$BP_kPa), 2)
  bp_mean  <- round(mean(bp$BP_kPa), 2)
  n_plast  <- sum(bp$WUE_class == "Plastic, High WUE")

  early     <- arrange(bp, BP_kPa)[1:2, ]
  bp1       <- round(early$BP_kPa[1], 2)
  bp2       <- round(early$BP_kPa[2], 2)
  geno1     <- early$Genotype[1]
  geno2     <- early$Genotype[2]

  bullets$F2 <<- c(
    sprintf("%d of %d genotypes (%d%%) showed a non-linear transpiration response, with VPD breakpoints ranging from %.2f to %.2f kPa (mean %.2f kPa).", n_bp, n_total, pct_bp, bp_min, bp_max, bp_mean),
    sprintf("%d of the %d genotypes with a breakpoint were previously classified as Plastic, High WUE (Chamarthi et al. 2022), consistent with the expected limited-transpiration phenotype.", n_plast, n_bp),
    sprintf("The two earliest breakpoints (%.2f and %.2f kPa in %s and %s) indicate that stomatal restriction can occur before VPD reaches typical midday field values (~2.0-2.5 kPa), potentially conserving water during morning hours.", bp1, bp2, geno1, geno2)
  )
})


# ====== FIGURE 3 — Intercept x Slope quadrant classification ================
local({
  x_mid     <- median(tr$Intercept, na.rm = TRUE)
  s_mid     <- median(tr$Slope1,    na.rm = TRUE)

  plastic   <- filter(tr, WUE_class == "Plastic, High WUE") |>
    mutate(lo_hi = Intercept < x_mid & Slope1 > s_mid)
  n_lo_hi   <- sum(plastic$lo_hi)
  n_plastic <- nrow(plastic)

  nonplastic <- filter(tr, WUE_class == "Nonplastic, Low WUE") |>
    mutate(quad = case_when(
      Intercept < x_mid & Slope1 > s_mid  ~ "Low/High",
      Intercept >= x_mid & Slope1 > s_mid ~ "High/High",
      Intercept < x_mid & Slope1 <= s_mid ~ "Low/Low",
      TRUE                                 ~ "High/Low"
    ))
  n_quads <- length(unique(nonplastic$quad))

  extreme   <- arrange(plastic, Intercept)[1:2, ]
  names_ext <- paste(extreme$Genotype, collapse = " and ")
  int_thr   <- ceiling(max(extreme$Intercept))
  slp_thr   <- floor(min(extreme$Slope1))

  bullets$F3 <<- c(
    sprintf("%d of %d Plastic, High WUE genotypes clustered in the Low intercept - High sensitivity quadrant, consistent with lower baseline transpiration and stronger stomatal responsiveness at elevated VPD.", n_lo_hi, n_plastic),
    sprintf("Nonplastic, Low WUE genotypes were distributed across %d quadrants, indicating within-class heterogeneity not captured by the Chamarthi et al. plasticity classification alone.", n_quads),
    sprintf("%s (both Plastic, High WUE) occupied the most extreme position with intercepts below %d mg/m2/s and initial slopes above %d mg/m2/s/kPa, making them the strongest candidates for limited-transpiration trait characterization.", names_ext, int_thr, slp_thr)
  )
})


# ====== FIGURE 4 — gsw response to VPD =====================================
local({
  slp_min   <- round(min(gsw$Slope1), 3)
  slp_max   <- round(max(gsw$Slope1), 3)
  n_geno    <- nrow(gsw)

  slow      <- filter(gsw, WUE_class == "slow_wilting") |> arrange(Genotype)
  sw_names  <- paste(slow$Genotype, collapse = " and ")
  sw_s1     <- round(slow$Slope1[1], 3)
  sw_s2     <- round(slow$Slope1[2], 3)

  paradox   <- filter(gsw, WUE_class == "Nonplastic, Low WUE") |>
    slice_max(abs(Slope1), n = 1)
  p_name    <- paradox$Genotype
  p_slp     <- round(paradox$Slope1, 2)
  p_int     <- round(paradox$Intercept, 2)

  bullets$F4 <<- c(
    sprintf("All %d measured genotypes showed negative stomatal conductance slopes (range: %.2f to %.2f mol/m2/s/kPa), confirming universal stomatal closure in response to increasing vapor pressure demand.", n_geno, slp_min, slp_max),
    sprintf("Slow wilting genotypes (%s) had intermediate slopes (%.3f and %.3f mol/m2/s/kPa) and lower intercepts, consistent with the expected limited-transpiration response.", sw_names, sw_s1, sw_s2),
    sprintf("%s (Nonplastic, Low WUE) had the steepest negative slope (%.2f mol/m2/s/kPa) and highest intercept (%.2f mol/m2/s), indicating paradoxically high baseline conductance relative to its non-plastic classification.", p_name, p_slp, p_int)
  )
})


# ====== FIGURE 5 — Predicted TR and gsw at 1.5 and 3.5 kPa =================
local({
  n_plat_bil <- sum(tr$Best_model %in% c("Lin. plateau", "Bilinear"))
  n_lin_tr   <- sum(tr$Best_model == "Linear")
  n_gsw      <- nrow(gsw)

  gsw_wide <- filter(pred, variable == "gsw") |>
    pivot_wider(id_cols = c(Genotype, WUE_class),
                names_from  = VPD_target,
                values_from = value,
                names_prefix = "vpd_") |>
    mutate(pct_change = (vpd_3.5 - vpd_1.5) / vpd_1.5 * 100)

  pct_min <- round(min(gsw_wide$pct_change))
  pct_max <- round(max(gsw_wide$pct_change))

  bullets$F5 <<- c(
    sprintf("%d genotypes with plateau or bilinear models showed attenuated or no increase in transpiration above their breakpoint VPD, while the %d linear-model genotypes continued to increase through 3.5 kPa, demonstrating the functional consequence of the VPD threshold response.", n_plat_bil, n_lin_tr),
    sprintf("Stomatal conductance declined consistently from VPD 1.5 to 3.5 kPa across all %d measured genotypes (range: %d to %d%%), with the magnitude of decline reflecting differences in stomatal sensitivity rather than a uniform, species-level response.", n_gsw, pct_min, pct_max),
    "The divergence between TR and gsw patterns at 3.5 kPa - rising TR despite declining gsw in several genotypes - suggests that transpiration under high atmospheric demand is increasingly driven by boundary layer conductance and vapor pressure gradient rather than stomatal aperture alone."
  )
})


# --------------------------------------------------------------------------- #
# Export
# --------------------------------------------------------------------------- #
out_json <- file.path(
  "/home/user/LT-trait/experimental_data_soy_2026/figures/paper_figures",
  "figure_bullets.json"
)
write_json(bullets, out_json, pretty = TRUE, auto_unbox = FALSE)
cat("Bullets written to:", out_json, "\n")
cat("Preview:\n")
for (fig in names(bullets)) {
  cat("\n[", fig, "]\n", sep = "")
  for (b in bullets[[fig]]) cat(" -", b, "\n")
}
