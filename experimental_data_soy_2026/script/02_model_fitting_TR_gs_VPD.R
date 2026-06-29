# =============================================================================
# 02_model_fitting_TR_gs_VPD.R
# Growth Chamber 2026 — Soybean Limited-Transpiration Experiment
#
# Purpose: Fit three model types (linear, bilinear, linear plateau) for
#          transpiration rate (TR) and stomatal conductance (gsw) vs VPD
#          using Day 1 data only.
#          Random effect: Rep (biological replicate).
#
# Model types:
#   (A) Linear           — TR = a + b * VPD
#   (B) Bilinear         — TR = a + b1 * VPD  for VPD <= BP (breakpoint)
#                               a + b1*BP + b2*(VPD-BP)  for VPD > BP
#   (C) Linear plateau   — TR = a + b * VPD  for VPD <= BP
#                               a + b * BP           for VPD > BP
#
# Packages:
#   lme4    — linear mixed-effects (Model A)
#   nlme    — nonlinear mixed-effects (Models B and C)
#   segmented — optional: breakpoint CI for Model B/C
#
# Run from project root (experimental_data_soy_2026/):
#   source("script/02_model_fitting_TR_gs_VPD.R")
# =============================================================================

library(readxl)
library(tidyverse)
library(lme4)       # lmer()
library(nlme)       # nlme(), lme()
library(patchwork)  # multi-panel figures

dir.create("figures", showWarnings = FALSE)


# --------------------------------------------------------------------------- #
# SECTION 1: LOAD DATA AND SUBSET TO DAY 1
# --------------------------------------------------------------------------- #

# Use the cleaned CSV if QC script has been run; otherwise reload from Excel
if (file.exists("data/data_clean.csv")) {
  df_all <- read_csv("data/data_clean.csv", show_col_types = FALSE) %>%
    mutate(Rep = factor(Rep), Day = factor(Day), Genotype = factor(Genotype))
} else {
  df_all <- read_excel("data/data_all.xlsx", sheet = "dataAll") %>%
    rename(Rep = RepNo, Day = day, TR = TR_mg_sqm_seg,
           RelTR = RelativeTransp, VPD = avg_vpd) %>%
    mutate(
      Rep      = factor(Rep),
      Day      = factor(Day, levels = c(1, 2), labels = c("Day 1", "Day 2")),
      Genotype = factor(Genotype),
      gsw      = rowMeans(cbind(gsw1, gsw2, gsw3), na.rm = FALSE)
    )
}

# ----- Day 1 only (per analysis plan) ----------------------------------------
day1 <- df_all %>%
  filter(Day == "Day 1")

cat("Day 1 dataset:\n")
cat("  Rows       :", nrow(day1), "\n")
cat("  Genotypes  :", nlevels(day1$Genotype), "\n")
cat("  Reps       :", nlevels(day1$Rep), "\n")
cat("  VPD range  :", round(min(day1$VPD), 2), "–",
                       round(max(day1$VPD), 2), "kPa\n\n")

# Stomatal conductance subset (LICOR data not always present on Day 1)
day1_gsw <- day1 %>% filter(!is.na(gsw))
cat("Day 1 subset with gsw data:", nrow(day1_gsw), "rows\n")
cat("  Genotypes with gsw:", paste(levels(droplevels(day1_gsw$Genotype)),
                                    collapse = ", "), "\n\n")


# --------------------------------------------------------------------------- #
# SECTION 2: HELPER FUNCTIONS
# --------------------------------------------------------------------------- #

# --- 2.1 Bilinear model function (two slopes, one breakpoint) ----------------
# y = a + b1*x                     for x <= bp
# y = a + b1*bp + b2*(x - bp)      for x > bp
bilinear <- function(x, a, b1, b2, bp) {
  ifelse(x <= bp,
         a + b1 * x,
         a + b1 * bp + b2 * (x - bp))
}

# --- 2.2 Linear plateau model function (slope then flat) ---------------------
# y = a + b*x    for x <= bp
# y = a + b*bp   for x > bp  (plateau value)
lin_plateau <- function(x, a, b, bp) {
  ifelse(x <= bp,
         a + b * x,
         a + b * bp)
}

# --- 2.3 Fit all three models for ONE genotype, ONE response variable --------
#   data    : data frame with columns VPD, y, Rep (factor)
#   response: label string for plots
#   geno    : genotype label for plots
fit_models_one <- function(data, response = "TR", geno = "unknown") {

  out <- list(geno = geno, response = response, n = nrow(data))

  # ---- (A) LINEAR mixed model ------------------------------------------------
  tryCatch({
    m_lin <- lme(y ~ VPD, random = ~1 | Rep, data = data,
                 method = "REML", na.action = na.omit)
    out$linear <- m_lin
  }, error = function(e) {
    message("Linear model failed for ", geno, ": ", conditionMessage(e))
    out$linear <<- NULL
  })

  # ---- (B) BILINEAR — fit via nls first to get starting values ---------------
  #   Starting values: midpoint of VPD range as breakpoint; b1 ≈ b2 ≈ slope of linear fit
  vrange <- range(data$VPD, na.rm = TRUE)
  bp0    <- mean(vrange)  # initial breakpoint guess

  # Rough OLS starting slope
  sl0 <- tryCatch(coef(lm(y ~ VPD, data = data))["VPD"], error = function(e) 1)
  a0  <- mean(data$y, na.rm = TRUE) - sl0 * mean(data$VPD, na.rm = TRUE)

  nls_bilin <- tryCatch(
    nls(y ~ bilinear(VPD, a, b1, b2, bp),
        data  = data,
        start = list(a = a0, b1 = sl0, b2 = sl0 * 0.1, bp = bp0),
        control = nls.control(maxiter = 200, tol = 1e-5)),
    error = function(e) { message("nls bilinear failed for ", geno, ": ", conditionMessage(e)); NULL }
  )

  if (!is.null(nls_bilin)) {
    sv <- as.list(coef(nls_bilin))
    tryCatch({
      m_bilin <- nlme(y ~ bilinear(VPD, a, b1, b2, bp),
                      fixed  = a + b1 + b2 + bp ~ 1,
                      random = a ~ 1 | Rep,
                      start  = c(a = sv$a, b1 = sv$b1, b2 = sv$b2, bp = sv$bp),
                      data   = data, na.action = na.omit,
                      control = nlmeControl(maxIter = 200, pnlsTol = 0.01))
      out$bilinear <- m_bilin
    }, error = function(e) {
      message("nlme bilinear failed for ", geno, "; storing nls result instead. ", conditionMessage(e))
      out$bilinear <<- nls_bilin
    })
  } else {
    out$bilinear <- NULL
  }

  # ---- (C) LINEAR PLATEAU ------------------------------------------------
  nls_lp <- tryCatch(
    nls(y ~ lin_plateau(VPD, a, b, bp),
        data  = data,
        start = list(a = a0, b = sl0, bp = bp0),
        lower = list(a = -Inf, b = 0, bp = vrange[1]),
        upper = list(a = Inf,  b = Inf, bp = vrange[2]),
        algorithm = "port",
        control = nls.control(maxiter = 200, tol = 1e-5)),
    error = function(e) { message("nls lin_plateau failed for ", geno, ": ", conditionMessage(e)); NULL }
  )

  if (!is.null(nls_lp)) {
    sv_lp <- as.list(coef(nls_lp))
    tryCatch({
      m_lp <- nlme(y ~ lin_plateau(VPD, a, b, bp),
                   fixed  = a + b + bp ~ 1,
                   random = a ~ 1 | Rep,
                   start  = c(a = sv_lp$a, b = sv_lp$b, bp = sv_lp$bp),
                   data   = data, na.action = na.omit,
                   control = nlmeControl(maxIter = 200, pnlsTol = 0.01))
      out$lin_plateau <- m_lp
    }, error = function(e) {
      message("nlme lin_plateau failed for ", geno, "; storing nls result instead. ", conditionMessage(e))
      out$lin_plateau <<- nls_lp
    })
  } else {
    out$lin_plateau <- NULL
  }

  out
}

# --- 2.4 Extract fixed-effect parameters from a fitted model -----------------
extract_params <- function(fit, model_type) {
  if (is.null(fit)) return(NULL)
  tryCatch({
    cf <- fixef(fit)  # works for both lme and nlme objects
    data.frame(model = model_type, param = names(cf), value = cf,
               row.names = NULL)
  }, error = function(e) {
    cf <- coef(fit)   # fallback for nls
    data.frame(model = model_type, param = names(cf), value = cf,
               row.names = NULL)
  })
}

# --- 2.5 Generate fitted curve for plotting ----------------------------------
fitted_curve <- function(fit, model_type, vpd_seq) {
  if (is.null(fit)) return(NULL)
  tryCatch({
    pars <- tryCatch(fixef(fit), error = function(e) coef(fit))
    if (model_type == "Linear") {
      y_hat <- pars["(Intercept)"] + pars["VPD"] * vpd_seq
    } else if (model_type == "Bilinear") {
      y_hat <- bilinear(vpd_seq, pars["a"], pars["b1"], pars["b2"], pars["bp"])
    } else if (model_type == "Lin. plateau") {
      y_hat <- lin_plateau(vpd_seq, pars["a"], pars["b"], pars["bp"])
    }
    data.frame(VPD = vpd_seq, y_hat = y_hat, model = model_type)
  }, error = function(e) NULL)
}


# --------------------------------------------------------------------------- #
# SECTION 3: FIT MODELS — TRANSPIRATION RATE (TR)
# --------------------------------------------------------------------------- #

cat("=== Fitting TR models (Day 1) ===\n")

genotypes <- levels(droplevels(day1$Genotype))

# Fit all three models for every genotype
tr_fits <- lapply(genotypes, function(g) {
  sub <- day1 %>% filter(Genotype == g) %>% select(VPD, y = TR, Rep)
  cat("  Genotype:", g, "(n =", nrow(sub), ")\n")
  fit_models_one(sub, response = "TR", geno = g)
})
names(tr_fits) <- genotypes


# --------------------------------------------------------------------------- #
# SECTION 4: FIT MODELS — STOMATAL CONDUCTANCE (gsw)
# --------------------------------------------------------------------------- #

cat("\n=== Fitting gsw models (Day 1, LICOR rows only) ===\n")

geno_gsw <- levels(droplevels(day1_gsw$Genotype))

gsw_fits <- lapply(geno_gsw, function(g) {
  sub <- day1_gsw %>% filter(Genotype == g) %>% select(VPD, y = gsw, Rep)
  cat("  Genotype:", g, "(n =", nrow(sub), ")\n")
  fit_models_one(sub, response = "gsw", geno = g)
})
names(gsw_fits) <- geno_gsw


# --------------------------------------------------------------------------- #
# SECTION 5: PARAMETER TABLES
# --------------------------------------------------------------------------- #

# --- 5.1 TR parameter table --------------------------------------------------
tr_params <- bind_rows(lapply(genotypes, function(g) {
  f <- tr_fits[[g]]
  bind_rows(
    extract_params(f$linear,      "Linear"),
    extract_params(f$bilinear,    "Bilinear"),
    extract_params(f$lin_plateau, "Lin. plateau")
  ) %>% mutate(Genotype = g)
})) %>%
  select(Genotype, model, param, value) %>%
  pivot_wider(names_from = param, values_from = value)

cat("\n--- TR model parameters (fixed effects) ---\n")
print(tr_params, n = 50)

write_csv(tr_params, "data/TR_model_parameters.csv")
cat("Saved: data/TR_model_parameters.csv\n")


# --- 5.2 gsw parameter table -------------------------------------------------
gsw_params <- bind_rows(lapply(geno_gsw, function(g) {
  f <- gsw_fits[[g]]
  bind_rows(
    extract_params(f$linear,      "Linear"),
    extract_params(f$bilinear,    "Bilinear"),
    extract_params(f$lin_plateau, "Lin. plateau")
  ) %>% mutate(Genotype = g)
})) %>%
  select(Genotype, model, param, value) %>%
  pivot_wider(names_from = param, values_from = value)

cat("\n--- gsw model parameters (fixed effects) ---\n")
print(gsw_params, n = 50)

write_csv(gsw_params, "data/gsw_model_parameters.csv")
cat("Saved: data/gsw_model_parameters.csv\n")


# --------------------------------------------------------------------------- #
# SECTION 6: AIC COMPARISON TABLE
# --------------------------------------------------------------------------- #

# Compare models by AIC (lower = better fit; valid for models fitted with ML)
# Note: lme() uses REML by default for final estimates but ML is needed for AIC
# Re-fit with method = "ML" for comparison

aic_tr <- bind_rows(lapply(genotypes, function(g) {
  f <- tr_fits[[g]]
  sub <- day1 %>% filter(Genotype == g) %>% select(VPD, y = TR, Rep)

  aic_lin <- tryCatch(AIC(update(f$linear, method = "ML")), error = function(e) NA)
  aic_bil <- tryCatch(AIC(f$bilinear), error = function(e) NA)
  aic_lp  <- tryCatch(AIC(f$lin_plateau), error = function(e) NA)

  data.frame(
    Genotype    = g,
    AIC_linear  = aic_lin,
    AIC_bilinear = aic_bil,
    AIC_linplat  = aic_lp,
    Best_model   = c("Linear", "Bilinear", "Lin. plateau")[
      which.min(c(aic_lin, aic_bil, aic_lp))]
  )
}))

cat("\n--- AIC comparison — TR models ---\n")
print(aic_tr)
write_csv(aic_tr, "data/TR_AIC_comparison.csv")


# --------------------------------------------------------------------------- #
# SECTION 7: FIGURES — TR MODELS
# --------------------------------------------------------------------------- #

vpd_seq <- seq(min(day1$VPD, na.rm = TRUE),
               max(day1$VPD, na.rm = TRUE), length.out = 200)

# WUE class colour palette
wue_cols <- c(
  "slow_wilting"        = "#1b7837",
  "high_wilting"        = "#d73027",
  "Plastic, High WUE"   = "#4393c3",
  "Nonplastic, Low WUE" = "#f4a582",
  "MAGIC parent"        = "#9970ab"
)

# --- Fig M1: TR vs VPD — all genotypes, Day 1 --------------------------------
p_tr_all <- ggplot(day1, aes(x = VPD, y = TR, colour = WUE_class, shape = Rep)) +
  geom_point(alpha = 0.7, size = 1.8) +
  facet_wrap(~Genotype, ncol = 5) +
  scale_colour_manual(values = wue_cols, name = "WUE class") +
  labs(
    title = "Fig M1: Transpiration rate vs VPD — Day 1 (all genotypes)",
    x     = "VPD (kPa)",
    y     = expression(TR~(mg~m^{-2}~s^{-1}))
  ) +
  theme_bw(base_size = 9) +
  theme(strip.text = element_text(size = 7), legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 2))

ggsave("figures/M1_TR_vs_VPD_Day1_allgeno.png",
       p_tr_all, width = 14, height = 8, dpi = 300)


# --- Fig M2: TR models overlaid per genotype ---------------------------------
# Build fitted curve data for all genotypes
curve_data_tr <- bind_rows(lapply(genotypes, function(g) {
  f <- tr_fits[[g]]
  bind_rows(
    fitted_curve(f$linear,      "Linear",      vpd_seq),
    fitted_curve(f$bilinear,    "Bilinear",    vpd_seq),
    fitted_curve(f$lin_plateau, "Lin. plateau", vpd_seq)
  ) %>% mutate(Genotype = g)
}))

model_cols <- c("Linear" = "#636363", "Bilinear" = "#2166ac", "Lin. plateau" = "#d95f02")
model_lty  <- c("Linear" = "dashed",  "Bilinear" = "solid",   "Lin. plateau" = "solid")

p_tr_models <- ggplot() +
  geom_point(data = day1,
             aes(x = VPD, y = TR, colour = WUE_class),
             alpha = 0.4, size = 1.2, show.legend = TRUE) +
  geom_line(data = curve_data_tr,
            aes(x = VPD, y = y_hat, linetype = model, colour = model),
            linewidth = 0.9, inherit.aes = FALSE) +
  scale_colour_manual(
    values = c(wue_cols, model_cols),
    breaks = names(model_cols),
    name   = "Model"
  ) +
  scale_linetype_manual(values = model_lty, name = "Model") +
  facet_wrap(~Genotype, ncol = 5, scales = "free_y") +
  labs(
    title = "Fig M2: TR models fitted to Day 1 data",
    x     = "VPD (kPa)",
    y     = expression(TR~(mg~m^{-2}~s^{-1}))
  ) +
  theme_bw(base_size = 9) +
  theme(strip.text = element_text(size = 7), legend.position = "bottom")

ggsave("figures/M2_TR_models_per_genotype.png",
       p_tr_models, width = 14, height = 8, dpi = 300)


# --- Fig M3: Breakpoint estimates — TR models --------------------------------
bp_tr <- tr_params %>%
  select(Genotype, model, bp) %>%
  filter(!is.na(bp), model %in% c("Bilinear", "Lin. plateau")) %>%
  left_join(day1 %>% distinct(Genotype, WUE_class), by = "Genotype")

if (nrow(bp_tr) > 0) {
  p_bp_tr <- ggplot(bp_tr, aes(x = Genotype, y = bp, colour = WUE_class,
                                shape = model)) +
    geom_hline(yintercept = 1.73, linetype = "dashed", colour = "grey40") +
    geom_point(size = 3, position = position_dodge(0.4)) +
    annotate("text", x = 0.5, y = 1.73, label = "BP = 1.73 kPa (lit.)",
             hjust = 0, vjust = -0.3, size = 3, colour = "grey30") +
    scale_colour_manual(values = wue_cols, name = "WUE class") +
    scale_shape_manual(values = c("Bilinear" = 16, "Lin. plateau" = 17),
                       name = "Model") +
    labs(
      title = "Fig M3: Breakpoint (VPD, kPa) from bilinear and linear plateau models — TR",
      x     = NULL,
      y     = "Breakpoint VPD (kPa)"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right")

  ggsave("figures/M3_TR_breakpoints.png",
         p_bp_tr, width = 10, height = 5, dpi = 300)
}


# --- Fig M4: Slopes (linear segment) — TR model B and C ----------------------
slope_tr <- tr_params %>%
  mutate(
    slope = case_when(
      model == "Linear"       ~ `VPD`,
      model == "Bilinear"     ~ b1,
      model == "Lin. plateau" ~ b
    )
  ) %>%
  filter(!is.na(slope)) %>%
  left_join(day1 %>% distinct(Genotype, WUE_class), by = "Genotype")

p_slope_tr <- ggplot(slope_tr, aes(x = reorder(Genotype, slope),
                                    y = slope, fill = WUE_class)) +
  geom_col() +
  facet_wrap(~model, ncol = 3) +
  scale_fill_manual(values = wue_cols, name = "WUE class") +
  coord_flip() +
  labs(
    title = "Fig M4: Transpiration sensitivity slope (mg m-2 s-1 per kPa) — Day 1",
    x     = NULL,
    y     = expression(Slope~(mg~m^{-2}~s^{-1}~kPa^{-1}))
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom")

ggsave("figures/M4_TR_slopes.png",
       p_slope_tr, width = 10, height = 6, dpi = 300)


# --------------------------------------------------------------------------- #
# SECTION 8: FIGURES — gsw MODELS (if data available)
# --------------------------------------------------------------------------- #

if (nrow(day1_gsw) > 5) {

  vpd_seq_gsw <- seq(min(day1_gsw$VPD, na.rm = TRUE),
                     max(day1_gsw$VPD, na.rm = TRUE), length.out = 200)

  # --- Fig M5: gsw vs VPD, Day 1 --------------------------------------------
  p_gsw_all <- ggplot(day1_gsw, aes(x = VPD, y = gsw, colour = WUE_class)) +
    geom_point(alpha = 0.8, size = 2) +
    facet_wrap(~Genotype, ncol = 5) +
    scale_colour_manual(values = wue_cols, name = "WUE class") +
    labs(
      title = "Fig M5: Stomatal conductance vs VPD — Day 1",
      x     = "VPD (kPa)",
      y     = expression(g[sw]~(mol~m^{-2}~s^{-1}))
    ) +
    theme_bw(base_size = 9) +
    theme(strip.text = element_text(size = 7), legend.position = "bottom")

  ggsave("figures/M5_gsw_vs_VPD_Day1_allgeno.png",
         p_gsw_all, width = 14, height = 8, dpi = 300)


  # --- Fig M6: gsw models per genotype (only if fits converged) -------------
  curve_data_gsw <- bind_rows(lapply(geno_gsw, function(g) {
    f <- gsw_fits[[g]]
    bind_rows(
      fitted_curve(f$linear,      "Linear",      vpd_seq_gsw),
      fitted_curve(f$bilinear,    "Bilinear",    vpd_seq_gsw),
      fitted_curve(f$lin_plateau, "Lin. plateau", vpd_seq_gsw)
    ) %>% mutate(Genotype = g)
  }))

  if (nrow(curve_data_gsw) > 0) {
    p_gsw_models <- ggplot() +
      geom_point(data = day1_gsw,
                 aes(x = VPD, y = gsw, colour = WUE_class),
                 alpha = 0.6, size = 2, show.legend = TRUE) +
      geom_line(data = curve_data_gsw,
                aes(x = VPD, y = y_hat, linetype = model, colour = model),
                linewidth = 1, inherit.aes = FALSE) +
      scale_colour_manual(
        values = c(wue_cols, model_cols),
        breaks = names(model_cols), name = "Model"
      ) +
      scale_linetype_manual(values = model_lty, name = "Model") +
      facet_wrap(~Genotype, ncol = 4, scales = "free_y") +
      labs(
        title = "Fig M6: gsw models fitted to Day 1 LICOR data",
        x     = "VPD (kPa)",
        y     = expression(g[sw]~(mol~m^{-2}~s^{-1}))
      ) +
      theme_bw(base_size = 9) +
      theme(strip.text = element_text(size = 7), legend.position = "bottom")

    ggsave("figures/M6_gsw_models_per_genotype.png",
           p_gsw_models, width = 12, height = 8, dpi = 300)
  }

  # --- Fig M7: Breakpoints — gsw models -------------------------------------
  bp_gsw <- gsw_params %>%
    select(Genotype, model, bp) %>%
    filter(!is.na(bp), model %in% c("Bilinear", "Lin. plateau")) %>%
    left_join(day1_gsw %>% distinct(Genotype, WUE_class), by = "Genotype")

  if (nrow(bp_gsw) > 0) {
    p_bp_gsw <- ggplot(bp_gsw, aes(x = Genotype, y = bp, colour = WUE_class,
                                    shape = model)) +
      geom_hline(yintercept = 1.73, linetype = "dashed", colour = "grey40") +
      geom_point(size = 3, position = position_dodge(0.4)) +
      scale_colour_manual(values = wue_cols, name = "WUE class") +
      scale_shape_manual(values = c("Bilinear" = 16, "Lin. plateau" = 17),
                         name = "Model") +
      labs(
        title = "Fig M7: Breakpoint (VPD, kPa) — gsw models",
        x     = NULL,
        y     = "Breakpoint VPD (kPa)"
      ) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggsave("figures/M7_gsw_breakpoints.png",
           p_bp_gsw, width = 8, height = 5, dpi = 300)
  }

} else {
  cat("\nSkipping gsw model figures: too few Day 1 rows with gsw data.\n")
}


# --------------------------------------------------------------------------- #
# SECTION 9: POOLED ANALYSIS BY WUE CLASS
# --------------------------------------------------------------------------- #

cat("\n=== Pooled analysis by WUE class ===\n")

# Linear mixed model: TR ~ VPD * WUE_class + (1|Rep) + (1|Genotype)
# Genotype is nested within WUE_class; both Rep and Genotype are random
pool_lin <- lme(TR ~ VPD * WUE_class, data = day1,
                random = list(Rep = ~1, Genotype = ~1),
                method = "REML", na.action = na.omit)

cat("\n--- Pooled linear model: TR ~ VPD * WUE_class ---\n")
print(summary(pool_lin))

# Pooled linear model: gsw
if (nrow(day1_gsw) > 10) {
  pool_lin_gsw <- lme(gsw ~ VPD * WUE_class, data = day1_gsw,
                      random = list(Rep = ~1, Genotype = ~1),
                      method = "REML", na.action = na.omit)
  cat("\n--- Pooled linear model: gsw ~ VPD * WUE_class ---\n")
  print(summary(pool_lin_gsw))
}


# --- Fig M8: Pooled TR vs VPD by WUE class ----------------------------------
p_pool <- ggplot(day1, aes(x = VPD, y = TR, colour = WUE_class)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
  scale_colour_manual(values = wue_cols, name = "WUE class") +
  labs(
    title = "Fig M8: TR vs VPD by WUE class — Day 1 (pooled across genotypes)",
    x     = "VPD (kPa)",
    y     = expression(TR~(mg~m^{-2}~s^{-1}))
  ) +
  theme_bw() +
  theme(legend.position = "right")

ggsave("figures/M8_TR_vs_VPD_by_WUEclass.png",
       p_pool, width = 8, height = 5, dpi = 300)

cat("\n=== Done ===\n")
cat("Figures saved to figures/\n")
cat("Parameter tables saved to data/\n")
