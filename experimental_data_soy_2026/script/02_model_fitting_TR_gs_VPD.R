# =============================================================================
# 02_model_fitting_TR_gs_VPD.R
# Growth Chamber 2026 — Soybean Limited-Transpiration Experiment
#
# Purpose : Fit three model types for TR and gsw vs VPD (Day 1 only):
#   (A) Linear           TR = a + b·VPD
#   (B) Bilinear         two slopes with a breakpoint (BP)
#   (C) Linear plateau   linear then flat at TR_max (classic LT response)
#
#   Random effect : Rep (replicate) via nlme::lme() / nlme::nlme()
#   Model selection : AIC (lower = better)
#
# Outputs (all written to figures/ and data/):
#   Fig1  TR models per genotype (panel grid, publication layout)
#   Fig2  TR by WUE class (pooled)
#   Fig3  Breakpoint bar chart
#   Fig4  Slope₁ ranking
#   Fig5  gsw models per genotype
#   Fig6  AIC heatmap
#   data/TR_best_model_summary.csv
#   data/gsw_best_model_summary.csv
#   data/TR_model_parameters.csv
#   data/gsw_model_parameters.csv
#
# Run from project root  (experimental_data_soy_2026/):
#   source("script/02_model_fitting_TR_gs_VPD.R")
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
  library(nlme)
  library(patchwork)
  library(scales)
})

dir.create("figures", showWarnings = FALSE)
dir.create("data",    showWarnings = FALSE)

# --------------------------------------------------------------------------- #
# COLOUR SYSTEM — University of Arkansas
# --------------------------------------------------------------------------- #
UA_RED   <- "#9D2235"
UA_WINE  <- "#6E1424"
UA_RED2  <- "#C8102E"
UA_BLUE  <- "#4393c3"
UA_ORG   <- "#d4824a"
UA_DARK  <- "#1A1A1A"
UA_GRAY  <- "#7A7A7A"

WUE_COLS <- c(
  "slow_wilting"        = UA_WINE,
  "high_wilting"        = UA_RED2,
  "Plastic, High WUE"   = UA_BLUE,
  "Nonplastic, Low WUE" = UA_ORG,
  "MAGIC parent"        = UA_DARK
)
WUE_SHAPES <- c(
  "slow_wilting"        = 21,
  "high_wilting"        = 22,
  "Plastic, High WUE"   = 24,
  "Nonplastic, Low WUE" = 23,
  "MAGIC parent"        = 25
)
MODEL_COLS <- c(
  "Linear"       = UA_GRAY,
  "Bilinear"     = UA_RED,
  "Lin. plateau" = UA_WINE
)
MODEL_LTY <- c(
  "Linear"       = "dashed",
  "Bilinear"     = "solid",
  "Lin. plateau" = "dotdash"
)

UA_THEME <- theme_bw(base_size = 9) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(linewidth = 0.25, colour = "grey90"),
    strip.background  = element_rect(fill = UA_RED, colour = NA),
    strip.text        = element_text(colour = "white", face = "bold", size = 7),
    axis.title        = element_text(size = 8),
    plot.title        = element_text(face = "bold", size = 9, colour = UA_DARK),
    legend.key.size   = unit(0.4, "cm"),
    legend.text       = element_text(size = 7.5),
    legend.title      = element_text(size = 8, face = "bold")
  )

WUE_ORDER <- c("slow_wilting", "high_wilting",
                "Plastic, High WUE", "Nonplastic, Low WUE", "MAGIC parent")

WUE_EXPECTED <- c(
  "slow_wilting"        = "Slow wilting / LT",
  "high_wilting"        = "High wilting / non-LT",
  "Plastic, High WUE"   = "Plastic – High WUE",
  "Nonplastic, Low WUE" = "Nonplastic – Low WUE",
  "MAGIC parent"        = "MAGIC parent"
)


# --------------------------------------------------------------------------- #
# 1. LOAD DATA — Day 1 only
# --------------------------------------------------------------------------- #
if (file.exists("data/data_clean.csv")) {
  df_all <- read_csv("data/data_clean.csv", show_col_types = FALSE) |>
    mutate(Rep = factor(Rep), Genotype = factor(Genotype),
           WUE_class = factor(WUE_class, levels = WUE_ORDER))
} else {
  df_all <- read_excel("data/data_all.xlsx", sheet = "dataAll") |>
    rename(Rep = RepNo, Day = day, TR = TR_mg_sqm_seg,
           RelTR = RelativeTransp, VPD = avg_vpd) |>
    mutate(
      Rep       = factor(Rep),
      Day_label = factor(ifelse(Day == 1, "Day 1", "Day 2"),
                         levels = c("Day 1","Day 2")),
      Genotype  = factor(Genotype),
      WUE_class = factor(WUE_class, levels = WUE_ORDER),
      gsw       = rowMeans(cbind(gsw1, gsw2, gsw3), na.rm = FALSE)
    )
}

day1     <- filter(df_all, Day_label == "Day 1")
day1_gsw <- filter(day1, !is.na(gsw))

genotypes <- levels(droplevels(day1$Genotype))
geno_gsw  <- levels(droplevels(day1_gsw$Genotype))

# Individual LICOR readings (gsw1/gsw2/gsw3) in long format for plotting
# Each row in day1_gsw yields 3 plotable points showing within-plant variability
day1_gsw_long <- day1_gsw |>
  select(Genotype, WUE_class, VPD, Rep, gsw1, gsw2, gsw3) |>
  pivot_longer(cols = c(gsw1, gsw2, gsw3),
               names_to  = "reading",
               values_to = "gsw_reading") |>
  filter(!is.na(gsw_reading))

# Lookup tables
wue_lookup <- day1 |>
  distinct(Genotype, WUE_class) |>
  deframe()   # named vector: Genotype -> WUE_class

cat("Day 1:", nrow(day1), "rows |", length(genotypes), "genotypes |",
    nrow(day1_gsw), "rows with gsw |",
    nrow(day1_gsw_long), "individual LICOR readings\n")


# --------------------------------------------------------------------------- #
# 2. MODEL FUNCTIONS
# --------------------------------------------------------------------------- #

# Bilinear: a + b1*x for x<=bp; a + b1*bp + b2*(x-bp) for x>bp
bilinear_fn <- function(x, a, b1, b2, bp) {
  ifelse(x <= bp, a + b1 * x, a + b1 * bp + b2 * (x - bp))
}

# Linear plateau: a + b*x for x<=bp; a + b*bp for x>bp
lin_plateau_fn <- function(x, a, b, bp) {
  ifelse(x <= bp, a + b * x, a + b * bp)
}


# --------------------------------------------------------------------------- #
# 3. FIT FUNCTION (per genotype)
# --------------------------------------------------------------------------- #
# Returns a list with elements: linear, bilinear, lin_plateau, best (name)
# Each element is the fitted model object (lme or nls fallback).

fit_one <- function(data, response_col = "TR", geno_label = "?") {
  y   <- data[[response_col]]
  vpd <- data$VPD
  rep <- data$Rep
  n   <- sum(!is.na(y) & !is.na(vpd))
  if (n < 4) return(NULL)

  vmin <- min(vpd, na.rm = TRUE)
  vmax <- max(vpd, na.rm = TRUE)
  bp0  <- (vmin + vmax) / 2
  # rough starting slope from OLS
  sl0  <- coef(lm(y ~ vpd, na.action = na.omit))["vpd"]
  if (is.na(sl0)) sl0 <- 1
  a0   <- mean(y, na.rm = TRUE) - sl0 * mean(vpd, na.rm = TRUE)

  df_fit <- data.frame(y = y, VPD = vpd, Rep = rep) |> na.omit()

  out <- list()

  # Keys use display names so helpers can look up by model name directly
  # ---- (A) LINEAR mixed model -----------------------------------------------
  out[["Linear"]] <- tryCatch(
    lme(y ~ VPD, random = ~1 | Rep, data = df_fit, method = "ML",
        control = lmeControl(opt = "optim")),
    error = function(e) {
      tryCatch(lm(y ~ VPD, data = df_fit), error = function(e2) NULL)
    }
  )

  # ---- (B) BILINEAR — nls starting values, then nlme ----------------------
  nls_bil <- tryCatch(
    nls(y ~ bilinear_fn(VPD, a, b1, b2, bp),
        data    = df_fit,
        start   = list(a = a0, b1 = sl0, b2 = sl0 * 0.1, bp = bp0),
        lower   = c(a = -Inf, b1 = 0,  b2 = -Inf, bp = vmin),
        upper   = c(a =  Inf, b1 = Inf, b2 =  Inf, bp = vmax),
        algorithm = "port",
        control = nls.control(maxiter = 300, tol = 1e-5)),
    error = function(e) NULL
  )

  if (!is.null(nls_bil)) {
    sv <- as.list(coef(nls_bil))
    out[["Bilinear"]] <- tryCatch(
      nlme(y ~ bilinear_fn(VPD, a, b1, b2, bp),
           fixed   = a + b1 + b2 + bp ~ 1,
           random  = a ~ 1 | Rep,
           start   = c(a = sv$a, b1 = sv$b1, b2 = sv$b2, bp = sv$bp),
           data    = df_fit,
           control = nlmeControl(maxIter = 200, pnlsTol = 0.05, opt = "nlminb")),
      error = function(e) nls_bil
    )
  }

  # ---- (C) LINEAR PLATEAU ---------------------------------------------------
  nls_lp <- tryCatch(
    nls(y ~ lin_plateau_fn(VPD, a, b, bp),
        data    = df_fit,
        start   = list(a = a0, b = sl0, bp = bp0),
        lower   = c(a = -Inf, b = 0,   bp = vmin),
        upper   = c(a =  Inf, b = Inf, bp = vmax),
        algorithm = "port",
        control = nls.control(maxiter = 300, tol = 1e-5)),
    error = function(e) NULL
  )

  if (!is.null(nls_lp)) {
    sv_lp <- as.list(coef(nls_lp))
    out[["Lin. plateau"]] <- tryCatch(
      nlme(y ~ lin_plateau_fn(VPD, a, b, bp),
           fixed   = a + b + bp ~ 1,
           random  = a ~ 1 | Rep,
           start   = c(a = sv_lp$a, b = sv_lp$b, bp = sv_lp$bp),
           data    = df_fit,
           control = nlmeControl(maxIter = 200, pnlsTol = 0.05, opt = "nlminb")),
      error = function(e) nls_lp
    )
  }

  # ---- Select best model by AIC --------------------------------------------
  get_aic <- function(m) {
    if (is.null(m)) return(NA_real_)
    tryCatch(AIC(m), error = function(e) NA_real_)
  }

  aics <- c(
    "Linear"       = get_aic(out[["Linear"]]),
    "Bilinear"     = get_aic(out[["Bilinear"]]),
    "Lin. plateau" = get_aic(out[["Lin. plateau"]])
  )
  out$AIC  <- aics
  out$best <- names(which.min(aics))
  out
}


# Null-coalescing helper (must be defined before use)
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a[[1]])) a else b

# --------------------------------------------------------------------------- #
# 4. EXTRACT PARAMETERS
# --------------------------------------------------------------------------- #
get_params <- function(fit_list, model_name) {
  m <- fit_list[[model_name]]
  if (is.null(m)) return(NULL)
  tryCatch(
    {
      cf <- if (inherits(m, "lme") || inherits(m, "nlme")) fixef(m) else coef(m)
      as.list(cf)
    },
    error = function(e) NULL
  )
}

# Build one-row summary for the best model
summarise_fit <- function(g, fit_list) {
  if (is.null(fit_list)) return(NULL)
  best <- fit_list$best
  if (is.null(best)) return(NULL)
  pars <- get_params(fit_list, best)
  if (is.null(pars)) return(NULL)
  wc <- as.character(wue_lookup[g])

  # Slope1: first slope; BP; Slope2 (bilinear only)
  slope1 <- switch(best,
    "Linear"       = pars[["VPD"]],
    "Bilinear"     = pars[["b1"]],
    "Lin. plateau" = pars[["b"]],
    NA_real_)
  bp     <- switch(best,
    "Bilinear"     = pars[["bp"]],
    "Lin. plateau" = pars[["bp"]],
    NA_real_)
  slope2 <- if (best == "Bilinear") pars[["b2"]] else NA_real_

  tibble(
    Genotype          = g,
    WUE_class         = wc,
    Expected_behavior = WUE_EXPECTED[wc],
    Best_model        = best,
    Slope1            = round(slope1, 4),
    BP_kPa            = round(bp, 3),
    Slope2            = round(slope2, 4),
    AIC_best          = round(fit_list$AIC[best], 1)
  )
}

# Build full parameter table (all three models)
full_param_table <- function(g, fit_list) {
  if (is.null(fit_list)) return(NULL)
  wc   <- as.character(wue_lookup[g])
  best <- fit_list$best
  aics <- fit_list$AIC

  bind_rows(lapply(c("Linear","Bilinear","Lin. plateau"), function(mn) {
    pars <- get_params(fit_list, mn)
    if (is.null(pars)) return(NULL)
    tibble(
      Genotype   = g,
      WUE_class  = wc,
      Model      = mn,
      Best_model = if (!is.null(best) && mn == best) "yes" else "",
      a          = round(pars[["a"]] %||% pars[["(Intercept)"]] %||% NA_real_, 3),
      Slope1     = round(switch(mn,
        "Linear"       = pars[["VPD"]],
        "Bilinear"     = pars[["b1"]],
        "Lin. plateau" = pars[["b"]], NA_real_), 4),
      Slope2     = round(if (mn == "Bilinear") pars[["b2"]] else NA_real_, 4),
      BP_kPa     = round(if (mn %in% c("Bilinear","Lin. plateau")) pars[["bp"]]
                         else NA_real_, 3),
      AIC        = round(aics[mn], 1),
      Delta_AIC  = round(aics[mn] - min(aics, na.rm = TRUE), 1)
    )
  }))
}


# --------------------------------------------------------------------------- #
# 5. CURVE PREDICTION HELPERS
# --------------------------------------------------------------------------- #
pred_curve <- function(fit_list, model_name, vpd_seq) {
  m <- fit_list[[model_name]]
  if (is.null(m)) return(NULL)
  pars <- get_params(fit_list, model_name)
  if (is.null(pars)) return(NULL)
  y_hat <- switch(model_name,
    "Linear"       = pars[["(Intercept)"]] + pars[["VPD"]] * vpd_seq,
    "Bilinear"     = bilinear_fn(vpd_seq, pars$a, pars$b1, pars$b2, pars$bp),
    "Lin. plateau" = lin_plateau_fn(vpd_seq, pars$a, pars$b, pars$bp),
    stop("Unknown model")
  )
  tibble(VPD = vpd_seq, y_hat = y_hat, Model = model_name)
}

# --------------------------------------------------------------------------- #
# 6. FIT ALL GENOTYPES
# --------------------------------------------------------------------------- #
cat("\nFitting TR models (Day 1) ...\n")
tr_fits <- setNames(
  lapply(genotypes, function(g) {
    sub <- filter(day1, Genotype == g) |> select(VPD, TR, Rep)
    fit <- fit_one(sub, response_col = "TR", geno_label = g)
    cat(sprintf("  %-15s best = %s\n", g,
                if (!is.null(fit)) fit$best else "FAILED"))
    fit
  }),
  genotypes
)

cat("\nFitting gsw models (Day 1, LICOR only) ...\n")
gsw_fits <- setNames(
  lapply(geno_gsw, function(g) {
    sub <- filter(day1_gsw, Genotype == g) |> select(VPD, gsw, Rep) |>
           rename(TR = gsw)
    fit <- fit_one(sub, response_col = "TR", geno_label = g)
    cat(sprintf("  %-15s best = %s\n", g,
                if (!is.null(fit)) fit$best else "FAILED"))
    fit
  }),
  geno_gsw
)


# --------------------------------------------------------------------------- #
# 7. PARAMETER TABLES
# --------------------------------------------------------------------------- #

# Best-model summary
tr_best <- bind_rows(lapply(genotypes, function(g) summarise_fit(g, tr_fits[[g]])))
gsw_best <- bind_rows(lapply(geno_gsw, function(g) summarise_fit(g, gsw_fits[[g]])))

# Full parameter table (all models)
tr_params  <- bind_rows(lapply(genotypes, function(g) full_param_table(g, tr_fits[[g]])))
gsw_params <- bind_rows(lapply(geno_gsw,  function(g) full_param_table(g, gsw_fits[[g]])))

# Rename gsw column in gsw fits (we used TR as proxy column name above)
# Parameter names are identical, column interpretation differs

write_csv(tr_best,   "data/TR_best_model_summary.csv")
write_csv(gsw_best,  "data/gsw_best_model_summary.csv")
write_csv(tr_params, "data/TR_model_parameters.csv")
write_csv(gsw_params,"data/gsw_model_parameters.csv")

cat("\n--- TR best-model summary ---\n")
print(tr_best, n = 20)
cat("\n--- gsw best-model summary ---\n")
print(gsw_best, n = 20)


# --------------------------------------------------------------------------- #
# 8. FIGURES
# --------------------------------------------------------------------------- #
VPD_SEQ <- seq(min(day1$VPD, na.rm = TRUE) * 0.95,
               max(day1$VPD, na.rm = TRUE) * 1.05, length.out = 300)

# Fixed y-axis limits shared across TR and gsw panels for fair comparison
TR_YLIM  <- c(0, max(day1$TR, na.rm = TRUE) * 1.05)
GSW_YLIM <- c(0, max(day1_gsw_long$gsw_reading, na.rm = TRUE) * 1.10)

# Shared legend elements
wue_legend <- scale_colour_manual(values = WUE_COLS, name = "WUE class",
                                   limits = names(WUE_COLS))
wue_fill   <- scale_fill_manual(values = WUE_COLS, name = "WUE class",
                                 limits = names(WUE_COLS))
mod_legend <- scale_colour_manual(values = MODEL_COLS, name = "Model",
                                   limits = names(MODEL_COLS))
mod_lty    <- scale_linetype_manual(values = MODEL_LTY, name = "Model",
                                    limits = names(MODEL_LTY))


# ---- FIG 1: TR per genotype with all three model curves -------------------
cat("Building Fig1 ...\n")

# Build predicted curve data frame for all genotypes
curves_tr <- bind_rows(lapply(genotypes, function(g) {
  fit <- tr_fits[[g]]
  if (is.null(fit)) return(NULL)
  bind_rows(lapply(c("Linear","Bilinear","Lin. plateau"), function(mn) {
    tryCatch(pred_curve(fit, mn, VPD_SEQ) |> mutate(Genotype = g),
             error = function(e) NULL)
  }))
}))

# Breakpoint data frame
bp_points_tr <- tr_best |>
  filter(!is.na(BP_kPa)) |>
  rowwise() |>
  mutate(y_bp = {
    fit <- tr_fits[[Genotype]]
    pars <- get_params(fit, Best_model)
    if (Best_model == "Bilinear")
      bilinear_fn(BP_kPa, pars$a, pars$b1, pars$b2, pars$bp)
    else
      lin_plateau_fn(BP_kPa, pars$a, pars$b, pars$bp)
  }) |>
  ungroup()

p1 <- ggplot() +
  # Data points (alpha varies with Rep for depth)
  geom_point(data = day1 |>
               mutate(Genotype = factor(Genotype, levels = genotypes)),
             aes(x = VPD, y = TR, fill = WUE_class, shape = WUE_class),
             colour = "white", size = 2.0, alpha = 0.75, stroke = 0.3) +
  scale_fill_manual(values = WUE_COLS, name = "WUE class") +
  scale_shape_manual(values = WUE_SHAPES, name = "WUE class") +
  # Model curves
  geom_line(data = curves_tr,
            aes(x = VPD, y = y_hat, colour = Model, linetype = Model,
                linewidth = Model),
            show.legend = TRUE) +
  mod_legend + mod_lty +
  scale_linewidth_manual(values = c("Linear" = 0.7, "Bilinear" = 1.4,
                                     "Lin. plateau" = 1.4),
                          name = "Model") +
  # Breakpoint markers
  geom_vline(data = bp_points_tr,
             aes(xintercept = BP_kPa),
             linetype = "dotted", colour = UA_DARK, linewidth = 0.5, alpha = 0.5) +
  geom_point(data = bp_points_tr,
             aes(x = BP_kPa, y = y_bp),
             shape = 25, size = 2.5, fill = UA_DARK, colour = "white",
             stroke = 0.4, inherit.aes = FALSE) +
  facet_wrap(~Genotype, ncol = 5, scales = "fixed") +
  coord_cartesian(ylim = TR_YLIM) +
  labs(
    title = "Fig 1  Transpiration rate vs VPD — Day 1  (linear, bilinear, linear-plateau models)",
    x = "VPD (kPa)",
    y = expression(TR~(mg~m^{-2}~s^{-1}))
  ) +
  UA_THEME +
  theme(legend.position = "bottom") +
  guides(
    fill      = guide_legend(nrow = 2, override.aes = list(size = 3, stroke = 0)),
    shape     = guide_legend(nrow = 2),
    colour    = guide_legend(nrow = 1, override.aes = list(linewidth = 1.2)),
    linetype  = guide_legend(nrow = 1),
    linewidth = guide_legend(nrow = 1)
  )

ggsave("figures/Fig1_TR_models_per_genotype.png", p1,
       width = 14, height = 9, dpi = 300)
cat("Saved Fig1\n")


# ---- FIG 2: TR by WUE class (pooled) ----------------------------------------
cat("Building Fig2 ...\n")

# Per-genotype best-curve data for thin background lines
curves_best_tr <- bind_rows(lapply(genotypes, function(g) {
  fit <- tr_fits[[g]]
  if (is.null(fit) || is.null(fit$best)) return(NULL)
  tryCatch(pred_curve(fit, fit$best, VPD_SEQ) |> mutate(Genotype = g),
           error = function(e) NULL)
})) |>
  left_join(day1 |> distinct(Genotype, WUE_class), by = "Genotype")

p2 <- ggplot() +
  # Thin per-genotype curves
  geom_line(data = curves_best_tr,
            aes(x = VPD, y = y_hat, group = Genotype, colour = WUE_class),
            linewidth = 0.6, alpha = 0.35) +
  # Data scatter
  geom_point(data = day1,
             aes(x = VPD, y = TR, colour = WUE_class, shape = WUE_class),
             size = 1.6, alpha = 0.50) +
  # Pooled OLS line per WUE class
  geom_smooth(data = day1,
              aes(x = VPD, y = TR, colour = WUE_class),
              method = "lm", se = TRUE, linewidth = 1.8,
              alpha = 0.15, formula = y ~ x) +
  wue_legend +
  scale_shape_manual(values = WUE_SHAPES, name = "WUE class") +
  facet_wrap(~WUE_class, ncol = 5, labeller = label_wrap_gen(15), scales = "fixed") +
  coord_cartesian(ylim = TR_YLIM) +
  labs(
    title = "Fig 2  TR vs VPD by WUE class — Day 1\n(thin = per-genotype best model; thick = pooled OLS +/- 95% CI)",
    x = "VPD (kPa)",
    y = expression(TR~(mg~m^{-2}~s^{-1}))
  ) +
  UA_THEME +
  theme(legend.position = "none",
        strip.text = element_text(size = 8))

ggsave("figures/Fig2_TR_by_WUEclass.png", p2,
       width = 14, height = 4, dpi = 300)
cat("Saved Fig2\n")


# ---- FIG 3: Breakpoints -------------------------------------------------------
cat("Building Fig3 ...\n")
bp_df <- tr_best |>
  filter(!is.na(BP_kPa)) |>
  mutate(Genotype = fct_reorder(Genotype, BP_kPa))

p3 <- ggplot(bp_df, aes(y = Genotype, x = BP_kPa, fill = WUE_class)) +
  geom_col(width = 0.65, colour = "white", linewidth = 0.3) +
  geom_vline(xintercept = 1.73, linetype = "dashed",
             colour = UA_DARK, linewidth = 0.7) +
  annotate("text", x = 1.73 + 0.04, y = 0.4, label = "lit. BP = 1.73 kPa",
           hjust = 0, size = 2.9, colour = UA_DARK) +
  wue_fill +
  labs(
    title = "Fig 3  Estimated breakpoint VPD (bilinear or linear-plateau models) — TR",
    x = "Breakpoint VPD (kPa)", y = NULL
  ) +
  UA_THEME +
  theme(legend.position = "right")

ggsave("figures/Fig3_TR_breakpoints.png", p3,
       width = 8, height = 4, dpi = 300)
cat("Saved Fig3\n")


# ---- FIG 4: Slope₁ ranking ---------------------------------------------------
cat("Building Fig4 ...\n")
slope_df <- tr_best |>
  filter(!is.na(Slope1)) |>
  mutate(Genotype = fct_reorder(Genotype, Slope1))

p4 <- ggplot(slope_df, aes(y = Genotype, x = Slope1, fill = WUE_class)) +
  geom_col(width = 0.65, colour = "white", linewidth = 0.3) +
  wue_fill +
  labs(
    title = "Fig 4  Transpiration sensitivity (Slope₁) by genotype — Day 1",
    x     = expression(Slope[1]~(mg~m^{-2}~s^{-1}~kPa^{-1})),
    y     = NULL
  ) +
  UA_THEME +
  theme(legend.position = "right")

ggsave("figures/Fig4_TR_slopes.png", p4,
       width = 8, height = 5, dpi = 300)
cat("Saved Fig4\n")


# ---- FIG 5: gsw per genotype ------------------------------------------------
if (length(geno_gsw) > 0) {
  cat("Building Fig5 ...\n")
  VPD_SEQ_GSW <- seq(min(day1_gsw$VPD, na.rm = TRUE) * 0.95,
                     max(day1_gsw$VPD, na.rm = TRUE) * 1.05, length.out = 300)

  curves_gsw <- bind_rows(lapply(geno_gsw, function(g) {
    fit <- gsw_fits[[g]]
    if (is.null(fit)) return(NULL)
    bind_rows(lapply(c("Linear","Bilinear","Lin. plateau"), function(mn) {
      tryCatch(pred_curve(fit, mn, VPD_SEQ_GSW) |> mutate(Genotype = g),
               error = function(e) NULL)
    }))
  }))

  bp_gsw_pts <- gsw_best |>
    filter(!is.na(BP_kPa))

  if (nrow(bp_gsw_pts) > 0) {
    bp_gsw_pts <- bp_gsw_pts |>
      rowwise() |>
      mutate(y_bp = {
        fit  <- gsw_fits[[Genotype]]
        pars <- get_params(fit, Best_model)
        if (Best_model == "Bilinear")
          bilinear_fn(BP_kPa, pars$a, pars$b1, pars$b2, pars$bp)
        else
          lin_plateau_fn(BP_kPa, pars$a, pars$b, pars$bp)
      }) |>
      ungroup()
  } else {
    bp_gsw_pts <- bp_gsw_pts |> mutate(y_bp = numeric(0))
  }

  p5 <- ggplot() +
    geom_point(data = day1_gsw_long |>
                 mutate(Genotype = factor(Genotype, levels = geno_gsw)),
               aes(x = VPD, y = gsw_reading, fill = WUE_class, shape = WUE_class),
               colour = "white", size = 2.5, alpha = 0.80, stroke = 0.3,
               position = position_jitter(width = 0.02, height = 0, seed = 42)) +
    scale_fill_manual(values = WUE_COLS, name = "WUE class") +
    scale_shape_manual(values = WUE_SHAPES, name = "WUE class") +
    geom_line(data = curves_gsw,
              aes(x = VPD, y = y_hat, colour = Model, linetype = Model,
                  linewidth = Model)) +
    mod_legend + mod_lty +
    scale_linewidth_manual(values = c("Linear" = 0.7,"Bilinear" = 1.4,
                                       "Lin. plateau" = 1.4), name = "Model") +
    {if (nrow(bp_gsw_pts) > 0)
       list(
         geom_vline(data = bp_gsw_pts, aes(xintercept = BP_kPa),
                    linetype = "dotted", colour = UA_DARK, linewidth = 0.5, alpha = 0.5),
         geom_point(data = bp_gsw_pts, aes(x = BP_kPa, y = y_bp),
                    shape = 25, size = 2.5, fill = UA_DARK, colour = "white",
                    stroke = 0.4, inherit.aes = FALSE)
       )
    } +
    facet_wrap(~Genotype, ncol = 4, scales = "fixed") +
    coord_cartesian(ylim = GSW_YLIM) +
    labs(
      title = "Fig 5  Stomatal conductance vs VPD — Day 1 (3 LICOR readings per plant; model fit to mean)",
      x     = "VPD (kPa)",
      y     = expression(g[sw]~(mol~m^{-2}~s^{-1}))
    ) +
    UA_THEME +
    theme(legend.position = "bottom") +
    guides(
      fill      = guide_legend(nrow = 2, override.aes = list(size = 3, stroke = 0)),
      shape     = guide_legend(nrow = 2),
      colour    = guide_legend(nrow = 1, override.aes = list(linewidth = 1.2)),
      linetype  = guide_legend(nrow = 1),
      linewidth = guide_legend(nrow = 1)
    )

  ggsave("figures/Fig5_gsw_models_per_genotype.png", p5,
         width = 12, height = 8, dpi = 300)
  cat("Saved Fig5\n")
}


# ---- FIG 6: AIC heatmap (TR) ------------------------------------------------
cat("Building Fig6 ...\n")
aic_heatmap_df <- tr_params |>
  filter(!is.na(Delta_AIC)) |>
  mutate(
    Model     = factor(Model, levels = c("Linear","Bilinear","Lin. plateau")),
    WUE_class = factor(WUE_class, levels = WUE_ORDER),
    Genotype  = fct_reorder(Genotype, as.integer(WUE_class))
  )

p6 <- ggplot(aic_heatmap_df, aes(x = Model, y = Genotype, fill = Delta_AIC)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f", Delta_AIC),
                colour = Delta_AIC > 6),
            size = 2.8, fontface = "bold") +
  scale_fill_gradient2(low = UA_RED, mid = "#F7C6CC", high = "white",
                       midpoint = 6, limits = c(0, NA),
                       name = "ΔAIC", na.value = "grey85") +
  scale_colour_manual(values = c("TRUE" = UA_DARK, "FALSE" = "white"),
                      guide = "none") +
  labs(
    title = "Fig 6  Model comparison — ΔAIC (TR, Day 1)\nΔAIC = 0 → best model for that genotype",
    x = NULL, y = NULL
  ) +
  UA_THEME +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 9, face = "bold")
  )

ggsave("figures/Fig6_TR_AIC_heatmap.png", p6,
       width = 7, height = 6, dpi = 300)
cat("Saved Fig6\n")


# ---- FIG 7: TR mean ± SD per genotype × VPD treatment ----------------------
cat("Building Fig7 ...\n")

tr_summary <- day1 |>
  group_by(Genotype, WUE_class, Treatment, VPD) |>
  summarise(
    mean_TR = mean(TR, na.rm = TRUE),
    sd_TR   = sd(TR, na.rm = TRUE),
    n_TR    = sum(!is.na(TR)),
    .groups = "drop"
  ) |>
  mutate(
    se_TR = sd_TR / sqrt(n_TR),
    WUE_class = factor(WUE_class, levels = WUE_ORDER)
  )

p7 <- ggplot(tr_summary, aes(x = VPD, y = mean_TR,
                               colour = WUE_class, fill = WUE_class,
                               shape = WUE_class, group = Genotype)) +
  geom_errorbar(aes(ymin = mean_TR - sd_TR, ymax = mean_TR + sd_TR),
                width = 0.05, linewidth = 0.55, alpha = 0.55) +
  geom_point(size = 3, alpha = 0.85, colour = "white", stroke = 0.4) +
  scale_fill_manual(values = WUE_COLS, name = "WUE class") +
  scale_colour_manual(values = WUE_COLS, name = "WUE class") +
  scale_shape_manual(values = WUE_SHAPES, name = "WUE class") +
  coord_cartesian(ylim = TR_YLIM) +
  labs(
    title = "Fig 7  Transpiration rate — genotype mean +/- SD across replicates",
    subtitle = "Day 1; one symbol per genotype per VPD treatment",
    x = "VPD (kPa)",
    y = expression(TR~(mg~m^{-2}~s^{-1}))
  ) +
  UA_THEME +
  theme(legend.position = "right")

ggsave("figures/Fig7_TR_mean_SD.png", p7,
       width = 10, height = 5.5, dpi = 300)
cat("Saved Fig7\n")


# ---- FIG 8: gsw mean ± SD (from 3 LICOR readings) per genotype × VPD ------
if (nrow(day1_gsw_long) > 0) {
  cat("Building Fig8 ...\n")

  gsw_summary <- day1_gsw_long |>
    group_by(Genotype, WUE_class, VPD) |>
    summarise(
      mean_gsw = mean(gsw_reading, na.rm = TRUE),
      sd_gsw   = sd(gsw_reading, na.rm = TRUE),
      n_gsw    = sum(!is.na(gsw_reading)),
      .groups  = "drop"
    ) |>
    filter(!is.na(VPD)) |>
    mutate(WUE_class = factor(WUE_class, levels = WUE_ORDER))

  p8 <- ggplot(gsw_summary, aes(x = VPD, y = mean_gsw,
                                  colour = WUE_class, fill = WUE_class,
                                  shape = WUE_class, group = Genotype)) +
    geom_errorbar(aes(ymin = mean_gsw - sd_gsw, ymax = mean_gsw + sd_gsw),
                  width = 0.05, linewidth = 0.55, alpha = 0.55) +
    geom_point(size = 3, alpha = 0.85, colour = "white", stroke = 0.4) +
    scale_fill_manual(values = WUE_COLS, name = "WUE class") +
    scale_colour_manual(values = WUE_COLS, name = "WUE class") +
    scale_shape_manual(values = WUE_SHAPES, name = "WUE class") +
    coord_cartesian(ylim = GSW_YLIM) +
    labs(
      title = "Fig 8  Stomatal conductance — genotype mean +/- SD across 3 LICOR readings",
      subtitle = "Day 1; one symbol per genotype per VPD treatment",
      x = "VPD (kPa)",
      y = expression(g[sw]~(mol~m^{-2}~s^{-1}))
    ) +
    UA_THEME +
    theme(legend.position = "right")

  ggsave("figures/Fig8_gsw_mean_SD.png", p8,
         width = 10, height = 5.5, dpi = 300)
  cat("Saved Fig8\n")
}


cat("\n=== Model fitting complete ===\n")
cat("Figures : figures/\n")
cat("Tables  : data/TR_best_model_summary.csv\n")
cat("          data/gsw_best_model_summary.csv\n")
cat("          data/TR_model_parameters.csv\n")
cat("          data/gsw_model_parameters.csv\n")
