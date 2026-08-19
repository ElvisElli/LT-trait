# ==============================================================================
# paper_figures.R
# Soybean Limited-Transpiration Experiment (GC 2026) — Paper Figures
#
# Run from this file's directory  (setwd to wherever paper_figures.R lives):
#   source("paper_figures.R")       # all 5 figures
#   # — OR — highlight one local({}) block and run it alone
#
# Each FIGURE section is completely self-contained:
#   1. Loads data_all.xlsx
#   2. Fits three candidate models per genotype with Rep as random effect
#   3. Selects the best model by AIC (lower = better) — printed to console
#   4. Builds and saves the TIFF
#
# Nothing is shared across sections: changing / running one figure does not
# affect any other.
#
# Data  : data_all.xlsx  (in this directory)
# Output: TIFF files     (in this directory)
# DPI   : set below — 600 for submission, 150 for a quick screen preview
# ==============================================================================

DPI    <- 600     # change to 150 for fast preview
OUTDIR <- "."     # TIFFs land next to this script

# Headless-safe TIFF device (Linux without X11)
options(bitmapType = "cairo")
tiff_cairo <- function(filename, width, height, units = "in", res, ...) {
  tiff(filename, width = width, height = height, units = units,
       res = res, compression = "lzw", type = "cairo")
}

# --------------------------------------------------------------------------- #
# Packages — loaded once at the top, available in all local() scopes
# --------------------------------------------------------------------------- #
suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(forcats)
  library(tibble)
  library(nlme)
  library(patchwork)
  library(ggrepel)
  library(ggh4x)
})

# --------------------------------------------------------------------------- #
# Colour system and shared ggplot theme
# --------------------------------------------------------------------------- #
UA_RED   <- "#9D2235"
UA_WINE  <- "#6E1424"
UA_BLUE  <- "#4393c3"
UA_GRAY  <- "#7A7A7A"
UA_DARK  <- "#1A1A1A"

# WUE class display order and colour mapping
WUE_ORDER <- c("slow_wilting", "high_wilting",
                "Plastic, High WUE", "Nonplastic, Low WUE", "MAGIC parent")

# Nature style: all dark green for consistency
WUE_COLS <- c(
  "slow_wilting"        = "#31a354",
  "high_wilting"        = "#31a354",
  "Plastic, High WUE"   = "#31a354",
  "Nonplastic, Low WUE" = "#31a354",
  "MAGIC parent"        = "#31a354"
)

UA_THEME <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 0.7, fill = NA),
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.7),
    strip.text       = element_text(colour = "black", face = "plain", size = 10),
    axis.title       = element_text(size = 9, face = "plain"),
    axis.text        = element_text(size = 8, face = "plain", colour = "black"),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    legend.key.size  = unit(0.45, "cm"),
    legend.text      = element_text(size = 8, face = "plain"),
    legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.5)
  )

# --------------------------------------------------------------------------- #
# Model functions — three candidate shapes for TR / gsw vs VPD
#
#   Linear       y = a + b1 * VPD               (stomata open linearly)
#   Bilinear     two-slope segmented with BP     (limited response above BP)
#   Lin. plateau flat after BP                   (classic LT trait response)
# --------------------------------------------------------------------------- #

bilinear_fn <- function(x, a, b1, b2, bp) {
  # Rises with slope b1 up to bp, then continues with slope b2
  ifelse(x <= bp,
         a + b1 * x,
         a + b1 * bp + b2 * (x - bp))
}

lin_plateau_fn <- function(x, a, b, bp) {
  # Rises with slope b up to bp, then stays flat (y_max = a + b * bp)
  ifelse(x <= bp,
         a + b * x,
         a + b * bp)
}

# --------------------------------------------------------------------------- #
# fit_all_models()
#
# Fits Linear, Bilinear, and Linear-plateau models to one genotype's data.
# Rep (replicate) enters as a random intercept via nlme.
# Returns a list with:
#   $Linear, $Bilinear, $"Lin. plateau"  — model objects
#   $AIC  — named numeric vector of AIC for each model
#   $best — name of the model with the lowest AIC
# --------------------------------------------------------------------------- #
fit_all_models <- function(data, response = "TR") {
  y   <- data[[response]]
  vpd <- data$VPD
  rep <- data$Rep
  n   <- sum(!is.na(y) & !is.na(vpd))
  if (n < 4) return(NULL)

  vmin <- min(vpd, na.rm = TRUE)
  vmax <- max(vpd, na.rm = TRUE)
  bp0  <- (vmin + vmax) / 2            # breakpoint starting value = midpoint
  sl0  <- coef(lm(y ~ vpd, na.action = na.omit))["vpd"]  # OLS slope as start
  if (is.na(sl0)) sl0 <- 1
  a0   <- mean(y, na.rm = TRUE) - sl0 * mean(vpd, na.rm = TRUE)

  df_fit <- data.frame(y = y, VPD = vpd, Rep = rep) |> na.omit()
  out    <- list()

  # ── (A) LINEAR — lme with random intercept per Rep ────────────────────────
  out[["Linear"]] <- tryCatch(
    lme(y ~ VPD, random = ~1 | Rep, data = df_fit, method = "ML",
        control = lmeControl(opt = "optim")),
    error = function(e) tryCatch(lm(y ~ VPD, data = df_fit),
                                 error = function(e2) NULL)
  )

  # ── (B) BILINEAR — nls starting values → nlme with random intercept ───────
  nls_bil <- tryCatch(
    nls(y ~ bilinear_fn(VPD, a, b1, b2, bp),
        data      = df_fit,
        start     = list(a = a0, b1 = sl0, b2 = sl0 * 0.1, bp = bp0),
        lower     = c(a = -Inf, b1 = 0,   b2 = -Inf, bp = vmin),
        upper     = c(a =  Inf, b1 = Inf, b2 =  Inf, bp = vmax),
        algorithm = "port",
        control   = nls.control(maxiter = 300, tol = 1e-5)),
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
      error = function(e) nls_bil   # fall back to nls if nlme fails
    )
  }

  # ── (C) LINEAR PLATEAU — nls starting values → nlme ─────────────────────
  nls_lp <- tryCatch(
    nls(y ~ lin_plateau_fn(VPD, a, b, bp),
        data      = df_fit,
        start     = list(a = a0, b = sl0, bp = bp0),
        lower     = c(a = -Inf, b = 0,   bp = vmin),
        upper     = c(a =  Inf, b = Inf, bp = vmax),
        algorithm = "port",
        control   = nls.control(maxiter = 300, tol = 1e-5)),
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

  # ── AIC-based model selection ─────────────────────────────────────────────
  get_aic <- function(m) tryCatch(AIC(m), error = function(e) NA_real_)
  aics <- c(
    "Linear"       = get_aic(out[["Linear"]]),
    "Bilinear"     = get_aic(out[["Bilinear"]]),
    "Lin. plateau" = get_aic(out[["Lin. plateau"]])
  )
  out$AIC  <- aics
  out$best <- names(which.min(aics))   # name of the best (lowest AIC) model
  out
}

# --------------------------------------------------------------------------- #
# extract_params()
#
# Pulls fixed-effect coefficients from a fitted model into a tidy list so the
# same downstream code works regardless of whether the model is lme or nlme.
# --------------------------------------------------------------------------- #
extract_params <- function(fit, model_name) {
  m <- fit[[model_name]]
  if (is.null(m)) return(NULL)
  tryCatch(
    as.list(if (inherits(m, "lme") || inherits(m, "nlme")) fixef(m) else coef(m)),
    error = function(e) NULL
  )
}

# --------------------------------------------------------------------------- #
# predict_curve()
#
# Evaluates the fitted curve for one genotype at a grid of VPD values.
# Returns a data frame with columns: VPD, y_hat.
# --------------------------------------------------------------------------- #
predict_curve <- function(fit, model_name, vpd_seq) {
  pars <- extract_params(fit, model_name)
  if (is.null(pars)) return(NULL)
  a  <- pars[["a"]] %||% pars[["(Intercept)"]] %||% NA_real_
  y  <- switch(model_name,
    "Linear"       = a + pars[["VPD"]] * vpd_seq,
    "Bilinear"     = bilinear_fn(vpd_seq, pars$a, pars$b1, pars$b2, pars$bp),
    "Lin. plateau" = lin_plateau_fn(vpd_seq, pars$a, pars$b, pars$bp),
    rep(NA_real_, length(vpd_seq))
  )
  data.frame(VPD = vpd_seq, y_hat = pmax(y, 0))
}

# Null-coalescing helper
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a[1])) a else b

cat("=== paper_figures.R loaded — building", DPI, "DPI TIFFs ===\n\n")


# ==============================================================================
# FIGURE 1
# TR response to VPD per genotype — best-fit model curve + raw Day 1 data
#
# For each genotype we fit three models and pick the winner by AIC.
# The panel strip shows the genotype name and the winning model type.
# ==============================================================================
local({
  cat("──────────────────────────────────────────────────\n")
  cat("FIGURE 1  TR response to VPD (best-fit per genotype)\n")
  cat("──────────────────────────────────────────────────\n")

  # ── 1. Load and prepare data ────────────────────────────────────────────────
  raw <- read_excel("data_all.xlsx", sheet = "dataAll")
  df  <- raw |>
    select(where(~ !all(is.na(.)))) |>
    rename(Rep = RepNo, Day = day, TR = TR_mg_sqm_seg, VPD = avg_vpd) |>
    mutate(
      Rep       = factor(Rep),
      Day_label = factor(ifelse(Day == 1, "Day 1", "Day 2"),
                         levels = c("Day 1", "Day 2")),
      Genotype  = factor(Genotype),
      WUE_class = factor(WUE_class, levels = WUE_ORDER)
    )
  day1      <- filter(df, Day_label == "Day 1")
  genotypes <- levels(droplevels(day1$Genotype))
  VPD_SEQ   <- seq(min(day1$VPD, na.rm = TRUE) * 0.95,
                   max(day1$VPD, na.rm = TRUE) * 1.05, length.out = 300)

  # ── 2. Fit models — AIC table printed to console ────────────────────────────
  cat(sprintf("\n  %-14s  %7s  %8s  %12s  %s\n",
              "Genotype", "Linear", "Bilinear", "Lin.plateau", "Best (AIC)"))
  cat("  ", strrep("-", 60), "\n", sep = "")

  tr_fits <- setNames(
    lapply(genotypes, function(g) {
      d   <- filter(day1, Genotype == g) |> select(TR, VPD, Rep)
      fit <- fit_all_models(d, response = "TR")
      if (!is.null(fit)) {
        cat(sprintf("  %-14s  %7.1f  %8.1f  %12.1f  %s\n",
                    g,
                    fit$AIC["Linear"], fit$AIC["Bilinear"],
                    fit$AIC["Lin. plateau"], fit$best))
      }
      fit
    }),
    genotypes
  )

  # ── 3. Build parameter summary (best model per genotype) ────────────────────
  best_params <- bind_rows(lapply(genotypes, function(g) {
    fit <- tr_fits[[g]]
    if (is.null(fit)) return(NULL)
    pars  <- extract_params(fit, fit$best)
    wue   <- as.character(day1 |> filter(Genotype == g) |> pull(WUE_class) |> unique())
    tibble(
      Genotype   = g,
      WUE_class  = wue,
      Best_model = fit$best,
      Intercept  = round(pars[["a"]] %||% pars[["(Intercept)"]], 2),
      Slope1     = round(pars[["VPD"]] %||% pars[["b1"]] %||% pars[["b"]], 3),
      BP_kPa     = if (!is.null(pars[["bp"]])) round(pars[["bp"]], 2) else NA_real_,
      Slope2     = if (!is.null(pars[["b2"]])) round(pars[["b2"]], 3) else NA_real_,
      AIC_best   = round(fit$AIC[fit$best], 1)
    )
  }))

  cat("\n  Best-model parameter summary:\n")
  print(best_params, n = Inf)

  # ── 4. Build prediction curves (best model only) ───────────────────────────
  curves <- bind_rows(lapply(genotypes, function(g) {
    fit <- tr_fits[[g]]
    if (is.null(fit)) return(NULL)
    d <- predict_curve(fit, fit$best, VPD_SEQ)
    if (is.null(d)) return(NULL)
    d |> mutate(Genotype = g)
  })) |>
    left_join(best_params |> select(Genotype, WUE_class, Best_model),
              by = "Genotype") |>
    mutate(WUE_class = factor(WUE_class, levels = WUE_ORDER))

  # Panel label: genotype + WUE class description
  panel_order <- best_params |>
    arrange(match(WUE_class, WUE_ORDER), Genotype) |>
    mutate(Panel_label = paste0(Genotype, "\n(", WUE_class, ")"))

  day1_lab   <- left_join(day1,   panel_order |> select(Genotype, Panel_label), by = "Genotype")
  curves_lab <- left_join(curves, panel_order |> select(Genotype, Panel_label), by = "Genotype")

  lvls       <- panel_order$Panel_label
  day1_lab   <- mutate(day1_lab,   Panel_label = factor(Panel_label, levels = lvls))
  curves_lab <- mutate(curves_lab, Panel_label = factor(Panel_label, levels = lvls))

  # ── 5. Plot ─────────────────────────────────────────────────────────────────
  p <- ggplot(day1_lab, aes(x = VPD, y = TR, colour = WUE_class)) +
    geom_point(size = 2.5, shape = 21, colour = "black", fill = "#31a354", stroke = 0.8, alpha = 0.9) +
    geom_line(data = curves_lab,
              aes(x = VPD, y = y_hat, colour = WUE_class),
              linewidth = 1.0, inherit.aes = FALSE) +
    scale_colour_manual(values = WUE_COLS, name = NULL) +
    facet_wrap(~Panel_label, ncol = 7, scales = "fixed") +
    coord_cartesian(ylim = c(0, 120)) +
    labs(x = "VPD (kPa)",
         y = expression(TR~(mg~m^{-2}~s^{-1}))) +
    UA_THEME +
    theme(legend.position = "none") +
    guides(colour = "none")

  ggsave(file.path(OUTDIR, "Figure1_TR_response_VPD.tiff"), p,
         width = 16, height = 6, dpi = DPI, device = tiff_cairo)
  cat("\n  Saved: Figure1_TR_response_VPD.tiff\n")

  # Also save as PNG for preview
  ggsave(file.path(OUTDIR, "Figure1_TR_response_VPD.png"), p,
         width = 16, height = 6, dpi = 150, device = "png")
  cat("  Saved: Figure1_TR_response_VPD.png (preview)\n\n")
})


# ==============================================================================
# FIGURE 2
# Genotypic variability in VPD breakpoints
#
# Only genotypes whose best model is Bilinear or Linear-plateau have a
# meaningful breakpoint.  Bars sorted by breakpoint value, coloured by WUE class.
# ==============================================================================
local({
  cat("──────────────────────────────────────────────────\n")
  cat("FIGURE 2  VPD breakpoints (non-linear models only)\n")
  cat("──────────────────────────────────────────────────\n")

  # ── 1. Load and prepare data ────────────────────────────────────────────────
  raw <- read_excel("data_all.xlsx", sheet = "dataAll")
  df  <- raw |>
    select(where(~ !all(is.na(.)))) |>
    rename(Rep = RepNo, Day = day, TR = TR_mg_sqm_seg, VPD = avg_vpd) |>
    mutate(
      Rep       = factor(Rep),
      Day_label = factor(ifelse(Day == 1, "Day 1", "Day 2"),
                         levels = c("Day 1", "Day 2")),
      Genotype  = factor(Genotype),
      WUE_class = factor(WUE_class, levels = WUE_ORDER)
    )
  day1      <- filter(df, Day_label == "Day 1")
  genotypes <- levels(droplevels(day1$Genotype))

  # ── 2. Fit models — same process as Figure 1 ────────────────────────────────
  cat(sprintf("\n  %-14s  %7s  %8s  %12s  %s\n",
              "Genotype", "Linear", "Bilinear", "Lin.plateau", "Best (AIC)"))
  cat("  ", strrep("-", 60), "\n", sep = "")

  tr_fits <- setNames(
    lapply(genotypes, function(g) {
      d   <- filter(day1, Genotype == g) |> select(TR, VPD, Rep)
      fit <- fit_all_models(d, response = "TR")
      if (!is.null(fit))
        cat(sprintf("  %-14s  %7.1f  %8.1f  %12.1f  %s\n",
                    g, fit$AIC["Linear"], fit$AIC["Bilinear"],
                    fit$AIC["Lin. plateau"], fit$best))
      fit
    }),
    genotypes
  )

  # ── 3. Extract breakpoints from non-linear best models ──────────────────────
  bp_df <- bind_rows(lapply(genotypes, function(g) {
    fit <- tr_fits[[g]]
    if (is.null(fit) || !fit$best %in% c("Bilinear", "Lin. plateau")) return(NULL)
    pars <- extract_params(fit, fit$best)
    wue  <- as.character(day1 |> filter(Genotype == g) |> pull(WUE_class) |> unique())
    tibble(
      Genotype  = g,
      WUE_class = wue,
      Best_model = fit$best,
      BP_kPa    = if (!is.null(pars[["bp"]])) round(pars[["bp"]], 3) else NA_real_
    )
  })) |>
    mutate(
      WUE_class = factor(WUE_class, levels = WUE_ORDER),
      Genotype  = fct_reorder(Genotype, BP_kPa)
    )

  cat(sprintf("\n  %d genotypes with a breakpoint (range %.2f – %.2f kPa, mean %.2f kPa)\n",
              nrow(bp_df), min(bp_df$BP_kPa), max(bp_df$BP_kPa), mean(bp_df$BP_kPa)))
  print(bp_df)

  # ── 4. Plot ─────────────────────────────────────────────────────────────────
  p <- ggplot(bp_df, aes(x = Genotype, y = BP_kPa, fill = WUE_class)) +
    geom_col(width = 0.65, colour = "white", linewidth = 0.4) +
    geom_hline(yintercept = mean(bp_df$BP_kPa), linetype = "dashed",
               colour = UA_GRAY, linewidth = 0.6) +
    annotate("text",
             x     = nrow(bp_df) * 0.97,
             y     = mean(bp_df$BP_kPa) + 0.04,
             label = sprintf("mean = %.2f kPa", mean(bp_df$BP_kPa)),
             hjust = 1, vjust = 0, size = 2.8, colour = UA_GRAY, fontface = "italic") +
    scale_fill_manual(values = WUE_COLS, name = "WUE class") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = NULL, y = "Breakpoint VPD (kPa)") +
    UA_THEME +
    theme(axis.text.x     = element_text(angle = 35, hjust = 1, size = 8),
          legend.position = "bottom")

  ggsave(file.path(OUTDIR, "Figure2_TR_breakpoints.tiff"), p,
         width = 6, height = 4.5, dpi = DPI, device = tiff_cairo)
  cat("\n  Saved: Figure2_TR_breakpoints.tiff\n\n")
})


# ==============================================================================
# FIGURE 3
# Genotype classification: intercept (baseline TR at low VPD) vs
# initial slope (VPD sensitivity) from the best-fit TR model
#
# Quadrant thresholds: medians across all genotypes.
# Plastic and Nonplastic WUE classes highlighted; others shown in grey.
# ==============================================================================
local({
  cat("──────────────────────────────────────────────────\n")
  cat("FIGURE 3  Intercept vs slope quadrant classification\n")
  cat("──────────────────────────────────────────────────\n")

  # ── 1. Load and prepare data ────────────────────────────────────────────────
  raw <- read_excel("data_all.xlsx", sheet = "dataAll")
  df  <- raw |>
    select(where(~ !all(is.na(.)))) |>
    rename(Rep = RepNo, Day = day, TR = TR_mg_sqm_seg, VPD = avg_vpd) |>
    mutate(
      Rep       = factor(Rep),
      Day_label = factor(ifelse(Day == 1, "Day 1", "Day 2"),
                         levels = c("Day 1", "Day 2")),
      Genotype  = factor(Genotype),
      WUE_class = factor(WUE_class, levels = WUE_ORDER)
    )
  day1      <- filter(df, Day_label == "Day 1")
  # Filter out: AG43XF2, PI471938, PI442012A, PI603171, PI603418D
  day1      <- filter(day1, !(Genotype %in% c("AG43XF2", "PI471938", "PI442012A", "PI603171", "PI603418D")))
  genotypes <- levels(droplevels(day1$Genotype))

  # ── 2. Fit models ────────────────────────────────────────────────────────────
  cat(sprintf("\n  %-14s  %7s  %8s  %12s  %s\n",
              "Genotype", "Linear", "Bilinear", "Lin.plateau", "Best (AIC)"))
  cat("  ", strrep("-", 60), "\n", sep = "")

  tr_fits <- setNames(
    lapply(genotypes, function(g) {
      d   <- filter(day1, Genotype == g) |> select(TR, VPD, Rep)
      fit <- fit_all_models(d, response = "TR")
      if (!is.null(fit))
        cat(sprintf("  %-14s  %7.1f  %8.1f  %12.1f  %s\n",
                    g, fit$AIC["Linear"], fit$AIC["Bilinear"],
                    fit$AIC["Lin. plateau"], fit$best))
      fit
    }),
    genotypes
  )

  # ── 3. Extract intercept and initial slope from best model ──────────────────
  fig_df <- bind_rows(lapply(genotypes, function(g) {
    fit  <- tr_fits[[g]]
    if (is.null(fit)) return(NULL)
    pars <- extract_params(fit, fit$best)
    wue  <- as.character(day1 |> filter(Genotype == g) |> pull(WUE_class) |> unique())
    tibble(
      Genotype   = g,
      WUE_class  = wue,
      Best_model = fit$best,
      AIC_best   = round(fit$AIC[fit$best], 1),
      Intercept  = round(pars[["a"]] %||% pars[["(Intercept)"]], 2),
      # Slope1 = initial slope regardless of model type
      Slope1     = round(pars[["VPD"]] %||% pars[["b1"]] %||% pars[["b"]], 3)
    )
  })) |>
    # Three-level grouping: highlight Plastic and Nonplastic, grey everything else
    mutate(
      WUE_group = case_when(
        WUE_class == "Plastic, High WUE"   ~ "Plastic, High WUE",
        WUE_class == "Nonplastic, Low WUE" ~ "Nonplastic, Low WUE",
        TRUE                                ~ "Other"
      )
    )

  # Quadrant thresholds (medians)
  x_mid     <- median(fig_df$Intercept, na.rm = TRUE)
  slope_mid <- median(fig_df$Slope1,    na.rm = TRUE)
  x_rng     <- range(fig_df$Intercept, na.rm = TRUE)

  cat(sprintf("\n  Quadrant thresholds: median intercept = %.1f  median slope = %.3f\n",
              x_mid, slope_mid))
  cat("\n  Full parameter table:\n")
  print(fig_df |> select(Genotype, WUE_class, Best_model, AIC_best, Intercept, Slope1))

  # Colour mapping for the three WUE groups (Nature style: all dark green)
  WUE_COLS_Q <- c(
    "Plastic, High WUE"   = "#31a354",
    "Nonplastic, Low WUE" = "#31a354",
    "Other"               = "#31a354"
  )

  # Shaded quadrant backgrounds (minimal, light)
  quad_shade <- data.frame(
    xmin = c(-Inf,  x_mid, -Inf,  x_mid),
    xmax = c(x_mid,  Inf,  x_mid,  Inf),
    ymin = c(slope_mid, slope_mid, -Inf,       -Inf),
    ymax = c(Inf,       Inf,       slope_mid,  slope_mid),
    fill = c("#f5f5f5", "#f5f5f5", "#ffffff",  "#ffffff")
  )

  # ── 4. Plot ─────────────────────────────────────────────────────────────────
  p <- ggplot(fig_df, aes(x = Intercept, y = Slope1, colour = WUE_group)) +
    # Shaded quadrant backgrounds
    geom_rect(data = quad_shade,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
              inherit.aes = FALSE, alpha = 0.55) +
    scale_fill_identity() +
    # Threshold lines
    geom_vline(xintercept = x_mid,     linetype = "dashed",
               colour = UA_GRAY, linewidth = 0.65) +
    geom_hline(yintercept = slope_mid, linetype = "dashed",
               colour = UA_GRAY, linewidth = 0.65) +
    # Quadrant text labels (corners)
    annotate("text", x = x_mid - diff(x_rng) * 0.03, y = Inf,
             label = "Low water use\nunder high VPD",
             hjust = 1, vjust = 2.0, size = 2.6,
             colour = "black", fontface = "plain", lineheight = 0.9) +
    annotate("text", x = x_rng[2] - diff(x_rng) * 0.05, y = Inf,
             label = "High water use\nunder high VPD",
             hjust = 1, vjust = 2.0, size = 2.6,
             colour = "black", fontface = "plain", lineheight = 0.9) +
    annotate("text", x = x_rng[1] + diff(x_rng) * 0.05, y = -Inf,
             label = "Low water use\nunder low VPD",
             hjust = 0, vjust = -1.0, size = 2.6,
             colour = "black", fontface = "plain", lineheight = 0.9) +
    annotate("text", x = x_rng[2] - diff(x_rng) * 0.05, y = -Inf,
             label = "High water use\nunder low VPD",
             hjust = 1, vjust = -1.0, size = 2.6,
             colour = "black", fontface = "plain", lineheight = 0.9) +
    # Points and labels (Nature style: dark green with black borders)
    geom_point(aes(shape = Best_model), size = 3.5, colour = "black", fill = "#31a354",
               stroke = 1.0, alpha = 0.9) +
    geom_label_repel(
      aes(label = Genotype),
      size = 2.5, box.padding = 0.4, point.padding = 0.25,
      max.overlaps = 20, label.size = 0.3, alpha = 0.85,
      show.legend = FALSE, seed = 42,
      fill = "white", colour = "black"
    ) +
    scale_colour_manual(values = WUE_COLS_Q, name = NULL,
                        breaks = c("Plastic, High WUE", "Nonplastic, Low WUE")) +
    scale_shape_manual(
      values = c("Linear" = 21, "Bilinear" = 22, "Lin. plateau" = 24),
      name = NULL
    ) +
    labs(
      x = expression(Intercept~~(mg~m^{-2}~s^{-1})),
      y = expression(Slope~~(mg~m^{-2}~s^{-1}~kPa^{-1}))
    ) +
    UA_THEME +
    theme(legend.position = "bottom",
          legend.box = "horizontal")

  ggsave(file.path(OUTDIR, "Figure3_intercept_slope_quadrants.tiff"), p,
         width = 9, height = 6.5, dpi = DPI, device = tiff_cairo)
  cat("\n  Saved: Figure3_intercept_slope_quadrants.tiff\n\n")
})


# ==============================================================================
# FIGURE 4
# Stomatal conductance (gsw) response to VPD per genotype
# Raw individual LICOR readings (3 per plant) + best-fit model curve
#
# Only genotypes with LICOR data (gsw1/gsw2/gsw3) are included.
# ==============================================================================
local({
  cat("──────────────────────────────────────────────────\n")
  cat("FIGURE 4  gsw response to VPD (best-fit per genotype)\n")
  cat("──────────────────────────────────────────────────\n")

  # ── 1. Load and prepare data ────────────────────────────────────────────────
  raw <- read_excel("data_all.xlsx", sheet = "dataAll")
  df  <- raw |>
    select(where(~ !all(is.na(.)))) |>
    rename(Rep = RepNo, Day = day, VPD = avg_vpd) |>
    mutate(
      Rep       = factor(Rep),
      Day_label = factor(ifelse(Day == 1, "Day 1", "Day 2"),
                         levels = c("Day 1", "Day 2")),
      Genotype  = factor(Genotype),
      WUE_class = factor(WUE_class, levels = WUE_ORDER),
      gsw       = rowMeans(cbind(gsw1, gsw2, gsw3), na.rm = FALSE)
    )
  # Only rows where all three LICOR readings exist
  all_na <- is.na(df$gsw1) & is.na(df$gsw2) & is.na(df$gsw3)
  df$gsw[all_na] <- NA

  day1     <- filter(df, Day_label == "Day 1")
  day1_gsw <- filter(day1, !is.na(gsw))
  geno_gsw <- levels(droplevels(day1_gsw$Genotype))

  # Long format: each plant contributes 3 individual LICOR readings
  day1_long <- day1_gsw |>
    pivot_longer(cols = c(gsw1, gsw2, gsw3),
                 names_to = "licor_rep", values_to = "gsw_reading") |>
    filter(!is.na(gsw_reading))

  VPD_SEQ <- seq(min(day1_gsw$VPD, na.rm = TRUE) * 0.95,
                 max(day1_gsw$VPD, na.rm = TRUE) * 1.05, length.out = 300)

  # ── 2. Fit models ────────────────────────────────────────────────────────────
  cat(sprintf("\n  %-14s  %7s  %8s  %12s  %s\n",
              "Genotype", "Linear", "Bilinear", "Lin.plateau", "Best (AIC)"))
  cat("  ", strrep("-", 60), "\n", sep = "")

  gsw_fits <- setNames(
    lapply(geno_gsw, function(g) {
      # Rename gsw → TR so fit_all_models() works without modification
      d   <- filter(day1_gsw, Genotype == g) |>
             select(gsw, VPD, Rep) |>
             rename(TR = gsw)
      fit <- fit_all_models(d, response = "TR")
      if (!is.null(fit))
        cat(sprintf("  %-14s  %7.1f  %8.1f  %12.1f  %s\n",
                    g, fit$AIC["Linear"], fit$AIC["Bilinear"],
                    fit$AIC["Lin. plateau"], fit$best))
      fit
    }),
    geno_gsw
  )

  # ── 3. Build parameter summary ───────────────────────────────────────────────
  best_params <- bind_rows(lapply(geno_gsw, function(g) {
    fit  <- gsw_fits[[g]]
    if (is.null(fit)) return(NULL)
    pars <- extract_params(fit, fit$best)
    wue  <- as.character(day1_gsw |> filter(Genotype == g) |> pull(WUE_class) |> unique())
    tibble(
      Genotype   = g,
      WUE_class  = wue,
      Best_model = fit$best,
      AIC_best   = round(fit$AIC[fit$best], 1),
      Intercept  = round(pars[["a"]] %||% pars[["(Intercept)"]], 4),
      Slope1     = round(pars[["VPD"]] %||% pars[["b1"]] %||% pars[["b"]], 4)
    )
  }))

  cat("\n  Best-model parameter summary (gsw):\n")
  print(best_params, n = Inf)

  # ── 4. Build prediction curves ───────────────────────────────────────────────
  curves <- bind_rows(lapply(geno_gsw, function(g) {
    fit <- gsw_fits[[g]]
    if (is.null(fit)) return(NULL)
    d <- predict_curve(fit, fit$best, VPD_SEQ)
    if (is.null(d)) return(NULL)
    d |> mutate(Genotype = g)
  })) |>
    left_join(best_params |> select(Genotype, WUE_class, Best_model),
              by = "Genotype") |>
    mutate(WUE_class = factor(WUE_class, levels = WUE_ORDER))

  # Panel labels
  panel_order <- best_params |>
    arrange(match(WUE_class, WUE_ORDER), Genotype) |>
    mutate(Panel_label = paste0(Genotype, "\n[", Best_model, "]"))

  day1_long_lab <- left_join(day1_long,
                              panel_order |> select(Genotype, Panel_label),
                              by = "Genotype") |>
    mutate(Panel_label = factor(Panel_label, levels = panel_order$Panel_label),
           WUE_class   = factor(WUE_class,   levels = WUE_ORDER))
  curves_lab <- left_join(curves,
                           panel_order |> select(Genotype, Panel_label),
                           by = "Genotype") |>
    mutate(Panel_label = factor(Panel_label, levels = panel_order$Panel_label))

  GSW_YLIM <- c(0, max(day1_long$gsw_reading, na.rm = TRUE) * 1.10)

  # ── 5. Plot ─────────────────────────────────────────────────────────────────
  p <- ggplot(day1_long_lab, aes(x = VPD, y = gsw_reading, colour = WUE_class)) +
    geom_point(size = 1.5, alpha = 0.65,
               position = position_jitter(width = 0.02, seed = 1)) +
    geom_line(data    = curves_lab,
              aes(x = VPD, y = y_hat, colour = WUE_class),
              linewidth = 0.9, inherit.aes = FALSE) +
    scale_colour_manual(values = WUE_COLS, name = "WUE class") +
    facet_wrap(~Panel_label, ncol = 5, scales = "fixed") +
    coord_cartesian(ylim = GSW_YLIM) +
    labs(x = "VPD (kPa)",
         y = expression(g[sw]~(mol~m^{-2}~s^{-1}))) +
    UA_THEME +
    theme(strip.text      = element_text(colour = "white", face = "bold",
                                         size = 6.5, lineheight = 1.0),
          legend.position = "bottom") +
    guides(colour = guide_legend(nrow = 2,
                                 override.aes = list(size = 2.5, linewidth = 0.8)))

  ggsave(file.path(OUTDIR, "Figure4_gsw_response_VPD.tiff"), p,
         width = 12, height = 6.5, dpi = DPI, device = tiff_cairo)
  cat("\n  Saved: Figure4_gsw_response_VPD.tiff\n\n")
})


# ==============================================================================
# FIGURE 5
# Predicted TR and gsw at VPD 1.5 kPa and 3.5 kPa per genotype
# Bar charts arranged by WUE class (coloured strips), one bar per VPD target.
#
# Predictions extracted from best-fit model at each target VPD.
# ==============================================================================
local({
  cat("──────────────────────────────────────────────────\n")
  cat("FIGURE 5  Predicted TR and gsw at 1.5 and 3.5 kPa\n")
  cat("──────────────────────────────────────────────────\n")

  VPD_TARGETS <- c(1.5, 3.5)           # kPa values to predict at
  VPD_COLS    <- c("1.5 kPa" = UA_BLUE, "3.5 kPa" = UA_WINE)

  # Short WUE labels and strip colours for facet strips
  wue_short <- c(
    "slow_wilting"        = "Slow wilting",
    "high_wilting"        = "High wilting",
    "Plastic, High WUE"   = "Plastic / High WUE",
    "Nonplastic, Low WUE" = "Nonplastic / Low WUE",
    "MAGIC parent"        = "MAGIC parent"
  )
  WUE_STRIP_COLS <- c(
    "Slow wilting"         = "#6E1424",
    "High wilting"         = "#C8102E",
    "Plastic / High WUE"   = "#4393c3",
    "Nonplastic / Low WUE" = "#d4824a",
    "MAGIC parent"         = "#555555"
  )

  # ── 1. Load and prepare data ────────────────────────────────────────────────
  raw <- read_excel("data_all.xlsx", sheet = "dataAll")
  df  <- raw |>
    select(where(~ !all(is.na(.)))) |>
    rename(Rep = RepNo, Day = day, TR = TR_mg_sqm_seg, VPD = avg_vpd) |>
    mutate(
      Rep       = factor(Rep),
      Day_label = factor(ifelse(Day == 1, "Day 1", "Day 2"),
                         levels = c("Day 1", "Day 2")),
      Genotype  = factor(Genotype),
      WUE_class = factor(WUE_class, levels = WUE_ORDER),
      gsw       = rowMeans(cbind(gsw1, gsw2, gsw3), na.rm = FALSE)
    )
  all_na <- is.na(df$gsw1) & is.na(df$gsw2) & is.na(df$gsw3)
  df$gsw[all_na] <- NA

  day1      <- filter(df, Day_label == "Day 1")
  day1_gsw  <- filter(day1, !is.na(gsw))
  genotypes <- levels(droplevels(day1$Genotype))
  geno_gsw  <- levels(droplevels(day1_gsw$Genotype))

  # ── 2a. Fit TR models ────────────────────────────────────────────────────────
  cat(sprintf("\n  TR fits:\n  %-14s  %7s  %8s  %12s  %s\n",
              "Genotype", "Linear", "Bilinear", "Lin.plateau", "Best (AIC)"))
  cat("  ", strrep("-", 60), "\n", sep = "")

  tr_fits <- setNames(
    lapply(genotypes, function(g) {
      d   <- filter(day1, Genotype == g) |> select(TR, VPD, Rep)
      fit <- fit_all_models(d, response = "TR")
      if (!is.null(fit))
        cat(sprintf("  %-14s  %7.1f  %8.1f  %12.1f  %s\n",
                    g, fit$AIC["Linear"], fit$AIC["Bilinear"],
                    fit$AIC["Lin. plateau"], fit$best))
      fit
    }),
    genotypes
  )

  # ── 2b. Fit gsw models ───────────────────────────────────────────────────────
  cat(sprintf("\n  gsw fits:\n  %-14s  %7s  %8s  %12s  %s\n",
              "Genotype", "Linear", "Bilinear", "Lin.plateau", "Best (AIC)"))
  cat("  ", strrep("-", 60), "\n", sep = "")

  gsw_fits <- setNames(
    lapply(geno_gsw, function(g) {
      d   <- filter(day1_gsw, Genotype == g) |>
             select(gsw, VPD, Rep) |> rename(TR = gsw)
      fit <- fit_all_models(d, response = "TR")
      if (!is.null(fit))
        cat(sprintf("  %-14s  %7.1f  %8.1f  %12.1f  %s\n",
                    g, fit$AIC["Linear"], fit$AIC["Bilinear"],
                    fit$AIC["Lin. plateau"], fit$best))
      fit
    }),
    geno_gsw
  )

  # ── 3. Predict at target VPDs ────────────────────────────────────────────────
  # Helper: evaluate fitted model at two VPD values
  point_pred <- function(fits_list, geno_list, vpd_targets) {
    bind_rows(lapply(geno_list, function(g) {
      fit  <- fits_list[[g]]
      if (is.null(fit)) return(NULL)
      pars <- extract_params(fit, fit$best)
      if (is.null(pars)) return(NULL)
      vals <- switch(fit$best,
        "Linear"       = (pars[["a"]] %||% pars[["(Intercept)"]]) + pars[["VPD"]] * vpd_targets,
        "Bilinear"     = bilinear_fn(vpd_targets, pars$a, pars$b1, pars$b2, pars$bp),
        "Lin. plateau" = lin_plateau_fn(vpd_targets, pars$a, pars$b, pars$bp),
        rep(NA_real_, length(vpd_targets))
      )
      tibble(Genotype = g, VPD_target = vpd_targets, value = pmax(vals, 0))
    }))
  }

  tr_pred  <- point_pred(tr_fits,  genotypes, VPD_TARGETS)
  gsw_pred <- point_pred(gsw_fits, geno_gsw,  VPD_TARGETS)

  cat("\n  Predicted TR and gsw:\n")
  print(left_join(tr_pred,  day1 |> distinct(Genotype, WUE_class), by = "Genotype"))

  # ── 4. Prepare panel data ────────────────────────────────────────────────────
  make_panel_df <- function(pred_df, day_df) {
    pred_df |>
      left_join(day_df |> distinct(Genotype, WUE_class), by = "Genotype") |>
      mutate(
        VPD_label   = factor(paste0(VPD_target, " kPa"),
                             levels = c("1.5 kPa", "3.5 kPa")),
        WUE_short   = wue_short[as.character(WUE_class)],
        Panel_label = paste0(Genotype, "\n[", WUE_short, "]"),
        wue_rank    = match(as.character(WUE_class), WUE_ORDER)
      ) |>
      mutate(Panel_label = factor(Panel_label,
                                  levels = unique(Panel_label[order(wue_rank, Genotype)])))
  }

  tr_plot  <- make_panel_df(tr_pred,  day1)
  gsw_plot <- make_panel_df(gsw_pred, day1_gsw)

  # Strip colours: look up from panel label (last word in brackets)
  strip_colours <- function(df) {
    lvls   <- levels(df$Panel_label)
    wshort <- gsub("^.*\\[(.*)\\]$", "\\1", lvls)
    setNames(WUE_STRIP_COLS[wshort], lvls)
  }

  # ── 5. Build upper (gsw) and lower (TR) panels ───────────────────────────────
  p_gsw <- ggplot(gsw_plot, aes(x = VPD_label, y = value, fill = VPD_label)) +
    geom_col(width = 0.6, colour = "white", linewidth = 0.4) +
    scale_fill_manual(values = VPD_COLS, name = "VPD") +
    facet_wrap2(~Panel_label, ncol = 5, scales = "fixed",
                strip = strip_themed(
                  background_x = elem_list_rect(fill = unname(strip_colours(gsw_plot))))) +
    coord_cartesian(ylim = c(0, max(gsw_plot$value, na.rm = TRUE) * 1.15)) +
    labs(x = NULL, y = expression(g[sw]~(mol~m^{-2}~s^{-1}))) +
    UA_THEME +
    theme(legend.position = "top",
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          strip.text   = element_text(colour = "white", face = "bold",
                                      size = 6.5, lineheight = 1.1))

  p_tr <- ggplot(tr_plot, aes(x = VPD_label, y = value, fill = VPD_label)) +
    geom_col(width = 0.6, colour = "white", linewidth = 0.4) +
    scale_fill_manual(values = VPD_COLS, name = "VPD") +
    facet_wrap2(~Panel_label, ncol = 5, scales = "fixed",
                strip = strip_themed(
                  background_x = elem_list_rect(fill = unname(strip_colours(tr_plot))))) +
    coord_cartesian(ylim = c(0, 110)) +
    labs(x = "VPD treatment",
         y = expression(TR~(mg~m^{-2}~s^{-1}))) +
    UA_THEME +
    theme(legend.position = "none",
          strip.text = element_text(colour = "white", face = "bold",
                                    size = 6.5, lineheight = 1.1))

  p5 <- p_gsw / p_tr + plot_layout(heights = c(1, 1.4))

  ggsave(file.path(OUTDIR, "Figure5_predicted_TR_gsw.tiff"), p5,
         width = 13, height = 9, dpi = DPI, device = tiff_cairo)
  cat("\n  Saved: Figure5_predicted_TR_gsw.tiff\n\n")
})

cat("=== All 5 figures complete. Files saved to:", normalizePath(OUTDIR), "===\n")
