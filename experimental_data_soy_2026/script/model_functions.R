## model_functions.R
## Shared model functions for the GC 2026 limited-transpiration analysis:
## three candidate TR / gsw ~ VPD models (linear, bilinear, linear-plateau),
## fitted per genotype with Rep as a random intercept (nlme) and compared by
## AIC. Extracted from figures/paper_figures/paper_figures.R so other scripts
## can `source("script/model_functions.R")`.
suppressPackageStartupMessages(library(nlme))

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
  # Multi-start nls: try several breakpoint starting values across the VPD
  # range and keep the converged fit with the lowest residual sum of squares
  # (a single midpoint start fails to converge for some genotypes).
  bp_starts <- unique(c(bp0, unname(quantile(vpd, c(0.25, 0.5, 0.75), na.rm = TRUE)), 1.5, 2.0, 2.5, 3.0))
  bp_starts <- bp_starts[bp_starts > vmin & bp_starts < vmax]
  best_nls  <- function(fits) {
    fits <- Filter(Negate(is.null), fits)
    if (!length(fits)) return(NULL)
    fits[[which.min(sapply(fits, deviance))]]
  }
  nls_bil <- best_nls(lapply(bp_starts, function(b0) tryCatch(
    nls(y ~ bilinear_fn(VPD, a, b1, b2, bp),
        data      = df_fit,
        start     = list(a = a0, b1 = sl0, b2 = sl0 * 0.1, bp = b0),
        lower     = c(a = -Inf, b1 = 0,   b2 = -Inf, bp = vmin),
        upper     = c(a =  Inf, b1 = Inf, b2 =  Inf, bp = vmax),
        algorithm = "port",
        control   = nls.control(maxiter = 300, tol = 1e-5)),
    error = function(e) NULL)))
  # Extra bilinear start taken from a linear-plateau fit (second slope = 0):
  # the bilinear surface is flat near the plateau solution, which is where the
  # port algorithm otherwise fails to converge for some genotypes.
  nls_lp0 <- best_nls(lapply(bp_starts, function(b0) tryCatch(
    nls(y ~ lin_plateau_fn(VPD, a, b, bp), data = df_fit,
        start = list(a = a0, b = sl0, bp = b0),
        lower = c(a = -Inf, b = 0, bp = vmin), upper = c(a = Inf, b = Inf, bp = vmax),
        algorithm = "port", control = nls.control(maxiter = 300, tol = 1e-5)),
    error = function(e) NULL)))
  if (!is.null(nls_lp0)) {
    sv0 <- as.list(coef(nls_lp0))
    nls_bil <- best_nls(list(nls_bil, tryCatch(
      nls(y ~ bilinear_fn(VPD, a, b1, b2, bp), data = df_fit,
          start = list(a = sv0$a, b1 = sv0$b, b2 = 0, bp = sv0$bp),
          lower = c(a = -Inf, b1 = 0, b2 = -Inf, bp = vmin),
          upper = c(a = Inf, b1 = Inf, b2 = Inf, bp = vmax),
          algorithm = "port", control = nls.control(maxiter = 300, tol = 1e-5)),
      error = function(e) NULL)))
  }
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
  nls_lp <- best_nls(lapply(bp_starts, function(b0) tryCatch(
    nls(y ~ lin_plateau_fn(VPD, a, b, bp),
        data      = df_fit,
        start     = list(a = a0, b = sl0, bp = b0),
        lower     = c(a = -Inf, b = 0,   bp = vmin),
        upper     = c(a =  Inf, b = Inf, bp = vmax),
        algorithm = "port",
        control   = nls.control(maxiter = 300, tol = 1e-5)),
    error = function(e) NULL)))
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
