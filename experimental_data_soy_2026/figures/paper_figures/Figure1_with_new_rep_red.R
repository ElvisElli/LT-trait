# ==============================================================================
# Figure 1 & Figure 1B: VPD Response with NEW REPLICATION in RED
# ==============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(forcats)
  library(tibble)
  library(nlme)
})

DPI <- 150  # Quick preview
OUTDIR <- "."

options(bitmapType = "cairo")
tiff_cairo <- function(filename, width, height, units = "in", res, ...) {
  tiff(filename, width = width, height = height, units = units,
       res = res, compression = "lzw", type = "cairo")
}

# Theme
UA_THEME <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 0.7, fill = NA),
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.7),
    strip.text = element_text(colour = "black", face = "plain", size = 14),
    axis.title = element_text(size = 14, face = "plain"),
    axis.text = element_text(size = 14, face = "plain", colour = "black"),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    legend.key.size = unit(0.45, "cm"),
    legend.text = element_text(size = 8, face = "plain"),
    legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.5)
  )

# Model functions
bilinear_fn <- function(x, a, b1, b2, bp) {
  ifelse(x <= bp, a + b1 * x, a + b1 * bp + b2 * (x - bp))
}

lin_plateau_fn <- function(x, a, b, bp) {
  ifelse(x <= bp, a + b * x, a + b * bp)
}

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a[1])) a else b

# Fit all models
fit_all_models <- function(data, response = "TR") {
  y <- data[[response]]
  vpd <- data$VPD
  rep <- data$Rep
  n <- sum(!is.na(y) & !is.na(vpd))
  if (n < 4) return(NULL)
  
  vmin <- min(vpd, na.rm = TRUE)
  vmax <- max(vpd, na.rm = TRUE)
  bp0 <- (vmin + vmax) / 2
  sl0 <- coef(lm(y ~ vpd, na.action = na.omit))["vpd"]
  if (is.na(sl0)) sl0 <- 1
  a0 <- mean(y, na.rm = TRUE) - sl0 * mean(vpd, na.rm = TRUE)
  
  df_fit <- data.frame(y = y, VPD = vpd, Rep = rep) |> na.omit()
  out <- list()
  
  out[["Linear"]] <- tryCatch(
    lme(y ~ VPD, random = ~1 | Rep, data = df_fit, method = "ML",
        control = lmeControl(opt = "optim")),
    error = function(e) tryCatch(lm(y ~ VPD, data = df_fit),
                                 error = function(e2) NULL)
  )
  
  nls_bil <- tryCatch(
    nls(y ~ bilinear_fn(VPD, a, b1, b2, bp),
        data = df_fit,
        start = list(a = a0, b1 = sl0, b2 = sl0 * 0.1, bp = bp0),
        lower = c(a = -Inf, b1 = 0, b2 = -Inf, bp = vmin),
        upper = c(a = Inf, b1 = Inf, b2 = Inf, bp = vmax),
        algorithm = "port",
        control = nls.control(maxiter = 300, tol = 1e-5)),
    error = function(e) NULL
  )
  if (!is.null(nls_bil)) {
    sv <- as.list(coef(nls_bil))
    out[["Bilinear"]] <- tryCatch(
      nlme(y ~ bilinear_fn(VPD, a, b1, b2, bp),
           fixed = a + b1 + b2 + bp ~ 1,
           random = a ~ 1 | Rep,
           start = c(a = sv$a, b1 = sv$b1, b2 = sv$b2, bp = sv$bp),
           data = df_fit,
           control = nlmeControl(maxIter = 200, pnlsTol = 0.05, opt = "nlminb")),
      error = function(e) nls_bil
    )
  }
  
  nls_lp <- tryCatch(
    nls(y ~ lin_plateau_fn(VPD, a, b, bp),
        data = df_fit,
        start = list(a = a0, b = sl0, bp = bp0),
        lower = c(a = -Inf, b = 0, bp = vmin),
        upper = c(a = Inf, b = Inf, bp = vmax),
        algorithm = "port",
        control = nls.control(maxiter = 300, tol = 1e-5)),
    error = function(e) NULL
  )
  if (!is.null(nls_lp)) {
    sv_lp <- as.list(coef(nls_lp))
    out[["Lin. plateau"]] <- tryCatch(
      nlme(y ~ lin_plateau_fn(VPD, a, b, bp),
           fixed = a + b + bp ~ 1,
           random = a ~ 1 | Rep,
           start = c(a = sv_lp$a, b = sv_lp$b, bp = sv_lp$bp),
           data = df_fit,
           control = nlmeControl(maxIter = 200, pnlsTol = 0.05, opt = "nlminb")),
      error = function(e) nls_lp
    )
  }
  
  get_aic <- function(m) tryCatch(AIC(m), error = function(e) NA_real_)
  aics <- c(
    "Linear" = get_aic(out[["Linear"]]),
    "Bilinear" = get_aic(out[["Bilinear"]]),
    "Lin. plateau" = get_aic(out[["Lin. plateau"]])
  )
  out$AIC <- aics
  out$best <- names(which.min(aics))
  out
}

extract_params <- function(fit, model_name) {
  m <- fit[[model_name]]
  if (is.null(m)) return(NULL)
  tryCatch(
    as.list(if (inherits(m, "lme") || inherits(m, "nlme")) fixef(m) else coef(m)),
    error = function(e) NULL
  )
}

predict_curve <- function(fit, model_name, vpd_seq) {
  pars <- extract_params(fit, model_name)
  if (is.null(pars)) return(NULL)
  a <- pars[["a"]] %||% pars[["(Intercept)"]] %||% NA_real_
  y <- switch(model_name,
    "Linear" = a + pars[["VPD"]] * vpd_seq,
    "Bilinear" = bilinear_fn(vpd_seq, pars$a, pars$b1, pars$b2, pars$bp),
    "Lin. plateau" = lin_plateau_fn(vpd_seq, pars$a, pars$b, pars$bp),
    rep(NA_real_, length(vpd_seq))
  )
  data.frame(VPD = vpd_seq, y_hat = pmax(y, 0))
}

# Load NEW data
cat("Loading updated data with new replication...\n")
raw <- read_excel("../data/data_all_9_10_2026.xlsx", sheet = "dataAll")
df <- raw |>
  select(where(~ !all(is.na(.)))) |>
  rename(Rep = RepNo, Day = day, TR = TR_mg_sqm_seg, VPD = avg_vpd) |>
  mutate(
    Rep = factor(Rep),
    Day_label = factor(ifelse(Day == 1, "Day 1", "Day 2"), levels = c("Day 1", "Day 2")),
    Genotype = factor(Genotype),
    WUE_class = factor(WUE_class),
    IsNewRep = (Rep == 5)  # Mark new replication
  )

day1 <- filter(df, Day_label == "Day 1")
day1 <- filter(day1, !(Genotype %in% c("PI398272", "PI404191", "PI442012A", "PI416937", "AG43XF2")))
genotypes <- levels(droplevels(day1$Genotype))
VPD_SEQ <- seq(min(day1$VPD, na.rm = TRUE) * 0.95, max(day1$VPD, na.rm = TRUE) * 1.05, length.out = 300)

# ========== FIGURE 1: 9-PANEL (SELECTED GENOTYPES) ==========
cat("\n=== FIGURE 1: 9-PANEL WITH NEW REP IN RED ===\n")

tr_fits <- setNames(
  lapply(genotypes, function(g) {
    d <- filter(day1, Genotype == g) |> select(TR, VPD, Rep)
    fit_all_models(d, response = "TR")
  }),
  genotypes
)

best_params <- bind_rows(lapply(genotypes, function(g) {
  fit <- tr_fits[[g]]
  if (is.null(fit)) return(NULL)
  pars <- extract_params(fit, fit$best)
  wue <- as.character(day1 |> filter(Genotype == g) |> pull(WUE_class) |> unique())
  tibble(
    Genotype = g,
    WUE_class = wue,
    Best_model = fit$best,
    BP_kPa = if (!is.null(pars[["bp"]])) round(pars[["bp"]], 2) else NA_real_
  )
}))

curves <- bind_rows(lapply(genotypes, function(g) {
  fit <- tr_fits[[g]]
  if (is.null(fit)) return(NULL)
  d <- predict_curve(fit, fit$best, VPD_SEQ)
  if (is.null(d)) return(NULL)
  d |> mutate(Genotype = g)
})) |>
  left_join(best_params |> select(Genotype, WUE_class), by = "Genotype")

# Custom order
custom_order <- c("PI398492", "PI398617", "PI424231", "PI603418D",
                  "PI377574", "PI603171", "PI398198", "A5959", "PI471938")

day1_panel <- day1 |>
  filter(Genotype %in% custom_order) |>
  mutate(
    sort_order = match(Genotype, custom_order),
    Panel_label = Genotype
  ) |>
  arrange(sort_order)

curves_panel <- curves |>
  filter(Genotype %in% custom_order) |>
  mutate(
    sort_order = match(Genotype, custom_order),
    Panel_label = Genotype
  ) |>
  arrange(sort_order)

lvls <- custom_order
day1_panel <- mutate(day1_panel, Panel_label = factor(Panel_label, levels = lvls))
curves_panel <- mutate(curves_panel, Panel_label = factor(Panel_label, levels = lvls))

# Plot with NEW points in RED
p1 <- ggplot(day1_panel, aes(x = VPD, y = TR)) +
  # OLD points (Rep 1-4) in green
  geom_point(data = filter(day1_panel, IsNewRep == FALSE),
             size = 2.5, shape = 21, colour = "black", fill = "#31a354", stroke = 0.8, alpha = 0.9) +
  # NEW points (Rep 5) in RED
  geom_point(data = filter(day1_panel, IsNewRep == TRUE),
             size = 2.5, shape = 21, colour = "black", fill = "#e41a1c", stroke = 0.8, alpha = 0.9) +
  # Curves (using old data)
  geom_line(data = curves_panel, aes(y = y_hat),
            linewidth = 1.0, colour = "#31a354", inherit.aes = FALSE) +
  facet_wrap(~Panel_label, ncol = 3, scales = "fixed") +
  coord_cartesian(ylim = c(0, 120)) +
  labs(x = "Vapor Pressure Deficit (kPa)",
       y = expression(Transpiration~(mg~m^{-2}~s^{-1})),
       title = "Figure 1: TR Response to VPD (9 selected genotypes)\nOld data (Rep 1-4) = Green | New data (Rep 5) = RED") +
  UA_THEME +
  theme(legend.position = "none")

ggsave("Figure1_9panel_NEW_REP_RED.png", p1, width = 12, height = 10, dpi = DPI)
cat("Saved: Figure1_9panel_NEW_REP_RED.png\n")

# ========== FIGURE 1B: ALL 14 GENOTYPES ==========
cat("\n=== FIGURE 1B: ALL 14 GENOTYPES ===\n")

day1_all <- day1 |>
  mutate(Panel_label = Genotype) |>
  mutate(Panel_label = factor(Panel_label, levels = levels(droplevels(Genotype))))

curves_all <- curves |>
  filter(Genotype %in% levels(droplevels(day1_all$Genotype))) |>
  mutate(Panel_label = Genotype) |>
  mutate(Panel_label = factor(Panel_label, levels = levels(droplevels(day1_all$Genotype))))

p2 <- ggplot(day1_all, aes(x = VPD, y = TR)) +
  geom_point(data = filter(day1_all, IsNewRep == FALSE),
             size = 2.5, shape = 21, colour = "black", fill = "#31a354", stroke = 0.8, alpha = 0.9) +
  geom_point(data = filter(day1_all, IsNewRep == TRUE),
             size = 2.5, shape = 21, colour = "black", fill = "#e41a1c", stroke = 0.8, alpha = 0.9) +
  geom_line(data = curves_all, aes(y = y_hat),
            linewidth = 1.0, colour = "#31a354", inherit.aes = FALSE) +
  facet_wrap(~Panel_label, ncol = 7, scales = "fixed") +
  coord_cartesian(ylim = c(0, 120)) +
  labs(x = "Vapor Pressure Deficit (kPa)",
       y = expression(Transpiration~(mg~m^{-2}~s^{-1})),
       title = "Figure 1B: TR Response to VPD (ALL 14 genotypes)\nOld data (Rep 1-4) = Green | New data (Rep 5) = RED") +
  UA_THEME +
  theme(legend.position = "none", strip.text = element_text(size = 10))

ggsave("Figure1B_14genotypes_NEW_REP_RED.png", p2, width = 16, height = 6, dpi = DPI)
cat("Saved: Figure1B_14genotypes_NEW_REP_RED.png\n")

cat("\n✓ Figures regenerated with new replication in RED\n")
quit(save = "no")
