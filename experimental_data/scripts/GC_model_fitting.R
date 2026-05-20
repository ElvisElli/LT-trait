## Growth chamber transpiration: linear, bilinear, and linear-plateau mixed effects models
## Fits all three model types per genotype × experiment using lme4 / nlme
## Output: model comparison table + fitted-curve figures

setwd("/home/user/LT-trait")
library(dplyr); library(tidyr); library(ggplot2); library(readxl)
library(lme4); library(lmerTest); library(nlme)
library(segmented)

dir.create("experimental_data/GC_analysis", recursive = TRUE, showWarnings = FALSE)

## ── Load combined clean data ──────────────────────────────────────────────────
d <- read.csv("experimental_data/GC_analysis/GC_combined_clean.csv") %>%
  filter(!is.na(RelTranspiration), !is.na(VPD)) %>%
  mutate(
    Genotype    = factor(Genotype),
    Experiment  = factor(Experiment),
    Replication = factor(Replication),
    Day         = factor(Day)
  )

## Use Day 2 only — cleaner signal (Day 1 shows higher variability / acclimation)
d2 <- filter(d, Day == 2)

cat("Observations (Day 2):", nrow(d2), "\n")
cat("Genotypes:", paste(levels(d2$Genotype), collapse = ", "), "\n\n")

## ── Helper: fit all three model types for one genotype subset ─────────────────
fit_models <- function(df, geno_name, exp_name) {

  ## ── 1. Linear mixed model ─────────────────────────────────────────────────
  lm_fit <- tryCatch(
    lmer(RelTranspiration ~ VPD + (1 | Replication), data = df, REML = FALSE),
    error = function(e) NULL
  )

  ## ── 2. Bilinear (piecewise) via segmented + lme via two-phase lm per rep ─
  ## Fit ordinary lm to get a segmented object, then extract breakpoint
  lm_base  <- lm(RelTranspiration ~ VPD, data = df)
  seg_fit  <- tryCatch(
    segmented(lm_base, seg.Z = ~VPD,
              psi = list(VPD = 2.0),
              control = seg.control(it.max = 100, tol = 1e-4)),
    error = function(e) NULL
  )

  ## Build bilinear predictor using estimated breakpoint (or fixed 2.0 if failed)
  bp <- if (!is.null(seg_fit)) seg_fit$psi[, "Est."] else 2.0
  bp <- max(min(bp, max(df$VPD) - 0.1), min(df$VPD) + 0.1)  # keep in range

  df <- df %>% mutate(
    VPD_lo  = pmin(VPD, bp),
    VPD_hi  = pmax(VPD - bp, 0)
  )
  bilin_fit <- tryCatch(
    lmer(RelTranspiration ~ VPD_lo + VPD_hi + (1 | Replication),
         data = df, REML = FALSE),
    error = function(e) NULL
  )

  ## ── 3. Linear plateau (nlme with self-start or grid search) ───────────────
  ## Model: T = a + b*VPD for VPD < plateau_vpd; T = a + b*plateau_vpd above
  ## Parameterised as nlme with fixed plateau_vpd searched by grid
  lp_grid <- tryCatch({
    bp_grid <- seq(min(df$VPD) + 0.2, max(df$VPD) - 0.2, by = 0.1)
    aic_vec <- sapply(bp_grid, function(bp_try) {
      df_g <- df %>% mutate(VPD_plateau = pmin(VPD, bp_try))
      m <- tryCatch(
        lmer(RelTranspiration ~ VPD_plateau + (1 | Replication),
             data = df_g, REML = FALSE),
        error = function(e) NULL
      )
      if (is.null(m)) Inf else AIC(m)
    })
    best_bp   <- bp_grid[which.min(aic_vec)]
    df_lp     <- df %>% mutate(VPD_plateau = pmin(VPD, best_bp))
    lmer(RelTranspiration ~ VPD_plateau + (1 | Replication),
         data = df_lp, REML = FALSE)
  }, error = function(e) NULL)

  ## ── Extract metrics ────────────────────────────────────────────────────────
  get_r2 <- function(m, df_sub) {
    if (is.null(m)) return(NA_real_)
    cor(fitted(m), df_sub$RelTranspiration[!is.na(df_sub$RelTranspiration)])^2
  }

  lp_bp <- if (!is.null(lp_grid)) {
    bp_grid2 <- seq(min(df$VPD) + 0.2, max(df$VPD) - 0.2, by = 0.1)
    aic_vec2 <- sapply(bp_grid2, function(bp_try) {
      df_g <- df %>% mutate(VPD_plateau = pmin(VPD, bp_try))
      m <- tryCatch(
        lmer(RelTranspiration ~ VPD_plateau + (1 | Replication),
             data = df_g, REML = FALSE),
        error = function(e) NULL
      )
      if (is.null(m)) Inf else AIC(m)
    })
    bp_grid2[which.min(aic_vec2)]
  } else NA_real_

  tibble(
    Experiment    = exp_name,
    Genotype      = geno_name,
    ## Linear
    Lin_slope     = if (!is.null(lm_fit)) fixef(lm_fit)["VPD"] else NA_real_,
    Lin_AIC       = if (!is.null(lm_fit)) AIC(lm_fit)          else NA_real_,
    Lin_R2        = get_r2(lm_fit, df),
    ## Bilinear
    Bilin_BP      = bp,
    Bilin_slope_lo = if (!is.null(bilin_fit)) fixef(bilin_fit)["VPD_lo"] else NA_real_,
    Bilin_slope_hi = if (!is.null(bilin_fit)) fixef(bilin_fit)["VPD_hi"] else NA_real_,
    Bilin_AIC     = if (!is.null(bilin_fit)) AIC(bilin_fit)    else NA_real_,
    Bilin_R2      = get_r2(bilin_fit, df),
    ## Linear plateau
    LP_BP         = lp_bp,
    LP_slope      = if (!is.null(lp_grid)) fixef(lp_grid)["VPD_plateau"] else NA_real_,
    LP_AIC        = if (!is.null(lp_grid)) AIC(lp_grid)                   else NA_real_,
    LP_R2         = get_r2(lp_grid, df),
    ## Best model by AIC
    Best_model    = {
      aics <- c(Linear = if (!is.null(lm_fit)) AIC(lm_fit) else Inf,
                Bilinear = if (!is.null(bilin_fit)) AIC(bilin_fit) else Inf,
                `Linear-Plateau` = if (!is.null(lp_grid)) AIC(lp_grid) else Inf)
      names(aics)[which.min(aics)]
    },
    dAIC_best     = {
      aics <- c(if (!is.null(lm_fit)) AIC(lm_fit) else Inf,
                if (!is.null(bilin_fit)) AIC(bilin_fit) else Inf,
                if (!is.null(lp_grid)) AIC(lp_grid) else Inf)
      min_a <- min(aics); sort(aics - min_a)[2]
    }
  )
}

## ── Run per genotype × experiment ─────────────────────────────────────────────
cat("Fitting models per genotype × experiment...\n")
results <- d2 %>%
  group_by(Experiment, Genotype) %>%
  group_modify(~ fit_models(.x, unique(.x$Genotype), unique(.x$Experiment))) %>%
  ungroup()

print(as.data.frame(results[, c("Experiment","Genotype","Best_model","Lin_R2","Bilin_R2",
                                 "LP_R2","Bilin_BP","LP_BP","dAIC_best")]), digits = 3)

write.csv(results, "experimental_data/GC_analysis/GC_model_comparison.csv",
          row.names = FALSE)
cat("\nSaved: GC_model_comparison.csv\n")

## ── Figure: fitted curves overlaid on data ───────────────────────────────────
## Generate smooth prediction lines for each model × genotype
pred_grid <- d2 %>%
  group_by(Experiment, Genotype, Expected_WUE) %>%
  summarise(VPD_min = min(VPD), VPD_max = max(VPD), .groups = "drop") %>%
  rowwise() %>%
  mutate(VPD_seq = list(seq(VPD_min, VPD_max, length.out = 60))) %>%
  unnest(VPD_seq) %>%
  rename(VPD = VPD_seq) %>%
  dplyr::select(-VPD_min, -VPD_max)

## Build prediction lines using base-R row lookup to avoid dplyr/MASS conflicts
preds_list <- lapply(seq_len(nrow(results)), function(i) {
  r      <- results[i, ]
  g      <- as.character(r$Genotype)
  e      <- as.character(r$Experiment)
  d2sub  <- d2[as.character(d2$Genotype) == g & as.character(d2$Experiment) == e, ]
  vpd_rng <- seq(min(d2sub$VPD), max(d2sub$VPD), length.out = 60)
  ewue   <- unique(d2sub$Expected_WUE)[1]

  ## Intercepts anchored at group centroid
  mu_t   <- mean(d2sub$RelTranspiration)
  mu_v   <- mean(d2sub$VPD)
  bp_b   <- r$Bilin_BP;  bp_lp <- r$LP_BP

  data.frame(
    Experiment    = e,
    Genotype      = g,
    Expected_WUE  = ewue,
    VPD           = vpd_rng,
    Best_model    = r$Best_model,
    pred_Linear   = mu_t + r$Lin_slope * (vpd_rng - mu_v),
    pred_Bilinear = {
      mu_lo <- mean(pmin(d2sub$VPD, bp_b))
      mu_hi <- mean(pmax(d2sub$VPD - bp_b, 0))
      int_b <- mu_t - r$Bilin_slope_lo * mu_lo - r$Bilin_slope_hi * mu_hi
      int_b + r$Bilin_slope_lo * pmin(vpd_rng, bp_b) +
              r$Bilin_slope_hi * pmax(vpd_rng - bp_b, 0)
    },
    pred_LP = {
      mu_lp <- mean(pmin(d2sub$VPD, bp_lp))
      int_lp <- mu_t - r$LP_slope * mu_lp
      int_lp + r$LP_slope * pmin(vpd_rng, bp_lp)
    }
  )
})

preds <- do.call(rbind, preds_list) %>%
  tidyr::pivot_longer(cols = starts_with("pred_"),
               names_to = "Model", values_to = "pred",
               names_prefix = "pred_") %>%
  mutate(Model = factor(Model, levels = c("Linear","Bilinear","LP"),
                        labels  = c("Linear","Bilinear","Linear-Plateau")))

## Facet label: Genotype + WUE quadrant
d2 <- d2 %>%
  mutate(facet_label = paste0(Genotype, "\n", Expected_WUE))
preds <- preds %>%
  mutate(facet_label = paste0(Genotype, "\n", Expected_WUE))

model_colours <- c("Linear" = "#636363",
                   "Bilinear" = "#2171b5",
                   "Linear-Plateau" = "#d95f02")

p_models <- ggplot() +
  geom_vline(xintercept = 2.0, linetype = "dashed",
             colour = "grey70", linewidth = 0.6) +
  geom_point(data = d2,
             aes(x = VPD, y = RelTranspiration),
             alpha = 0.45, size = 1.8, colour = "grey30") +
  geom_line(data = preds,
            aes(x = VPD, y = pred, colour = Model, linewidth = Model)) +
  facet_wrap(Experiment ~ facet_label, nrow = 2, scales = "free_x") +
  scale_colour_manual(values = model_colours, name = "Model") +
  scale_linewidth_manual(values = c("Linear" = 0.7, "Bilinear" = 1.1,
                                     "Linear-Plateau" = 1.1), guide = "none") +
  scale_y_continuous(limits = c(0.2, 1.1), breaks = seq(0.2, 1.0, 0.2)) +
  labs(x     = "Vapour pressure deficit (kPa)",
       y     = "Relative transpiration",
       title = "Model comparison: linear, bilinear, and linear-plateau (Day 2, mixed effects)",
       subtitle = "Dashed line = VPD 2.0 kPa reference  ·  Points = individual observations") +
  theme_bw(base_size = 11) +
  theme(legend.position    = "bottom",
        strip.background   = element_rect(fill = "grey92", colour = NA),
        strip.text         = element_text(face = "bold", size = 8),
        panel.grid.minor   = element_blank(),
        plot.title         = element_text(face = "bold", size = 11),
        plot.subtitle      = element_text(colour = "grey40", size = 8.5))

ggsave("experimental_data/GC_analysis/GC_Fig_model_fits.png",
       plot = p_models, width = 26, height = 14, units = "cm", dpi = 150, bg = "white")
ggsave("experimental_data/GC_analysis/GC_Fig_model_fits.tiff",
       plot = p_models, width = 26, height = 14, units = "cm", dpi = 300,
       compression = "lzw", bg = "white")

## ── AIC comparison heatmap ───────────────────────────────────────────────────
wue_key <- d2 %>% dplyr::select(Experiment, Genotype, Expected_WUE) %>% distinct()
aic_long <- results %>%
  left_join(wue_key, by = c("Experiment","Genotype")) %>%
  dplyr::select(Experiment, Genotype, Expected_WUE, Lin_AIC, Bilin_AIC, LP_AIC) %>%
  pivot_longer(cols = ends_with("_AIC"), names_to = "Model", values_to = "AIC") %>%
  mutate(Model = recode(Model,
    "Lin_AIC"   = "Linear",
    "Bilin_AIC" = "Bilinear",
    "LP_AIC"    = "Linear-Plateau")) %>%
  group_by(Experiment, Genotype) %>%
  mutate(delta_AIC = AIC - min(AIC, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(label = paste0(Genotype, "\n(", sub(",.*", "", Expected_WUE), ")"))

p_aic <- ggplot(aic_long,
                aes(x = Model, y = label, fill = delta_AIC)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = round(delta_AIC, 1)), size = 3, fontface = "bold") +
  facet_wrap(~ Experiment, nrow = 1) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#d73027",
                       midpoint = 2, name = "dAIC from\nbest model",
                       limits = c(0, NA)) +
  labs(x     = NULL, y = NULL,
       title = "Model fit comparison: dAIC (lower = better)",
       subtitle = "Blue = best fit · Red = poorest fit · Values = dAIC from best model per genotype") +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text       = element_text(face = "bold"),
        panel.grid       = element_blank(),
        plot.title       = element_text(face = "bold"),
        plot.subtitle    = element_text(colour = "grey40", size = 8.5),
        axis.text.y      = element_text(size = 8))

ggsave("experimental_data/GC_analysis/GC_Fig_AIC_heatmap.png",
       plot = p_aic, width = 20, height = 14, units = "cm", dpi = 150, bg = "white")
ggsave("experimental_data/GC_analysis/GC_Fig_AIC_heatmap.tiff",
       plot = p_aic, width = 20, height = 14, units = "cm", dpi = 300,
       compression = "lzw", bg = "white")

cat("Saved: GC_Fig_model_fits, GC_Fig_AIC_heatmap\n")

## ── Print final summary table for report ──────────────────────────────────────
cat("\n====== FINAL MODEL SUMMARY TABLE ======\n")
tab_final <- results %>% mutate(across(where(is.numeric), ~ round(.x, 3)))
print(as.data.frame(tab_final[, c("Experiment","Genotype","Best_model","Lin_slope","Lin_R2",
                                   "Bilin_BP","Bilin_slope_lo","Bilin_slope_hi","Bilin_R2",
                                   "LP_BP","LP_slope","LP_R2","dAIC_best")]), width = 120)
