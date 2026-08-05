# ==============================================================================
# paper_figures.R
# Soybean Limited-Transpiration Experiment (GC 2026) — Paper Figures
#
# Run from this file's directory:
#   setwd("<path>/figures/paper_figures")
#   source("paper_figures.R")
#
# Each figure section is fully self-contained (local scope via local({})).
# Changing one figure does not affect the others.
#
# Output : TIFF files in this directory
# Data   : data_all.xlsx (included in this directory)
# DPI    : set below — 6000 for submission, 600 for quick preview
# ==============================================================================

DPI    <- 600           # 600 DPI for submission
OUTDIR <- "."           # TIFFs written to this directory

# Headless-safe bitmap device (needed on Linux without X11)
options(bitmapType = "cairo")

# Cairo-based TIFF device function for ggsave (headless-safe, LZW compressed).
# ggplot2 4.x maps ggsave's `dpi` → `res` only when the formals contain "res".
tiff_cairo <- function(filename, width, height, units = "in", res, ...) {
  tiff(filename, width = width, height = height, units = units,
       res = res, compression = "lzw", type = "cairo")
}

# --------------------------------------------------------------------------- #
# Shared packages — loaded once, available to all local() scopes
# --------------------------------------------------------------------------- #
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
  library(ggrepel)
  library(ggh4x)
})

# --------------------------------------------------------------------------- #
# Shared colour / theme constants — visible inside every local() scope
# --------------------------------------------------------------------------- #
UA_RED   <- "#9D2235"
UA_WINE  <- "#6E1424"
UA_BLUE  <- "#4393c3"
UA_GRAY  <- "#7A7A7A"
UA_DARK  <- "#1A1A1A"

WUE_ORDER <- c("slow_wilting", "high_wilting",
                "Plastic, High WUE", "Nonplastic, Low WUE", "MAGIC parent")

WUE_COLS <- c(
  "slow_wilting"        = UA_WINE,
  "high_wilting"        = "#C8102E",
  "Plastic, High WUE"   = UA_BLUE,
  "Nonplastic, Low WUE" = "#d4824a",
  "MAGIC parent"        = "#333333"
)

UA_THEME <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(linewidth = 0.3, colour = "grey90"),
    strip.background  = element_rect(fill = UA_RED, colour = NA),
    strip.text        = element_text(colour = "white", face = "bold", size = 7),
    axis.title        = element_text(size = 9),
    plot.title        = element_text(face = "bold", size = 10, colour = UA_DARK),
    legend.key.size   = unit(0.45, "cm"),
    legend.text       = element_text(size = 8)
  )

# Curve prediction from stored parameters (no re-fitting needed)
pred_from_params <- function(vpd_seq, best_model, a, b1, bp, b2) {
  if (best_model == "Linear") {
    a + b1 * vpd_seq
  } else if (best_model == "Lin. plateau") {
    ifelse(vpd_seq <= bp, a + b1 * vpd_seq, a + b1 * bp)
  } else if (best_model == "Bilinear") {
    ifelse(vpd_seq <= bp, a + b1 * vpd_seq, a + b1 * bp + b2 * (vpd_seq - bp))
  } else {
    rep(NA_real_, length(vpd_seq))
  }
}

cat("paper_figures.R loaded — building 5 figures at", DPI, "DPI\n")


# ==============================================================================
# FIGURE 1
# Transpiration rate response to VPD per genotype — raw data + best-fit model
# ==============================================================================
local({
  cat("Building Figure 1 ...\n")

  # --- Data ---
  raw <- read_excel("data_all.xlsx", sheet = "dataAll")
  df <- raw |>
    select(where(~ !all(is.na(.)))) |>
    rename(Rep = RepNo, Day = day, TR = TR_mg_sqm_seg,
           VPD = avg_vpd, LeafArea = LeafArea_squarecm) |>
    mutate(
      Rep       = factor(Rep),
      Day_label = factor(ifelse(Day == 1, "Day 1", "Day 2"),
                         levels = c("Day 1", "Day 2")),
      Genotype  = factor(Genotype),
      WUE_class = factor(WUE_class, levels = WUE_ORDER)
    )
  day1 <- filter(df, Day_label == "Day 1")

  params <- read_csv("../../data/TR_best_model_summary.csv",
                     show_col_types = FALSE)

  TR_YLIM <- c(0, 110)

  # Generate prediction curves per genotype
  vpd_range <- range(day1$VPD, na.rm = TRUE)
  VPD_SEQ   <- seq(vpd_range[1] * 0.95, vpd_range[2] * 1.05, length.out = 300)

  curves <- params |>
    rowwise() |>
    mutate(
      vpd = list(VPD_SEQ),
      y   = list(pred_from_params(VPD_SEQ, Best_model,
                                  Intercept, Slope1, BP_kPa, Slope2))
    ) |>
    ungroup() |>
    select(Genotype, WUE_class, Best_model, vpd, y) |>
    unnest(cols = c(vpd, y)) |>
    filter(!is.na(y), y >= 0) |>
    mutate(WUE_class = factor(WUE_class, levels = WUE_ORDER))

  # Build model-label for strip
  model_label <- params |>
    select(Genotype, Best_model) |>
    mutate(Panel_label = paste0(Genotype, "\n[", Best_model, "]"),
           Panel_label = factor(Panel_label))

  day1_lab <- left_join(day1, model_label, by = "Genotype") |>
    mutate(Panel_label = factor(Panel_label,
                                levels = levels(model_label$Panel_label)))
  curves_lab <- left_join(curves, model_label, by = "Genotype") |>
    mutate(Panel_label = factor(Panel_label,
                                levels = levels(model_label$Panel_label)))

  p <- ggplot(day1_lab, aes(x = VPD, y = TR, colour = WUE_class)) +
    geom_point(size = 1.5, alpha = 0.65) +
    geom_line(data = curves_lab,
              aes(x = vpd, y = y, colour = WUE_class),
              linewidth = 0.9, inherit.aes = FALSE) +
    scale_colour_manual(values = WUE_COLS, name = "WUE class") +
    facet_wrap(~Panel_label, ncol = 5, scales = "fixed") +
    coord_cartesian(ylim = TR_YLIM) +
    labs(x = "VPD (kPa)",
         y = expression(TR~(mg~m^{-2}~s^{-1}))) +
    UA_THEME +
    theme(
      strip.text       = element_text(colour = "white", face = "bold",
                                      size = 6.5, lineheight = 1.0),
      legend.position  = "bottom"
    ) +
    guides(colour = guide_legend(nrow = 2,
                                 override.aes = list(size = 2.5, linewidth = 0.8)))

  ggsave(file.path(OUTDIR, "Figure1_TR_response_VPD.tiff"), p,
         width = 13, height = 8, dpi = DPI, device = tiff_cairo)
  cat("  Saved Figure 1\n")
})


# ==============================================================================
# FIGURE 2
# Genotypic variability in VPD breakpoints
# ==============================================================================
local({
  cat("Building Figure 2 ...\n")

  params <- read_csv("../../data/TR_best_model_summary.csv",
                     show_col_types = FALSE)

  bp_df <- params |>
    filter(Best_model %in% c("Bilinear", "Lin. plateau"),
           !is.na(BP_kPa)) |>
    mutate(
      WUE_class = factor(WUE_class, levels = WUE_ORDER),
      Genotype  = fct_reorder(Genotype, BP_kPa)
    )

  p <- ggplot(bp_df, aes(x = Genotype, y = BP_kPa, fill = WUE_class)) +
    geom_col(width = 0.65, colour = "white", linewidth = 0.4) +
    geom_hline(yintercept = mean(bp_df$BP_kPa), linetype = "dashed",
               colour = UA_GRAY, linewidth = 0.6) +
    annotate("text",
             x = nrow(bp_df) * 0.97, y = mean(bp_df$BP_kPa) + 0.04,
             label = sprintf("mean = %.2f kPa", mean(bp_df$BP_kPa)),
             hjust = 1, vjust = 0, size = 2.8, colour = UA_GRAY,
             fontface = "italic") +
    scale_fill_manual(values = WUE_COLS, name = "WUE class") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = NULL,
         y = "Breakpoint VPD (kPa)") +
    UA_THEME +
    theme(
      axis.text.x     = element_text(angle = 35, hjust = 1, size = 8),
      legend.position = "bottom"
    )

  ggsave(file.path(OUTDIR, "Figure2_TR_breakpoints.tiff"), p,
         width = 6, height = 4.5, dpi = DPI, device = tiff_cairo)
  cat("  Saved Figure 2\n")
})


# ==============================================================================
# FIGURE 3
# Genotype classification: intercept vs slope quadrant scatter
# ==============================================================================
local({
  cat("Building Figure 3 ...\n")

  params <- read_csv("../../data/TR_best_model_summary.csv",
                     show_col_types = FALSE)

  # Highlight only Plastic / Nonplastic; grey the rest
  fig_df <- params |>
    mutate(WUE_group = case_when(
      WUE_class == "Plastic, High WUE"   ~ "Plastic, High WUE",
      WUE_class == "Nonplastic, Low WUE" ~ "Nonplastic, Low WUE",
      TRUE                                ~ "Other"
    ))

  WUE_COLS_Q <- c(
    "Plastic, High WUE"   = UA_BLUE,
    "Nonplastic, Low WUE" = "#d4824a",
    "Other"               = "#AAAAAA"
  )

  x_mid     <- median(fig_df$Intercept, na.rm = TRUE)
  slope_mid <- median(fig_df$Slope1,    na.rm = TRUE)
  x_rng     <- range(fig_df$Intercept, na.rm = TRUE)
  y_rng     <- range(fig_df$Slope1,    na.rm = TRUE)

  quad_shade <- data.frame(
    xmin = c(-Inf,  x_mid, -Inf,  x_mid),
    xmax = c(x_mid, Inf,   x_mid, Inf),
    ymin = c(slope_mid, slope_mid, -Inf,    -Inf),
    ymax = c(Inf,       Inf,       slope_mid, slope_mid),
    fill = c("#e8f0fb", "#fde8e8", "#eafaea", "#fef6e0")
  )

  p <- ggplot(fig_df, aes(x = Intercept, y = Slope1, colour = WUE_group)) +
    geom_rect(data = quad_shade,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
              inherit.aes = FALSE, alpha = 0.55) +
    scale_fill_identity() +
    geom_vline(xintercept = x_mid,     linetype = "dashed",
               colour = UA_GRAY, linewidth = 0.65) +
    geom_hline(yintercept = slope_mid, linetype = "dashed",
               colour = UA_GRAY, linewidth = 0.65) +
    # Quadrant labels — blue quadrant at its upper-right corner
    annotate("text",
             x = x_mid - diff(x_rng) * 0.03, y = Inf,
             label = "Low intercept\nHigh sensitivity",
             hjust = 1, vjust = 2.0, size = 2.6,
             colour = "#555555", fontface = "italic", lineheight = 0.9) +
    annotate("text",
             x = x_rng[2] - diff(x_rng) * 0.05, y = Inf,
             label = "High intercept\nHigh sensitivity",
             hjust = 1, vjust = 2.0, size = 2.6,
             colour = "#555555", fontface = "italic", lineheight = 0.9) +
    annotate("text",
             x = x_rng[1] + diff(x_rng) * 0.05, y = -Inf,
             label = "Low intercept\nLow sensitivity",
             hjust = 0, vjust = -1.0, size = 2.6,
             colour = "#555555", fontface = "italic", lineheight = 0.9) +
    annotate("text",
             x = x_rng[2] - diff(x_rng) * 0.05, y = -Inf,
             label = "High intercept\nLow sensitivity",
             hjust = 1, vjust = -1.0, size = 2.6,
             colour = "#555555", fontface = "italic", lineheight = 0.9) +
    geom_point(aes(shape = Best_model), size = 4, alpha = 0.92) +
    geom_label_repel(
      aes(label = Genotype),
      size = 2.7, box.padding = 0.5, point.padding = 0.3,
      max.overlaps = 20, label.size = 0.2, alpha = 0.92,
      show.legend = FALSE, seed = 42
    ) +
    scale_colour_manual(values = WUE_COLS_Q, name = NULL,
                        breaks = c("Plastic, High WUE", "Nonplastic, Low WUE")) +
    scale_shape_manual(
      values = c("Linear" = 16, "Bilinear" = 17, "Lin. plateau" = 15),
      name = NULL
    ) +
    labs(
      x = expression(Intercept~~(mg~m^{-2}~s^{-1})),
      y = expression(Slope~~(mg~m^{-2}~s^{-1}~kPa^{-1}))
    ) +
    UA_THEME +
    theme(legend.position = "right")

  ggsave(file.path(OUTDIR, "Figure3_intercept_slope_quadrants.tiff"), p,
         width = 9, height = 6.5, dpi = DPI, device = tiff_cairo)
  cat("  Saved Figure 3\n")
})


# ==============================================================================
# FIGURE 4
# Stomatal conductance response to VPD per genotype — raw LICOR readings + model
# ==============================================================================
local({
  cat("Building Figure 4 ...\n")

  # --- Data ---
  raw <- read_excel("data_all.xlsx", sheet = "dataAll")
  df <- raw |>
    select(where(~ !all(is.na(.)))) |>
    rename(Rep = RepNo, Day = day, VPD = avg_vpd) |>
    mutate(
      Rep       = factor(Rep),
      Day_label = factor(ifelse(Day == 1, "Day 1", "Day 2"),
                         levels = c("Day 1", "Day 2")),
      Genotype  = factor(Genotype),
      WUE_class = factor(WUE_class, levels = WUE_ORDER),
      gsw_mean  = rowMeans(cbind(gsw1, gsw2, gsw3), na.rm = FALSE)
    )
  all_na <- is.na(df$gsw1) & is.na(df$gsw2) & is.na(df$gsw3)
  df$gsw_mean[all_na] <- NA

  day1     <- filter(df, Day_label == "Day 1")
  day1_gsw <- filter(day1, !is.na(gsw_mean))

  # Individual LICOR readings (3 per plant)
  day1_gsw_long <- day1_gsw |>
    pivot_longer(cols = c(gsw1, gsw2, gsw3),
                 names_to = "licor_rep", values_to = "gsw_reading") |>
    filter(!is.na(gsw_reading))

  params <- read_csv("../../data/gsw_best_model_summary.csv",
                     show_col_types = FALSE)
  geno_gsw <- unique(day1_gsw$Genotype)

  GSW_YLIM <- c(0, max(day1_gsw_long$gsw_reading, na.rm = TRUE) * 1.10)

  vpd_range <- range(day1_gsw$VPD, na.rm = TRUE)
  VPD_SEQ   <- seq(vpd_range[1] * 0.95, vpd_range[2] * 1.05, length.out = 300)

  curves <- params |>
    filter(Genotype %in% geno_gsw) |>
    rowwise() |>
    mutate(
      vpd = list(VPD_SEQ),
      y   = list(pred_from_params(VPD_SEQ, Best_model,
                                  Intercept, Slope1, BP_kPa, Slope2))
    ) |>
    ungroup() |>
    select(Genotype, WUE_class, Best_model, vpd, y) |>
    unnest(cols = c(vpd, y)) |>
    filter(!is.na(y), y >= 0) |>
    mutate(WUE_class = factor(WUE_class, levels = WUE_ORDER))

  model_label <- params |>
    filter(Genotype %in% geno_gsw) |>
    select(Genotype, Best_model) |>
    mutate(Panel_label = paste0(Genotype, "\n[", Best_model, "]"),
           Panel_label = factor(Panel_label))

  day1_gsw_long_lab <- left_join(day1_gsw_long, model_label, by = "Genotype") |>
    mutate(Panel_label = factor(Panel_label,
                                levels = levels(model_label$Panel_label)),
           WUE_class   = factor(WUE_class, levels = WUE_ORDER))
  curves_lab <- left_join(curves, model_label, by = "Genotype") |>
    mutate(Panel_label = factor(Panel_label,
                                levels = levels(model_label$Panel_label)))

  p <- ggplot(day1_gsw_long_lab,
              aes(x = VPD, y = gsw_reading, colour = WUE_class)) +
    geom_point(size = 1.5, alpha = 0.65,
               position = position_jitter(width = 0.02, seed = 1)) +
    geom_line(data = curves_lab,
              aes(x = vpd, y = y, colour = WUE_class),
              linewidth = 0.9, inherit.aes = FALSE) +
    scale_colour_manual(values = WUE_COLS, name = "WUE class") +
    facet_wrap(~Panel_label, ncol = 5, scales = "fixed") +
    coord_cartesian(ylim = GSW_YLIM) +
    labs(x = "VPD (kPa)",
         y = expression(g[sw]~(mol~m^{-2}~s^{-1}))) +
    UA_THEME +
    theme(
      strip.text      = element_text(colour = "white", face = "bold",
                                     size = 6.5, lineheight = 1.0),
      legend.position = "bottom"
    ) +
    guides(colour = guide_legend(nrow = 2,
                                 override.aes = list(size = 2.5, linewidth = 0.8)))

  ggsave(file.path(OUTDIR, "Figure4_gsw_response_VPD.tiff"), p,
         width = 12, height = 6.5, dpi = DPI, device = tiff_cairo)
  cat("  Saved Figure 4\n")
})


# ==============================================================================
# FIGURE 5
# Predicted TR and gsw at VPD 1.5 and 3.5 kPa — per genotype, WUE-class strips
# ==============================================================================
local({
  cat("Building Figure 5 ...\n")

  tr_params  <- read_csv("../../data/TR_best_model_summary.csv",
                          show_col_types = FALSE)
  gsw_params <- read_csv("../../data/gsw_best_model_summary.csv",
                          show_col_types = FALSE)

  VPD_PRED <- c(1.5, 3.5)
  VPD_COLS <- c("1.5 kPa" = UA_BLUE, "3.5 kPa" = UA_WINE)

  WUE_ORDER_F <- c("slow_wilting", "high_wilting",
                   "Plastic, High WUE", "Nonplastic, Low WUE", "MAGIC parent")

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

  TR_YLIM  <- c(0, 110)
  GSW_MAX  <- 0.45   # approximate max from data; adjust if needed

  # Compute predictions from stored parameters
  compute_preds <- function(params_df, vpd_targets) {
    params_df |>
      rowwise() |>
      mutate(
        preds = list(
          tibble(VPD_target = vpd_targets,
                 value = pred_from_params(vpd_targets, Best_model,
                                          Intercept, Slope1, BP_kPa, Slope2))
        )
      ) |>
      ungroup() |>
      select(Genotype, WUE_class, Best_model, preds) |>
      unnest(preds)
  }

  tr_pred  <- compute_preds(tr_params, VPD_PRED)
  gsw_pred <- compute_preds(gsw_params, VPD_PRED)

  make_panel_df <- function(df) {
    df |>
      mutate(
        VPD_label  = factor(paste0(VPD_target, " kPa"),
                            levels = c("1.5 kPa", "3.5 kPa")),
        WUE_short  = wue_short[WUE_class],
        Panel_label = paste0(Genotype, "\n[", WUE_short, "]")
      ) |>
      mutate(
        wue_rank    = match(WUE_class, WUE_ORDER_F),
        Panel_label = factor(Panel_label,
                             levels = unique(Panel_label[order(wue_rank, Genotype)]))
      )
  }

  tr_plot  <- make_panel_df(tr_pred)
  gsw_plot <- make_panel_df(gsw_pred)

  strip_fill_for <- function(df) {
    lvls   <- levels(df$Panel_label)
    wshort <- gsub("^.*\\[(.*)\\]$", "\\1", lvls)
    setNames(WUE_STRIP_COLS[wshort], lvls)
  }

  gsw_strip_cols <- strip_fill_for(gsw_plot)
  tr_strip_cols  <- strip_fill_for(tr_plot)

  p_gsw <- ggplot(gsw_plot, aes(x = VPD_label, y = value, fill = VPD_label)) +
    geom_col(width = 0.6, colour = "white", linewidth = 0.4) +
    scale_fill_manual(values = VPD_COLS, name = "VPD") +
    facet_wrap2(~Panel_label, ncol = 5, scales = "fixed",
                strip = strip_themed(
                  background_x = elem_list_rect(
                    fill = unname(gsw_strip_cols)))) +
    coord_cartesian(ylim = c(0, GSW_MAX * 1.10)) +
    labs(x = NULL,
         y = expression(g[sw]~(mol~m^{-2}~s^{-1}))) +
    UA_THEME +
    theme(legend.position = "top",
          axis.text.x     = element_blank(),
          axis.ticks.x    = element_blank(),
          strip.text      = element_text(colour = "white", face = "bold",
                                         size = 6.5, lineheight = 1.1))

  p_tr <- ggplot(tr_plot, aes(x = VPD_label, y = value, fill = VPD_label)) +
    geom_col(width = 0.6, colour = "white", linewidth = 0.4) +
    scale_fill_manual(values = VPD_COLS, name = "VPD") +
    facet_wrap2(~Panel_label, ncol = 5, scales = "fixed",
                strip = strip_themed(
                  background_x = elem_list_rect(
                    fill = unname(tr_strip_cols)))) +
    coord_cartesian(ylim = TR_YLIM) +
    labs(x = "VPD treatment",
         y = expression(TR~(mg~m^{-2}~s^{-1}))) +
    UA_THEME +
    theme(legend.position = "none",
          strip.text      = element_text(colour = "white", face = "bold",
                                         size = 6.5, lineheight = 1.1))

  p5 <- p_gsw / p_tr +
    plot_layout(heights = c(1, 1.4))

  ggsave(file.path(OUTDIR, "Figure5_predicted_TR_gsw.tiff"), p5,
         width = 13, height = 9, dpi = DPI, device = tiff_cairo)
  cat("  Saved Figure 5\n")
})

cat("\nAll 5 paper figures saved to:", normalizePath(OUTDIR), "\n")
