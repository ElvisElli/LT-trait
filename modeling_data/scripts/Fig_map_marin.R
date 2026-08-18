# ==============================================================================
# Fig_map_marin.R  (v4)
#
# Marin et al. (2022)-style yield map — Nature publishing aesthetic.
#
# Layout (panels aligned, shared legend at top):
#   (A) Map × 2 facets (LT | Non-LT) — pie charts at 6 sites, no map text
#   (B) Bar chart — mean ± SE across 6 locations; stacked rainfed + gap
#
# Colours  (consistent between map and bar):
#   Water-limited yield  (rainfed)  = #0571b0  blue
#   Potential yield      (irrigated gap)  = #ca0020  red
#
# Data:   modeling_data/data/SoybeanSSM_sims2.xlsx
# Output: modeling_data/figures/Fig_map_marin.tiff (600 DPI, LZW)
#         modeling_data/figures/Fig_map_marin.png  (300 DPI, git-tracked)
#
# Run from the LT-trait project root.
# Windows: set project_root below; Linux/Mac: default used automatically.
#
# Required packages:
#   install.packages(c("readxl","dplyr","tidyr","ggplot2","patchwork","maps"))
#
# Author: Elvis F. Elli | eelli@uark.edu
# ==============================================================================


# ── 1.  Packages  ─────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(readxl)    # read SSM Excel data
  library(dplyr)     # data wrangling
  library(tidyr)     # pivot_wider / pivot_longer
  library(ggplot2)   # all plotting
  library(cowplot)   # plot_grid + get_legend to combine panels
  library(maps)      # US state boundary polygons
})


# ── 2.  Design constants — edit these to tune the figure  ─────────────────────

# Colours
COL_RAIN <- "#0571b0"   # water-limited yield   (rainfed)   — blue
COL_GAP  <- "#ca0020"   # potential yield        (irr gap)   — red

FILL_VALS   <- c("Rainfed" = COL_RAIN, "Gap" = COL_GAP)
FILL_LABELS <- c("Rainfed" = "Water-limited yield",
                 "Gap"     = "Yield gap by water deficit")
FILL_BREAKS <- c("Rainfed", "Gap")

# Facet strip labels
GENO_LABELS <- c("LT" = "Limited transpiration",
                 "Non-LT" = "Non-limited transpiration")

# Map aspect: cos-corrected for study-area mean latitude (≈ 38.5°N)
# MAP_RATIO: 1 degree latitude takes MAP_RATIO× more pixels than 1 degree longitude.
# MAP_ASPECT: panel height/width that gives the same visual scale as coord_fixed(MAP_RATIO).
#   = MAP_RATIO × (y_range / x_range)
# Using theme(aspect.ratio = MAP_ASPECT) instead of coord_fixed() lets cowplot
# resize panels for strip alignment while preserving round circles.
MAP_LAT    <- 38.5
MAP_RATIO  <- 1 / cos(MAP_LAT * pi / 180)   # ≈ 1.277

# Map extent — zoom in on AR / MO / NE with a small margin around the pies
MAP_XLIM <- c(-106, -86)    # longitude window
MAP_YLIM <- c(31.5,  44.5)  # latitude  window

MAP_ASPECT <- MAP_RATIO * diff(MAP_YLIM) / diff(MAP_XLIM)  # ≈ 0.830

# Pie size — increase PIE_MAX_R to make all pies larger (still proportional)
PIE_MAX_R <- 1.5    # x-radius (lon degrees) for the site with highest potential yield
PIE_N_PTS <- 100    # polygon smoothness (points per full circle)

# Border appearance
BORDER_COL <- "black"
BORDER_LW  <- 0.7   # linewidth for panel borders, strips, and pie outlines


# ── 3.  Site coordinates  ─────────────────────────────────────────────────────

sites <- tibble(
  Location = c("Jonesboro", "Marianna", "Columbia",
               "Norborne",  "Lincoln",  "North Platte"),
  State    = c("AR", "AR", "MO", "MO", "NE", "NE"),
  lat      = c( 35.84,  34.78,  38.95,  39.30,  40.81,  41.12),
  lon      = c(-90.70, -90.76, -92.33, -93.67, -96.70, -100.77)
)


# ── 4.  Load and summarise SSM data  ─────────────────────────────────────────

ssm_raw <- read_excel("modeling_data/data/SoybeanSSM_sims2.xlsx")

ssm <- ssm_raw |>
  filter(Crop %in% c("SOYBEAN-MG4",       "SOYBEAN-LT1.5-MG4",
                     "SOYBEAN-LT1.5-MG3", "SOYBEAN-MG3")) |>
  mutate(
    Genotype   = factor(
      ifelse(Crop %in% c("SOYBEAN-MG4", "SOYBEAN-MG3"), "Non-LT", "LT"),
      levels = c("LT", "Non-LT")        # LT on left panel
    ),
    water_mgmt = if_else(water_mgmt == "Rainfed", "Rainfed", "Irrigated")
  )

# ── 4a. 30-yr mean per site × genotype × water management  (for pie charts) ──
yield_site <- ssm |>
  group_by(Location, Genotype, water_mgmt) |>
  summarise(Yield_mean = mean(Yield_kgha, na.rm = TRUE), .groups = "drop") |>
  left_join(sites, by = "Location")

yield_wide <- yield_site |>
  pivot_wider(names_from = water_mgmt, values_from = Yield_mean) |>
  mutate(
    Gap   = Irrigated - Rainfed,
    frac  = Rainfed / Irrigated,    # rainfed fraction of potential (= pie blue fraction)
    r_pie = sqrt(Irrigated) /       # pie x-radius ∝ √(potential yield)
              sqrt(max(Irrigated, na.rm = TRUE)) * PIE_MAX_R
  )

# ── 4b. Location-level means → cross-location mean ± SE  (for bar chart) ──────
#
# Step 1: average each site over its 30 years × 3 planting dates
# Step 2: treat each site as one observation; compute mean and SE (n = 6) per genotype

bar_loc <- ssm |>
  group_by(Location, Genotype, water_mgmt) |>
  summarise(Yield_mean = mean(Yield_kgha, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = water_mgmt, values_from = Yield_mean) |>
  mutate(Gap = Irrigated - Rainfed)

bar_summary <- bar_loc |>
  group_by(Genotype) |>
  summarise(
    Rain_mean = mean(Rainfed,   na.rm = TRUE),
    Gap_mean  = mean(Gap,       na.rm = TRUE),
    Irr_mean  = mean(Irrigated, na.rm = TRUE),
    Irr_se    = sd(Irrigated,   na.rm = TRUE) / sqrt(n()),
    n_loc     = n(),
    .groups   = "drop"
  ) |>
  mutate(Genotype = factor(Genotype, levels = c("LT", "Non-LT")))


# ── 5.  Pie polygon builder  ──────────────────────────────────────────────────
#
# Slice "Rainfed"  (blue): angle = frac × 2π  (water-limited proportion)
# Slice "Gap"      (red) : remaining angle    (potential gain from irrigation)
# y-radius corrected by MAP_RATIO → visually round circles under coord_fixed()

make_pie_polys <- function(row, ratio = MAP_RATIO, n_pts = PIE_N_PTS) {
  x0   <- row[["lon"]]
  y0   <- row[["lat"]]
  r    <- row[["r_pie"]]
  ry   <- r / ratio          # corrected y-radius
  frac <- row[["frac"]]
  geno <- as.character(row[["Genotype"]])
  loc  <- row[["Location"]]

  n1 <- max(3L, ceiling(n_pts * frac))
  n2 <- max(3L, ceiling(n_pts * (1 - frac)))

  t_rain <- seq(-pi / 2, -pi / 2 + frac * 2 * pi,       length.out = n1)
  t_gap  <- seq(-pi / 2 + frac * 2 * pi, -pi / 2 + 2 * pi, length.out = n2)

  rbind(
    data.frame(
      Genotype = geno, Location = loc, slice = "Rainfed",
      px = c(x0, x0 + r * cos(t_rain), x0),
      py = c(y0, y0 + ry * sin(t_rain), y0),
      grp = paste(geno, loc, "rain", sep = "_"),
      stringsAsFactors = FALSE
    ),
    data.frame(
      Genotype = geno, Location = loc, slice = "Gap",
      px = c(x0, x0 + r * cos(t_gap), x0),
      py = c(y0, y0 + ry * sin(t_gap), y0),
      grp = paste(geno, loc, "gap", sep = "_"),
      stringsAsFactors = FALSE
    )
  )
}

pie_polys <- do.call(rbind,
  lapply(seq_len(nrow(yield_wide)), function(i) make_pie_polys(yield_wide[i, ]))
) |>
  mutate(
    slice    = factor(slice, levels = c("Rainfed", "Gap")),
    Genotype = factor(Genotype, levels = c("LT", "Non-LT"))
  )


# ── 6.  US base-map layers  ───────────────────────────────────────────────────

us_states    <- map_data("state")
focus_states <- c("arkansas", "missouri", "nebraska")
state_bg     <- filter(us_states, region %in% focus_states)


# ── 7.  Shared Nature-style themes  ───────────────────────────────────────────

theme_nature_map <- function(base_size = 8) {
  theme_bw(base_size = base_size, base_family = "sans") +
  theme(
    axis.text         = element_blank(),
    axis.ticks        = element_blank(),
    axis.title        = element_blank(),
    panel.grid        = element_blank(),
    # Thick black panel border
    panel.border      = element_rect(colour = BORDER_COL, linewidth = BORDER_LW),
    # White strip with thick black border and dark text
    strip.background  = element_rect(fill = "white", colour = BORDER_COL,
                                     linewidth = BORDER_LW),
    strip.text        = element_text(colour = "grey10", face = "bold",
                                     size = base_size + 0.5),
    # Legend collected at top by patchwork — suppress here
    legend.position   = "top",
    legend.direction  = "horizontal",
    legend.title      = element_blank(),
    legend.text       = element_text(size = base_size),
    legend.key.size   = unit(0.55, "lines"),
    legend.key        = element_rect(colour = "grey40", linewidth = 0.35),
    legend.spacing.x  = unit(0.35, "cm"),
    legend.margin     = margin(0, 0, 2, 0),
    plot.margin       = margin(1, 2, 1, 2, "mm"),
    panel.spacing     = unit(0.4, "lines")
  )
}

theme_nature_bar <- function(base_size = 8) {
  theme_nature_map(base_size) +
  theme(
    axis.text.x        = element_text(size = base_size + 0.5, colour = "grey10",
                                      face = "bold"),
    axis.text.y        = element_text(size = base_size - 0.5, colour = "grey20"),
    axis.title.y       = element_text(size = base_size, colour = "grey10",
                                      angle = 90, margin = margin(r = 3)),
    axis.ticks         = element_line(linewidth = 0.4, colour = "grey30"),
    axis.ticks.length  = unit(1.5, "mm"),
    # Subtle horizontal grid only
    panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey85"),
    # Strip sized to match height of map panel strip
    strip.background   = element_rect(fill = "white", colour = BORDER_COL,
                                      linewidth = BORDER_LW),
    strip.text         = element_text(colour = "grey10", face = "bold",
                                      size = base_size + 0.5),
    legend.position    = "none"   # legend collected by patchwork into map slot
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# PANELS (A) and (B) — MAP  (one ggplot per genotype)
#
# Splitting into two separate ggplot objects lets cowplot align their strip
# titles at exactly the same height as panel (C)'s strip.
# ══════════════════════════════════════════════════════════════════════════════

# Helper — builds one map panel for a single genotype level
make_map_panel <- function(geno, show_legend = FALSE) {
  ggplot() +

    # All US states — light grey background
    geom_polygon(
      data = us_states,
      aes(x = long, y = lat, group = group),
      fill = "grey92", colour = "grey75", linewidth = 0.15
    ) +

    # Study states (AR / MO / NE) — darker fill
    geom_polygon(
      data = state_bg,
      aes(x = long, y = lat, group = group),
      fill = "grey65", colour = BORDER_COL, linewidth = 0.45
    ) +

    # Pie slices for this genotype only
    geom_polygon(
      data   = filter(pie_polys, Genotype == geno),
      aes(x = px, y = py, group = grp, fill = slice),
      colour = BORDER_COL, linewidth = 0.28
    ) +

    # Single-level facet → strip title matches bar panel strip height
    facet_wrap(~Genotype, labeller = labeller(Genotype = GENO_LABELS)) +

    scale_fill_manual(
      values = FILL_VALS,
      labels = FILL_LABELS,
      breaks = FILL_BREAKS,
      name   = NULL
    ) +

    coord_cartesian(xlim = MAP_XLIM, ylim = MAP_YLIM) +

    theme_nature_map() +

    # aspect.ratio = MAP_ASPECT is mathematically equivalent to coord_fixed(MAP_RATIO)
    # for this xlim/ylim window, but allows cowplot to align panel tops across plots.
    theme(
      aspect.ratio   = MAP_ASPECT,
      legend.position = if (show_legend) "top" else "none"
    ) +

    guides(fill = guide_legend(nrow = 1,
                               override.aes = list(colour = "grey30",
                                                   linewidth = 0.3)))
}

# Build one panel per genotype (no legend — legend extracted separately)
map_lt    <- make_map_panel("LT")
map_nonlt <- make_map_panel("Non-LT")


# ══════════════════════════════════════════════════════════════════════════════
# PANEL (B) — BAR CHART
#
# X-axis: LT and Non-LT (mean across 6 locations).
# Bars stacked: water-limited (blue) on bottom, potential gain (red) on top.
# Error bar placed at the irrigated total (Irr_mean ± Irr_se).
# No horizontal cap whiskers; y-axis in kg ha⁻¹.
# ══════════════════════════════════════════════════════════════════════════════

# Long format for stacked geom_col
bar_long <- bar_summary |>
  select(Genotype, Rain_mean, Gap_mean) |>
  pivot_longer(
    cols      = c(Rain_mean, Gap_mean),
    names_to  = "slice",
    values_to = "Yield_kgha"
  ) |>
  mutate(
    slice    = recode(slice, Rain_mean = "Rainfed", Gap_mean = "Gap"),
    slice    = factor(slice, levels = c("Rainfed", "Gap")),   # Rainfed on bottom
    Genotype = factor(Genotype, levels = c("LT", "Non-LT"))
  )

bar_panel <- ggplot() +

  # Stacked bars (Rainfed bottom = blue, Gap top = red)
  geom_col(
    data  = bar_long,
    aes(x = Genotype, y = Yield_kgha, fill = slice),
    width   = 0.5,
    colour  = BORDER_COL,
    linewidth = 0.4
  ) +

  # Error bar at the top of the irrigated total (Irr_mean ± Irr_se)
  # width = 0 → no horizontal cap whiskers
  # The bar top equals Irr_mean, so ymin/ymax straddle the bar top correctly
  geom_errorbar(
    data  = bar_summary,
    aes(x    = Genotype,
        ymin  = Irr_mean - Irr_se,
        ymax  = Irr_mean + Irr_se),
    width     = 0,
    linewidth = 0.7,
    colour    = "grey10"
  ) +

  scale_fill_manual(values = FILL_VALS, guide = "none") +

  scale_y_continuous(
    expand = expansion(mult = c(0, 0.12)),
    name   = expression(Yield~(kg~ha^{-1})),
    labels = scales::comma,
    breaks = seq(0, 10000, by = 2000)
  ) +

  labs(x = NULL) +

  # Dummy facet label so strip height matches map panel — keeps alignment
  facet_wrap(~"All locations", strip.position = "top") +

  theme_nature_bar()


# ══════════════════════════════════════════════════════════════════════════════
# COMBINE — ggdraw with fixed coordinates
#
# cowplot's align= cannot equalize strips when aspect.ratio is set on some
# panels: it adjusts plot.margin ABOVE the strip, which shifts strips down.
# ggdraw places each panel at the same y and height → strip tops are at
# exactly the same position in all three panels regardless of aspect ratio.
#
# Width fractions: rel_widths 1:1:0.5 → total = 2.5
#   map_lt    0.00–0.40   (40%)
#   map_nonlt 0.40–0.80   (40%)
#   bar       0.80–1.00   (20%)
# Panels occupy y = 0.07–1.00 (93%), legend y = 0.00–0.07 (7%).
# ══════════════════════════════════════════════════════════════════════════════

# Shared legend as a grob (horizontal, placed in the bottom strip)
legend_grob <- get_legend(
  make_map_panel("LT", show_legend = TRUE) +
  theme(legend.position  = "bottom",
        legend.direction = "horizontal",
        legend.margin    = margin(0, 0, 0, 0))
)

combined <- ggdraw() +
  draw_plot(map_lt,      x = 0.00, y = 0.07, width = 0.40, height = 0.93) +
  draw_plot(map_nonlt,   x = 0.40, y = 0.07, width = 0.40, height = 0.93) +
  draw_plot(bar_panel,   x = 0.80, y = 0.07, width = 0.20, height = 0.93) +
  draw_grob(legend_grob, x = 0.50, y = 0.035, hjust = 0.5, vjust = 0.5) +
  draw_label("A", x = 0.005, y = 0.995, hjust = 0, vjust = 1,
             fontface = "bold", size = 9.5, fontfamily = "sans") +
  draw_label("B", x = 0.405, y = 0.995, hjust = 0, vjust = 1,
             fontface = "bold", size = 9.5, fontfamily = "sans") +
  draw_label("C", x = 0.805, y = 0.995, hjust = 0, vjust = 1,
             fontface = "bold", size = 9.5, fontfamily = "sans")


# ══════════════════════════════════════════════════════════════════════════════
# SAVE
# TIFF 600 DPI — submission quality; *.tiff in .gitignore (not committed)
# PNG  300 DPI — git-tracked for record and sharing
# ══════════════════════════════════════════════════════════════════════════════

tiff_type <- if (.Platform$OS.type == "windows") "windows" else "cairo"

tiff_file <- "modeling_data/figures/Fig_map_marin.tiff"
tiff(tiff_file, width = 20, height = 8, units = "cm",
     res = 600, compression = "lzw", type = tiff_type)
print(combined)
dev.off()
message("Saved TIFF: ", tiff_file)

png_file <- "modeling_data/figures/Fig_map_marin.png"
ggsave(png_file, combined, width = 20, height = 8, units = "cm",
       dpi = 300, bg = "white")
message("Saved PNG:  ", png_file)


