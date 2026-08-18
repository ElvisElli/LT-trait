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
PIE_MAX_R <- 1.05   # x-radius (lon degrees) — 30% smaller than 1.5
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

  # Start at top (π/2), sweep counterclockwise:
  #   blue (water-limited) opens to the LEFT, red (yield gap) to the RIGHT
  t_rain <- seq(pi / 2, pi / 2 + frac * 2 * pi,       length.out = n1)
  t_gap  <- seq(pi / 2 + frac * 2 * pi, pi / 2 + 2 * pi, length.out = n2)

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
    strip.text        = element_text(colour = "grey10", face = "plain",
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
# PANEL (A) — MAP  (both LT and Non-LT as one figure)
#
# Two facets side-by-side: LT on left, Non-LT on right.
# This is treated as a single panel in the combined layout (paired with C).
# ══════════════════════════════════════════════════════════════════════════════

map_panel_full <- ggplot() +

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

  # Pie slices for all genotypes
  geom_polygon(
    data = pie_polys,
    aes(x = px, y = py, group = grp, fill = slice),
    colour = BORDER_COL, linewidth = 0.28
  ) +

  # Facet: two panels (LT | Non-LT)
  facet_wrap(~Genotype, ncol = 2,
             labeller = labeller(Genotype = GENO_LABELS)) +

  scale_fill_manual(
    values = FILL_VALS,
    labels = FILL_LABELS,
    breaks = FILL_BREAKS,
    name   = NULL
  ) +

  coord_cartesian(xlim = MAP_XLIM, ylim = MAP_YLIM, expand = FALSE) +

  theme_nature_map() +

  guides(fill = guide_legend(nrow = 1,
                             override.aes = list(colour = "grey30",
                                                 linewidth = 0.3)))


# ══════════════════════════════════════════════════════════════════════════════
# PANEL (C) — BAR CHART  (stacked: water-limited bottom, yield gap top)
#
# X-axis: LT and Non-LT (mean across 6 locations).
# Stacking order mirrors the pie: blue (water-limited) on bottom,
# red (yield gap) on top. Error bar at the irrigated total ± SE.
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
    slice    = factor(slice, levels = c("Gap", "Rainfed")),  # Gap on top, Rainfed (blue) on bottom
    Genotype = factor(Genotype, levels = c("LT", "Non-LT"))
  )

bar_panel <- ggplot() +

  # Stacked: blue (water-limited) on bottom, red (yield gap) on top
  geom_col(
    data  = bar_long,
    aes(x = Genotype, y = Yield_kgha, fill = slice),
    width     = 0.5,
    colour    = BORDER_COL,
    linewidth = 0.4
  ) +

  # Error bar at the irrigated total (Irr_mean ± Irr_se); no cap whiskers
  geom_errorbar(
    data  = bar_summary,
    aes(x = Genotype, ymin = Irr_mean - Irr_se, ymax = Irr_mean + Irr_se),
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

  # Dummy facet strip so height matches map panels
  facet_wrap(~"All locations", strip.position = "top") +

  theme_nature_bar()


# ══════════════════════════════════════════════════════════════════════════════
# COMBINE — plot_grid: legend on top, then A (map) and B (bar) below
#
# Step 1: extract shared legend from map
# Step 2: remove legend from both panels
# Step 3: arrange map and bar side-by-side with labels A and B
# Step 4: stack legend above
# ══════════════════════════════════════════════════════════════════════════════

# 1. Extract legend
shared_legend <- get_legend(map_panel_full)

# 2. Panels with legend removed
map_no_legend <- map_panel_full + theme(legend.position = "none")
bar_no_legend <- bar_panel + theme(legend.position = "none")

# 3. Two panels side-by-side: A (map, 2:1) and B (bar, 1:1)
panels <- plot_grid(
  map_no_legend, bar_no_legend,
  ncol        = 2,
  rel_widths  = c(2, 1),
  labels      = c("A", "B"),
  label_size  = 9.5,
  label_fontface = "plain",
  label_fontfamily = "sans",
  align       = "hv",
  axis        = "tb"
)

# 4. Stack: legend on top, panels below
combined <- plot_grid(
  shared_legend, panels,
  ncol        = 1,
  rel_heights = c(0.08, 1)
)


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


