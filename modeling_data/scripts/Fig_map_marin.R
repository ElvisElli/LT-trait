# ==============================================================================
# Fig_map_marin.R  (v2)
#
# Marin et al. (2022)-style yield map — Nature publishing aesthetic.
#
# Layout:
#   [shared legend at top — LT vs Non-LT, rainfed vs yield gap]
#   (A) Map × 2 panels (LT | Non-LT) — pie charts at 6 sites
#   (B) Bar chart — all locations, no state faceting, no cap whiskers
#
# Key design choices:
#   • Pie circles computed as polygons with aspect-ratio correction so
#     they are visually round under coord_fixed()
#   • Nature style: clean white panels, no axis ticks/labels, grey base map,
#     dark strip labels, thin borders, sans-serif font
#   • Shared 4-level fill legend (genotype × slice) at top
#   • Error bars: SD across 6 site-level means — no horizontal caps
#
# Data:   modeling_data/data/SoybeanSSM_sims2.xlsx
# Output: modeling_data/figures/Fig_map_marin.tiff (600 DPI)
#         modeling_data/figures/Fig_map_marin.png  (300 DPI, git-tracked)
#
# Run from the LT-trait project root.
# Windows: edit project_root below; Linux/Mac: default path used automatically.
#
# Packages (install if missing):
#   install.packages(c("readxl","dplyr","tidyr","ggplot2",
#                      "patchwork","maps","ggrepel"))
# (ggforce no longer needed — pies are built with geom_polygon)
#
# Author: Elvis F. Elli | eelli@uark.edu
# ==============================================================================


# ── 0.  Working directory  ────────────────────────────────────────────────────

if (.Platform$OS.type == "windows") {
  project_root <- "C:/Users/YourName/LT-trait"   # <-- adjust for your machine
  setwd(project_root)
} else {
  setwd("/home/user/LT-trait")
}


# ── 1.  Packages  ─────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(readxl)    # read SSM Excel data
  library(dplyr)     # data wrangling
  library(tidyr)     # pivot_wider / pivot_longer
  library(ggplot2)   # all plotting
  library(patchwork) # combine panels
  library(maps)      # US state boundary polygons
  library(ggrepel)   # repelled site labels
})


# ── 2.  Design constants  ─────────────────────────────────────────────────────

# Project colour convention: LT = orange, Non-LT = green
# Dark tint = water-limited (rainfed);  light tint = yield gap (to potential)
COL_LT_RAIN  <- "#d95f02"   # LT   rainfed  (dark orange)
COL_LT_GAP   <- "#fdbe85"   # LT   gap      (light orange)
COL_NLT_RAIN <- "#31a354"   # Non-LT rainfed (dark green)
COL_NLT_GAP  <- "#a1d99b"   # Non-LT gap    (light green)

FILL_VALS <- c(
  "LT|Rainfed"     = COL_LT_RAIN,
  "LT|Gap"         = COL_LT_GAP,
  "Non-LT|Rainfed" = COL_NLT_RAIN,
  "Non-LT|Gap"     = COL_NLT_GAP
)

FILL_LABELS <- c(
  "LT|Rainfed"     = "LT \u2014 water-limited (rainfed)",
  "LT|Gap"         = "LT \u2014 yield gap (to potential)",
  "Non-LT|Rainfed" = "Non-LT \u2014 water-limited (rainfed)",
  "Non-LT|Gap"     = "Non-LT \u2014 yield gap (to potential)"
)

FILL_BREAKS <- c("LT|Rainfed", "LT|Gap", "Non-LT|Rainfed", "Non-LT|Gap")

# Strip labels for the two map panels
GENO_LABELS <- c(
  "LT"     = "Limited transpiration",
  "Non-LT" = "Non-limited transpiration"
)

# Latitude correction factor so pies look visually round.
# A lon/lat coordinate system with coord_fixed(MAP_RATIO) stretches y by
# MAP_RATIO relative to x.  Compensate pie y-radius by dividing by MAP_RATIO.
MAP_LAT      <- 38.5                           # mean latitude of study area
MAP_RATIO    <- 1 / cos(MAP_LAT * pi / 180)   # ≈ 1.277
PIE_MAX_R    <- 1.55   # maximum pie x-radius (degrees lon) for highest-yield site
PIE_N_PTS    <- 90     # polygon points per full circle (smoothness)


# ── 3.  Site coordinates  ─────────────────────────────────────────────────────

sites <- tibble(
  Location = c("Jonesboro", "Marianna", "Columbia",
               "Norborne",  "Lincoln",  "North Platte"),
  State    = c("AR", "AR", "MO", "MO", "NE", "NE"),
  lat      = c( 35.84,  34.78,  38.95,  39.30,  40.81,  41.12),
  lon      = c(-90.70, -90.76, -92.33, -93.67, -96.70, -100.77),
  code     = c("JO",   "MA",   "CO",   "NO",   "LI",   "NP")
)


# ── 4.  Load and summarise SSM data  ─────────────────────────────────────────

ssm_raw <- read_excel("modeling_data/data/SoybeanSSM_sims2.xlsx")

ssm <- ssm_raw |>
  filter(Crop %in% c("SOYBEAN-MG4",       "SOYBEAN-LT1.5-MG4",
                     "SOYBEAN-LT1.5-MG3", "SOYBEAN-MG3")) |>
  mutate(
    # LT panel on the left → set as first factor level
    Genotype   = factor(
      ifelse(Crop %in% c("SOYBEAN-MG4", "SOYBEAN-MG3"), "Non-LT", "LT"),
      levels = c("LT", "Non-LT")
    ),
    water_mgmt = if_else(water_mgmt == "Rainfed", "Rainfed", "Irrigated")
  )

# 30-year × planting-date mean yield per site × genotype × water management
yield_site <- ssm |>
  group_by(Location, Genotype, water_mgmt) |>
  summarise(Yield_mean = mean(Yield_kgha, na.rm = TRUE), .groups = "drop") |>
  left_join(sites, by = "Location")

# Wide format: Rainfed and Irrigated columns + pie geometry
yield_wide <- yield_site |>
  pivot_wider(names_from = water_mgmt, values_from = Yield_mean) |>
  mutate(
    Gap   = Irrigated - Rainfed,
    frac  = Rainfed / Irrigated,          # share of potential achieved by rainfed
    r_pie = sqrt(Irrigated) /             # pie x-radius ∝ √(potential yield)
              sqrt(max(Irrigated, na.rm = TRUE)) * PIE_MAX_R
  )


# ── 5.  Pie polygon builder  ──────────────────────────────────────────────────
#
# Computes pie-slice polygon coordinates in lon/lat data space.
# y-radius = r / MAP_RATIO so the slice looks visually circular
# when drawn with coord_fixed(MAP_RATIO).

make_pie_polys <- function(row, ratio = MAP_RATIO, n_pts = PIE_N_PTS) {
  x0   <- row[["lon"]]
  y0   <- row[["lat"]]
  r    <- row[["r_pie"]]
  ry   <- r / ratio            # corrected y-radius for visual roundness
  frac <- row[["frac"]]
  geno <- as.character(row[["Genotype"]])
  loc  <- row[["Location"]]

  n1 <- max(3L, ceiling(n_pts * frac))
  n2 <- max(3L, ceiling(n_pts * (1 - frac)))

  t_rain <- seq(-pi / 2, -pi / 2 + frac * 2 * pi, length.out = n1)
  t_gap  <- seq(-pi / 2 + frac * 2 * pi, -pi / 2 + 2 * pi, length.out = n2)

  rbind(
    data.frame(
      Genotype = geno, Location = loc,
      slice    = "Rainfed",
      lon0     = x0, lat0 = y0, r_pie = r, frac = frac,
      px       = c(x0, x0 + r  * cos(t_rain), x0),
      py       = c(y0, y0 + ry * sin(t_rain), y0),
      grp      = paste(geno, loc, "rain", sep = "_"),
      stringsAsFactors = FALSE
    ),
    data.frame(
      Genotype = geno, Location = loc,
      slice    = "Gap",
      lon0     = x0, lat0 = y0, r_pie = r, frac = frac,
      px       = c(x0, x0 + r  * cos(t_gap), x0),
      py       = c(y0, y0 + ry * sin(t_gap), y0),
      grp      = paste(geno, loc, "gap", sep = "_"),
      stringsAsFactors = FALSE
    )
  )
}

pie_polys <- do.call(rbind, lapply(seq_len(nrow(yield_wide)), function(i)
  make_pie_polys(yield_wide[i, ])
)) |>
  mutate(
    fill_key = paste(Genotype, slice, sep = "|"),
    Genotype = factor(Genotype, levels = c("LT", "Non-LT"))
  )


# ── 6.  Pie label positions  ──────────────────────────────────────────────────
#
# rainfed yield (white text) inside the rainfed slice at ~42 % of radius
# irrigated yield (dark text) above each pie

pie_labels <- yield_wide |>
  mutate(
    # mid-angle of the rainfed arc
    ang_mid  = -pi / 2 + frac * pi,
    lbl_x    = lon + r_pie * 0.42 * cos(ang_mid),
    lbl_y    = lat + (r_pie / MAP_RATIO) * 0.42 * sin(ang_mid),
    # potential yield label above pie top
    lbl_top_y = lat + r_pie / MAP_RATIO + 0.38,
    lbl_rain  = paste0(round(Rainfed   / 1000, 1), " t"),
    lbl_irr   = paste0(round(Irrigated / 1000, 1), " t")
  )


# ── 7.  US map base layers  ──────────────────────────────────────────────────

us_states    <- map_data("state")
focus_states <- c("arkansas", "missouri", "nebraska")
state_bg     <- filter(us_states, region %in% focus_states)

# State name positions (inside each state)
state_ann <- tibble(
  lon   = c(-91.5,  -93.5,  -99.3),
  lat   = c( 33.6,   36.8,   43.1),
  label = c("Arkansas", "Missouri", "Nebraska")
)


# ── 8.  Nature-style shared theme  ────────────────────────────────────────────

theme_nature <- function(base_size = 8) {
  theme_bw(base_size = base_size, base_family = "sans") +
  theme(
    # Map panels: no axes
    axis.text         = element_blank(),
    axis.ticks        = element_blank(),
    axis.title        = element_blank(),
    panel.grid        = element_blank(),
    # Thin clean border
    panel.border      = element_rect(colour = "grey35", linewidth = 0.4),
    # Dark strip with white bold text
    strip.background  = element_rect(fill = "grey20", colour = "grey20"),
    strip.text        = element_text(colour = "white", face = "bold",
                                     size = base_size),
    # Legend: compact, at top, horizontal, no title
    legend.position   = "top",
    legend.direction  = "horizontal",
    legend.title      = element_blank(),
    legend.text       = element_text(size = base_size - 1),
    legend.key.size   = unit(0.55, "lines"),
    legend.key        = element_rect(colour = "grey50", linewidth = 0.3),
    legend.spacing.x  = unit(0.25, "cm"),
    legend.margin     = margin(0, 0, 2, 0),
    # Tight plot margins
    plot.margin       = margin(1, 2, 1, 2, "mm"),
    panel.spacing     = unit(0.35, "lines")
  )
}

# Bar-panel theme adds back axis text and ticks
theme_nature_bar <- function(base_size = 8) {
  theme_nature(base_size) +
  theme(
    axis.text.x  = element_text(size = base_size - 1, colour = "grey20"),
    axis.text.y  = element_text(size = base_size - 1, colour = "grey20"),
    axis.title.y = element_text(size = base_size, colour = "grey15", angle = 90,
                                margin = margin(r = 3)),
    axis.ticks   = element_line(linewidth = 0.3, colour = "grey40"),
    axis.ticks.length = unit(1.5, "mm"),
    panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey88"),
    legend.position    = "none"   # legend carried by map panel only
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# PANEL (A) — MAP  (LT | Non-LT), faceted
# ══════════════════════════════════════════════════════════════════════════════

map_panel <- ggplot() +

  # ── US background (all states, light grey) ────────────────────────────────
  geom_polygon(
    data = us_states,
    aes(x = long, y = lat, group = group),
    fill = "grey93", colour = "grey78", linewidth = 0.18
  ) +

  # ── Study states: AR / MO / NE (white, with visible borders) ─────────────
  geom_polygon(
    data = state_bg,
    aes(x = long, y = lat, group = group),
    fill = "white", colour = "grey30", linewidth = 0.4
  ) +

  # ── State name labels ─────────────────────────────────────────────────────
  # No Genotype column → appears in BOTH facet panels automatically
  geom_text(
    data = state_ann,
    aes(x = lon, y = lat, label = label),
    fontface = "bold.italic", size = 2.3, colour = "grey45",
    inherit.aes = FALSE
  ) +

  # ── Pie slice polygons  ───────────────────────────────────────────────────
  # pie_polys has Genotype → each panel draws only its own genotype's pies
  geom_polygon(
    data = pie_polys,
    aes(x = px, y = py, group = grp, fill = fill_key),
    colour = "white", linewidth = 0.35
  ) +

  # ── Rainfed yield: white text inside each pie slice ───────────────────────
  geom_text(
    data = pie_labels,
    aes(x = lbl_x, y = lbl_y, label = lbl_rain),
    size = 2.0, colour = "white", fontface = "bold", inherit.aes = FALSE
  ) +

  # ── Potential (irrigated) yield: dark text above pie ─────────────────────
  geom_text(
    data = pie_labels,
    aes(x = lon, y = lbl_top_y, label = lbl_irr),
    size = 1.9, colour = "grey15", inherit.aes = FALSE
  ) +

  # ── Site code below each pie (JO / MA / CO / NO / LI / NP) ──────────────
  # sites has no Genotype → code appears in both panels
  geom_text(
    data = yield_wide |> distinct(Location, lon, lat, r_pie, .keep_all = FALSE) |>
      left_join(sites |> select(Location, code), by = "Location"),
    aes(x = lon, y = lat - r_pie / MAP_RATIO - 0.45, label = code),
    size = 2.2, colour = "grey20", fontface = "bold", inherit.aes = FALSE
  ) +

  # ── Facets: left = LT, right = Non-LT ────────────────────────────────────
  facet_wrap(~Genotype, ncol = 2,
             labeller = labeller(Genotype = GENO_LABELS)) +

  # ── Shared fill scale (4 levels) ─────────────────────────────────────────
  scale_fill_manual(
    values = FILL_VALS,
    labels = FILL_LABELS,
    breaks = FILL_BREAKS,
    name   = NULL
  ) +

  # Correct aspect ratio for US ~38.5°N latitude
  coord_fixed(MAP_RATIO, xlim = c(-106, -84), ylim = c(31.0, 46.5)) +

  theme_nature() +
  guides(fill = guide_legend(nrow = 2, override.aes = list(colour = NA)))


# ══════════════════════════════════════════════════════════════════════════════
# PANEL (B) — BAR CHART  (all locations, no state facets)
#
# For each of 6 sites, two bars side by side (LT and Non-LT),
# each stacked rainfed (dark) + gap (light).
# Error bar = SD across 30 simulation years × 3 planting dates at that site.
# No horizontal whisker caps (width = 0).
# ══════════════════════════════════════════════════════════════════════════════

# ── Location code lookup (named vector: Location → code) ──────────────────
site_codes <- setNames(sites$code, sites$Location)
CODE_LEVELS <- c("JO", "MA", "CO", "NO", "LI", "NP")

# ── Site-level bars (rainfed + gap stacked, per genotype) ─────────────────
bar_long <- yield_wide |>
  pivot_longer(cols = c(Rainfed, Gap),
               names_to  = "slice",
               values_to = "Yield_kgha") |>
  mutate(
    fill_key = paste(as.character(Genotype), slice, sep = "|"),
    slice    = factor(slice, levels = c("Rainfed", "Gap")),
    code     = factor(site_codes[Location], levels = CODE_LEVELS)
  )

# ── Error bar: SD across 30 years × 3 planting dates at each site ─────────
err_data <- ssm |>
  filter(water_mgmt == "Irrigated") |>
  group_by(Location, Genotype) |>
  summarise(
    Irr_mean = mean(Yield_kgha, na.rm = TRUE),
    Irr_sd   = sd(Yield_kgha,   na.rm = TRUE),
    .groups  = "drop"
  ) |>
  left_join(sites |> select(Location, code), by = "Location") |>
  mutate(
    code     = factor(code, levels = CODE_LEVELS),
    Genotype = factor(Genotype, levels = c("LT", "Non-LT"))
  )

bar_panel <- ggplot() +

  # ── Stacked bars (rainfed dark, gap light) ────────────────────────────────
  geom_col(
    data     = bar_long,
    aes(x    = code, y = Yield_kgha / 1000, fill = fill_key,
        group = Genotype),
    position = position_dodge(width = 0.7),
    width    = 0.62,
    colour   = NA
  ) +

  # ── Error bars — no horizontal caps ───────────────────────────────────────
  geom_errorbar(
    data     = err_data,
    aes(x    = code,
        ymin  = (Irr_mean - Irr_sd) / 1000,
        ymax  = (Irr_mean + Irr_sd) / 1000,
        group = Genotype),
    position = position_dodge(width = 0.7),
    width    = 0,             # no horizontal cap whiskers
    linewidth = 0.5,
    colour   = "grey20"
  ) +

  # ── State separator lines (thin vertical guides between state groups) ──────
  geom_vline(xintercept = c(2.5, 4.5),
             linewidth = 0.3, colour = "grey70", linetype = "dashed") +

  # ── State labels along the top ────────────────────────────────────────────
  annotate("text",
           x = c(1.5, 3.5, 5.5), y = Inf,
           label = c("AR", "MO", "NE"),
           vjust = -0.3, size = 2.3, fontface = "bold", colour = "grey30") +

  scale_fill_manual(values = FILL_VALS, guide = "none") +

  scale_y_continuous(
    expand = expansion(mult = c(0, 0.08)),
    name   = expression(Yield~(t~ha^{-1})),
    breaks = seq(0, 8, by = 2)
  ) +

  coord_cartesian(clip = "off") +   # allow state labels above plot area
  labs(x = NULL) +

  theme_nature_bar()


# ══════════════════════════════════════════════════════════════════════════════
# COMBINE with patchwork
# Map (left, 3/4 width) + bar chart (right, 1/4 width)
# tag_levels = "A" adds bold (A) / (B) panel labels (no figure title)
# ══════════════════════════════════════════════════════════════════════════════

combined <- (map_panel | bar_panel) +
  plot_layout(widths = c(3, 1)) +
  plot_annotation(
    tag_levels = "A",
    theme      = theme(
      plot.tag = element_text(face = "bold", size = 9, family = "sans")
    )
  )


# ══════════════════════════════════════════════════════════════════════════════
# SAVE
# TIFF at 600 DPI (LZW) — full-quality, not git-tracked (*.tiff in .gitignore)
# PNG  at 300 DPI       — git-tracked
# ══════════════════════════════════════════════════════════════════════════════

tiff_type <- if (.Platform$OS.type == "windows") "windows" else "cairo"

tiff_file <- "modeling_data/figures/Fig_map_marin.tiff"
tiff(tiff_file, width = 30, height = 14, units = "cm",
     res = 600, compression = "lzw", type = tiff_type)
print(combined)
dev.off()
message("Saved TIFF: ", tiff_file)

png_file <- "modeling_data/figures/Fig_map_marin.png"
ggsave(png_file, combined,
       width = 30, height = 14, units = "cm",
       dpi = 300, bg = "white")
message("Saved PNG:  ", png_file)
