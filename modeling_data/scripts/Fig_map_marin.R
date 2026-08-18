# ==============================================================================
# Fig_map_marin.R  (v3)
#
# Marin et al. (2022)-style yield map — Nature publishing aesthetic.
#
# Layout:
#   [shared legend at top: Water-limited yield | Potential yield]
#   (A) Map × 2 panels (LT | Non-LT) — pie charts at 6 sites, no text
#   (B) Bar chart — mean across 6 locations; SE across location averages
#
# Colours (shared between map and bar):
#   Water-limited (rainfed) = #ca0020   (red)
#   Potential gain  (gap)   = #0571b0   (blue)
#
# Data:   modeling_data/data/SoybeanSSM_sims2.xlsx
# Output: modeling_data/figures/Fig_map_marin.tiff (600 DPI, LZW)
#         modeling_data/figures/Fig_map_marin.png  (300 DPI, git-tracked)
#
# Run from the LT-trait project root.
# Windows: set project_root below; Linux/Mac: default path used automatically.
#
# Packages needed (install if missing):
#   install.packages(c("readxl","dplyr","tidyr","ggplot2","patchwork","maps"))
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
})


# ── 2.  Design constants  ─────────────────────────────────────────────────────

# Two-colour scheme: same in both map and bar panels
COL_RAIN <- "#ca0020"   # water-limited (rainfed) yield
COL_GAP  <- "#0571b0"   # yield gap (potential above rainfed)

# Named fill for scale_fill_manual — only two levels
FILL_VALS   <- c("Rainfed" = COL_RAIN, "Gap" = COL_GAP)
FILL_LABELS <- c("Rainfed" = "Water-limited yield",
                 "Gap"     = "Potential yield")
FILL_BREAKS <- c("Rainfed", "Gap")

# Strip labels for the two map panels
GENO_LABELS <- c(
  "LT"     = "Limited transpiration",
  "Non-LT" = "Non-limited transpiration"
)

# Latitude correction so pie circles look visually round.
# coord_fixed(MAP_RATIO) stretches y by MAP_RATIO relative to x.
# Compensate by dividing the pie y-radius by MAP_RATIO.
MAP_LAT   <- 38.5
MAP_RATIO <- 1 / cos(MAP_LAT * pi / 180)   # ≈ 1.277

PIE_MAX_R <- 1.55   # max pie x-radius in lon degrees (for highest-yield site)
PIE_N_PTS <- 90     # polygon points per full circle (controls smoothness)


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

# ── 4a. Site-level 30-yr means (for pie charts) ───────────────────────────────
yield_site <- ssm |>
  group_by(Location, Genotype, water_mgmt) |>
  summarise(Yield_mean = mean(Yield_kgha, na.rm = TRUE), .groups = "drop") |>
  left_join(sites, by = "Location")

yield_wide <- yield_site |>
  pivot_wider(names_from = water_mgmt, values_from = Yield_mean) |>
  mutate(
    Gap   = Irrigated - Rainfed,
    frac  = Rainfed / Irrigated,              # rainfed fraction of potential
    r_pie = sqrt(Irrigated) /                 # pie x-radius ∝ √(potential yield)
              sqrt(max(Irrigated, na.rm = TRUE)) * PIE_MAX_R
  )

# ── 4b. Location-level means then average + SE across locations (for bar) ─────
#
# Step 1: collapse years × planting-dates into one mean per location
# Step 2: compute mean and SE of those 6 location-level means per genotype
#
# SE = SD / sqrt(n), with n = 6 locations per genotype

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
    Irr_se    = sd(Irrigated,   na.rm = TRUE) / sqrt(n()),   # SE of irrigated
    n_loc     = n(),
    .groups   = "drop"
  )


# ── 5.  Pie polygon builder  ──────────────────────────────────────────────────
#
# Builds two polygon slices per site × genotype:
#   "Rainfed" slice — angle proportional to rainfed / irrigated
#   "Gap"     slice — remaining angle
#
# Both use the SAME two colours regardless of genotype.
# y-radius is corrected so circles are visually round with coord_fixed(MAP_RATIO).

make_pie_polys <- function(row, ratio = MAP_RATIO, n_pts = PIE_N_PTS) {
  x0   <- row[["lon"]]
  y0   <- row[["lat"]]
  r    <- row[["r_pie"]]
  ry   <- r / ratio          # corrected y-radius → visually round circle
  frac <- row[["frac"]]
  geno <- as.character(row[["Genotype"]])
  loc  <- row[["Location"]]

  n1 <- max(3L, ceiling(n_pts * frac))
  n2 <- max(3L, ceiling(n_pts * (1 - frac)))

  t_rain <- seq(-pi / 2, -pi / 2 + frac * 2 * pi,       length.out = n1)
  t_gap  <- seq(-pi / 2 + frac * 2 * pi, -pi / 2 + 2 * pi, length.out = n2)

  rbind(
    data.frame(
      Genotype = geno, Location = loc,
      slice = "Rainfed",
      px    = c(x0, x0 + r * cos(t_rain), x0),
      py    = c(y0, y0 + ry * sin(t_rain), y0),
      grp   = paste(geno, loc, "rain", sep = "_"),
      stringsAsFactors = FALSE
    ),
    data.frame(
      Genotype = geno, Location = loc,
      slice = "Gap",
      px    = c(x0, x0 + r * cos(t_gap), x0),
      py    = c(y0, y0 + ry * sin(t_gap), y0),
      grp   = paste(geno, loc, "gap", sep = "_"),
      stringsAsFactors = FALSE
    )
  )
}

pie_polys <- do.call(rbind,
  lapply(seq_len(nrow(yield_wide)), function(i) make_pie_polys(yield_wide[i, ]))
) |>
  mutate(Genotype = factor(Genotype, levels = c("LT", "Non-LT")))


# ── 6.  US base map layers  ───────────────────────────────────────────────────

us_states    <- map_data("state")
focus_states <- c("arkansas", "missouri", "nebraska")
state_bg     <- filter(us_states, region %in% focus_states)


# ── 7.  Nature-style themes  ──────────────────────────────────────────────────

theme_nature_map <- function(base_size = 8) {
  theme_bw(base_size = base_size, base_family = "sans") +
  theme(
    axis.text         = element_blank(),
    axis.ticks        = element_blank(),
    axis.title        = element_blank(),
    panel.grid        = element_blank(),
    panel.border      = element_rect(colour = "grey40", linewidth = 0.35),
    # White strip with dark text (instead of dark strip)
    strip.background  = element_rect(fill = "white", colour = "grey40",
                                     linewidth = 0.35),
    strip.text        = element_text(colour = "grey10", face = "bold",
                                     size = base_size),
    legend.position   = "top",
    legend.direction  = "horizontal",
    legend.title      = element_blank(),
    legend.text       = element_text(size = base_size),
    legend.key.size   = unit(0.55, "lines"),
    legend.key        = element_rect(colour = "grey40", linewidth = 0.3),
    legend.spacing.x  = unit(0.3, "cm"),
    legend.margin     = margin(0, 0, 3, 0),
    plot.margin       = margin(1, 2, 1, 2, "mm"),
    panel.spacing     = unit(0.4, "lines")
  )
}

theme_nature_bar <- function(base_size = 8) {
  theme_nature_map(base_size) +
  theme(
    axis.text.x        = element_text(size = base_size, colour = "grey15",
                                      face = "bold"),
    axis.text.y        = element_text(size = base_size - 1, colour = "grey20"),
    axis.title.y       = element_text(size = base_size, colour = "grey10",
                                      angle = 90, margin = margin(r = 3)),
    axis.ticks         = element_line(linewidth = 0.3, colour = "grey40"),
    axis.ticks.length  = unit(1.5, "mm"),
    panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey88"),
    legend.position    = "none"    # legend lives on the map panel
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# PANEL (A) — MAP
#
# Two facets (LT | Non-LT).
# Surrounding states: light grey background.
# Study states (AR, MO, NE): slightly darker grey — stands out from surroundings.
# Pie charts: no text labels, black borders.
# ══════════════════════════════════════════════════════════════════════════════

map_panel <- ggplot() +

  # ── Background: all US states (very light grey) ───────────────────────────
  geom_polygon(
    data = us_states,
    aes(x = long, y = lat, group = group),
    fill = "grey90", colour = "grey72", linewidth = 0.18
  ) +

  # ── Study states: AR / MO / NE — darker than surroundings ────────────────
  geom_polygon(
    data = state_bg,
    aes(x = long, y = lat, group = group),
    fill = "grey68", colour = "grey25", linewidth = 0.45
  ) +

  # ── Pie slice polygons (black border, shared 2-colour fill) ───────────────
  geom_polygon(
    data = pie_polys,
    aes(x = px, y = py, group = grp, fill = slice),
    colour = "black", linewidth = 0.25
  ) +

  # ── Facets ─────────────────────────────────────────────────────────────────
  facet_wrap(~Genotype, ncol = 2,
             labeller = labeller(Genotype = GENO_LABELS)) +

  # ── Shared 2-level fill scale ──────────────────────────────────────────────
  scale_fill_manual(
    values = FILL_VALS,
    labels = FILL_LABELS,
    breaks = FILL_BREAKS,
    name   = NULL
  ) +

  coord_fixed(MAP_RATIO, xlim = c(-106, -84), ylim = c(31.0, 46.5)) +

  theme_nature_map() +
  guides(fill = guide_legend(nrow = 1, override.aes = list(colour = "grey30")))


# ══════════════════════════════════════════════════════════════════════════════
# PANEL (B) — BAR CHART
#
# Mean yield across 6 locations per genotype.
# Bars stacked: water-limited (rainfed) + gap.
# Error bar on the total bar (irrigated mean ± SE across 6 locations).
# No horizontal cap whiskers (width = 0).
# Y-axis in kg ha⁻¹.
# ══════════════════════════════════════════════════════════════════════════════

# Stacked bar data (two rows per genotype: Rainfed and Gap)
bar_long <- bar_summary |>
  pivot_longer(
    cols      = c(Rain_mean, Gap_mean),
    names_to  = "slice",
    values_to = "Yield_kgha"
  ) |>
  mutate(
    slice = recode(slice, Rain_mean = "Rainfed", Gap_mean = "Gap"),
    slice = factor(slice, levels = c("Rainfed", "Gap"))   # Rainfed on bottom
  )

bar_panel <- ggplot() +

  # ── Stacked bars ──────────────────────────────────────────────────────────
  geom_col(
    data  = bar_long,
    aes(x = Genotype, y = Yield_kgha, fill = slice),
    width = 0.55, colour = NA
  ) +

  # ── Error bar at the total bar top (irrigated mean ± SE) ─────────────────
  # Positioned at Irr_mean so the bar touches the middle of the error bar.
  # width = 0 removes the horizontal cap whiskers.
  geom_errorbar(
    data  = bar_summary,
    aes(x    = Genotype,
        ymin  = Irr_mean - Irr_se,
        ymax  = Irr_mean + Irr_se),
    width     = 0,          # no horizontal whisker caps
    linewidth = 0.6,
    colour    = "grey10"
  ) +

  scale_fill_manual(values = FILL_VALS, guide = "none") +

  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1)),
    name   = expression(Yield~(kg~ha^{-1})),
    labels = scales::comma,       # add thousands separator: 1,000 / 4,000 etc.
    breaks = seq(0, 8000, by = 2000)
  ) +

  labs(x = NULL) +

  theme_nature_bar()


# ══════════════════════════════════════════════════════════════════════════════
# COMBINE with patchwork
# Map (3/4 width) + bar chart (1/4 width); panel tags (A) and (B)
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
# TIFF at 600 DPI (LZW) — submission quality; *.tiff in .gitignore (not tracked)
# PNG  at 300 DPI       — git-tracked for sharing
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
