# ==============================================================================
# Fig_map_marin.R
#
# Purpose:
#   Marin et al. (2022)-style yield map for the LT-trait soybean analysis.
#   Two side-by-side US maps (LT | Non-LT), each with a pie chart at every
#   simulation site showing how much of the potential (irrigated) yield is
#   achieved under rainfed (water-limited) conditions.
#   A side bar panel shows mean yield ± SD by state for both genotypes.
#
# Figure layout:
#   ┌─────────────────────────────────────┬────────────┐
#   │  Map panel (faceted: LT | Non-LT)  │  Bar panel │
#   │  Pie at each of 6 sites:           │  by state  │
#   │    dark slice  = rainfed yield     │  (AR/MO/NE)│
#   │    light slice = yield gap         │            │
#   │    pie area ∝ potential yield      │            │
#   └─────────────────────────────────────┴────────────┘
#
# Data:   modeling_data/data/SoybeanSSM_sims2.xlsx
# Output: modeling_data/figures/Fig_map_marin.tiff  (600 DPI, LZW)
#
# Run from the LT-trait project root directory.
# On Windows: set project_root below before running.
#
# Required packages (install if missing):
#   install.packages(c("readxl","dplyr","tidyr","ggplot2",
#                      "ggforce","patchwork","maps","ggrepel"))
#
# Author: Elvis F. Elli | eelli@uark.edu
# ==============================================================================


# ── 0.  Working directory  ────────────────────────────────────────────────────
#
# Set project_root to the folder that CONTAINS 'modeling_data/'.
#   Windows example: project_root <- "C:/Users/YourName/Documents/LT-trait"
#   Linux / Mac:     the default below is used automatically.

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
  library(ggforce)   # geom_arc_bar  (pie slices in map coordinates)
  library(patchwork) # combine map panel + bar panel
  library(maps)      # US state boundary polygons
  library(ggrepel)   # non-overlapping site-name labels
})


# ── 2.  Design constants  ─────────────────────────────────────────────────────
#
# Colour convention follows the project standard:
#   LT      (Limited Transpiration)  =  orange family
#   Non-LT  (Non-limited)            =  green  family
#
# Dark colour = water-limited (rainfed, actual) yield
# Light colour = yield gap (potential − actual)

COL_LT_RAIN  <- "#d95f02"   # LT   rainfed slice  (dark orange)
COL_LT_GAP   <- "#fdbe85"   # LT   gap slice      (light orange)
COL_NLT_RAIN <- "#31a354"   # Non-LT rainfed slice (dark green)
COL_NLT_GAP  <- "#a1d99b"   # Non-LT gap slice    (light green)

# Named fill vector — used by BOTH map and bar panels so colours are identical
FILL_VALS <- c(
  "LT|Rainfed"      = COL_LT_RAIN,
  "LT|Gap"          = COL_LT_GAP,
  "Non-LT|Rainfed"  = COL_NLT_RAIN,
  "Non-LT|Gap"      = COL_NLT_GAP
)

# Strip labels for the two map panels
GENO_LABELS <- c(
  "LT"     = "Limited transpiration (LT)",
  "Non-LT" = "Non-limited transpiration"
)

PIE_MAX_R <- 1.55   # maximum pie radius in map degrees (for the highest-yield site)


# ── 3.  Site coordinates  ─────────────────────────────────────────────────────

sites <- tibble(
  Location = c("Jonesboro", "Marianna", "Columbia",
               "Norborne",  "Lincoln",  "North Platte"),
  State    = c("AR", "AR", "MO", "MO", "NE", "NE"),
  lat      = c( 35.84,  34.78,  38.95,  39.30,  40.81,  41.12),
  lon      = c(-90.70, -90.76, -92.33, -93.67, -96.70, -100.77)
)

# Short site codes for labelling (printed below each pie)
site_codes <- c(
  "Jonesboro"   = "JO",  "Marianna"    = "MA",
  "Columbia"    = "CO",  "Norborne"    = "NO",
  "Lincoln"     = "LI",  "North Platte" = "NP"
)


# ── 4.  Load and summarise SSM data  ─────────────────────────────────────────
#
# We use only MG4 and MG3 crops (the two maturity groups in the study).
# LT genotypes: SOYBEAN-LT1.5-MG4 / SOYBEAN-LT1.5-MG3
# Non-LT:       SOYBEAN-MG4        / SOYBEAN-MG3
#
# Summary: 30-year mean yield per Location × Genotype × water_mgmt
# (averaged over simulation years and planting dates)

ssm_raw <- read_excel("modeling_data/data/SoybeanSSM_sims2.xlsx")

ssm <- ssm_raw |>
  filter(Crop %in% c("SOYBEAN-MG4",      "SOYBEAN-LT1.5-MG4",
                     "SOYBEAN-LT1.5-MG3", "SOYBEAN-MG3")) |>
  mutate(
    Genotype   = factor(
      ifelse(Crop %in% c("SOYBEAN-MG4", "SOYBEAN-MG3"), "Non-LT", "LT"),
      levels = c("LT", "Non-LT")   # LT first so LT panel appears on the left
    ),
    water_mgmt = if_else(water_mgmt == "Rainfed", "Rainfed", "Irrigated")
  )

# ── 4a. Site-level mean yield (used to build pie charts and map labels)
yield_site <- ssm |>
  group_by(Location, Genotype, water_mgmt) |>
  summarise(Yield_mean = mean(Yield_kgha, na.rm = TRUE), .groups = "drop") |>
  left_join(sites, by = "Location")

# Wide: one row per Location × Genotype  →  columns Rainfed and Irrigated
yield_wide <- yield_site |>
  pivot_wider(names_from = water_mgmt, values_from = Yield_mean) |>
  mutate(
    Gap    = Irrigated - Rainfed,            # exploitable yield gap
    frac   = Rainfed / Irrigated,            # fraction of potential achieved
    r_pie  = sqrt(Irrigated) /               # pie radius ∝ √(potential yield)
               sqrt(max(Irrigated, na.rm = TRUE)) * PIE_MAX_R,
    code   = site_codes[Location]            # short code for label
  )

# ── 4b. Year-level SD for bar-chart error bars
#        (variability across all years × planting dates within state)
err_irr <- ssm |>
  filter(water_mgmt == "Irrigated") |>
  left_join(sites |> select(Location, State), by = "Location") |>
  group_by(State, Genotype) |>
  summarise(
    Irr_mean = mean(Yield_kgha, na.rm = TRUE),
    Irr_sd   = sd(Yield_kgha,   na.rm = TRUE),
    .groups  = "drop"
  )


# ── 5.  Build arc data for geom_arc_bar (pie slices)  ─────────────────────────
#
# Each pie = two arcs drawn in lon/lat coordinates:
#   Slice "Rainfed" : from top (−π/2) clockwise by frac × 2π
#   Slice "Gap"     : the remaining arc to complete the circle
#
# Faceting by Genotype means the LT panel shows only LT arcs (and vice versa).

arcs <- bind_rows(
  yield_wide |>
    mutate(slice = "Rainfed",
           start = -pi / 2,
           end   = -pi / 2 + frac * 2 * pi),
  yield_wide |>
    mutate(slice = "Gap",
           start = -pi / 2 + frac * 2 * pi,
           end   = -pi / 2 + 2 * pi)
) |>
  mutate(fill_key = paste(as.character(Genotype), slice, sep = "|"))


# ── 6.  Yield value labels printed around each pie  ──────────────────────────
#
# White text inside rainfed slice  = rainfed yield (t ha⁻¹)
# Grey text above pie              = irrigated (potential) yield (t ha⁻¹)

pie_labels <- yield_wide |>
  mutate(
    lbl_rain = paste0(round(Rainfed   / 1000, 1), " t"),
    lbl_irr  = paste0(round(Irrigated / 1000, 1), " t")
  )


# ── 7.  US base-map layers  ──────────────────────────────────────────────────

us_states    <- map_data("state")
focus_states <- c("arkansas", "missouri", "nebraska")
state_bg     <- filter(us_states, region %in% focus_states)

# State name positions (placed inside each state's territory)
state_labels <- tibble(
  lon   = c(-91.5,  -93.5,  -99.3),
  lat   = c( 33.5,   36.7,   43.1),
  label = c("Arkansas", "Missouri", "Nebraska")
)


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE PANEL A — MAP with pie charts
# facet_wrap(~Genotype):  left panel = LT,  right panel = Non-LT
# ══════════════════════════════════════════════════════════════════════════════

map_panel <- ggplot() +

  # ── Base map: all US states (light grey) ─────────────────────────────────
  geom_polygon(
    data = us_states,
    aes(x = long, y = lat, group = group),
    fill = "grey97", colour = "grey80", linewidth = 0.2
  ) +

  # ── Focus states: AR, MO, NE (slightly tinted) ───────────────────────────
  geom_polygon(
    data = state_bg,
    aes(x = long, y = lat, group = group),
    fill = "#edf3fb", colour = "grey55", linewidth = 0.4
  ) +

  # ── State name labels ─────────────────────────────────────────────────────
  # These have no Genotype column → they appear in both facet panels
  geom_text(
    data = state_labels,
    aes(x = lon, y = lat, label = label),
    fontface = "bold.italic", size = 2.7, colour = "grey35",
    inherit.aes = FALSE
  ) +

  # ── Pie slices (geom_arc_bar from ggforce) ────────────────────────────────
  # arcs has a Genotype column → each panel shows only its genotype's pies
  geom_arc_bar(
    data = arcs,
    aes(x0 = lon, y0 = lat,
        r0 = 0, r = r_pie,
        start = start, end = end,
        fill  = fill_key),
    colour = "white", linewidth = 0.35
  ) +

  # ── Rainfed yield label — inside the rainfed slice (white, bold) ──────────
  # Placed at ~42 % of the radius along the mid-angle of the rainfed arc
  geom_text(
    data = pie_labels,
    aes(
      x     = lon + r_pie * 0.42 * cos(-pi / 2 + frac * pi),
      y     = lat + r_pie * 0.42 * sin(-pi / 2 + frac * pi),
      label = lbl_rain
    ),
    size = 1.9, colour = "white", fontface = "bold", inherit.aes = FALSE
  ) +

  # ── Potential (irrigated) yield label — above each pie ────────────────────
  geom_text(
    data = pie_labels,
    aes(x = lon, y = lat + r_pie + 0.4, label = lbl_irr),
    size = 1.8, colour = "grey15", inherit.aes = FALSE
  ) +

  # ── Short site-code below each pie (e.g. "JO", "MA") ─────────────────────
  # site_codes data has no Genotype → the code appears in both panels,
  # which is intentional (same 6 sites in each panel)
  geom_text(
    data = yield_wide |> distinct(Location, lon, lat, r_pie, code),
    aes(x = lon, y = lat - r_pie - 0.5, label = code),
    size = 2.2, colour = "grey25", fontface = "bold", inherit.aes = FALSE
  ) +

  # ── Facet: left = LT, right = Non-LT ─────────────────────────────────────
  facet_wrap(~Genotype, ncol = 2,
             labeller = labeller(Genotype = GENO_LABELS)) +

  # ── Fill scale (4 levels: 2 genotypes × 2 slices) ────────────────────────
  scale_fill_manual(
    values = FILL_VALS,
    labels = c(
      "LT|Rainfed"     = "LT — water-limited (rainfed)",
      "LT|Gap"         = "LT — yield gap to potential",
      "Non-LT|Rainfed" = "Non-LT — water-limited (rainfed)",
      "Non-LT|Gap"     = "Non-LT — yield gap to potential"
    ),
    breaks = c("LT|Rainfed", "LT|Gap",
               "Non-LT|Rainfed", "Non-LT|Gap"),
    name   = NULL
  ) +

  coord_fixed(1.3, xlim = c(-106, -84), ylim = c(31.5, 46.0)) +
  labs(x = NULL, y = NULL) +

  theme_bw(base_size = 10) +
  theme(
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    panel.grid       = element_blank(),
    # Strip background: dark with white bold text
    strip.background = element_rect(fill = "grey20", colour = "grey20"),
    strip.text       = element_text(colour = "white", face = "bold", size = 9.5),
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.title     = element_blank(),
    legend.text      = element_text(size = 8.5),
    legend.key.size  = unit(0.7, "lines"),
    panel.border     = element_rect(colour = "grey60"),
    panel.spacing    = unit(0.5, "lines")
  ) +
  guides(
    fill = guide_legend(nrow = 2, override.aes = list(colour = NA))
  )


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE PANEL B — BAR CHART  (yield magnitude by state)
#
# Stacked bars: rainfed (bottom, dark) + yield gap (top, light)
# Error bar: ± SD of irrigated yield across all years × planting dates
#            at sites within that state
# ══════════════════════════════════════════════════════════════════════════════

# ── State-level means for stacked bars ───────────────────────────────────────
bar_means <- yield_site |>
  # yield_site already has State from the sites join in section 4a
  group_by(State, Genotype, water_mgmt) |>
  summarise(Yield_mean = mean(Yield_mean, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = water_mgmt, values_from = Yield_mean) |>
  mutate(Gap = Irrigated - Rainfed)

# Long format for stacked geom_col (two rows per State × Genotype)
bar_long <- bar_means |>
  pivot_longer(cols = c(Rainfed, Gap),
               names_to  = "slice",
               values_to = "Yield_kgha") |>
  mutate(
    fill_key = paste(as.character(Genotype), slice, sep = "|"),
    State    = factor(State, levels = c("AR", "MO", "NE")),
    # Rainfed stacks first (bottom), Gap on top
    slice    = factor(slice, levels = c("Rainfed", "Gap"))
  )

# ── Error bars on the irrigated (potential) total  ────────────────────────────
# Use state-level SD from err_irr (computed in section 4b)
err_bar <- err_irr |>
  mutate(
    State    = factor(State, levels = c("AR", "MO", "NE")),
    Genotype = factor(Genotype, levels = c("LT", "Non-LT"))
  )

bar_panel <- ggplot() +

  # ── Stacked bars (rainfed dark, gap light) ────────────────────────────────
  geom_col(
    data  = bar_long,
    aes(x = Genotype, y = Yield_kgha / 1000, fill = fill_key),
    width = 0.60, colour = "white", linewidth = 0.3
  ) +

  # ── Error bars (SD of irrigated yield: across years × sites in state) ─────
  geom_errorbar(
    data  = err_bar,
    aes(x    = Genotype,
        ymin  = (Irr_mean - Irr_sd) / 1000,
        ymax  = (Irr_mean + Irr_sd) / 1000),
    width     = 0.22,
    linewidth = 0.55,
    colour    = "grey20"
  ) +

  # ── Rainfed yield value printed inside the bar (lower half) ───────────────
  geom_text(
    data = bar_means,
    aes(x     = Genotype,
        y     = Rainfed / 1000 / 2,
        label = paste0(round(Rainfed / 1000, 1), " t")),
    size = 2.3, colour = "white", fontface = "bold"
  ) +

  scale_fill_manual(values = FILL_VALS, guide = "none") +

  scale_y_continuous(
    expand = expansion(mult = c(0, 0.14)),
    name   = expression(Yield~(t~ha^{-1})),
    breaks = seq(0, 8, by = 2)
  ) +

  # One panel per state (AR / MO / NE) stacked vertically
  facet_wrap(~State, ncol = 1, strip.position = "right") +

  labs(x = NULL, caption = "Error bars: ± SD across years and sites") +

  theme_bw(base_size = 9) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    # Right-side strip (state names)
    strip.background   = element_rect(fill = "grey20", colour = "grey20"),
    strip.text         = element_text(colour = "white", face = "bold",
                                      size = 8.5, angle = 0),
    axis.text.x        = element_text(size = 8),
    panel.border       = element_rect(colour = "grey60"),
    panel.spacing      = unit(0.35, "lines"),
    plot.caption       = element_text(size = 6.5, colour = "grey45",
                                      hjust = 0.5)
  )


# ══════════════════════════════════════════════════════════════════════════════
# COMBINE panels with patchwork
# Map takes 3/4 of total width; bar chart takes 1/4
# ══════════════════════════════════════════════════════════════════════════════

combined <- map_panel + bar_panel +
  plot_layout(widths = c(3, 1)) +
  plot_annotation(
    title    = "Water-limited and potential soybean yield at six US locations",
    subtitle = paste0(
      "Pie area \u221d irrigated (potential) yield  \u00b7  ",
      "Dark slice = rainfed yield (t ha\u207b\u00b9, labelled in white)  \u00b7  ",
      "Light slice = yield gap (irrigated \u2212 rainfed)  \u00b7  ",
      "Bar: state mean (n = sites \u00d7 30 years)"
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 12, hjust = 0),
      plot.subtitle = element_text(colour = "grey40", size = 8.5, hjust = 0)
    )
  )


# ══════════════════════════════════════════════════════════════════════════════
# SAVE
#
# (a) TIFF at 600 DPI — full-quality submission file
#     Note: *.tiff is listed in .gitignore for this repository (modeling_data
#     follows a PNG-only convention for git tracking).  The TIFF is written to
#     disk and usable locally but will not be committed.
#
# (b) PNG at 300 DPI — tracked in git for sharing / record-keeping
#
# tiff() type = "cairo" on Linux/Mac (headless, no X11 needed)
#              "windows" on Windows
# ══════════════════════════════════════════════════════════════════════════════

tiff_type <- if (.Platform$OS.type == "windows") "windows" else "cairo"

# (a) TIFF — high-resolution, not git-tracked
tiff_file <- "modeling_data/figures/Fig_map_marin.tiff"
tiff(tiff_file,
     width   = 30, height = 15, units = "cm",
     res     = 600, compression = "lzw",
     type    = tiff_type)
print(combined)
dev.off()
message("Saved TIFF: ", tiff_file)

# (b) PNG — git-tracked
png_file <- "modeling_data/figures/Fig_map_marin.png"
ggsave(png_file, combined,
       width = 30, height = 15, units = "cm",
       dpi   = 300, bg = "white")
message("Saved PNG:  ", png_file)
