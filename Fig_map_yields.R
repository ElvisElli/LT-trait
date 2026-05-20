## Map figure: water-limited (rainfed) vs potential (irrigated) yield
## as pie charts at AR, MO, NE sites — LT vs non-LT soybean

setwd("/home/user/LT-trait")

library(dplyr)
library(tidyr)
library(ggplot2)
library(maps)
library(readxl)

## ── Helper: generate pie polygon coordinates ──────────────────────────────────

make_pie <- function(cx, cy, r, values, labels, n_pts = 100) {
  total  <- sum(values)
  angles <- cumsum(c(0, values / total * 2 * pi)) - pi / 2  # start from top
  purrr::map_dfr(seq_along(values), function(i) {
    theta <- seq(angles[i], angles[i + 1], length.out = n_pts)
    tibble(
      x     = c(cx, cx + r * cos(theta), cx),
      y     = c(cy, cy + r * sin(theta), cy),
      slice = labels[i]
    )
  })
}

## ── Site coordinates ──────────────────────────────────────────────────────────

sites <- tibble(
  State = c("AR",   "MO",   "NE"),
  lat   = c(36.08,  38.95,  40.82),
  lon   = c(-94.17, -92.33, -96.70),
  label = c("Fayetteville\nAR", "Columbia\nMO", "Lincoln\nNE")
)

## ── Yield data from SSM ───────────────────────────────────────────────────────

ssm.model <- read_excel("SoybeanSSM_sims2.xlsx")

yield_data <- ssm.model %>%
  filter(Crop %in% c("SOYBEAN-MG4", "SOYBEAN-LT1.5-MG4",
                     "SOYBEAN-LT1.5-MG3", "SOYBEAN-MG3")) %>%
  mutate(
    Genotype = ifelse(Crop %in% c("SOYBEAN-MG4", "SOYBEAN-MG3"),
                      "Non-LT", "LT"),
    State = case_when(
      grepl("AR", Manag, ignore.case = TRUE) ~ "AR",
      grepl("MO", Manag, ignore.case = TRUE) ~ "MO",
      grepl("NE", Manag, ignore.case = TRUE) ~ "NE",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(State)) %>%
  group_by(State, Genotype, water_mgmt) %>%
  summarise(Yield = mean(Yield_kgha, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = water_mgmt, values_from = Yield) %>%
  left_join(sites, by = "State") %>%
  ## offset LT west, non-LT east so pies don't overlap
  mutate(lon_pie = ifelse(Genotype == "LT", lon - 1.4, lon + 1.4),
         ## radius scaled to potential (irrigated) yield
         radius  = Irrigated / max(Irrigated) * 1.6)

## ── Build pie polygons ────────────────────────────────────────────────────────

pie_polys <- yield_data %>%
  rowwise() %>%
  mutate(poly = list(
    make_pie(
      cx     = lon_pie,
      cy     = lat,
      r      = radius,
      values = c(Rainfed, Irrigated - Rainfed),   # rainfed + extra for irrigated
      labels = c("Water-limited\n(Rainfed)", "Extra potential\n(Irrigated)")
    ) %>% mutate(Genotype = Genotype, State = State)
  )) %>%
  ungroup() %>%
  select(poly) %>%
  unnest(poly)

## Assign fill colours per genotype × slice
pie_polys <- pie_polys %>%
  mutate(fill_key = paste(Genotype, slice, sep = " — "))

fill_vals <- c(
  "Non-LT — Water-limited\n(Rainfed)"     = "#74c476",   # light green
  "Non-LT — Extra potential\n(Irrigated)"  = "#006d2c",   # dark green
  "LT — Water-limited\n(Rainfed)"         = "#fdae6b",   # light orange
  "LT — Extra potential\n(Irrigated)"     = "#d95f02"    # dark orange
)

## ── Base map ─────────────────────────────────────────────────────────────────

us_states <- map_data("state")
focus_states <- c("arkansas", "missouri", "nebraska")

## crop to relevant extent
xlim <- c(-105, -84)
ylim <- c(32,   46)

p_map <- ggplot() +

  ## all CONUS states — light grey background
  geom_polygon(data = us_states,
               aes(x = long, y = lat, group = group),
               fill = "grey93", colour = "white", linewidth = 0.4) +

  ## highlight the three study states
  geom_polygon(data = filter(us_states, region %in% focus_states),
               aes(x = long, y = lat, group = group),
               fill = "grey82", colour = "white", linewidth = 0.4) +

  ## pie polygons
  geom_polygon(data = pie_polys,
               aes(x = x, y = y, group = interaction(Genotype, State, slice),
                   fill = fill_key),
               colour = "grey20", linewidth = 0.25) +

  ## site labels
  geom_text(data = sites,
            aes(x = lon, y = lat - 2.4, label = label),
            size = 3, fontface = "bold", colour = "grey20") +

  ## genotype labels above each pie
  geom_text(data = yield_data,
            aes(x = lon_pie, y = lat + radius + 0.35,
                label = ifelse(Genotype == "LT",
                               "Limited\ntranspiration",
                               "Non-limited\ntranspiration")),
            size = 2.4, colour = "grey20", linewidth = 0.3) +

  scale_fill_manual(
    values = fill_vals,
    labels = c("Non-LT — Water-limited\n(Rainfed)"    = "Non-LT  |  Water-limited (rainfed)",
               "Non-LT — Extra potential\n(Irrigated)" = "Non-LT  |  Potential gain (irrigated)",
               "LT — Water-limited\n(Rainfed)"        = "LT  |  Water-limited (rainfed)",
               "LT — Extra potential\n(Irrigated)"    = "LT  |  Potential gain (irrigated)"),
    breaks = c("Non-LT — Water-limited\n(Rainfed)",
               "Non-LT — Extra potential\n(Irrigated)",
               "LT — Water-limited\n(Rainfed)",
               "LT — Extra potential\n(Irrigated)")
  ) +

  coord_fixed(1.3, xlim = xlim, ylim = ylim) +

  labs(
    x    = NULL, y = NULL,
    fill = NULL,
    title    = "Water-limited vs. potential soybean yield",
    subtitle = "Pie area ∝ potential (irrigated) yield  |  Inner slice = water-limited (rainfed) yield"
  ) +

  theme_bw(base_size = 11) +
  theme(
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    panel.grid       = element_blank(),
    legend.position  = "bottom",
    legend.key.size  = unit(0.7, "lines"),
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(colour = "grey40", size = 9)
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

ggsave("Fig_map_yields.png",  plot = p_map,
       width = 20, height = 16, units = "cm", dpi = 150, bg = "white")
ggsave("Fig_map_yields.tiff", plot = p_map,
       width = 20, height = 16, units = "cm", dpi = 300,
       compression = "lzw", bg = "white")
message("Saved: Fig_map_yields (.png + .tiff)")
