## Map figure: water-limited (rainfed) vs potential (irrigated) yield
## Pie charts via ggforce::geom_arc_bar at 6 locations — LT vs non-LT soybean

setwd("/home/user/LT-trait")

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggforce)
library(maps)
library(readxl)

## ── Site coordinates (2 per state) ───────────────────────────────────────────

sites <- tibble(
  Location = c("Jonesboro", "Marianna",
                "Columbia",  "Norborne",
                "Lincoln",   "North Platte"),
  State    = c("AR", "AR", "MO", "MO", "NE", "NE"),
  lat      = c(35.84, 34.78, 38.95, 39.30, 40.81, 41.12),
  lon      = c(-90.70, -90.76, -92.33, -93.67, -96.70, -100.77)
)

## ── Yield data ────────────────────────────────────────────────────────────────

ssm.model <- read_excel("SoybeanSSM_sims2.xlsx")

yield_data <- ssm.model %>%
  filter(Crop %in% c("SOYBEAN-MG4", "SOYBEAN-LT1.5-MG4",
                     "SOYBEAN-LT1.5-MG3", "SOYBEAN-MG3")) %>%
  mutate(Genotype = ifelse(Crop %in% c("SOYBEAN-MG4", "SOYBEAN-MG3"),
                           "Non-limited transpiration",
                           "Limited transpiration"),
         Genotype = factor(Genotype,
                           levels = c("Non-limited transpiration",
                                      "Limited transpiration"))) %>%
  group_by(Location, Genotype, water_mgmt) %>%
  summarise(Yield = mean(Yield_kgha, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = water_mgmt, values_from = Yield) %>%
  left_join(sites, by = "Location") %>%
  ## offset LT west, non-LT east so pies don't overlap
  mutate(lon_pie = ifelse(Genotype == "Limited transpiration",
                          lon - 1.2, lon + 1.2),
         ## radius in map-degree units, scaled to potential yield
         r = Irrigated / max(Irrigated) * 1.4)

## ── Build geom_arc_bar pie data ───────────────────────────────────────────────
## geom_arc_bar needs: x0, y0, r, r0, start, end, fill

make_arcs <- function(df) {
  df %>%
    rowwise() %>%
    mutate(
      total       = Irrigated,
      frac_rain   = Rainfed  / total,
      frac_irr    = 1 - frac_rain,
      ## water-limited slice (rainfed)
      start_rain  = -pi / 2,
      end_rain    = start_rain + frac_rain * 2 * pi,
      ## potential-gain slice (irrigated beyond rainfed)
      start_irr   = end_rain,
      end_irr     = start_rain + 2 * pi
    ) %>%
    ungroup() %>%
    ## pivot to one row per slice
    select(Location, Genotype, State, lat, lon, lon_pie, r,
           Rainfed, Irrigated,
           start_rain, end_rain, start_irr, end_irr) %>%
    pivot_longer(cols = c(start_rain, start_irr),
                 names_to  = "slice_start_col",
                 values_to = "start") %>%
    mutate(
      end   = ifelse(slice_start_col == "start_rain", end_rain, end_irr),
      slice = ifelse(slice_start_col == "start_rain",
                     "Water-limited (Rainfed)",
                     "Potential gain (Irrigated)")
    ) %>%
    select(Location, Genotype, State, lat, lon, lon_pie, r,
           Rainfed, Irrigated, start, end, slice)
}

arc_data <- make_arcs(yield_data)

## ── Colour palette ────────────────────────────────────────────────────────────

fill_colours <- c(
  "Non-limited transpiration|Water-limited (Rainfed)"   = "#a1d99b",
  "Non-limited transpiration|Potential gain (Irrigated)" = "#238b45",
  "Limited transpiration|Water-limited (Rainfed)"       = "#fdae6b",
  "Limited transpiration|Potential gain (Irrigated)"    = "#d95f02"
)

arc_data <- arc_data %>%
  mutate(fill_key = paste(Genotype, slice, sep = "|"))

## ── Bubble summary for geom_point legend ─────────────────────────────────────

bubble_ref <- tibble(
  lon_pie = -104, lat = c(34.5, 36.5, 38.5),
  r_raw   = yield_data$Irrigated[c(1,3,5)],   # example yields
  r       = yield_data$r[c(1,3,5)],
  label   = paste0(round(yield_data$Irrigated[c(1,3,5)]), " kg/ha")
)

## ── Map layers ────────────────────────────────────────────────────────────────

us_states   <- map_data("state")
focus_states <- c("arkansas", "missouri", "nebraska")

p_map <- ggplot() +

  ## base CONUS
  geom_polygon(data = us_states,
               aes(x = long, y = lat, group = group),
               fill = "grey95", colour = "white", linewidth = 0.35) +

  ## highlight study states
  geom_polygon(data = filter(us_states, region %in% focus_states),
               aes(x = long, y = lat, group = group),
               fill = "#deebf7", colour = "white", linewidth = 0.35) +

  ## site dot (geom_point) — size proportional to potential (irrigated) yield,
  ## drawn beneath pies as a subtle anchor
  geom_point(data = yield_data,
             aes(x = lon_pie, y = lat,
                 size  = Irrigated,
                 shape = Genotype),
             colour = "grey40", fill = NA, stroke = 0.8,
             show.legend = TRUE) +

  ## pie charts via geom_arc_bar
  geom_arc_bar(data = arc_data,
               aes(x0 = lon_pie, y0 = lat,
                   r0 = 0, r = r,
                   start = start, end = end,
                   fill  = fill_key),
               colour = "grey25", linewidth = 0.2) +

  ## state labels
  annotate("text",
           x = c(-92.5, -93.5, -99.5),
           y = c(34.1,  37.0,  42.8),
           label = c("Arkansas", "Missouri", "Nebraska"),
           fontface = "bold", size = 3.2, colour = "grey30") +

  ## location labels
  geom_text(data = sites,
            aes(x = lon, y = lat, label = Location),
            nudge_y = -1.05, size = 2.5, colour = "grey20", fontface = "italic") +

  ## scales
  scale_fill_manual(
    values = fill_colours,
    labels = c(
      "Non-limited transpiration|Water-limited (Rainfed)"    = "Non-LT  —  Water-limited (rainfed)",
      "Non-limited transpiration|Potential gain (Irrigated)"  = "Non-LT  —  Potential gain (irrigated)",
      "Limited transpiration|Water-limited (Rainfed)"        = "LT  —  Water-limited (rainfed)",
      "Limited transpiration|Potential gain (Irrigated)"     = "LT  —  Potential gain (irrigated)"
    ),
    breaks = c(
      "Non-limited transpiration|Water-limited (Rainfed)",
      "Non-limited transpiration|Potential gain (Irrigated)",
      "Limited transpiration|Water-limited (Rainfed)",
      "Limited transpiration|Potential gain (Irrigated)"
    ),
    name = "Yield component"
  ) +

  scale_size_continuous(
    name   = "Potential yield (kg ha⁻¹)",
    range  = c(2, 8),
    breaks = c(2000, 3000, 4000),
    labels = c("2 000", "3 000", "4 000")
  ) +

  scale_shape_manual(
    name   = "Genotype",
    values = c("Non-limited transpiration" = 21,
               "Limited transpiration"      = 24)
  ) +

  coord_fixed(1.3, xlim = c(-105, -84), ylim = c(32, 46)) +

  labs(
    x        = NULL, y = NULL,
    title    = "Water-limited vs. potential soybean yield across six locations",
    subtitle = "Pie area ∝ potential (irrigated) yield  ·  Inner slice = water-limited (rainfed) yield"
  ) +

  theme_bw(base_size = 11) +
  theme(
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    panel.grid       = element_blank(),
    legend.position  = "bottom",
    legend.box       = "vertical",
    legend.key.size  = unit(0.75, "lines"),
    legend.title     = element_text(face = "bold", size = 9),
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(colour = "grey40", size = 9),
    panel.border     = element_rect(colour = "grey70")
  ) +
  guides(
    fill  = guide_legend(order = 1, nrow = 2, byrow = TRUE,
                         override.aes = list(colour = NA)),
    size  = guide_legend(order = 2),
    shape = guide_legend(order = 3,
                         override.aes = list(size = 3, stroke = 1))
  )

ggsave("Fig_map_yields.png",  plot = p_map,
       width = 22, height = 17, units = "cm", dpi = 150, bg = "white")
ggsave("Fig_map_yields.tiff", plot = p_map,
       width = 22, height = 17, units = "cm", dpi = 300,
       compression = "lzw", bg = "white")
message("Saved: Fig_map_yields (.png + .tiff)")
