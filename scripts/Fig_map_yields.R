## Fig — Map: water-limited vs potential yield at 6 locations
## Part 1 (p):  pie-chart map (preserved, path fixed)
## Part 2 (p2): cleaner facet_grid map — 2×2 panels (primary output)

setwd("/home/user/LT-trait")
library(dplyr); library(tidyr); library(ggplot2); library(ggforce)
library(maps); library(readxl); library(ggrepel)

## ── Site coordinates ──────────────────────────────────────────────────────────
sites <- tibble(
  Location = c("Jonesboro","Marianna","Columbia","Norborne","Lincoln","North Platte"),
  State    = c("AR","AR","MO","MO","NE","NE"),
  lat      = c(35.84, 34.78, 38.95, 39.30, 40.81, 41.12),
  lon      = c(-90.70,-90.76,-92.33,-93.67,-96.70,-100.77)
)

## ── Yield data (shared) ───────────────────────────────────────────────────────
ssm <- read_excel("data/SoybeanSSM_sims2.xlsx")

yield_long <- ssm %>%
  filter(Crop %in% c("SOYBEAN-MG4","SOYBEAN-LT1.5-MG4",
                     "SOYBEAN-LT1.5-MG3","SOYBEAN-MG3")) %>%
  mutate(Genotype = factor(
    ifelse(Crop %in% c("SOYBEAN-MG4","SOYBEAN-MG3"),
           "Non-LT", "LT"),
    levels = c("Non-LT","LT")),
    water_mgmt = factor(
      ifelse(water_mgmt == "Rainfed", "Rainfed", "Irrigated"),
      levels = c("Rainfed","Irrigated"),
      labels = c("Rainfed (water-limited)","Irrigated (potential)")
    )
  ) %>%
  group_by(Location, Genotype, water_mgmt) %>%
  summarise(Yield = mean(Yield_kgha, na.rm = TRUE), .groups = "drop") %>%
  left_join(sites, by = "Location")

## ── Map base layers ───────────────────────────────────────────────────────────
us_states    <- map_data("state")
focus_states <- c("arkansas","missouri","nebraska")
state_bg     <- filter(us_states, region %in% focus_states)

## ── State label positions ─────────────────────────────────────────────────────
state_labels <- tibble(
  lon   = c(-91.5, -93.5, -99.0),
  lat   = c(33.6,   36.7,  43.3),
  label = c("Arkansas", "Missouri", "Nebraska")
)

## ══════════════════════════════════════════════════════════════════════════════
## PART 1 — PIE-CHART MAP  (p)
## ══════════════════════════════════════════════════════════════════════════════

yield_wide <- ssm %>%
  filter(Crop %in% c("SOYBEAN-MG4","SOYBEAN-LT1.5-MG4",
                     "SOYBEAN-LT1.5-MG3","SOYBEAN-MG3")) %>%
  mutate(Genotype = factor(
    ifelse(Crop %in% c("SOYBEAN-MG4","SOYBEAN-MG3"),
           "Non-limited transpiration","Limited transpiration"),
    levels = c("Non-limited transpiration","Limited transpiration"))) %>%
  group_by(Location, Genotype, water_mgmt) %>%
  summarise(Yield = mean(Yield_kgha, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = water_mgmt, values_from = Yield) %>%
  left_join(sites, by = "Location") %>%
  mutate(
    lon_pie    = ifelse(Genotype == "Limited transpiration",
                        lon - 1.35, lon + 1.35),
    r          = sqrt(Irrigated) / sqrt(max(Irrigated)) * 1.5,
    frac_rain  = Rainfed / Irrigated,
    start_rain = -pi/2,
    end_rain   = -pi/2 + frac_rain * 2 * pi,
    start_irr  = end_rain,
    end_irr    = -pi/2 + 2 * pi
  )

arcs <- bind_rows(
  yield_wide %>% mutate(slice = "Water-limited (Rainfed)",
                        start = start_rain, end = end_rain),
  yield_wide %>% mutate(slice = "Potential gain (Irrigated)",
                        start = start_irr,  end = end_irr)
) %>%
  mutate(fill_key = paste(as.character(Genotype), slice, sep = "|"))

fill_colours <- c(
  "Non-limited transpiration|Water-limited (Rainfed)"    = "#a1d99b",
  "Non-limited transpiration|Potential gain (Irrigated)" = "#238b45",
  "Limited transpiration|Water-limited (Rainfed)"        = "#fdbe85",
  "Limited transpiration|Potential gain (Irrigated)"     = "#d95f02"
)

pie_labels <- yield_wide %>%
  mutate(
    label_rain = paste0(round(Rainfed  / 1000, 1), "t"),
    label_irr  = paste0(round(Irrigated/ 1000, 1), "t")
  )

p <- ggplot() +
  geom_polygon(data = us_states,
               aes(x = long, y = lat, group = group),
               fill = "grey97", colour = "grey75", linewidth = 0.3) +
  geom_polygon(data = state_bg,
               aes(x = long, y = lat, group = group),
               fill = "#dbeaf5", colour = "grey65", linewidth = 0.4) +
  geom_point(data = yield_wide,
             aes(x = lon_pie, y = lat, size = Irrigated, colour = Genotype),
             shape = 1, stroke = 1.0, alpha = 0.5, show.legend = FALSE) +
  geom_arc_bar(data = arcs,
               aes(x0 = lon_pie, y0 = lat, r0 = 0, r = r,
                   start = start, end = end, fill = fill_key),
               colour = "grey25", linewidth = 0.2) +
  geom_text(data = pie_labels,
            aes(x = lon_pie + r * 0.38 * cos(-pi/2 + frac_rain * pi),
                y = lat     + r * 0.38 * sin(-pi/2 + frac_rain * pi),
                label = label_rain),
            size = 2.0, colour = "grey10", fontface = "bold") +
  geom_text(data = pie_labels,
            aes(x = lon_pie, y = lat + r + 0.42, label = label_irr),
            size = 2.2, colour = "grey10") +
  geom_text(data = yield_wide,
            aes(x = lon_pie, y = lat - r - 0.38,
                label = ifelse(Genotype == "Limited transpiration",
                               "LT", "Non-LT"),
                colour = Genotype),
            size = 2.4, fontface = "bold", show.legend = FALSE) +
  geom_text(data = state_labels,
            aes(x = lon, y = lat, label = label),
            fontface = "bold.italic", size = 3.0, colour = "grey30") +
  scale_fill_manual(
    values = fill_colours,
    labels = c(
      "Non-limited transpiration|Water-limited (Rainfed)"    = "Non-LT  ·  Rainfed (water-limited)",
      "Non-limited transpiration|Potential gain (Irrigated)" = "Non-LT  ·  Irrigated gain (potential)",
      "Limited transpiration|Water-limited (Rainfed)"        = "LT  ·  Rainfed (water-limited)",
      "Limited transpiration|Potential gain (Irrigated)"     = "LT  ·  Irrigated gain (potential)"
    ),
    breaks = c(
      "Non-limited transpiration|Water-limited (Rainfed)",
      "Non-limited transpiration|Potential gain (Irrigated)",
      "Limited transpiration|Water-limited (Rainfed)",
      "Limited transpiration|Potential gain (Irrigated)"
    ),
    name = "Yield component"
  ) +
  scale_colour_manual(
    values = c("Non-limited transpiration" = "#238b45",
               "Limited transpiration"      = "#d95f02"),
    guide  = "none"
  ) +
  scale_size_continuous(
    range  = c(2, 9),
    name   = expression(Potential~yield~(kg~ha^{-1})),
    breaks = c(3000, 5000, 7000),
    labels = c("3 000","5 000","7 000")
  ) +
  coord_fixed(1.3, xlim = c(-106, -84), ylim = c(32.5, 45.5)) +
  labs(
    x        = NULL, y = NULL,
    title    = "Water-limited vs. potential soybean yield across six locations",
    subtitle = paste0("Pie area ∝ irrigated (potential) yield  ·  ",
                      "Rainfed yield labelled inside  ·  ",
                      "Irrigated yield above each pie (t ha⁻¹)")
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    panel.grid       = element_blank(),
    legend.position  = "bottom",
    legend.box       = "horizontal",
    legend.title     = element_text(face = "bold", size = 9),
    legend.text      = element_text(size = 8),
    legend.key.size  = unit(0.75, "lines"),
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(colour = "grey40", size = 8.5),
    panel.border     = element_rect(colour = "grey70")
  ) +
  guides(
    fill = guide_legend(order = 1, nrow = 2,
                        override.aes = list(colour = NA)),
    size = guide_legend(order = 2,
                        override.aes = list(shape = 21,
                                            colour = "grey50",
                                            fill   = NA))
  )

ggsave("figures/Fig_map_yields_pie.png", plot = p,
       width = 22, height = 17, units = "cm", dpi = 150, bg = "white")
message("Saved: figures/Fig_map_yields_pie.png")

## ══════════════════════════════════════════════════════════════════════════════
## PART 2 — FACET MAP  (p2)  — PRIMARY OUTPUT
## Rows = water management (Rainfed / Irrigated)
## Cols = Genotype (Non-LT / LT)
## ══════════════════════════════════════════════════════════════════════════════

## Colour palette matching spec
geno_colours <- c("Non-LT" = "#31a354", "LT" = "#d95f02")

## Scale point size: normalise within the full data range so circles are
## comparable across all four panels
yield_long2 <- yield_long   # already has both factor levels as labels

## Mean yield per panel per state — annotated as text on the map
state_means <- yield_long2 %>%
  group_by(Genotype, water_mgmt, State) %>%
  summarise(mean_yield = mean(Yield, na.rm = TRUE), .groups = "drop") %>%
  left_join(
    tibble(
      State = c("AR","MO","NE"),
      s_lon = c(-91.5, -92.8, -98.5),
      s_lat = c(33.4,  36.0,  42.8)
    ),
    by = "State"
  ) %>%
  mutate(label = paste0(round(mean_yield / 1000, 1), " t/ha"))

## Nudge duplicate AR/MO/NE site points (Jonesboro & Marianna are very close)
## — use the original lon/lat, no horizontal shift needed for facet approach
p2 <- ggplot() +

  ## ── Basemap: all states (grey background) ─────────────────────────────────
  geom_polygon(
    data = us_states,
    aes(x = long, y = lat, group = group),
    fill = "grey96", colour = "grey80", linewidth = 0.25
  ) +

  ## ── Focus states: AR / MO / NE highlighted ────────────────────────────────
  geom_polygon(
    data = state_bg,
    aes(x = long, y = lat, group = group),
    fill = "#eef4fb", colour = "grey55", linewidth = 0.55
  ) +

  ## ── State name labels ─────────────────────────────────────────────────────
  geom_text(
    data = state_labels,
    aes(x = lon, y = lat, label = label),
    fontface = "bold.italic", size = 2.8, colour = "grey35",
    inherit.aes = FALSE
  ) +

  ## ── Site circles: filled, sized by yield, coloured by genotype ───────────
  geom_point(
    data = yield_long2,
    aes(x = lon, y = lat, size = Yield, fill = Genotype),
    shape = 21, colour = "white", stroke = 0.5, alpha = 0.90
  ) +

  ## ── Yield value label at each site ────────────────────────────────────────
  geom_label(
    data = yield_long2,
    aes(x = lon, y = lat,
        label = paste0(round(Yield / 1000, 1), "t"),
        colour = Genotype),
    size = 2.4, fontface = "bold",
    label.size = 0.25, fill = "white", alpha = 0.80,
    nudge_y = -0.90, label.padding = unit(0.12, "lines")
  ) +

  ## ── Location name labels (repelled to avoid overlap) ─────────────────────
  geom_text_repel(
    data = yield_long2 %>% distinct(Location, lon, lat),
    aes(x = lon, y = lat, label = Location),
    size = 2.2, colour = "grey20",
    nudge_x = 0.5, nudge_y = 0.65,
    segment.size = 0.3, segment.colour = "grey60",
    min.segment.length = 0.2,
    inherit.aes = FALSE,
    max.overlaps = 20
  ) +

  ## ── Mean-yield-per-state annotation ───────────────────────────────────────
  geom_label(
    data = state_means,
    aes(x = s_lon, y = s_lat, label = paste0("State mean\n", label)),
    size = 2.2, colour = "grey15", fill = "white",
    label.size = 0.3, alpha = 0.85,
    fontface = "italic", inherit.aes = FALSE
  ) +

  ## ── Facet: rows = water_mgmt, cols = Genotype ─────────────────────────────
  facet_grid(water_mgmt ~ Genotype) +

  ## ── Scales ────────────────────────────────────────────────────────────────
  scale_fill_manual(values = geno_colours, name = "Genotype") +
  scale_colour_manual(values = geno_colours, guide = "none") +
  scale_size_area(
    max_size = 14,
    name     = expression(Yield~(kg~ha^{-1})),
    breaks   = c(2000, 4000, 6000),
    labels   = c("2 000", "4 000", "6 000")
  ) +

  coord_fixed(1.3, xlim = c(-106, -84), ylim = c(32.0, 45.8)) +

  labs(
    x        = NULL, y = NULL,
    title    = "Soybean yield across six US locations: Non-LT vs LT genotypes",
    subtitle = paste0("Circle area ∝ mean simulated yield  ·  ",
                      "Values in t ha⁻¹  ·  ",
                      "AR = Arkansas, MO = Missouri, NE = Nebraska")
  ) +

  theme_bw(base_size = 11) +
  theme(
    axis.text         = element_blank(),
    axis.ticks        = element_blank(),
    panel.grid        = element_blank(),
    strip.background  = element_rect(fill = "grey15", colour = "grey15"),
    strip.text        = element_text(colour = "white", face = "bold", size = 10),
    legend.position   = "bottom",
    legend.box        = "horizontal",
    legend.title      = element_text(face = "bold", size = 9),
    legend.text       = element_text(size = 8),
    legend.key.size   = unit(0.8, "lines"),
    plot.title        = element_text(face = "bold", size = 12),
    plot.subtitle     = element_text(colour = "grey40", size = 8.5),
    panel.border      = element_rect(colour = "grey60"),
    panel.spacing     = unit(0.6, "lines")
  ) +
  guides(
    fill = guide_legend(order = 1,
                        override.aes = list(size = 5, shape = 21,
                                            colour = "white", stroke = 0.5)),
    size = guide_legend(order = 2,
                        override.aes = list(shape = 21,
                                            fill   = "grey60",
                                            colour = "white"))
  )

ggsave("figures/Fig_map_yields_facet.png", plot = p2,
       width = 24, height = 18, units = "cm", dpi = 150, bg = "white")
message("Saved: figures/Fig_map_yields_facet.png")
