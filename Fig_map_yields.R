## Fig 5 — Map: water-limited vs potential yield at 6 locations
## ggforce::geom_arc_bar pies + geom_point scaled by yield

setwd("/home/user/LT-trait")
library(dplyr); library(tidyr); library(ggplot2); library(ggforce); library(maps); library(readxl)

## ── Site coordinates ──────────────────────────────────────────────────────────
sites <- tibble(
  Location = c("Jonesboro","Marianna","Columbia","Norborne","Lincoln","North Platte"),
  State    = c("AR","AR","MO","MO","NE","NE"),
  lat      = c(35.84, 34.78, 38.95, 39.30, 40.81, 41.12),
  lon      = c(-90.70,-90.76,-92.33,-93.67,-96.70,-100.77)
)

## ── Yield data ────────────────────────────────────────────────────────────────
ssm <- read_excel("SoybeanSSM_sims2.xlsx")

yield_data <- ssm %>%
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
    r          = sqrt(Irrigated) / sqrt(max(Irrigated)) * 1.5,   # area ∝ yield
    frac_rain  = Rainfed / Irrigated,
    start_rain = -pi/2,
    end_rain   = -pi/2 + frac_rain * 2 * pi,
    start_irr  = end_rain,
    end_irr    = -pi/2 + 2 * pi
  )

## ── Arc data (two slices per pie) ────────────────────────────────────────────
arcs <- bind_rows(
  yield_data %>% mutate(slice = "Water-limited (Rainfed)",
                        start = start_rain, end = end_rain),
  yield_data %>% mutate(slice = "Potential gain (Irrigated)",
                        start = start_irr,  end = end_irr)
) %>%
  mutate(fill_key = paste(as.character(Genotype), slice, sep = "|"))

fill_colours <- c(
  "Non-limited transpiration|Water-limited (Rainfed)"    = "#a1d99b",
  "Non-limited transpiration|Potential gain (Irrigated)" = "#238b45",
  "Limited transpiration|Water-limited (Rainfed)"        = "#fdbe85",
  "Limited transpiration|Potential gain (Irrigated)"     = "#d95f02"
)

## ── Yield labels at pie centre ────────────────────────────────────────────────
pie_labels <- yield_data %>%
  mutate(
    label_rain = paste0(round(Rainfed / 1000, 1), "t"),
    label_irr  = paste0(round(Irrigated / 1000, 1), "t")
  )

## ── State label positions ─────────────────────────────────────────────────────
state_labels <- tibble(
  lon   = c(-91.5, -93.5, -99.0),
  lat   = c(33.6,   36.7,  43.0),
  label = c("Arkansas", "Missouri", "Nebraska")
)

## ── Genotype legend key positions ─────────────────────────────────────────────
geno_legend <- tibble(
  lon_pie = c(-104.5, -104.5),
  lat     = c(34.2,   35.5),
  r       = 0.55,
  Genotype = factor(c("Non-limited transpiration","Limited transpiration"),
                    levels = c("Non-limited transpiration","Limited transpiration")),
  label   = c("Non-LT", "LT")
)
geno_arcs <- bind_rows(
  geno_legend %>% mutate(slice="full", start=-pi/2, end=3*pi/2,
                         fill_key=paste(as.character(Genotype),"Water-limited (Rainfed)",sep="|"))
)

## ── Map ───────────────────────────────────────────────────────────────────────
us_states    <- map_data("state")
focus_states <- c("arkansas","missouri","nebraska")

p <- ggplot() +
  ## basemap
  geom_polygon(data = us_states,
               aes(x = long, y = lat, group = group),
               fill = "grey97", colour = "grey75", linewidth = 0.3) +
  geom_polygon(data = filter(us_states, region %in% focus_states),
               aes(x = long, y = lat, group = group),
               fill = "#dbeaf5", colour = "grey65", linewidth = 0.4) +

  ## geom_point — circle scaled by potential (irrigated) yield as background anchor
  geom_point(data = yield_data,
             aes(x = lon_pie, y = lat, size = Irrigated, colour = Genotype),
             shape = 1, stroke = 1.0, alpha = 0.5, show.legend = FALSE) +

  ## pie charts
  geom_arc_bar(data = arcs,
               aes(x0 = lon_pie, y0 = lat, r0 = 0, r = r,
                   start = start, end = end, fill = fill_key),
               colour = "grey25", linewidth = 0.2) +

  ## rainfed yield label (inner slice)
  geom_text(data = pie_labels,
            aes(x = lon_pie + r * 0.38 * cos(-pi/2 + frac_rain * pi),
                y = lat     + r * 0.38 * sin(-pi/2 + frac_rain * pi),
                label = label_rain),
            size = 2.0, colour = "grey10", fontface = "bold") +

  ## potential yield label (outer ring, top)
  geom_text(data = pie_labels,
            aes(x = lon_pie, y = lat + r + 0.42, label = label_irr),
            size = 2.2, colour = "grey10") +

  ## genotype label below each pie
  geom_text(data = yield_data,
            aes(x = lon_pie, y = lat - r - 0.38,
                label = ifelse(Genotype == "Limited transpiration",
                               "LT", "Non-LT"),
                colour = Genotype),
            size = 2.4, fontface = "bold", show.legend = FALSE) +

  ## state labels
  geom_text(data = state_labels,
            aes(x = lon, y = lat, label = label),
            fontface = "bold.italic", size = 3.0, colour = "grey30") +

  ## scales
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

ggsave("Fig_map_yields.png",  plot = p, width = 22, height = 17,
       units = "cm", dpi = 150, bg = "white")
ggsave("Fig_map_yields.tiff", plot = p, width = 22, height = 17,
       units = "cm", dpi = 300, compression = "lzw", bg = "white")
message("Saved: Fig_map_yields")
