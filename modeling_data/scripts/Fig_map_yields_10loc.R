## Fig_map_yields_10loc.R
## Water-limited yield vs potential yield at the 10 updated SSM locations.
##
## Panel A: one pie per site (green = rainfed / water-limited yield as a share
##          of the irrigated potential, yellow = yield gap by water deficit),
##          one map for the Non-LT genotype and one for the LT1.5 genotype.
##          Pies are drawn next to the sites (leader line to the true position)
##          so the four Arkansas sites and the two western Nebraska sites do
##          not overlap; the water-limited share is printed under each pie.
## Panel B: mean attainable yield + yield gap across all locations (SE among sites).
##
## Input : modeling_data/data/updated-outputs/all_yearly.csv
## Output: modeling_data/figures/Fig_map_yields_10loc.png  (A + B)
##         modeling_data/figures/Fig_map_yields_10loc_A.png, _B.png
##         modeling_data/data/yield_gap_10loc_summary.csv
##
## Run from the repo root:  Rscript modeling_data/scripts/Fig_map_yields_10loc.R
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(ggforce)
  library(maps); library(readr); library(patchwork)
})

GREEN  <- "#31a354"
YELLOW <- "#ffffb2"

## true site position and the position where its pie is drawn
sites <- tribble(
  ~Location,      ~lat,  ~lon,     ~pie_lat, ~pie_lon,
  ## pies sit on an evenly spaced ring around the three states
  "North Platte", 41.12, -100.77,   43.9,    -101.8,
  "Lincoln",      40.81, -96.70,    43.9,    -97.6,
  "Albany",       40.25, -94.33,    43.9,    -93.4,
  "Novelty",      40.01, -92.21,    43.9,    -89.2,
  "Jonesboro",    35.84, -90.70,    39.6,    -87.0,
  "Keiser",       35.67, -90.09,    35.5,    -87.0,
  "Marianna",     34.77, -90.76,    31.2,    -89.6,
  "Rohwer",       33.80, -91.27,    31.2,    -94.2,
  "Mount Vernon", 37.10, -93.82,    35.5,    -99.8,
  "Eustis",       40.66, -100.03,   38.6,    -102.8
)

## ── Data ──────────────────────────────────────────────────────────────────────
yr <- read_csv("modeling_data/data/updated-outputs/all_yearly.csv",
               show_col_types = FALSE) |>
  mutate(Location   = trimws(gsub("­", "", Location)),
         Yield_kgha = as.numeric(Yield_kgha)) |>
  filter(Genotype_class %in% c("Non-LT", "LT1.5")) |>
  mutate(Genotype = factor(Genotype_class,
                           levels = c("Non-LT", "LT1.5"),
                           labels = c("Non-limited transpiration (Non-LT)",
                                      "Limited transpiration (LT)")))

site <- yr |>
  group_by(Location, Genotype, water_mgmt) |>
  summarise(Yield = mean(Yield_kgha, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = water_mgmt, values_from = Yield) |>
  mutate(WL_share = Rainfed / Irrigated, Gap = Irrigated - Rainfed) |>
  left_join(sites, by = "Location")
write_csv(site |> select(-pie_lat, -pie_lon), "modeling_data/data/yield_gap_10loc_summary.csv")

summ <- site |>
  group_by(Genotype) |>
  summarise(Rainfed_m = mean(Rainfed), Irrigated_m = mean(Irrigated),
            Rainfed_se = sd(Rainfed) / sqrt(n()),
            Irrigated_se = sd(Irrigated) / sqrt(n()), .groups = "drop") |>
  mutate(Gap_pct = 100 * (Irrigated_m - Rainfed_m) / Irrigated_m)
print(summ)
nl <- filter(summ, grepl("Non-LT", Genotype)); lt <- filter(summ, !grepl("Non-LT", Genotype))
cat(sprintf("LT water-limited yield vs Non-LT: %+.1f%%\n", 100 * (lt$Rainfed_m / nl$Rainfed_m - 1)))
cat(sprintf("Non-LT potential yield vs LT:     %+.1f%%\n", 100 * (nl$Irrigated_m / lt$Irrigated_m - 1)))

## ── Base map ──────────────────────────────────────────────────────────────────
us_states <- map_data("state")
focus     <- filter(us_states, region %in% c("arkansas", "missouri", "nebraska"))
counties  <- filter(map_data("county"), region %in% c("arkansas", "missouri", "nebraska"))
state_lab <- tibble(label = c("Arkansas", "Missouri", "Nebraska"),
                    lon = c(-92.9, -92.4, -98.4), lat = c(35.0, 38.4, 41.7))

## pie geometry: constant radius (potential yield is shown in panel B),
## green slice = water-limited share
pies <- site |>
  mutate(r = 1.05,
         end_green = 2 * pi * WL_share,
         pct = sprintf("%.0f%%", 100 * WL_share))
arcs <- bind_rows(
  pies |> transmute(Genotype, pie_lon, pie_lat, r, start = 0, end = end_green,
                    slice = "Water-limited yield"),
  pies |> transmute(Genotype, pie_lon, pie_lat, r, start = end_green, end = 2 * pi,
                    slice = "Yield gap by water deficit")
) |> mutate(slice = factor(slice, levels = c("Water-limited yield", "Yield gap by water deficit")))

pA <- ggplot() +
  geom_polygon(data = us_states, aes(long, lat, group = group),
               fill = "#ececec", colour = "#c4c4c4", linewidth = 0.25) +
  geom_polygon(data = focus, aes(long, lat, group = group),
               fill = "#c9c9c9", colour = NA) +
  ## county limits in the background of the three states
  geom_polygon(data = counties, aes(long, lat, group = group),
               fill = NA, colour = "#9a9a9a", linewidth = 0.15) +
  geom_polygon(data = focus, aes(long, lat, group = group),
               fill = NA, colour = "black", linewidth = 0.5) +
  geom_text(data = state_lab, aes(lon, lat, label = label),
            colour = "#3c3c3c", size = 3.6, fontface = "italic") +
  ## leader line from the site to its pie, and the site itself
  geom_segment(data = pies, aes(x = lon, y = lat, xend = pie_lon, yend = pie_lat),
               colour = "#444444", linewidth = 0.4) +
  geom_point(data = pies, aes(lon, lat), size = 1.8, colour = "black") +
  geom_arc_bar(data = arcs,
               aes(x0 = pie_lon, y0 = pie_lat, r0 = 0, r = r, start = start, end = end, fill = slice),
               colour = "black", linewidth = 0.3) +
  ## site name above and water-limited share below each pie
  geom_text(data = pies, aes(pie_lon, pie_lat + r + 0.5, label = Location), size = 3.3) +
  geom_text(data = pies, aes(pie_lon, pie_lat - r - 0.5, label = pct), size = 3.3, fontface = "bold") +
  scale_fill_manual(values = c("Water-limited yield" = GREEN,
                               "Yield gap by water deficit" = YELLOW), name = NULL) +
  facet_wrap(~Genotype, ncol = 2) +
  coord_fixed(1, xlim = c(-104.6, -85.0), ylim = c(29.4, 45.8), expand = FALSE) +
  theme_bw(base_size = 12) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(), panel.grid = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", linewidth = 0.8),
        strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.8),
        strip.text = element_text(size = 12),
        legend.position = "top", legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 11),
        plot.background = element_rect(fill = "white", colour = NA))

## ── Panel B ───────────────────────────────────────────────────────────────────
bars <- summ |>
  mutate(G = factor(ifelse(grepl("Non-LT", Genotype), "Non-LT", "LT"), levels = c("Non-LT", "LT")))
bars_long <- bind_rows(
  bars |> transmute(G, slice = "Water-limited yield", value = Rainfed_m),
  bars |> transmute(G, slice = "Yield gap by water deficit", value = Irrigated_m - Rainfed_m)
) |> mutate(slice = factor(slice, levels = c("Yield gap by water deficit", "Water-limited yield")))

pB <- ggplot() +
  geom_col(data = bars_long, aes(G, value, fill = slice), width = 0.6,
           colour = "black", linewidth = 0.4) +
  geom_errorbar(data = bars, aes(G, ymin = Rainfed_m - Rainfed_se, ymax = Rainfed_m + Rainfed_se),
                width = 0, linewidth = 0.6) +
  geom_errorbar(data = bars, aes(G, ymin = Irrigated_m - Irrigated_se, ymax = Irrigated_m + Irrigated_se),
                width = 0, linewidth = 0.6) +
  geom_text(data = bars, aes(G, Irrigated_m + 450, label = sprintf("%.0f%% gap", Gap_pct)), size = 3.4) +
  scale_fill_manual(values = c("Water-limited yield" = GREEN,
                               "Yield gap by water deficit" = YELLOW), guide = "none") +
  scale_y_continuous(labels = scales::comma, limits = c(0, 8000),
                     breaks = c(0, 2000, 4000, 6000, 8000), expand = c(0, 0)) +
  facet_wrap(~"All locations") +
  labs(x = NULL, y = expression(Yield~(kg~ha^{-1}))) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 0.8),
        strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.8),
        strip.text = element_text(size = 12), axis.text = element_text(colour = "black"),
        plot.background = element_rect(fill = "white", colour = NA))

ggsave("modeling_data/figures/Fig_map_yields_10loc_A.png", pA, width = 12.5, height = 6.2, dpi = 300, bg = "white")
ggsave("modeling_data/figures/Fig_map_yields_10loc_B.png", pB, width = 4.2, height = 6.2, dpi = 300, bg = "white")
p <- (pA + pB) + plot_layout(widths = c(2.4, 0.55)) + plot_annotation(tag_levels = "A")
ggsave("modeling_data/figures/Fig_map_yields_10loc.png", p, width = 17, height = 6.2, dpi = 300, bg = "white")
cat("Saved Fig_map_yields_10loc*.png\n")
