## Fig_map_study_sites.R
## Materials & Methods map: county soybean harvested area in the study region
## (AR, MO, NE, plus KS and OK for context) with the 10 SSM simulation sites.
##
## Acreage input (optional): a USDA NASS Quick Stats CSV export saved as
##   modeling_data/data/nass_soybean_acres_harvested_county.csv
## Quick Stats query (https://quickstats.nass.usda.gov):
##   Program = SURVEY, Sector = CROPS, Group = FIELD CROPS, Commodity = SOYBEANS,
##   Data Item = "SOYBEANS - ACRES HARVESTED", Domain = TOTAL,
##   Geographic Level = COUNTY, State = ARKANSAS, KANSAS, MISSOURI, NEBRASKA,
##   OKLAHOMA, Year = e.g. 2020-2024  ->  "Get Data" -> "Spreadsheet" (CSV)
## The mean across the exported years is mapped. Without the file the map is
## drawn with county outlines only, so the layout can be checked.
##
## Output: modeling_data/figures/Fig_map_study_sites.png / .tiff (600 dpi)
## Run from the repo root:  Rscript modeling_data/scripts/Fig_map_study_sites.R
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(maps); library(readr); library(ggrepel) })

NASS_FILE <- "modeling_data/data/nass_soybean_acres_harvested_county.csv"
STATES    <- c("arkansas", "kansas", "missouri", "nebraska", "oklahoma")
STUDY     <- c("arkansas", "missouri", "nebraska")

sites <- tribble(
  ~Location,      ~lat,  ~lon,
  "Jonesboro",    35.84, -90.70,
  "Keiser",       35.67, -90.09,
  "Marianna",     34.77, -90.76,
  "Rohwer",       33.80, -91.27,
  "Albany",       40.25, -94.33,
  "Mount Vernon", 37.10, -93.82,
  "Novelty",      40.01, -92.21,
  "Eustis",       40.66, -100.03,
  "Lincoln",      40.81, -96.70,
  "North Platte", 41.12, -100.77
)

## ── County polygons with FIPS codes ─────────────────────────────────────────────
cty <- map_data("county") |> filter(region %in% STATES) |>
  mutate(polyname = paste(region, subregion, sep = ","))
fips <- maps::county.fips |> mutate(polyname = sub(":.*", "", polyname)) |> distinct(polyname, .keep_all = TRUE)
cty <- cty |> left_join(fips, by = "polyname")
st  <- map_data("state") |> filter(region %in% STATES)

## ── Acreage (mean across exported years), if available ───────────────────────
have_acres <- file.exists(NASS_FILE)
if (have_acres) {
  nass <- read_csv(NASS_FILE, show_col_types = FALSE)
  acres <- nass |>
    filter(!is.na(`County ANSI`)) |>
    mutate(fips = as.integer(`State ANSI`) * 1000L + as.integer(`County ANSI`),
           Value = suppressWarnings(as.numeric(gsub(",", "", Value)))) |>
    filter(!is.na(Value)) |>
    group_by(fips) |> summarise(acres = mean(Value), n_years = n(), .groups = "drop")
  cat("Acreage: ", nrow(acres), "counties, years:", paste(range(nass$Year), collapse = "-"), "\n")
  cty <- cty |> left_join(acres, by = "fips")
} else {
  cat("No NASS file at", NASS_FILE, "- drawing county outlines only\n")
  cty$acres <- NA_real_
}

state_lab <- tibble(label = c("AR", "KS", "MO", "NE", "OK"),
                    lon = c(-92.4, -98.4, -92.6, -99.6, -97.4), lat = c(34.6, 38.4, 38.2, 41.6, 35.4))

p <- ggplot() +
  geom_polygon(data = cty, aes(long, lat, group = group, fill = acres),
               colour = "#bdbdbd", linewidth = 0.12) +
  geom_polygon(data = st, aes(long, lat, group = group), fill = NA, colour = "black",
               linewidth = ifelse(st$region %in% STUDY, 0.9, 0.45)) +
  geom_text(data = state_lab, aes(lon, lat, label = label), size = 4.2, colour = "#222222") +
  geom_point(data = sites, aes(lon, lat), shape = 21, size = 4.2, fill = "#ffe600",
             colour = "black", stroke = 0.8) +
  geom_text_repel(data = sites, aes(lon, lat, label = Location), size = 3.3,
                  box.padding = 0.45, point.padding = 0.35, min.segment.length = 0,
                  segment.colour = "#555555", segment.size = 0.3, seed = 3,
                  bg.colour = "white", bg.r = 0.12) +
  scale_fill_gradient(low = "#f7fcb9", high = "#006d2c", na.value = "white",
                      name = "Soybean harvested\narea (acres)",
                      labels = function(x) paste0(x / 1000, "k"),
                      limits = if (have_acres) c(0, NA) else c(0, 1)) +
  coord_quickmap(xlim = c(-104.3, -88.9), ylim = c(32.9, 43.2), expand = FALSE) +
  scale_x_continuous(breaks = seq(-104, -90, 2)) +
  scale_y_continuous(breaks = seq(34, 42, 2)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 0.7),
        axis.text = element_text(colour = "black"),
        legend.position = if (have_acres) c(0.995, 0.995) else "none",
        legend.justification = c(1, 1), legend.key.height = unit(0.32, "cm"),
        legend.margin = margin(3, 4, 3, 4),
        legend.background = element_rect(fill = "white", colour = NA),
        legend.title = element_text(size = 10), legend.text = element_text(size = 9),
        plot.background = element_rect(fill = "white", colour = NA))

ggsave("modeling_data/figures/Fig_map_study_sites.png", p, width = 7.5, height = 6.2, dpi = 300, bg = "white")
ggsave("modeling_data/figures/Fig_map_study_sites.tiff", p, width = 7.5, height = 6.2, dpi = 600, bg = "white",
       device = "tiff", compression = "lzw")
cat("Saved modeling_data/figures/Fig_map_study_sites.png / .tiff",
    if (have_acres) "(with acreage)\n" else "(outlines only; acreage file missing)\n")
