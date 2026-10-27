## Fig_soil_water_season.R
## Seasonal course of root-zone soil water under rainfed conditions, limited
## transpiration (LT1.5) vs non-limited transpiration (Non-LT), one panel per
## simulation location. Line = mean across the 30 weather-years x 3 sowing
## dates, shaded band = mean +/- 1 SD across those seasons. Each season is cut
## at physiological maturity (R7.PM from the yearly output).
##
## Input : modeling_data/data/updated-outputs/<Location>_output.xlsx (sheet DailyO)
## Output: modeling_data/figures/Fig_soil_water_season_ATSW.png   (mm)
##         modeling_data/figures/Fig_soil_water_season_FTSW.png   (fraction)
##         modeling_data/data/soil_water_season_summary.csv
## Run from the repo root:  Rscript modeling_data/scripts/Fig_soil_water_season.R
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(readxl); library(readr) })

COLS <- c("Non-limited transpiration" = "#31a354", "Limited transpiration" = "#d95f02")
LOC_ORDER <- c("Jonesboro", "Keiser", "Marianna", "Rohwer",
               "Albany", "Mount Vernon", "Novelty",
               "Eustis", "Lincoln", "North Platte")

files <- list.files("modeling_data/data/updated-outputs", pattern = "_output\\.xlsx$", full.names = TRUE)
cat("Reading", length(files), "location files ...\n")

daily <- bind_rows(lapply(files, function(f) {
  cat("  ", basename(f), "\n")
  maturity <- read_excel(f, sheet = "YearlyO", guess_max = 5000) |>
    select(sName, Pyear, R7.PM) |> mutate(R7.PM = as.numeric(R7.PM))
  read_excel(f, sheet = "DailyO", col_types = "guess", guess_max = 5000) |>
    select(sName, Location, Manag, Crop, Pyear, DAP, ATSWRZ, FTSWRZ) |>
    left_join(maturity, by = c("sName", "Pyear")) |>
    filter(!is.na(R7.PM), DAP <= R7.PM) |>
    filter(grepl("RFD", Manag), Crop %in% c("SOYBEAN-MG3", "SOYBEAN-MG4",
                                            "SOYBEAN-LT1.5-MG3", "SOYBEAN-LT1.5-MG4"),
           !is.na(DAP), DAP >= 1) |>
    mutate(Location = trimws(gsub("­", "", Location)),
           Genotype = ifelse(grepl("LT1.5", Crop), "Limited transpiration", "Non-limited transpiration"),
           Sowing   = sub(".*-(ELY|MID|LAT).*", "\\1", Manag),
           season   = paste(Pyear, Sowing),
           ATSWRZ = as.numeric(ATSWRZ), FTSWRZ = as.numeric(FTSWRZ))
}))
cat("Rows:", nrow(daily), "| locations:", paste(sort(unique(daily$Location)), collapse = ", "), "\n")
cat("Seasons per location x genotype:\n"); print(daily |> distinct(Location, Genotype, season) |> count(Location, Genotype) |> head(4))

## summarise across seasons at each day after planting; keep days where at
## least half of the seasons are still growing so the tail is not driven by
## a few long seasons
summ <- daily |>
  group_by(Location, Genotype) |> mutate(n_seasons = n_distinct(season)) |> ungroup() |>
  group_by(Location, Genotype, DAP, n_seasons) |>
  summarise(n = n(), across(c(ATSWRZ, FTSWRZ), list(mean = \(v) mean(v, na.rm = TRUE), sd = \(v) sd(v, na.rm = TRUE))), .groups = "drop") |>
  filter(n >= 0.5 * n_seasons) |>
  mutate(Location = factor(Location, levels = LOC_ORDER),
         Genotype = factor(Genotype, levels = names(COLS)))
write_csv(summ, "modeling_data/data/soil_water_season_summary.csv")

thm <- theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#ececec"),
        panel.border = element_rect(colour = "black", linewidth = 0.7),
        strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.7),
        strip.text = element_text(size = 12, colour = "black"),
        axis.text = element_text(colour = "black"), legend.position = "top",
        legend.title = element_blank(), plot.background = element_rect(fill = "white", colour = NA))

make_plot <- function(var, ylab) {
  m <- summ[[paste0(var, "_mean")]]; s <- summ[[paste0(var, "_sd")]]
  d <- summ |> mutate(mean = m, lo = m - s, hi = m + s)
  if (var == "FTSWRZ") d <- d |> mutate(lo = pmax(lo, 0), hi = pmin(hi, 1))
  ggplot(d, aes(DAP, mean, colour = Genotype, fill = Genotype)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.22, colour = NA) +
    geom_line(linewidth = 1.0) +
    scale_colour_manual(values = COLS) + scale_fill_manual(values = COLS) +
    facet_wrap(~Location, ncol = 5) +
    labs(x = "Days after planting", y = ylab) + thm
}

pA <- make_plot("ATSWRZ", "Root-zone transpirable soil water (mm)")
ggsave("modeling_data/figures/Fig_soil_water_season_ATSW.png", pA, width = 14, height = 6.5, dpi = 300, bg = "white")
pF <- make_plot("FTSWRZ", "Fraction of transpirable soil water")
ggsave("modeling_data/figures/Fig_soil_water_season_FTSW.png", pF, width = 14, height = 6.5, dpi = 300, bg = "white")
cat("Saved Fig_soil_water_season_ATSW.png and Fig_soil_water_season_FTSW.png\n")

## headline numbers: mean difference LT - Non-LT over the season, per location
diff <- summ |> select(Location, Genotype, DAP, ATSWRZ_mean, FTSWRZ_mean) |>
  pivot_wider(names_from = Genotype, values_from = c(ATSWRZ_mean, FTSWRZ_mean)) |>
  group_by(Location) |>
  summarise(ATSW_diff_mm = mean(`ATSWRZ_mean_Limited transpiration` - `ATSWRZ_mean_Non-limited transpiration`, na.rm = TRUE),
            FTSW_diff = mean(`FTSWRZ_mean_Limited transpiration` - `FTSWRZ_mean_Non-limited transpiration`, na.rm = TRUE),
            ATSW_diff_mm_max = max(`ATSWRZ_mean_Limited transpiration` - `ATSWRZ_mean_Non-limited transpiration`, na.rm = TRUE)) |>
  mutate(across(where(is.numeric), ~ round(., 2)))
print(diff, n = Inf)
