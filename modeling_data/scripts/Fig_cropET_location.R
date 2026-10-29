## Fig_cropET_location.R
## Cumulative crop ET under rainfed conditions, per simulation location.
##
##   Fig_cropET_cumulative_location.png : cumulative crop ET (mm), LT1.5 vs
##       Non-LT, mean across seasons, one panel per location
##   Fig_cropET_difference_location.png : LT minus Non-LT cumulative crop ET
##       (mm), computed season by season (same site, year, sowing date and
##       maturity group, so the weather is identical), mean +/- 1 SD across
##       seasons; below zero = LT has used less water so far
##
## x axis = cumulative thermal time from sowing (SSM daily temperature units,
## DTU), as in Fig_soil_water_season_state.R. Each pair is cut at the earlier
## of the two R7 dates. Crop ET = soil evaporation + transpiration (CE + CTR).
##
## Input : modeling_data/data/updated-outputs/<Location>_output.xlsx
## Output: modeling_data/figures/Fig_cropET_cumulative_location.png / .tiff
##         modeling_data/figures/Fig_cropET_difference_location.png / .tiff
##         modeling_data/data/cropET_location_summary.csv
## Run from the repo root:  Rscript modeling_data/scripts/Fig_cropET_location.R
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(readxl); library(readr) })
invisible(Sys.setlocale("LC_CTYPE", "C.UTF-8"))

COLS <- c("Non-limited transpiration" = "#31a354", "Limited transpiration" = "#d95f02")
LOC_ORDER <- c("Jonesboro", "Keiser", "Marianna", "Rohwer",
               "Albany", "Mount Vernon", "Novelty",
               "Eustis", "Lincoln", "North Platte")
CROPS <- c("SOYBEAN-MG3", "SOYBEAN-MG4", "SOYBEAN-LT1.5-MG3", "SOYBEAN-LT1.5-MG4")
BIN <- 50   # deg C d

files <- list.files("modeling_data/data/updated-outputs", pattern = "_output\\.xlsx$", full.names = TRUE)
cat("Reading", length(files), "location files ...\n")
read_one <- function(f) {
  yearly <- read_excel(f, sheet = "YearlyO", guess_max = 5000) |>
    select(sName, Pyear, R7 = R7.PM) |> mutate(R7 = as.numeric(R7))
  read_excel(f, sheet = "DailyO", guess_max = 5000) |>
    select(sName, Location, Manag, Crop, Pyear, DAP, DTU, CE, CTR) |>
    filter(grepl("RFD", Manag), Crop %in% CROPS) |>
    left_join(yearly, by = c("sName", "Pyear")) |>
    mutate(across(c(DAP, DTU, CE, CTR), as.numeric))
}
daily <- bind_rows(lapply(files, read_one)) |>
  mutate(Location = trimws(gsub("­", "", Location)),
         Genotype = ifelse(grepl("LT1.5", Crop), "LT", "NonLT"),
         MG       = sub(".*-(MG[34])$", "\\1", Crop),
         pair     = paste(Location, Manag, MG, Pyear),
         cumET    = CE + CTR) |>
  arrange(pair, Genotype, DAP) |>
  group_by(pair, Genotype) |> mutate(GDD = cumsum(coalesce(DTU, 0))) |> ungroup()

## pair the two genotypes day by day (same weather -> same GDD at the same DAP)
paired <- daily |>
  select(Location, pair, Genotype, DAP, GDD, R7, cumET) |>
  pivot_wider(names_from = Genotype, values_from = c(R7, cumET, GDD)) |>
  filter(!is.na(cumET_LT), !is.na(cumET_NonLT), DAP >= 1, DAP <= pmin(R7_LT, R7_NonLT)) |>
  mutate(GDD = GDD_NonLT, diff = cumET_LT - cumET_NonLT,
         GDD_bin = (floor(GDD / BIN) + 0.5) * BIN)
cat("Paired seasons:", n_distinct(paired$pair), "\n")

## summaries per location and GDD bin (season mean within bin, then across seasons)
by_season <- paired |>
  group_by(Location, pair, GDD_bin) |>
  summarise(LT = mean(cumET_LT), NonLT = mean(cumET_NonLT), diff = mean(diff), .groups = "drop")
summ <- by_season |>
  group_by(Location) |> mutate(n_pairs = n_distinct(pair)) |> ungroup() |>
  group_by(Location, GDD_bin, n_pairs) |>
  summarise(n = n(), LT = mean(LT), NonLT = mean(NonLT),
            diff_sd = sd(diff), diff = mean(diff), .groups = "drop") |>
  filter(n >= 0.5 * n_pairs) |>
  mutate(Location = factor(Location, levels = LOC_ORDER))
write_csv(summ, "modeling_data/data/cropET_location_summary.csv")

thm <- theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#ececec"),
        panel.border = element_rect(colour = "black", linewidth = 0.7),
        strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.7),
        strip.text = element_text(size = 12, colour = "black"),
        axis.text = element_text(colour = "black"), legend.position = "top",
        legend.title = element_blank(), legend.text = element_text(size = 12),
        plot.background = element_rect(fill = "white", colour = NA))
xlab <- "Cumulative thermal time (°C d)"

## 1) cumulative ET, both genotypes
long <- summ |> select(Location, GDD_bin, LT, NonLT) |>
  pivot_longer(c(LT, NonLT), names_to = "g", values_to = "cumET") |>
  mutate(Genotype = factor(ifelse(g == "LT", "Limited transpiration", "Non-limited transpiration"), levels = names(COLS)))
p1 <- ggplot(long, aes(GDD_bin, cumET, colour = Genotype)) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = COLS) +
  facet_wrap(~Location, ncol = 5) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = xlab, y = "Cumulative crop ET (mm)") + thm
ggsave("modeling_data/figures/Fig_cropET_cumulative_location.png", p1, width = 14, height = 6.2, dpi = 300, bg = "white")
ggsave("modeling_data/figures/Fig_cropET_cumulative_location.tiff", p1, width = 14, height = 6.2, dpi = 600, bg = "white",
       device = "tiff", compression = "lzw")

## 2) difference LT - Non-LT
p2 <- ggplot(summ, aes(GDD_bin, diff)) +
  geom_hline(yintercept = 0, colour = "#555555", linewidth = 0.5) +
  geom_ribbon(aes(ymin = diff - diff_sd, ymax = diff + diff_sd), fill = "#d95f02", alpha = 0.15) +
  geom_line(colour = "#d95f02", linewidth = 1.2) +
  facet_wrap(~Location, ncol = 5) +
  labs(x = xlab, y = "Cumulative crop ET, LT minus Non-LT (mm)") + thm
ggsave("modeling_data/figures/Fig_cropET_difference_location.png", p2, width = 14, height = 6.2, dpi = 300, bg = "white")
ggsave("modeling_data/figures/Fig_cropET_difference_location.tiff", p2, width = 14, height = 6.2, dpi = 600, bg = "white",
       device = "tiff", compression = "lzw")
cat("Saved Fig_cropET_cumulative_location and Fig_cropET_difference_location (.png, .tiff)\n")

## headline numbers: largest deficit of LT and the value at the end
cat("\nLT minus Non-LT cumulative crop ET (mm), mean across paired seasons:\n")
summ |> group_by(Location) |>
  summarise(max_saving = round(min(diff), 1), at_GDD = GDD_bin[which.min(diff)],
            at_end = round(last(diff), 1), .groups = "drop") |> print()
