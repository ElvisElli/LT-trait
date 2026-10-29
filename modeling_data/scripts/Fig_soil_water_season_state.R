## Fig_soil_water_season_state.R
## Seasonal course of root-zone soil water under rainfed conditions, limited
## transpiration (LT1.5) vs non-limited transpiration (Non-LT), one panel per
## state (Arkansas, Missouri, Nebraska), on a thermal-time axis.
##
##   x axis : cumulative growing degree days from planting, the running sum of
##            the SSM daily thermal units (DTU, deg C d)
##   lines  : mean across all seasons (30 weather-years x 3 sowing dates) of all
##            locations in the state; each season is averaged within 25 deg C d
##            bins first, then across seasons
##   band   : mean +/- 1 SD across those seasons (light)
##   arrows : mean thermal time of R1 (flowering), R5 (beginning seed filling)
##            and R7 (physiological maturity) in that state, from the yearly
##            output; each season is cut at R7
##
## Input : modeling_data/data/updated-outputs/<Location>_output.xlsx
##         (sheets DailyO and YearlyO)
## Output: modeling_data/figures/Fig_soil_water_season_state_FTSW.png/.tiff
##         modeling_data/figures/Fig_soil_water_season_state_ATSW.png/.tiff
##         modeling_data/data/soil_water_season_state_summary.csv
##         modeling_data/data/soil_water_season_state_stages.csv
## Run from the repo root:  Rscript modeling_data/scripts/Fig_soil_water_season_state.R
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(readxl); library(readr) })
invisible(Sys.setlocale("LC_CTYPE", "C.UTF-8"))   # degree sign in the axis label

COLS  <- c("Non-limited transpiration" = "#31a354", "Limited transpiration" = "#d95f02")
STATE <- c(Jonesboro = "Arkansas", Keiser = "Arkansas", Marianna = "Arkansas", Rohwer = "Arkansas",
           Albany = "Missouri", `Mount Vernon` = "Missouri", Novelty = "Missouri",
           Eustis = "Nebraska", Lincoln = "Nebraska", `North Platte` = "Nebraska")
CROPS <- c("SOYBEAN-MG3", "SOYBEAN-MG4", "SOYBEAN-LT1.5-MG3", "SOYBEAN-LT1.5-MG4")
BIN   <- 25    # deg C d

## ── Read daily and yearly output, rainfed, Non-LT and LT1.5 ───────────────────
files <- list.files("modeling_data/data/updated-outputs", pattern = "_output\\.xlsx$", full.names = TRUE)
cat("Reading", length(files), "location files ...\n")
read_one <- function(f) {
  yearly <- read_excel(f, sheet = "YearlyO", guess_max = 5000) |>
    select(sName, Pyear, R1 = R1.SEL.SIL, R5 = R5.BSG, R7 = R7.PM) |>
    mutate(across(c(R1, R5, R7), as.numeric))
  read_excel(f, sheet = "DailyO", guess_max = 5000) |>
    select(sName, Location, Manag, Crop, Pyear, DAP, DTU, ATSWRZ, FTSWRZ) |>
    filter(grepl("RFD", Manag), Crop %in% CROPS) |>
    left_join(yearly, by = c("sName", "Pyear")) |>
    mutate(across(c(DAP, DTU, ATSWRZ, FTSWRZ), as.numeric))
}
daily <- bind_rows(lapply(files, read_one)) |>
  mutate(Location = trimws(gsub("­", "", Location)),
         State    = unname(STATE[Location]),
         Genotype = ifelse(grepl("LT1.5", Crop), "Limited transpiration", "Non-limited transpiration"),
         season   = paste(Location, sName, Pyear)) |>
  arrange(season, DAP) |>
  group_by(season) |>
  mutate(GDD = cumsum(coalesce(DTU, 0))) |>          # cumulative thermal time from planting
  ungroup()
cat("Rows:", nrow(daily), "| seasons:", n_distinct(daily$season), "\n")

## ── Phenology on the thermal-time axis ─────────────────────────────────────────
stage_gdd <- daily |>
  group_by(State, Genotype, season) |>
  summarise(R1 = GDD[match(first(R1), DAP)], R5 = GDD[match(first(R5), DAP)],
            R7 = GDD[match(first(R7), DAP)], .groups = "drop")
stages <- stage_gdd |>
  pivot_longer(c(R1, R5, R7), names_to = "Stage", values_to = "GDD") |>
  group_by(State, Stage) |> summarise(GDD = mean(GDD, na.rm = TRUE), .groups = "drop")
cat("\nMean thermal time of each stage (deg C d), per state and genotype:\n")
print(stage_gdd |> group_by(State, Genotype) |>
        summarise(across(c(R1, R5, R7), \(v) round(mean(v, na.rm = TRUE))), .groups = "drop"))
write_csv(stages, "modeling_data/data/soil_water_season_state_stages.csv")

## ── Bin each season by thermal time (up to R7), then average across seasons ───
st <- daily |>
  filter(DAP >= 1, DAP <= R7) |>
  mutate(GDD_bin = (floor(GDD / BIN) + 0.5) * BIN) |>
  group_by(State, Genotype, season, GDD_bin) |>
  summarise(ATSWRZ = mean(ATSWRZ, na.rm = TRUE), FTSWRZ = mean(FTSWRZ, na.rm = TRUE), .groups = "drop") |>
  group_by(State, Genotype) |> mutate(n_seasons = n_distinct(season)) |> ungroup() |>
  group_by(State, Genotype, GDD_bin, n_seasons) |>
  summarise(n = n(),
            ATSWRZ_sd = sd(ATSWRZ), ATSWRZ = mean(ATSWRZ),
            FTSWRZ_sd = sd(FTSWRZ), FTSWRZ = mean(FTSWRZ), .groups = "drop") |>
  filter(n >= 0.5 * n_seasons) |>                     # bins reached by at least half the seasons
  mutate(State = factor(State, levels = c("Arkansas", "Missouri", "Nebraska")),
         Genotype = factor(Genotype, levels = names(COLS)))
write_csv(st, "modeling_data/data/soil_water_season_state_summary.csv")
stages <- stages |> mutate(State = factor(State, levels = levels(st$State)))

## ── Plot ───────────────────────────────────────────────────────────────────────
thm <- theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#ececec"),
        panel.border = element_rect(colour = "black", linewidth = 0.7),
        strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.7),
        strip.text = element_text(size = 14, colour = "black"),
        axis.text = element_text(colour = "black"), legend.position = "top",
        legend.title = element_blank(), legend.text = element_text(size = 13),
        plot.caption = element_text(size = 11, colour = "#444444", hjust = 0),
        plot.background = element_rect(fill = "white", colour = NA))

plot_state <- function(var, ylab, ymax) {
  d <- st |> mutate(m = .data[[var]], s = .data[[paste0(var, "_sd")]],
                    lo = pmax(m - s, 0), hi = pmin(m + s, ymax))
  ## arrows sit in a strip just under the data, pointing at the x axis
  a_tail <- -0.11 * ymax; a_head <- -0.025 * ymax; lab_y <- -0.15 * ymax
  ggplot(d, aes(GDD_bin, m, colour = Genotype, fill = Genotype)) +
    geom_hline(yintercept = 0, colour = "#9a9a9a", linewidth = 0.4) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.12, colour = NA) +
    geom_line(linewidth = 1.3) +
    geom_segment(data = stages, aes(x = GDD, xend = GDD, y = a_tail, yend = a_head),
                 inherit.aes = FALSE, colour = "black", linewidth = 0.6,
                 arrow = arrow(length = unit(0.16, "cm"), type = "closed")) +
    geom_text(data = stages, aes(x = GDD, y = lab_y, label = Stage),
              inherit.aes = FALSE, size = 3.8, colour = "black") +
    scale_colour_manual(values = COLS) + scale_fill_manual(values = COLS) +
    facet_wrap(~State, ncol = 3) +
    scale_y_continuous(limits = c(-0.19 * ymax, ymax), breaks = pretty(c(0, ymax), 4),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Cumulative growing degree days (\u00b0C d)", y = ylab,
         caption = "Arrows: R1 flowering, R5 beginning seed filling, R7 physiological maturity (mean thermal time in each state).") +
    thm
}

pF <- plot_state("FTSWRZ", "Fraction of transpirable soil water", ymax = 1)
ggsave("modeling_data/figures/Fig_soil_water_season_state_FTSW.png", pF, width = 13, height = 5.2, dpi = 300, bg = "white")
ggsave("modeling_data/figures/Fig_soil_water_season_state_FTSW.tiff", pF, width = 13, height = 5.2, dpi = 600, bg = "white",
       device = "tiff", compression = "lzw", type = "cairo")
ymax_mm <- ceiling(max(pmin(st$ATSWRZ + st$ATSWRZ_sd, Inf), na.rm = TRUE) / 25) * 25
pA <- plot_state("ATSWRZ", "Root-zone transpirable soil water (mm)", ymax = ymax_mm)
ggsave("modeling_data/figures/Fig_soil_water_season_state_ATSW.png", pA, width = 13, height = 5.2, dpi = 300, bg = "white")
ggsave("modeling_data/figures/Fig_soil_water_season_state_ATSW.tiff", pA, width = 13, height = 5.2, dpi = 600, bg = "white",
       device = "tiff", compression = "lzw", type = "cairo")
cat("Saved Fig_soil_water_season_state_FTSW and _ATSW (.png, .tiff)\n")

## headline numbers per state
st |> select(State, Genotype, GDD_bin, ATSWRZ, FTSWRZ) |>
  pivot_wider(names_from = Genotype, values_from = c(ATSWRZ, FTSWRZ)) |>
  group_by(State) |>
  summarise(ATSW_diff_mm = mean(`ATSWRZ_Limited transpiration` - `ATSWRZ_Non-limited transpiration`, na.rm = TRUE),
            FTSW_diff    = mean(`FTSWRZ_Limited transpiration` - `FTSWRZ_Non-limited transpiration`, na.rm = TRUE),
            FTSW_min_NonLT = min(`FTSWRZ_Non-limited transpiration`, na.rm = TRUE),
            FTSW_min_LT    = min(`FTSWRZ_Limited transpiration`, na.rm = TRUE)) |>
  mutate(across(where(is.numeric), ~ round(., 2))) |> print()
