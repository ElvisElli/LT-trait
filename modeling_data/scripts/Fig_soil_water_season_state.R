## Fig_soil_water_season_state.R
## Seasonal course of root-zone soil water under rainfed conditions, limited
## transpiration (LT1.5) vs non-limited transpiration (Non-LT), one panel per
## state (Arkansas, Missouri, Nebraska). Lines are the mean across all
## seasons (30 weather-years x 3 sowing dates) of all locations in the state,
## no variability band. Each season is cut at physiological maturity.
##
## Input : modeling_data/data/soil_water_season_summary.csv
##         (per-location daily means and n, from Fig_soil_water_season.R)
## Output: modeling_data/figures/Fig_soil_water_season_state_ATSW.png  (mm)
##         modeling_data/figures/Fig_soil_water_season_state_FTSW.png  (fraction)
##         modeling_data/data/soil_water_season_state_summary.csv
## Run from the repo root:  Rscript modeling_data/scripts/Fig_soil_water_season_state.R
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(readr); library(patchwork) })

COLS  <- c("Non-limited transpiration" = "#31a354", "Limited transpiration" = "#d95f02")
STATE <- c(Jonesboro = "Arkansas", Keiser = "Arkansas", Marianna = "Arkansas", Rohwer = "Arkansas",
           Albany = "Missouri", `Mount Vernon` = "Missouri", Novelty = "Missouri",
           Eustis = "Nebraska", Lincoln = "Nebraska", `North Platte` = "Nebraska")

loc <- read_csv("modeling_data/data/soil_water_season_summary.csv", show_col_types = FALSE) |>
  mutate(State = STATE[Location])

## pooled mean across the seasons of every location in the state
## (location means weighted by their number of seasons still growing that day)
st <- loc |>
  group_by(State, Genotype, DAP) |>
  summarise(ATSWRZ = weighted.mean(ATSWRZ_mean, n),
            FTSWRZ = weighted.mean(FTSWRZ_mean, n),
            n_loc = n_distinct(Location), n = sum(n), .groups = "drop") |>
  group_by(State, Genotype) |> filter(n >= 0.5 * max(n)) |> ungroup() |>
  mutate(State = factor(State, levels = c("Arkansas", "Missouri", "Nebraska")),
         Genotype = factor(Genotype, levels = names(COLS)))
write_csv(st, "modeling_data/data/soil_water_season_state_summary.csv")

thm <- theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#ececec"),
        panel.border = element_rect(colour = "black", linewidth = 0.7),
        strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.7),
        strip.text = element_text(size = 14, colour = "black"),
        axis.text = element_text(colour = "black"), legend.position = "top",
        legend.title = element_blank(), legend.text = element_text(size = 13),
        plot.background = element_rect(fill = "white", colour = NA))

plot_state <- function(var, ylab) {
  ggplot(st, aes(DAP, .data[[var]], colour = Genotype)) +
    geom_line(linewidth = 1.3) +
    scale_colour_manual(values = COLS) +
    facet_wrap(~State, ncol = 3) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
    labs(x = "Days after planting", y = ylab) + thm
}

pA <- plot_state("ATSWRZ", "Root-zone transpirable soil water (mm)")
ggsave("modeling_data/figures/Fig_soil_water_season_state_ATSW.png", pA, width = 13, height = 4.8, dpi = 300, bg = "white")
pF <- plot_state("FTSWRZ", "Fraction of transpirable soil water")
ggsave("modeling_data/figures/Fig_soil_water_season_state_FTSW.png", pF, width = 13, height = 4.8, dpi = 300, bg = "white")
cat("Saved Fig_soil_water_season_state_ATSW.png and _FTSW.png\n")

## headline numbers per state
st |> select(State, Genotype, DAP, ATSWRZ, FTSWRZ) |>
  pivot_wider(names_from = Genotype, values_from = c(ATSWRZ, FTSWRZ)) |>
  group_by(State) |>
  summarise(ATSW_diff_mm = mean(`ATSWRZ_Limited transpiration` - `ATSWRZ_Non-limited transpiration`, na.rm = TRUE),
            ATSW_diff_max = max(`ATSWRZ_Limited transpiration` - `ATSWRZ_Non-limited transpiration`, na.rm = TRUE),
            FTSW_diff = mean(`FTSWRZ_Limited transpiration` - `FTSWRZ_Non-limited transpiration`, na.rm = TRUE),
            FTSW_min_NonLT = min(`FTSWRZ_Non-limited transpiration`, na.rm = TRUE),
            FTSW_min_LT = min(`FTSWRZ_Limited transpiration`, na.rm = TRUE)) |>
  mutate(across(where(is.numeric), ~ round(., 2))) |> print()
