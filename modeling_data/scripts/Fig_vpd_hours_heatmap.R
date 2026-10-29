## Fig_vpd_hours_heatmap.R
## Hours per day with VPD above 2 kPa at the 10 SSM simulation sites (AR, MO, NE),
## mid-March to September, averaged over the simulation years for each calendar day
## and binned (< 1, 1-2, 2-4, > 4 h).
##
## METHOD = "met" (default): model-estimated hourly VPD from the site weather files
##   used for the simulations (modeling_data/outputs/*.met, daily TMAX/TMIN, no humidity):
##   hourly T = sine curve between TMIN and TMAX with the peak at 14:00;
##   air vapour pressure = es(TMIN) all day (dew point = TMIN, as in SSM);
##   VPD(h) = es(T_h) - es(TMIN), es = 0.6108 * exp(17.27 T / (T + 237.3)) (FAO-56 eq. 11)
## METHOD = "power": hourly VPD from NASA POWER hourly T2M and RH2M (local solar time),
##   VPD = es(T) * (1 - RH/100); cached in modeling_data/data/nasapower_hourly_T2M_RH2M_10loc.csv
##
## Output: modeling_data/figures/Fig_vpd_hours_heatmap.png / .tiff
## and a copy for the deck at presentation/img/vpd_heatmap.png.
## Run from the repo root:  Rscript modeling_data/scripts/Fig_vpd_hours_heatmap.R
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(readr) })

METHOD     <- "met"
YEARS      <- 1990:2019   # SSM simulation period
VPD_THRESH <- 2
CACHE      <- "modeling_data/data/nasapower_hourly_T2M_RH2M_10loc.csv"
OUT        <- "modeling_data/figures/Fig_vpd_hours_heatmap"

sites <- tribble(
  ~Location,      ~State, ~lat,  ~lon,    ~met,
  "Rohwer",       "AR",   33.80, -91.27,  "Rowher_AR.met",
  "Marianna",     "AR",   34.77, -90.76,  "Marianna_AR.met",
  "Keiser",       "AR",   35.67, -90.09,  "Keiser_AR.met",
  "Jonesboro",    "AR",   35.84, -90.70,  "Jonesboro_AR.met",
  "Mount Vernon", "MO",   37.10, -93.82,  "MountVernon_MO.met",
  "Novelty",      "MO",   40.01, -92.21,  "Novelty_MO.met",
  "Albany",       "MO",   40.25, -94.33,  "Albany_MO.met",
  "Eustis",       "NE",   40.66, -100.03, "Eustis_NE.met",
  "Lincoln",      "NE",   40.81, -96.70,  "Lincoln_NE.met",
  "North Platte", "NE",   41.12, -100.77, "NorthPlatte_NE.met"
)
es <- function(t) 0.6108 * exp(17.27 * t / (t + 237.3))

## ── Hours per day above threshold ───────────────────────────────────────────────
if (METHOD == "met") {
  hrs <- 0:23
  daily <- bind_rows(lapply(seq_len(nrow(sites)), function(i) {
    read_table(file.path("modeling_data/outputs", sites$met[i]), skip = 6, show_col_types = FALSE,
               col_names = c("year", "day", "radn", "maxt", "mint", "rain")) |>
      filter(year %in% YEARS) |>
      mutate(Location = sites$Location[i], State = sites$State[i],
             date = as.Date(day - 1, origin = paste0(year, "-01-01")))
  })) |>
    rowwise() |>
    mutate(hours = sum(es((maxt + mint) / 2 + (maxt - mint) / 2 * cos(2 * pi * (hrs - 14) / 24)) -
                       es(mint) > VPD_THRESH)) |>
    ungroup()
} else {
  library(jsonlite)
  get_power <- function(lat, lon, y0, y1) {
    url <- sprintf(paste0("https://power.larc.nasa.gov/api/temporal/hourly/point?parameters=T2M,RH2M",
                          "&community=AG&longitude=%.2f&latitude=%.2f&start=%d0101&end=%d1231",
                          "&format=JSON&time-standard=LST"), lon, lat, y0, y1)
    p <- fromJSON(url)$properties$parameter
    tibble(time = names(p$T2M), T2M = unlist(p$T2M), RH2M = unlist(p$RH2M)[names(p$T2M)])
  }
  hourly <- if (file.exists(CACHE)) read_csv(CACHE, show_col_types = FALSE, col_types = "cciddd") else NULL
  if (is.null(hourly) || !all(YEARS %in% unique(substr(hourly$time, 1, 4)))) {
    hourly <- bind_rows(lapply(seq_len(nrow(sites)), function(i) {
      cat("Downloading", sites$Location[i], "\n")
      get_power(sites$lat[i], sites$lon[i], min(YEARS), max(YEARS)) |>
        mutate(Location = sites$Location[i], State = sites$State[i], .before = 1)
    })) |> mutate(hour = as.integer(substr(time, 9, 10)), .after = time)
    write_csv(hourly, CACHE)
  }
  daily <- hourly |>
    filter(T2M > -900, RH2M > -900, as.integer(substr(time, 1, 4)) %in% YEARS) |>
    mutate(date = as.Date(substr(time, 1, 8), "%Y%m%d"),
           vpd  = es(T2M) * (1 - pmin(RH2M, 100) / 100)) |>
    group_by(Location, State, date) |>
    summarise(hours = sum(vpd > VPD_THRESH), n = n(), .groups = "drop") |>
    filter(n >= 20)
}

## ── Mean across years per calendar day ──────────────────────────────────────────
daily <- daily |>
  mutate(md = as.Date(format(date, "2001-%m-%d"))) |>   # common non-leap year for the x axis
  filter(!is.na(md), md >= as.Date("2001-03-15"), md <= as.Date("2001-09-30")) |>
  group_by(Location, State, md) |>
  summarise(hours = mean(hours), .groups = "drop") |>
  mutate(cls = cut(hours, c(-Inf, 1, 2, 4, Inf), right = FALSE,
                   labels = c("< 1 h", "1 to 2 h", "2 to 4 h", "> 4 h")),
         site = factor(paste0(Location, ", ", State),
                       levels = rev(paste0(sites$Location, ", ", sites$State))))

cat("Method:", METHOD, " years:", paste(range(YEARS), collapse = "-"), " days per site:",
    paste(range(count(daily, site)$n), collapse = "-"), "\n")
print(daily |> group_by(State) |> summarise(mean_hours = round(mean(hours), 2),
                                            days_over_4h = round(sum(hours >= 4) / n_distinct(Location))))

## ── Plot ────────────────────────────────────────────────────────────────────────
p <- ggplot(daily, aes(md, site, fill = cls)) +
  geom_tile(height = 0.92) +
  scale_fill_manual(values = c("< 1 h" = "#d9d9d9", "1 to 2 h" = "#b2e2e2",
                               "2 to 4 h" = "#66c2a4", "> 4 h" = "#238b45"),
                    drop = FALSE,
                    name = sprintf("Hours/day\nabove %g kPa\n(%d-%d mean)", VPD_THRESH, min(YEARS), max(YEARS))) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b", expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 0.6),
        axis.text = element_text(colour = "black"),
        axis.ticks.y = element_blank(),
        legend.title = element_text(size = 10), legend.text = element_text(size = 9),
        plot.background = element_rect(fill = "white", colour = NA))

ggsave(paste0(OUT, ".png"), p, width = 7.5, height = 3.8, dpi = 300)
ggsave(paste0(OUT, ".tiff"), p, width = 7.5, height = 3.8, dpi = 600, compression = "lzw")
file.copy(paste0(OUT, ".png"), "presentation/img/vpd_heatmap.png", overwrite = TRUE)
cat("Saved", paste0(OUT, ".png / .tiff"), "and presentation/img/vpd_heatmap.png\n")
