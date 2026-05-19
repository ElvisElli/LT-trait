
## Script name: nasapower_PM_ETo.R ====

## Script objective: Download NASA POWER daily climate data and compute
## FAO-56 Penman-Monteith reference evapotranspiration (ETo) ====

## Author: Elvis F. Elli ====

## Script was created on: 2026-05-19 ====

## Cleaning up environment ====
rm(list=ls())

## Libraries ====
library(nasapower)
library(tidyverse)
library(lubridate)

## ── 1. Site parameters ────────────────────────────────────────────────────────

sites <- tibble(
  site     = c("Fayetteville_AR", "Columbia_MO", "Lincoln_NE"),
  lon      = c(-94.17,            -92.33,         -96.70),
  lat      = c( 36.08,             38.95,          40.82),
  elev_m   = c( 381,               274,            357)   # metres above sea level
)

years <- 1990:2023   # NASA POWER coverage

## ── 2. Download NASA POWER daily data ─────────────────────────────────────────

pull_nasapower <- function(lon, lat, years) {
  get_power(
    community  = "ag",
    lonlat     = c(lon, lat),
    pars       = c("T2M_MAX", "T2M_MIN", "RH2M", "WS2M",
                   "ALLSKY_SFC_SW_DWN", "ALLSKY_SFC_LW_DWN"),
    dates      = c(paste0(min(years), "-01-01"), paste0(max(years), "-12-31")),
    temporal_api = "daily"
  )
}

raw_list <- sites %>%
  rowwise() %>%
  mutate(data = list(pull_nasapower(lon, lat, years))) %>%
  ungroup()

## ── 3. FAO-56 Penman-Monteith ETo ─────────────────────────────────────────────
# ETo [mm/day] = (0.408·Δ·(Rn-G) + γ·(900/(T+273))·u2·(es-ea))
#                / (Δ + γ·(1+0.34·u2))
# Reference: Allen et al. (1998), FAO Irrigation & Drainage Paper 56.

calc_eto <- function(df, elev_m) {

  df %>%
    mutate(
      ## Basic met variables
      Tmax   = T2M_MAX,
      Tmin   = T2M_MIN,
      Tmean  = (Tmax + Tmin) / 2,
      RH     = RH2M,
      u2     = WS2M,
      Rs     = ALLSKY_SFC_SW_DWN,   # MJ m-2 day-1
      Rl_down = ALLSKY_SFC_LW_DWN,  # MJ m-2 day-1

      ## Atmospheric pressure (kPa) – FAO eq. 7
      P      = 101.3 * ((293 - 0.0065 * elev_m) / 293) ^ 5.26,

      ## Psychrometric constant γ (kPa °C-1) – FAO eq. 8
      gamma  = 0.000665 * P,

      ## Saturation vapour pressure (kPa) – FAO eq. 11 & 12
      es_max = 0.6108 * exp(17.27 * Tmax / (Tmax + 237.3)),
      es_min = 0.6108 * exp(17.27 * Tmin / (Tmin + 237.3)),
      es     = (es_max + es_min) / 2,

      ## Actual vapour pressure from mean RH – FAO eq. 19
      ea     = RH / 100 * es,

      ## Slope of saturation vapour pressure curve Δ (kPa °C-1) – FAO eq. 13
      Delta  = 4098 * (0.6108 * exp(17.27 * Tmean / (Tmean + 237.3))) /
               (Tmean + 237.3) ^ 2,

      ## Net solar (shortwave) radiation Rns (MJ m-2 day-1) – FAO eq. 38
      ## albedo α = 0.23 for grass reference surface
      Rns    = (1 - 0.23) * Rs,

      ## Net longwave radiation Rnl (MJ m-2 day-1) – FAO eq. 39
      ## Stefan-Boltzmann σ = 4.903e-9 MJ K-4 m-2 day-1
      sigma  = 4.903e-9,
      Rnl    = sigma *
               ((( Tmax + 273.16) ^ 4 + (Tmin + 273.16) ^ 4) / 2) *
               (0.34 - 0.14 * sqrt(ea)) *
               pmax(0, (1.35 * (Rs / (0.75 * Rs + 0.00001)) - 0.35)),

      ## Net radiation Rn (MJ m-2 day-1) – FAO eq. 40
      Rn     = Rns - Rnl,

      ## Soil heat flux G ≈ 0 for daily time-step – FAO eq. 42
      G      = 0,

      ## FAO-56 PM ETo (mm day-1)
      ETo    = (0.408 * Delta * (Rn - G) +
                gamma * (900 / (Tmean + 273)) * u2 * (es - ea)) /
               (Delta + gamma * (1 + 0.34 * u2))
    ) %>%
    select(YEAR, DOY, YYYYMMDD, Tmax, Tmin, Tmean, RH, u2,
           Rs, Rn, es, ea, ETo)
}

## Apply ETo calculation to each site
eto_list <- raw_list %>%
  rowwise() %>%
  mutate(eto = list(calc_eto(data, elev_m))) %>%
  ungroup()

## Combine into a single tidy data frame
eto_all <- eto_list %>%
  select(site, lat, lon, eto) %>%
  unnest(eto) %>%
  mutate(Date = as.Date(YYYYMMDD),
         Month = month(Date),
         Year  = year(Date))

## ── 4. Summary statistics ──────────────────────────────────────────────────────

## Monthly mean ETo by site
eto_monthly <- eto_all %>%
  group_by(site, Year, Month) %>%
  summarise(ETo_mm_month = sum(ETo, na.rm = TRUE), .groups = "drop")

## Long-term monthly climatology
eto_clim <- eto_monthly %>%
  group_by(site, Month) %>%
  summarise(ETo_mean = mean(ETo_mm_month, na.rm = TRUE),
            ETo_sd   = sd(ETo_mm_month,   na.rm = TRUE),
            .groups  = "drop")

## Annual ETo
eto_annual <- eto_all %>%
  group_by(site, Year) %>%
  summarise(ETo_mm_year = sum(ETo, na.rm = TRUE), .groups = "drop")

## ── 5. Plots ───────────────────────────────────────────────────────────────────

## Monthly climatology
eto_clim %>%
  ggplot(aes(x = Month, y = ETo_mean, colour = site, fill = site)) +
  geom_ribbon(aes(ymin = ETo_mean - ETo_sd,
                  ymax = ETo_mean + ETo_sd), alpha = 0.2, colour = NA) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 1:12,
                     labels = month.abb) +
  labs(x = NULL,
       y = "ETo (mm month⁻¹)",
       title = "FAO-56 Penman-Monteith ETo – monthly climatology",
       subtitle = paste0(min(years), "–", max(years)),
       colour = "Site", fill = "Site") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

ggsave("ETo_monthly_climatology.tiff",
       width = 20, height = 14, units = "cm",
       dpi = 600, compression = "lzw", bg = "white")

## Annual ETo time series
eto_annual %>%
  ggplot(aes(x = Year, y = ETo_mm_year, colour = site)) +
  geom_line(size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", size = 0.6) +
  labs(x = NULL,
       y = "ETo (mm year⁻¹)",
       title = "Annual reference evapotranspiration (FAO-56 PM)",
       colour = "Site") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

ggsave("ETo_annual_trend.tiff",
       width = 20, height = 14, units = "cm",
       dpi = 600, compression = "lzw", bg = "white")

## ── 6. Export ──────────────────────────────────────────────────────────────────

write_csv(eto_all,     "nasapower_ETo_daily.csv")
write_csv(eto_monthly, "nasapower_ETo_monthly.csv")
write_csv(eto_annual,  "nasapower_ETo_annual.csv")
