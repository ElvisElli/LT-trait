
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
                   "ALLSKY_SFC_SW_DWN"),
    dates      = c(paste0(min(years), "-01-01"), paste0(max(years), "-12-31")),
    temporal_api = "daily"
  )
}

raw_list <- sites %>%
  rowwise() %>%
  mutate(data = list(pull_nasapower(lon, lat, years))) %>%
  ungroup()

## ── 3. FAO-56 Penman-Monteith ETo ─────────────────────────────────────────────
# Reference: Allen et al. (1998), FAO Irrigation & Drainage Paper 56.
# WS2M from NASA POWER is already in m/s; no unit conversion applied.

calc_eto <- function(df, lat, elev_m) {

  LAT <- lat
  ALT <- elev_m

  df %>%
    mutate(
      Tmax     = T2M_MAX,
      Tmin     = T2M_MIN,
      RH_avg   = RH2M,
      Rad_down = ALLSKY_SFC_SW_DWN,          # MJ m-2 day-1
      WS       = WS2M,                        # m/s from NASA POWER

      ## solar geometry
      decl     = 23.45 * sin(((360 * (DOY - 80)) / 365) * (pi / 180)),
      hn       = (acos(-1 * (tan(LAT * (pi / 180)) * tan(decl * (pi / 180))))) *
                 (180 / pi),
      N        = (2 / 15) * hn,
      dD       = 1 + 0.033 * cos(((DOY * 360) / 365) * (pi / 180)),
      Qo       = 37.6 * dD *
                 ((pi / 180) * hn * sin(LAT * (pi / 180)) * sin(decl * (pi / 180)) +
                  cos(LAT * (pi / 180)) * cos(decl * (pi / 180)) * sin(hn * (pi / 180))),

      ## vapour pressures
      tmed     = (Tmax + Tmin) / 2,
      es       = 0.6108 * 10 ^ (7.5 * tmed / (237.3 + tmed)),
      ea       = es * (RH_avg / 100),
      vpd      = es - ea,

      ## net radiation
      BOC      = Rad_down * (1 - 0.23),
      Qgo      = (0.75 + 2e-5 * ALT) * Qo,
      BOL      = -4.903e-9 * ((tmed + 273) ^ 4 *
                  (0.34 - 0.14 * ea ^ 0.5) *
                  (1.35 * (Rad_down / Qgo) - 0.35)),
      Rn       = BOC + BOL,

      ## slope of vapour pressure curve
      s.vapor  = (4098 * es) / (tmed + 237.3) ^ 2,

      ## FAO-56 PM ETo (mm day-1); γ = 0.062 kPa °C-1
      ETo      = (0.408 * Rn * s.vapor +
                  0.062 * 900 * WS * (es - ea) / (tmed + 273)) /
                 (s.vapor + 0.062 * (1 + 0.34 * WS))
    ) %>%
    select(YEAR, DOY, YYYYMMDD, Tmax, Tmin, tmed, RH_avg, WS,
           Rad_down, Qo, Rn, es, ea, vpd, ETo)
}

## Apply ETo calculation to each site
eto_list <- raw_list %>%
  rowwise() %>%
  mutate(eto = list(calc_eto(data, lat, elev_m))) %>%
  ungroup()

## Combine into a single tidy data frame
eto_all <- eto_list %>%
  select(site, lat, lon, eto) %>%
  unnest(eto) %>%
  mutate(Date  = as.Date(YYYYMMDD),
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
