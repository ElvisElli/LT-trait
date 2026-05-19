## Headless runner — sources the ETo pipeline without RStudio dependencies.
## Skips SSM merge (SoybeanSSM_sims2.xlsx not available in this environment).

setwd("/home/user/LT-trait")

library(nasapower)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(scales)
library(purrr)

## ── 1. Site parameters ────────────────────────────────────────────────────────

sites <- tibble(
  State  = c("AR",              "MO",           "NE"),
  site   = c("Fayetteville_AR", "Columbia_MO",  "Lincoln_NE"),
  lon    = c(-94.17,            -92.33,          -96.70),
  lat    = c( 36.08,             38.95,           40.82),
  elev_m = c( 381,               274,             357)
)

years <- 1990:2023

## ── 2. Download NASA POWER daily data ─────────────────────────────────────────

pull_nasapower <- function(lon, lat, years) {
  get_power(
    community    = "ag",
    lonlat       = c(lon, lat),
    pars         = c("T2M_MAX", "T2M_MIN", "RH2M", "WS2M", "ALLSKY_SFC_SW_DWN"),
    dates        = c(paste0(min(years), "-01-01"), paste0(max(years), "-12-31")),
    temporal_api = "daily"
  )
}

message("Downloading NASA POWER data for all sites...")
raw_list <- sites %>%
  rowwise() %>%
  mutate(data = list(pull_nasapower(lon, lat, years))) %>%
  ungroup()
message("Download complete.")

## ── 3. FAO-56 PM ETo ──────────────────────────────────────────────────────────

calc_eto <- function(df, lat, elev_m) {
  LAT <- lat
  ALT <- elev_m
  df %>%
    mutate(
      Tmax     = T2M_MAX,
      Tmin     = T2M_MIN,
      RH_avg   = RH2M,
      Rad_down = ALLSKY_SFC_SW_DWN,
      WS       = WS2M,
      decl     = 23.45 * sin(((360 * (DOY - 80)) / 365) * (pi / 180)),
      hn       = acos(-tan(LAT * (pi/180)) * tan(decl * (pi/180))) * (180/pi),
      N        = (2/15) * hn,
      dD       = 1 + 0.033 * cos(((DOY * 360) / 365) * (pi / 180)),
      Qo       = 37.6 * dD *
                 ((pi/180) * hn * sin(LAT*(pi/180)) * sin(decl*(pi/180)) +
                  cos(LAT*(pi/180)) * cos(decl*(pi/180)) * sin(hn*(pi/180))),
      tmed     = (Tmax + Tmin) / 2,
      es       = 0.6108 * exp(17.27 * tmed / (tmed + 237.3)),
      ea       = es * (RH_avg / 100),
      vpd      = es - ea,
      BOC      = Rad_down * (1 - 0.23),
      Qgo      = (0.75 + 2e-5 * ALT) * Qo,
      BOL      = -4.903e-9 * ((tmed + 273)^4 *
                  (0.34 - 0.14 * sqrt(ea)) *
                  (1.35 * (Rad_down / Qgo) - 0.35)),
      Rn       = BOC + BOL,
      s.vapor  = (4098 * es) / (tmed + 237.3)^2,
      ETo      = pmax(0,
                   (0.408 * Rn * s.vapor +
                    0.062 * 900 * WS * vpd / (tmed + 273)) /
                   (s.vapor + 0.062 * (1 + 0.34 * WS)))
    ) %>%
    select(YEAR, DOY, YYYYMMDD, Tmax, Tmin, tmed, RH_avg, WS,
           Rad_down, Qo, Rn, es, ea, vpd, ETo)
}

eto_list <- raw_list %>%
  rowwise() %>%
  mutate(eto = list(calc_eto(data, lat, elev_m))) %>%
  ungroup()

eto_daily <- eto_list %>%
  select(State, site, lat, lon, eto) %>%
  unnest(eto) %>%
  mutate(Date  = as.Date(YYYYMMDD),
         Month = month(Date),
         Year  = year(Date))

## ── 4. Seasonal ETo ───────────────────────────────────────────────────────────

eto_seasonal <- eto_daily %>%
  filter(DOY >= 120, DOY <= 300) %>%
  group_by(State, site, Year) %>%
  summarise(ETo_season_mm = sum(ETo, na.rm = TRUE), .groups = "drop")

## ── 5. Climatology plot ───────────────────────────────────────────────────────

eto_clim <- eto_daily %>%
  filter(DOY >= 120, DOY <= 300) %>%
  group_by(site, Year, Month) %>%
  summarise(ETo_mm = sum(ETo, na.rm = TRUE), .groups = "drop") %>%
  group_by(site, Month) %>%
  summarise(ETo_mean = mean(ETo_mm),
            ETo_sd   = sd(ETo_mm),
            .groups  = "drop")

p_clim <- eto_clim %>%
  ggplot(aes(x = Month, y = ETo_mean, colour = site, fill = site)) +
  geom_ribbon(aes(ymin = pmax(0, ETo_mean - ETo_sd),
                  ymax = ETo_mean + ETo_sd),
              alpha = 0.18, colour = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.2) +
  scale_x_continuous(breaks = 5:10, labels = month.abb[5:10]) +
  labs(x = NULL, y = "ETo (mm month⁻¹)",
       title = "Growing-season ETo climatology (DOY 120–300)",
       subtitle = paste0(min(years), "–", max(years)),
       colour = "Site", fill = "Site") +
  theme_bw(base_size = 12) +
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold"))

ggsave("ETo_climatology_growingseason.tiff", plot = p_clim,
       width = 20, height = 13, units = "cm",
       dpi = 300, compression = "lzw", bg = "white")
message("Saved: ETo_climatology_growingseason.tiff")

## ── 6. Export CSVs ────────────────────────────────────────────────────────────

write.csv(eto_daily,    "nasapower_ETo_daily.csv",    row.names = FALSE)
write.csv(eto_seasonal, "nasapower_ETo_seasonal.csv", row.names = FALSE)
message("Saved: nasapower_ETo_daily.csv, nasapower_ETo_seasonal.csv")
message("Done. Run nasapower_PM_ETo.R in RStudio with SoybeanSSM_sims2.xlsx to generate the merged Kc plot.")
