## Headless runner — ETo pipeline + merged Kc plot (requires SoybeanSSM_sims2.xlsx).
## Saves PNG (for GitHub preview) alongside TIFF (for publication).

setwd("/home/user/LT-trait")

library(nasapower)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(scales)
library(purrr)
library(readxl)

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
       width = 20, height = 13, units = "cm", dpi = 300,
       compression = "lzw", bg = "white")
ggsave("ETo_climatology_growingseason.png", plot = p_clim,
       width = 20, height = 13, units = "cm", dpi = 150, bg = "white")
message("Saved: ETo_climatology_growingseason (.tiff + .png)")

## ── 6. Load and process SSM data ─────────────────────────────────────────────

ssm.model <- read_excel("SoybeanSSM_sims2.xlsx")

ssm_proc <- ssm.model %>%
  filter(Crop %in% c("SOYBEAN-MG4", "SOYBEAN-LT1.5-MG4",
                     "SOYBEAN-LT1.5-MG3", "SOYBEAN-MG3")) %>%
  mutate(
    Genotype = ifelse(Crop %in% c("SOYBEAN-MG4", "SOYBEAN-MG3"),
                      "Non-limited transpiration", "Limited transpiration"),
    Genotype = factor(Genotype,
                      levels = c("Non-limited transpiration", "Limited transpiration")),
    State = case_when(
      grepl("AR", Manag, ignore.case = TRUE) ~ "AR",
      grepl("MO", Manag, ignore.case = TRUE) ~ "MO",
      grepl("NE", Manag, ignore.case = TRUE) ~ "NE",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(State)) %>%
  rename(Year = Pyear)

## ── 7. Merge ETo with SSM data ────────────────────────────────────────────────

merged <- ssm_proc %>%
  left_join(eto_seasonal, by = c("State", "Year")) %>%
  filter(!is.na(ETo_season_mm), !is.na(SUMETT))

## ── 8. Main plot: crop ET vs ETo — LT vs non-LT ──────────────────────────────

kc_slopes <- merged %>%
  group_by(Genotype, water_mgmt) %>%
  summarise(
    slope = coef(lm(SUMETT ~ ETo_season_mm))[["ETo_season_mm"]],
    x_pos = quantile(ETo_season_mm, 0.75, na.rm = TRUE),
    y_pos = quantile(SUMETT,        0.10, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("Kc = ", round(slope, 2)))

p_main <- merged %>%
  ggplot(aes(x = ETo_season_mm, y = SUMETT,
             colour = Genotype, fill = Genotype)) +
  geom_point(alpha = 0.45, shape = 21, size = 2.2,
             stroke = 0.3, colour = "grey20") +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.1, alpha = 0.15) +
  geom_text(data = kc_slopes,
            aes(x = x_pos, y = y_pos, label = label, colour = Genotype),
            size = 3.3, fontface = "italic", show.legend = FALSE) +
  facet_grid(water_mgmt ~ State) +
  scale_colour_manual(values = c("Non-limited transpiration" = "#31a354",
                                 "Limited transpiration"      = "#d95f02")) +
  scale_fill_manual(values   = c("Non-limited transpiration" = "#31a354",
                                 "Limited transpiration"      = "#d95f02")) +
  scale_x_continuous(breaks = breaks_pretty(4)) +
  scale_y_continuous(breaks = breaks_pretty(4)) +
  labs(
    x       = "Seasonal ETo (mm, DOY 120–300)",
    y       = "Season crop evapotranspiration — SUMETT (mm)",
    colour  = NULL, fill = NULL,
    title   = "Crop water demand relative to atmospheric evaporative demand",
    subtitle = "Kc slope: how much crop ET rises per unit of ETo — lower slope = more water conservative"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "bottom",
    legend.key.size  = unit(0.8, "lines"),
    strip.background = element_rect(fill = "grey92", colour = NA),
    strip.text       = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(colour = "grey40", size = 9),
    axis.title       = element_text(size = 10)
  )

ggsave("ETo_vs_cropET_LT_nonLT.tiff", plot = p_main,
       width = 26, height = 18, units = "cm", dpi = 300,
       compression = "lzw", bg = "white")
ggsave("ETo_vs_cropET_LT_nonLT.png", plot = p_main,
       width = 26, height = 18, units = "cm", dpi = 150, bg = "white")
message("Saved: ETo_vs_cropET_LT_nonLT (.tiff + .png)")

## ── 9. Export CSVs ────────────────────────────────────────────────────────────

write.csv(eto_daily,    "nasapower_ETo_daily.csv",    row.names = FALSE)
write.csv(eto_seasonal, "nasapower_ETo_seasonal.csv", row.names = FALSE)
write.csv(merged,       "merged_ETo_cropET.csv",      row.names = FALSE)
message("Saved CSVs: daily, seasonal, merged")
message("All done.")
