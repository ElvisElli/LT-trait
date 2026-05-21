
## Script name: nasapower_PM_ETo.R ====

## Script objective: Download NASA POWER daily climate data, compute FAO-56
## Penman-Monteith ETo, and merge with DSSAT SSM output to compare seasonal
## crop water demand (SUMETT) against atmospheric evaporative demand (ETo)
## for limited-transpiration (LT) vs. non-limited-transpiration soybean ====

## Author: Elvis F. Elli ====

## Script was created on: 2026-05-19 ====

## ── Path configuration ────────────────────────────────────────────────────────
## Cloud / Linux: runs as-is, no changes needed
## Windows (local): edit the two paths below and uncomment them
# project_root <- "C:/Users/YourName/Projects/LT-trait"   # ← git repo location
# box_data     <- "C:/Users/YourName/Box/LT-trait-data"   # ← large files in Box
if (.Platform$OS.type == "windows") {
  # setwd(project_root)                                    # ← uncomment on Windows
} else {
  setwd("/home/user/LT-trait")
}
## ─────────────────────────────────────────────────────────────────────────────

## Cleaning up environment ====
rm(list=ls())

## Libraries ====
library(nasapower)
library(tidyverse)
library(lubridate)
library(readxl)

## Set working directory ====
setwd("/home/user/LT-trait")


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

raw_list <- sites %>%
  rowwise() %>%
  mutate(data = list(pull_nasapower(lon, lat, years))) %>%
  ungroup()


## ── 3. FAO-56 Penman-Monteith ETo ─────────────────────────────────────────────
# es: FAO-56 standard exp form (Allen et al. 1998, eq. 11) — input in Celsius.
# The base-10 Tetens form (10^[7.5·T/(237.3+T)]) is also for Celsius and agrees
# within ~0.1%; the exp form is used here to match the FAO-56 reference exactly.
#
# WS2M from NASA POWER is already in m/s — no unit conversion needed.
# Negative daily ETo (possible in winter when Rn < 0) is clamped to 0;
# agronomically meaningless during the non-growing season.

calc_eto <- function(df, lat, elev_m) {
  LAT <- lat
  ALT <- elev_m
  df %>%
    mutate(
      Tmax     = T2M_MAX,
      Tmin     = T2M_MIN,
      RH_avg   = RH2M,
      Rad_down = ALLSKY_SFC_SW_DWN,   # MJ m-2 day-1
      WS       = WS2M,                 # m/s

      ## solar geometry
      decl  = 23.45 * sin(((360 * (DOY - 80)) / 365) * (pi / 180)),
      hn    = acos(-tan(LAT * (pi/180)) * tan(decl * (pi/180))) * (180/pi),
      N     = (2/15) * hn,
      dD    = 1 + 0.033 * cos(((DOY * 360) / 365) * (pi / 180)),
      Qo    = 37.6 * dD *
              ((pi/180) * hn * sin(LAT*(pi/180)) * sin(decl*(pi/180)) +
               cos(LAT*(pi/180)) * cos(decl*(pi/180)) * sin(hn*(pi/180))),

      ## vapour pressures
      tmed    = (Tmax + Tmin) / 2,
      es      = 0.6108 * exp(17.27 * tmed / (tmed + 237.3)),  # kPa, FAO-56 eq.11
      ea      = es * (RH_avg / 100),
      vpd     = es - ea,

      ## net radiation
      BOC  = Rad_down * (1 - 0.23),
      Qgo  = (0.75 + 2e-5 * ALT) * Qo,
      BOL  = -4.903e-9 * ((tmed + 273)^4 *
              (0.34 - 0.14 * sqrt(ea)) *
              (1.35 * (Rad_down / Qgo) - 0.35)),
      Rn   = BOC + BOL,

      ## slope of vapour pressure curve & ETo
      s.vapor = (4098 * es) / (tmed + 237.3)^2,

      ## γ = 0.062 kPa °C-1 (fixed); WS in m/s
      ETo = pmax(0,
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


## ── 4. Seasonal ETo (growing season DOY 120–300 ≈ May–Oct) ───────────────────

eto_seasonal <- eto_daily %>%
  filter(DOY >= 120, DOY <= 300) %>%
  group_by(State, site, Year) %>%
  summarise(ETo_season_mm = sum(ETo, na.rm = TRUE), .groups = "drop")


## ── 5. Load and process SSM simulation data ───────────────────────────────────

ssm.model <- read_excel("modeling_data/data/SoybeanSSM_sims2.xlsx")

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


## ── 6. Merge ETo with SSM data ────────────────────────────────────────────────

merged <- ssm_proc %>%
  left_join(eto_seasonal, by = c("State", "Year")) %>%
  filter(!is.na(ETo_season_mm), !is.na(SUMETT))


## ── 7. Main plot: crop ET vs ETo — LT vs non-LT ──────────────────────────────

## Fit Kc slope per genotype x water management
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
  facet_grid(Genotype ~ State) +
  scale_colour_manual(values = c("Non-limited transpiration" = "#31a354",
                                 "Limited transpiration"      = "#d95f02")) +
  scale_fill_manual(values   = c("Non-limited transpiration" = "#31a354",
                                 "Limited transpiration"      = "#d95f02")) +
  scale_x_continuous(breaks  = scales::breaks_pretty(4)) +
  scale_y_continuous(breaks  = scales::breaks_pretty(4)) +
  labs(
    x        = "Seasonal ETo (mm, DOY 120–300)",
    y        = "Season crop evapotranspiration — SUMETT (mm)",
    colour   = NULL, fill = NULL,
    title    = "Crop water demand relative to atmospheric evaporative demand",
    subtitle = "Kc slope: how much crop ET rises per unit of ETo — lower slope = more water conservative\n(cf. Marin et al. 2016, 2019: Kc decreases at high ETo environments)"
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

p_main

ggsave("modeling_data/figures/Fig_kc_seasonal_scatter.png",
       plot = p_main,
       width = 26, height = 18, units = "cm",
       dpi = 150, bg = "white")


## ── 7b. Marin-style daily ETc vs ETo scatter (following Marin et al. 2016, 2019) ──
## Convert seasonal totals to mean daily values:
## ETc_daily = SUMETT / R7.PM  (mm d-1, over the crop cycle)
## ETo_daily = ETo_season_mm / 181  (mm d-1, DOY 120-300 mean)
## The regression slope = seasonal Kc; Marin et al. found Kc inversely related to ETo

daily_etc <- merged %>%
  filter(!is.na(R7.PM), R7.PM > 0, !is.na(SUMETT)) %>%
  mutate(
    ETc_daily = SUMETT / R7.PM,
    ETo_daily = ETo_season_mm / 181,
    Genotype  = factor(Genotype, levels = c("Non-limited transpiration","Limited transpiration"))
  )

kc_daily <- daily_etc %>%
  group_by(Genotype, water_mgmt) %>%
  summarise(
    Kc_slope = coef(lm(ETc_daily ~ ETo_daily))[["ETo_daily"]],
    x_pos    = quantile(ETo_daily, 0.75),
    y_pos    = quantile(ETc_daily, 0.15),
    .groups  = "drop"
  ) %>%
  mutate(label = paste0("Kc = ", round(Kc_slope, 2)))

p_daily <- daily_etc %>%
  ggplot(aes(x = ETo_daily, y = ETc_daily,
             colour = Genotype, fill = Genotype)) +
  geom_point(alpha = 0.35, shape = 21, size = 1.8,
             stroke = 0.3, colour = "grey20") +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.1, alpha = 0.15) +
  geom_text(data = kc_daily,
            aes(x = x_pos, y = y_pos, label = label, colour = Genotype),
            size = 3.3, fontface = "italic", show.legend = FALSE) +
  facet_wrap(water_mgmt ~ Genotype, nrow = 2, labeller = label_wrap_gen(20)) +
  scale_colour_manual(values = c("Non-limited transpiration" = "#31a354",
                                 "Limited transpiration"      = "#d95f02")) +
  scale_fill_manual(values   = c("Non-limited transpiration" = "#31a354",
                                 "Limited transpiration"      = "#d95f02")) +
  labs(
    x        = expression(Mean~daily~ETo~(mm~d^{-1})),
    y        = expression(Mean~daily~crop~ET~(mm~d^{-1})),
    colour   = NULL, fill = NULL,
    title    = "Daily crop ET vs reference ET: LT limits water use at high evaporative demand",
    subtitle = "Following Marin et al. (2016, 2019): Kc slope from lm(ETc ~ ETo); each point = one simulation year"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position   = "bottom",
    strip.background  = element_rect(fill = "grey92", colour = NA),
    strip.text        = element_text(face = "bold", size = 9),
    panel.grid.minor  = element_blank(),
    plot.title        = element_text(face = "bold", size = 12),
    plot.subtitle     = element_text(colour = "grey40", size = 9)
  )

ggsave("modeling_data/figures/Fig_marin_daily_etc_eto.png",
       plot = p_daily, width = 24, height = 18, units = "cm", dpi = 150, bg = "white")


## ── 8. Growing-season ETo climatology by site ─────────────────────────────────

eto_clim <- eto_daily %>%
  filter(DOY >= 120, DOY <= 300) %>%
  group_by(site, Year, Month) %>%
  summarise(ETo_mm = sum(ETo, na.rm = TRUE), .groups = "drop") %>%
  group_by(site, Month) %>%
  summarise(ETo_mean = mean(ETo_mm),
            ETo_sd   = sd(ETo_mm),
            .groups  = "drop")

eto_clim %>%
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

ggsave("modeling_data/figures/ETo_climatology_growingseason.png",
       width = 20, height = 13, units = "cm",
       dpi = 150, bg = "white")


## ── 9. Export ──────────────────────────────────────────────────────────────────

write_csv(eto_daily,    "modeling_data/data/nasapower_ETo_daily.csv")
write_csv(eto_seasonal, "modeling_data/data/nasapower_ETo_seasonal.csv")
write_csv(merged,       "modeling_data/data/merged_ETo_cropET.csv")
