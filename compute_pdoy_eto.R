## Compute planting-date-specific seasonal ETo using actual Pdoy + R7.PM windows
## Output: nasapower_ETo_pdoy.csv — one row per State x Year x Planting_date

setwd("/home/user/LT-trait")
library(dplyr); library(readxl)

## ── Daily ETo (3 representative sites, one per state) ────────────────────────
eto_daily <- read.csv("nasapower_ETo_daily.csv") %>%
  mutate(Date = as.Date(YYYYMMDD))

## State abbreviation from site name
eto_daily <- eto_daily %>%
  mutate(State = case_when(
    grepl("AR", site) ~ "AR",
    grepl("MO", site) ~ "MO",
    grepl("NE", site) ~ "NE"
  ))

## ── Growing-season windows from SSM (Pdoy + R7.PM, capped at DOY 330) ────────
ssm <- read_excel("SoybeanSSM_sims2.xlsx") %>%
  mutate(State = case_when(
    grepl("AR", Manag) ~ "AR",
    grepl("MO", Manag) ~ "MO",
    grepl("NE", Manag) ~ "NE"
  )) %>%
  rename(Year = Pyear) %>%
  mutate(harvest_doy = pmin(Pdoy + R7.PM, 330))

## Mean Pdoy and harvest_doy by State x Year x Planting_date
windows <- ssm %>%
  group_by(State, Year, Planting_date) %>%
  summarise(
    pdoy_start  = round(mean(Pdoy,        na.rm = TRUE)),
    pdoy_end    = round(mean(harvest_doy, na.rm = TRUE)),
    .groups     = "drop"
  )

## ── Sum ETo over each window ──────────────────────────────────────────────────
## For each State x Year x Planting_date, filter daily ETo and sum
compute_window_eto <- function(state, year, doy_start, doy_end) {
  eto_daily %>%
    filter(State == state, YEAR == year,
           DOY >= doy_start, DOY <= doy_end) %>%
    summarise(ETo_pdoy_mm = sum(ETo, na.rm = TRUE)) %>%
    pull(ETo_pdoy_mm)
}

message("Computing planting-date-specific ETo windows...")
pdoy_eto <- windows %>%
  rowwise() %>%
  mutate(ETo_pdoy_mm = compute_window_eto(State, Year,
                                           pdoy_start, pdoy_end)) %>%
  ungroup()

write.csv(pdoy_eto, "nasapower_ETo_pdoy.csv", row.names = FALSE)
message("Saved: nasapower_ETo_pdoy.csv (", nrow(pdoy_eto), " rows)")

## Quick summary
cat("\n=== Mean ETo by State x Planting_date (planting-specific window) ===\n")
pdoy_eto %>%
  group_by(State, Planting_date) %>%
  summarise(ETo_mean = round(mean(ETo_pdoy_mm), 1),
            ETo_sd   = round(sd(ETo_pdoy_mm),   1),
            .groups  = "drop") %>%
  print()
