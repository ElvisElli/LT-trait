# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo does

Agronomic analysis for **limited-transpiration (LT) soybean** trait evaluation. Two main pipelines:

1. **`nasapower_PM_ETo.R`** — Downloads daily NASA POWER climate data for three US sites (Fayetteville AR, Columbia MO, Lincoln NE), computes FAO-56 Penman-Monteith ETo, then merges with DSSAT SSM simulation output to compare seasonal crop water demand (SUMETT) against atmospheric evaporative demand — i.e. the growing-season crop coefficient (Kc) for LT vs. non-LT genotypes.

2. **`soyeban-ssm-11-9.R`** — Reads DSSAT SSM simulation output (`SoybeanSSM_sims2.xlsx`) and produces yield, ET, and transpiration figures comparing LT and non-LT soybean under rainfed and irrigated conditions across AR, MO, NE.

## Running scripts

**In RStudio (full pipeline including merged Kc plot):**
```r
# Open nasapower_PM_ETo.R or soyeban-ssm-11-9.R and source it.
# setwd() is handled automatically via rstudioapi.
# SoybeanSSM_sims2.xlsx must be in the same directory.
```

**Headless / cloud environment (ETo only, no Excel dependency):**
```bash
Rscript run_eto.R
```

## Key data file

`SoybeanSSM_sims2.xlsx` — DSSAT SSM simulation output. Required columns used across both scripts:

| Column | Description |
|---|---|
| `Crop` | Genotype ID (e.g. `SOYBEAN-MG4`, `SOYBEAN-LT1.5-MG4`) |
| `Manag` | Management label — contains state code (AR / MO / NE) |
| `Pyear` | Simulation year |
| `water_mgmt` | `"Rainfed"` or `"Irrigated"` |
| `SUMETT` | Seasonal total crop evapotranspiration (mm) |
| `SUMET2` / `SUMET3` | ET split by growth stage (sowing→R5, R5→harvest) |
| `Yield_kgha` | Grain yield (kg/ha) |
| `SRAINT` | Seasonal rainfall (mm) |
| `MTMAXT` / `MTMINT` | Mean seasonal max/min temperature (°C) |

## ETo formula

FAO-56 Penman-Monteith with these specifics:
- **es**: `0.6108 * exp(17.27 * tmed / (tmed + 237.3))` (Celsius input, FAO-56 eq. 11)
- **Net radiation**: computed from solar geometry (`decl/hn/dD/Qo`) + `BOC + BOL` (not from a separate longwave download)
- **γ**: fixed at `0.062 kPa °C⁻¹`
- **Wind**: NASA POWER `WS2M` is already in m/s — no unit conversion
- **Negative ETo**: clamped to 0 via `pmax(0, ETo)` (can occur in winter when Rn < 0)
- Growing season window: **DOY 120–300** (≈ May–October)

## Genotype classification

Across both scripts, crops are collapsed to two levels:
```r
Genotype = ifelse(Crop %in% c("SOYBEAN-MG4", "SOYBEAN-MG3"),
                  "Non-limited transpiration", "Limited transpiration")
```
Colour convention: non-LT = `#31a354` (green), LT = `#d95f02` (orange).

## Outputs

| File | Tracked in git | Notes |
|---|---|---|
| `nasapower_ETo_seasonal.csv` | Yes | Annual growing-season ETo by site |
| `nasapower_ETo_daily.csv` | No (.gitignore) | ~7 MB, regenerate with `run_eto.R` |
| `merged_ETo_cropET.csv` | Yes (after SSM run) | Joined ETo + SSM data |
| `ETo_climatology_growingseason.tiff` | Yes | Monthly ETo climatology plot |
| `ETo_vs_cropET_LT_nonLT.tiff` | Yes (after SSM run) | Main Kc scatter plot |
