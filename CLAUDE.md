# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo does

Agronomic analysis for **limited-transpiration (LT) soybean** trait evaluation, targeting a manuscript for *Agricultural Water Management*. Core argument: the LT trait produces a genotype-dependent crop coefficient (Kc) that challenges the FAO-56 assumption of a fixed, species-level Kc.

Two main pipelines:

1. **`nasapower_PM_ETo.R`** — Downloads daily NASA POWER climate data for three US sites (Fayetteville AR, Columbia MO, Lincoln NE), computes FAO-56 Penman-Monteith ETo, then merges with DSSAT SSM simulation output to compare seasonal crop water demand (SUMETT) against atmospheric evaporative demand.

2. **`soyeban-ssm-11-9.R`** — Reads DSSAT SSM simulation output (`SoybeanSSM_sims2.xlsx`) and produces yield, ET, and transpiration figures comparing LT and non-LT soybean under rainfed and irrigated conditions across AR, MO, NE.

## Running scripts

**Full pipeline in RStudio (requires `SoybeanSSM_sims2.xlsx`):**
```r
# Source nasapower_PM_ETo.R — setwd() handled via rstudioapi
```

**Headless / cloud (ETo + all figures, no RStudio needed):**
```bash
Rscript run_eto.R              # download NASA POWER + compute ETo + all figures
Rscript compute_pdoy_eto.R     # planting-date-specific ETo (requires nasapower_ETo_daily.csv)
Rscript Fig_sinclair_lt_trait.R
Rscript Fig_kc_scatter.R
Rscript Fig_water_productivity.R
Rscript Fig_map_yields.R
Rscript Fig_kc_summary.R
Rscript build_bridge_fig.R
```

## Key data file

`SoybeanSSM_sims2.xlsx` — DSSAT SSM simulation output (1990–2019, 6 locations, 2 160 rows). Required columns:

| Column | Description |
|---|---|
| `Crop` | Genotype ID (`SOYBEAN-MG4`, `SOYBEAN-LT1.5-MG4`, etc.) |
| `Manag` | Management label — contains state code (AR / MO / NE) and water regime |
| `Location` | Site name (Jonesboro, Marianna, Columbia, Norborne, Lincoln, North Platte) |
| `Pyear` | Simulation year |
| `Pdoy` | Planting day-of-year |
| `R7.PM` | Days from planting to physiological maturity — used for ETo window |
| `Planting_date` | "Early" / "Medium" / "Late" |
| `water_mgmt` | `"Rainfed"` or `"Irrigated"` |
| `SUMETT` | Seasonal total crop evapotranspiration (mm) |
| `SUMET2` / `SUMET3` | ET split by growth stage (sowing→R5, R5→harvest) |
| `Yield_kgha` | Grain yield (kg/ha) |

## ETo formula

FAO-56 Penman-Monteith with these specifics:
- **es**: `0.6108 * exp(17.27 * tmed / (tmed + 237.3))` (Celsius, FAO-56 eq. 11)
- **Net radiation**: computed from solar geometry (`decl/hn/dD/Qo`) + `BOC + BOL`
- **γ**: fixed at `0.062 kPa °C⁻¹`
- **Wind**: NASA POWER `WS2M` is already in m/s — no conversion
- **Negative ETo**: clamped to 0 via `pmax(0, ETo)`
- **Seasonal window**: DOY 120–300 for climatology; `Pdoy` to `Pdoy + R7.PM` (capped at DOY 330) for Kc computation

## Genotype classification

```r
Genotype = ifelse(Crop %in% c("SOYBEAN-MG4", "SOYBEAN-MG3"),
                  "Non-limited transpiration", "Limited transpiration")
```
Colour convention: non-LT = `#31a354` (green), LT = `#d95f02` (orange).

## Key results (for manuscript)

| Metric | LT | Non-LT |
|---|---|---|
| Seasonal Kc (irrigated, pdoy-specific) | 1.27 ± 0.20 | 1.78 ± 0.29 |
| Water savings vs Non-LT (irrigated) | 246–344 mm season⁻¹ | — |
| Rainfed yield vs Non-LT | +3.6 to +14.7% | — |
| Water productivity (rainfed) | 6.28 kg mm⁻¹ | 5.63 kg mm⁻¹ |
| ANCOVA genotype × ETo interaction | p < 0.001 | — |
| FAO-56 reference Kc (mid-season) | 1.15 | 1.15 |

## Figure inventory

| File | Figure | Status |
|---|---|---|
| `Fig_sinclair_lt_trait.R/png/tiff` | Fig 1 — Sinclair leaf transpiration vs VPD | Done |
| `ETo_climatology_growingseason.R/png/tiff` | Fig 2 — Seasonal ETo climatology by site | Done |
| `Fig_kc_scatter.R/png/tiff` | Fig 3 — Kc scatter irrigated-only (pdoy ETo) | Done |
| `Fig_water_productivity.R/png/tiff` | Fig 4 — Water productivity violin+boxplot | Done |
| `Fig_map_yields.R/png/tiff` | Fig 5 — Map with yield pie charts, 6 locations | Done |
| `Fig_kc_summary.R/png/tiff` | Fig 6 — Kc by state vs FAO-56 reference | Done |
| `build_bridge_fig.R/png/tiff` | Supplementary — Sinclair + Kc bridge panel | Done |

## Outputs

| File | Tracked in git | Notes |
|---|---|---|
| `nasapower_ETo_seasonal.csv` | Yes | Annual growing-season ETo (DOY 120–300) by site |
| `nasapower_ETo_pdoy.csv` | Yes | Planting-date-specific ETo by State × Year × Planting_date |
| `nasapower_ETo_daily.csv` | No (.gitignore) | ~7 MB — regenerate with `run_eto.R` |
| `merged_ETo_cropET.csv` | Yes | Joined ETo (DOY 120–300) + SSM data |
| `MANUSCRIPT_DRAFT.md` | Yes | Full first draft targeting Agric Water Management |
| `PLAN.md` | Yes | Figure plan, analysis TODOs, and roadmap |
| All `Fig_*.png/tiff` | Yes | Publication figures (PNG for GitHub, TIFF for journal) |
