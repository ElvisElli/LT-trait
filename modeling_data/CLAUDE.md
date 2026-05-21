# modeling_data — CLAUDE.md

## Project context
Sinclair Simplified Soybean Model (SSM) simulations paired with NASA POWER
climate data to evaluate the limited-transpiration (LT) trait effect on crop
water use (ETc), crop coefficient (Kc), yield, and water productivity.
PI: Elvis F. Elli · RIG project 04/2024–02/2026

## Folder structure
```
modeling_data/
├── data/           — input and processed data files
├── figures/        — all PNG figures (no TIFF)
├── references/     — key papers (Chamarthi 2022, etc.)
├── reports/        — SSM_analysis_report.Rmd + rendered HTML
├── scripts/        — all R analysis scripts
└── stats/          — model output tables (future use)
```

## Key data files

`data/SoybeanSSM_sims2.xlsx` — Sinclair SSM output (1990–2019, 6 locations, 2 160 rows)

| Column | Description |
|---|---|
| `Crop` | Genotype ID (`SOYBEAN-MG4`, `SOYBEAN-LT1.5-MG4`, etc.) |
| `Location` | Site name (Jonesboro, Marianna, Columbia, Norborne, Lincoln, North Platte) |
| `Pyear` / `Year` | Simulation year |
| `Pdoy` | Planting day-of-year |
| `R7.PM` | Days from planting to physiological maturity |
| `Planting_date` | "Early" / "Medium" / "Late" |
| `water_mgmt` | `"Rainfed"` or `"Irrigated"` |
| `SUMETT` | Seasonal total crop evapotranspiration (mm) |
| `Yield_kgha` | Grain yield (kg/ha) |

`data/merged_ETo_cropET.csv` — primary analysis file (SSM joined with NASA POWER ETo)

`data/nasapower_ETo_pdoy.csv` — planting-date-specific ETo by State × Year × Pdoy

`data/nasapower_ETo_seasonal.csv` — annual growing-season ETo (DOY 120–300)

`data/nasapower_ETo_daily.csv` — **NOT tracked in git** (~7 MB); regenerate with `scripts/run_eto.R`

## Running the pipeline
```r
setwd("/home/user/LT-trait")
source("modeling_data/scripts/nasapower_PM_ETo.R")       # ETo download + ETc/Kc plots
source("modeling_data/scripts/Fig_map_yields.R")          # yield maps (pie + facet)
source("modeling_data/scripts/Fig_kc_scatter.R")          # Kc scatter
source("modeling_data/scripts/Fig_water_productivity.R")  # water productivity
source("modeling_data/scripts/Fig_kc_summary.R")          # Kc summary vs FAO-56
```

**HTML report:**
```r
setwd("/home/user/LT-trait")
rmarkdown::render("modeling_data/reports/SSM_analysis_report.Rmd",
                  output_file = "modeling_data/reports/SSM_analysis_report.html")
```

## Genotype classification
```r
Genotype = ifelse(Crop %in% c("SOYBEAN-MG4", "SOYBEAN-MG3"),
                  "Non-limited transpiration", "Limited transpiration")
```
Color convention: Non-LT = `#31a354` (green), LT = `#d95f02` (orange)

## Key results

| Metric | LT | Non-LT |
|---|---|---|
| Seasonal Kc (irrigated, pdoy-specific) | 1.27 ± 0.20 | 1.78 ± 0.29 |
| Water savings vs Non-LT (irrigated) | 246–344 mm season⁻¹ | — |
| Rainfed yield vs Non-LT | +3.6 to +14.7% | — |
| Water productivity (rainfed) | 6.28 kg mm⁻¹ | 5.63 kg mm⁻¹ |
| ANCOVA genotype × ETo interaction | p < 0.001 | — |
| FAO-56 reference Kc (mid-season) | 1.15 | 1.15 |
