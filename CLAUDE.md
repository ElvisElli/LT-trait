# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo does

Agronomic analysis for **limited-transpiration (LT) soybean** trait evaluation, targeting a manuscript for *Agricultural Water Management*. Core argument: the LT trait produces a genotype-dependent crop coefficient (Kc) that challenges the FAO-56 assumption of a fixed, species-level Kc. The Marin et al. (2016, 2019) empirical finding — Kc decreases as ETo increases for coupled crops including soybean — is mechanistically explained by the LT trait, which genetically encodes this atmospheric-coupling-driven Kc reduction.

## Folder structure

```
LT-trait/
├── CLAUDE.md                   ← this file
├── MANUSCRIPT_DRAFT.md         ← first draft targeting Agric Water Management
├── RIG_Final_Report.md         ← RIG 2024-2026 final report
├── PLAN.md                     ← figure plan and analysis roadmap
├── scripts/                    ← all R analysis scripts
│   ├── nasapower_PM_ETo.R      ← main pipeline: download ETo + ETo vs cropET plots
│   ├── run_eto.R               ← headless ETo download (no RStudio needed)
│   ├── compute_pdoy_eto.R      ← planting-date-specific ETo windows
│   ├── Fig_sinclair_lt_trait.R ← Fig 1 — Sinclair leaf VPD response
│   ├── Fig_kc_scatter.R        ← Fig 3 — Kc scatter (irrigated, pdoy ETo)
│   ├── Fig_kc_summary.R        ← Fig 6 — Kc by state vs FAO-56 reference
│   ├── Fig_map_yields.R        ← Fig 5 — map: pie + facet_wrap versions
│   ├── Fig_water_productivity.R← Fig 4 — water productivity violin+box
│   ├── build_bridge_fig.R      ← Supplementary — Sinclair + Kc bridge
│   └── soybean_ssm_analysis.R  ← SSM plots: yield, ET, temperature
├── figures/                    ← all PNG figures (no TIFF)
│   ├── Fig_sinclair_lt_trait.png
│   ├── ETo_climatology_growingseason.png
│   ├── Fig_kc_scatter.png
│   ├── Fig_kc_seasonal_scatter.png  ← new (LT/Non-LT rows × State cols)
│   ├── Fig_marin_daily_etc_eto.png  ← new (Marin-style daily ETc vs ETo)
│   ├── Fig_water_productivity.png
│   ├── Fig_map_yields_pie.png       ← pie-chart map
│   ├── Fig_map_yields_facet.png     ← new: facet_grid(water_mgmt ~ Genotype)
│   ├── Fig_kc_summary.png
│   └── Fig_bridge_sinclair_vs_kc.png
├── data/                       ← input data (tracked in git except large files)
│   ├── SoybeanSSM_sims2.xlsx   ← Sinclair SSM output (2 160 rows)
│   ├── merged_ETo_cropET.csv   ← joined ETo + SSM data (primary analysis file)
│   ├── nasapower_ETo_pdoy.csv  ← planting-date-specific ETo by State×Year×Pdoy
│   └── nasapower_ETo_seasonal.csv ← annual growing-season ETo (DOY 120-300)
├── reports/
│   ├── SSM_analysis_report.Rmd ← RMarkdown source
│   └── SSM_analysis_report.html← rendered HTML report
└── experimental_data/          ← growth chamber experiments (see its own CLAUDE.md)
    ├── CLAUDE.md
    ├── GC_analysis_report.html
    ├── AUDIT_REPORT.md
    ├── references/             ← key PDFs (Marin 2016, 2019, DeMelo 2026, etc.)
    ├── scripts/
    ├── cleaned_data/
    ├── figures/                ← PNG only (no TIFF)
    └── stats/
```

## Running scripts

**Full pipeline (requires `data/SoybeanSSM_sims2.xlsx`):**
```r
setwd("/home/user/LT-trait")
source("scripts/nasapower_PM_ETo.R")   # ETo download + all ETo/Kc plots
source("scripts/Fig_map_yields.R")     # yield maps (pie + facet)
source("scripts/Fig_kc_scatter.R")     # Kc scatter
source("scripts/Fig_water_productivity.R")
source("scripts/Fig_kc_summary.R")
```

**HTML report:**
```r
setwd("/home/user/LT-trait")
rmarkdown::render("reports/SSM_analysis_report.Rmd",
                  output_file = "reports/SSM_analysis_report.html")
```

**Headless / cloud (ETo only, no RStudio):**
```bash
Rscript scripts/run_eto.R              # download NASA POWER + compute ETo
Rscript scripts/compute_pdoy_eto.R    # requires nasapower_ETo_daily.csv
```

## Key data file

`data/SoybeanSSM_sims2.xlsx` — Sinclair SSM simulation output (1990–2019, 6 locations, 2 160 rows).

| Column | Description |
|---|---|
| `Crop` | Genotype ID (`SOYBEAN-MG4`, `SOYBEAN-LT1.5-MG4`, etc.) |
| `Manag` | Management label — contains state code (AR / MO / NE) and water regime |
| `Location` | Site name (Jonesboro, Marianna, Columbia, Norborne, Lincoln, North Platte) |
| `Pyear` / `Year` | Simulation year |
| `Pdoy` | Planting day-of-year |
| `R7.PM` | Days from planting to physiological maturity — used for ETo window |
| `Planting_date` | "Early" / "Medium" / "Late" |
| `water_mgmt` | `"Rainfed"` or `"Irrigated"` |
| `SUMETT` | Seasonal total crop evapotranspiration (mm) |
| `SUMET2` / `SUMET3` | ET split by growth stage (sowing→R5, R5→harvest) |
| `Yield_kgha` | Grain yield (kg/ha) |

`data/merged_ETo_cropET.csv` — primary analysis file (SSM joined with NASA POWER ETo):

Additional columns: `Genotype` ("Non-limited transpiration" / "Limited transpiration"), `State` (AR/MO/NE), `site` (weather station), `ETo_season_mm` (DOY 120-300 seasonal ETo)

## ETo formula (FAO-56 Penman-Monteith)

- **es**: `0.6108 * exp(17.27 * tmed / (tmed + 237.3))` (°C, FAO-56 eq. 11)
- **Net radiation**: solar geometry (decl/hn/dD/Qo) + BOC + BOL
- **γ**: fixed at `0.062 kPa °C⁻¹`
- **Wind**: NASA POWER `WS2M` is already in m/s — no conversion needed
- **Negative ETo**: clamped to 0 via `pmax(0, ETo)`
- **Seasonal window**: DOY 120–300 for climatology; `Pdoy` to `Pdoy + R7.PM` (capped at DOY 330) for Kc

## Genotype classification

```r
Genotype = ifelse(Crop %in% c("SOYBEAN-MG4", "SOYBEAN-MG3"),
                  "Non-limited transpiration", "Limited transpiration")
```
**Colour convention:** Non-LT = `#31a354` (green), LT = `#d95f02` (orange).

## Marin et al. key insights

Marin et al. (2016, 2019) showed that **Kc decreases as ETo increases** for coupled crops (including soybean). The mechanism: high atmospheric demand induces stomatal closure → reduces transpiration below fixed-Kc prediction. The LT trait genetically encodes this response. Our Marin-style plot: ETc_daily = SUMETT/R7.PM vs ETo_daily = ETo_season_mm/181; regression slope = Kc. Expected: LT has lower slope than Non-LT.

De Melo et al. (2026): standard empirical Kc derivations overlook soil/plant hydraulics; process-based models better represent water dynamics.

## Key results (for manuscript)

| Metric | LT | Non-LT |
|---|---|---|
| Seasonal Kc (irrigated, pdoy-specific) | 1.27 ± 0.20 | 1.78 ± 0.29 |
| Water savings vs Non-LT (irrigated) | 246–344 mm season⁻¹ | — |
| Rainfed yield vs Non-LT | +3.6 to +14.7% | — |
| Water productivity (rainfed) | 6.28 kg mm⁻¹ | 5.63 kg mm⁻¹ |
| ANCOVA genotype × ETo interaction | p < 0.001 | — |
| FAO-56 reference Kc (mid-season) | 1.15 | 1.15 |
| Marin et al. (2019) whole-cycle soybean Kc | 0.83 ± 0.39 | (Brazilian field measurement) |

## Outputs tracked in git

| File | Notes |
|---|---|
| `data/merged_ETo_cropET.csv` | Primary analysis file |
| `data/nasapower_ETo_pdoy.csv` | Planting-date-specific ETo |
| `data/nasapower_ETo_seasonal.csv` | Annual growing-season ETo |
| `data/SoybeanSSM_sims2.xlsx` | SSM simulation output |
| `figures/*.png` | All figures (PNG only) |
| `reports/SSM_analysis_report.html` | HTML analysis report |
| `nasapower_ETo_daily.csv` | **NOT tracked** (~7 MB, regenerate with `scripts/run_eto.R`) |
