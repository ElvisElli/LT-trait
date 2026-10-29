# Handoff: ASA 2026 presentation and 2026 analyses

Last updated: 2026-10-29. Read this first when picking the work back up.

## Where things stand

**Deck:** `ASA-2026-Elli-LT-trait.pptx` (repo root, 21 slides), built on the 2025
ASA deck as a template. Rebuild it with `python3 presentation/build_deck.py`; all
slide text and layout live in that script (`build_slides()`), static images in
`presentation/img/`, the 2025 template in `presentation/template/`.

| Slide | Content | Data basis |
|---|---|---|
| 1-5 | Title, drought adaptation, LT trait (kept from 2025) | - |
| 6-9 | Hypothesis, objectives, methods | - |
| 10 | TR vs VPD, 9 accessions, 3 x 3 | `experimental_data_soy_2026/script/08_fig_TR_presentation.R` |
| 11-14 | Breakpoints, quadrants, gsw, candidate table | 14 accessions, 5 reps (paper figures); **not yet aligned with the 9 on slide 10** |
| 16 | Yield-gap map, 10 SSM sites | `modeling_data/scripts/Fig_map_yields_10loc.R` |
| 17 | Yield-potential trade-off bar | same script, panel B |
| 18-21 | Conclusions, next steps, team, thank you | - |

**Figures made but not in the deck yet** (PNG + 600 dpi TIFF):

- `experimental_data_soy_2026/figures/presentation/Fig_TR_bars_1.5_3.5kPa.png`: predicted TR at 1.5 and 3.5 kPa with 95% CI (palette options `_BrBG`, `_PuOr`, `_RdBu`, `_Greens`, `_BluesTeal`; BrBG is the default, user has not picked yet)
- `modeling_data/figures/Fig_soil_water_season_state_FTSW.png` / `_ATSW.png`: rainfed soil water by state on cumulative thermal time, R1/R5/R7 arrows
- `modeling_data/figures/Fig_cropET_season_state_daily.png` / `_cumulative.png`: rainfed crop ET by state
- `modeling_data/figures/Fig_cropET_difference_location.png` / `Fig_cropET_cumulative_location.png`: rainfed crop ET per location, paired LT minus Non-LT
- `modeling_data/figures/Fig_map_study_sites.png`: M&M map, **layout only, acreage shading missing** (see blocked items)
- `experimental_data_soy_2026/figures/paper_figures/Supplementary_Tables_model_selection.docx`: Tables S1-S2 (candidate models, AIC, fit type)

## Open items and decisions waiting on the user

1. **Missouri and Nebraska simulations mostly did not mature.** In the 10-site run,
   471 of 540 MO and 450 of 540 NE rainfed seasons ended with `MATYP = 4` (stop day
   reached), with "R7" at 155-185 days. This affects every 10-site result (map,
   slide 17, soil water, ET). Options: rerun SSM for MO/NE with a later stop day or
   earlier MG, or drop `MATYP != 1` seasons. Not changed yet.
2. **Thermal-time axis.** Soil water / ET plots use the SSM daily temperature
   units (DTU); they sit 4-6 deg C below the daily mean, so R1 comes at about
   980 deg C d at Rohwer versus about 790 with conventional base-10 GDD. User asked about
   this; offered to switch all thermal-time plots to base-10 GDD from `TMP`.
3. **Slides 11-14** still use the 14-accession numbers; slide 10 shows 9. Propagate
   the 9-accession basis if the user wants consistency (breakpoints then 1.7-3.4 kPa,
   4 of 9 nonlinear).
4. **Bar-plot palette** for `Fig_TR_bars_1.5_3.5kPa` (default BrBG).
5. **Why LT has a longer rainfed cycle:** SSM accelerates development under water
   deficit (`WSFD`/`WSSD` in `modeling_data/data/updated-outputs/readme.xlsx`); LT is
   less stressed. Identical phenology under irrigation. The SSM crop parameter file
   is needed to report the actual threshold/coefficient.

## Blocked: needs environment settings (user action)

- Other weather sources blocked as of 2026-10-29: `daymet.ornl.gov`,
  `thredds.northwestknowledge.net` (gridMET), `services.nacse.org` (PRISM),
  `www.ncei.noaa.gov`, `mesonet.agron.iastate.edu`. `power.larc.nasa.gov` works.

## Data notes

- **VPD heatmap** (slides 5 and 8, `presentation/img/vpd_heatmap.png`) is now reproducible:
  `Rscript modeling_data/scripts/Fig_vpd_hours_heatmap.R`. NASA POWER hourly T2M/RH2M at the
  10 SSM sites, hours/day with VPD > 2 kPa averaged over 2016-2025 per calendar day (`YEARS`
  in the script). The hourly cache (`nasapower_hourly_T2M_RH2M_10loc.csv`, 32 MB) is gitignored
  and re-downloads in about 40 s. The old image had no script, and its legend had a 2-3 h gap ("3 to 4 h").

- **Study-sites map complete** (`modeling_data/figures/Fig_map_study_sites.png`; TIFF is gitignored).
  County acreage: `modeling_data/data/nass_soybean_acres_harvested_county.csv`, USDA NASS
  Quick Stats API (SURVEY, "SOYBEANS - ACRES HARVESTED", TOTAL, COUNTY level; AR, KS, MO,
  NE, OK; 2020-2024; downloaded 2026-10-29). The map shows the 2020-2024 mean for 333 counties.
  "OTHER COUNTIES" rows (no County ANSI) are kept in the CSV and filtered out by the script.
  `NASS_API_KEY` is set in the environment, and `quickstats.nass.usda.gov` is now reachable.

- Growth chamber: the analysis uses `experimental_data_soy_2026/data/data_all_9_10_2026.xlsx`
  (5 replications). Run QC with `GC_INPUT=data_all_9_10_2026.xlsx Rscript script/01_data_inspection_QC.R`
  from `experimental_data_soy_2026/`.
- One impossible gsw value (A5959, rep 5, treatment 2, 29.3 mol m-2 s-1) is set to
  NA in the QC script.
- Presentation figures drop PI398272, PI404191, PI442012A (2 reps), AG43XF2,
  PI377574 (high scatter) and rep 5 of PI424231 (offset). Paper figures keep all 14.
- PI416937: bilinear `nls` never converges in R, so it is reported as linear.
- `modeling_data/data/updated-outputs/North_Platte_output.xlsx` was renamed (a
  soft hyphen in the old file name broke R).

## Setting up R in a fresh cloud session

R is not preinstalled. What worked:

```bash
apt-get update -qq && apt-get install -y -qq r-base r-base-dev r-cran-ggplot2 r-cran-dplyr \
  r-cran-tidyr r-cran-purrr r-cran-forcats r-cran-tibble r-cran-readxl r-cran-readr \
  r-cran-ggrepel r-cran-maps r-cran-ggforce r-cran-patchwork r-cran-scales r-cran-stringr \
  r-cran-jsonlite r-cran-base64enc r-cran-ggbeeswarm r-cran-officer r-cran-flextable
# ggh4x has no apt/CRAN build for R 4.3 here; install the 0.2.8 tag from GitHub
git clone -q --depth 1 --branch v0.2.8 https://github.com/teunbrand/ggh4x /tmp/ggh4x && R CMD INSTALL /tmp/ggh4x
```

Run scripts that write the degree sign or other Unicode with `LC_ALL=C.UTF-8`.
