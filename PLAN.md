# Project Plan — Agricultural Water Management Paper

## Target journal
**Agricultural Water Management** (Elsevier) — IF ~6.7
Direct fit: Kc methodology, irrigation scheduling, water-limited vs potential yield.
Word limit: ~8 000 words main text. Max 6 figures.

---

## Paper title (working)
"Rethinking crop coefficients for climate-smart soybean: limited-transpiration
genotypes exhibit genotype-dependent water use across the US Mid-South and Midwest"

---

## Core argument
FAO-56 assigns a fixed seasonal Kc per crop, independent of genotype. Using 30 years
of Sinclair SSM simulations at 6 locations and NASA POWER ETo, we show that:
1. LT genotypes have Kc = 1.06 vs Non-LT Kc = 1.48 under potential (irrigated) conditions — a 28% difference
2. LT genotypes save 246–344 mm of irrigation water per season
3. Despite lower water use, LT genotypes yield 3.6–14.7% MORE under rainfed conditions
4. LT water productivity (6.28 kg mm⁻¹) exceeds Non-LT under rainfed (5.63 kg mm⁻¹)
5. Genotype × ETo interaction is highly significant (ANCOVA p < 0.001)

---

## Figure plan (6 figures max)

| Fig | Content | Script | Status |
|-----|---------|--------|--------|
| 1 | Sinclair et al. leaf transpiration vs VPD — 4 genotypes | Fig_sinclair_lt_trait.R | Done |
| 2 | Seasonal ETo climatology — 3 sites 1990–2023 | run_eto.R | Done |
| 3 | ETo vs SUMETT scatter — LT vs Non-LT, irrigated only | Fig_kc_scatter.R | TODO |
| 4 | Water productivity (yield/ET) — boxplot by genotype x site x water_mgmt | Fig_water_productivity.R | TODO |
| 5 | Map — water-limited vs potential yield at 6 locations | Fig_map_yields.R | Done (needs polish) |
| 6 | Bridge/synthesis — Kc comparison LT vs Non-LT vs FAO-56 standard | Fig_kc_summary.R | TODO |

---

## Analysis TODO

- [ ] Fig 3: Kc scatter irrigated-only — clean ETo vs SUMETT by genotype, annotate Kc ratio
- [ ] Fig 4: Water productivity boxplots — yield/SUMETT by genotype, site, water_mgmt
- [ ] Fig 6: Bar/point plot comparing Kc values: LT (this study), Non-LT (this study), FAO-56 soybean standard (~1.15), Sinclair genotypes (transpiration slope as proxy)
- [ ] Stats: ANCOVA results table (slope, SE, p-values) for manuscript Table 1
- [ ] Stats: Compute year-to-year CV of Kc by genotype — does LT show lower variability?
- [ ] Map: polish Fig 5 — add yield numbers inside pies, improve legend
- [ ] Update CLAUDE.md with new scripts and outputs

## Code improvements needed
- [ ] Separate ETo computation by planting date × location for a tighter Kc window
- [ ] Add run_eto.R to .gitignore outputs section (only seasonal CSV needed)
- [ ] Consolidate all figure scripts into a single `make_figures.R` master runner
