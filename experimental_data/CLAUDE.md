# experimental_data — CLAUDE.md

## Project context
Growth chamber (GC) transpiration experiments testing VPD responses in soybean
accessions classified by WUE plasticity (Chamarthi et al., 2022 framework).
PI: Elvis F. Elli · RIG project 04/2024–02/2026

## Folder structure
raw_data/GC1/   — GC Experiment 1: PI471938, PI507408, PI548431, PI603457A (6 VPD levels, has SC data)
raw_data/GC2/   — GC Experiment 2: PI398245, PI567488B, PI377574, PI424397 (6 VPD levels, Final sheet in Data analysis 2.xlsx)
raw_data/GC3/   — GC Experiment 3: PI398617, PI070229, PI423802, PI603465C (5 VPD levels)
raw_data/LI600/ — LI-600 porometer stomatal conductance (concurrent with GC1)
raw_data/field_SC/ — Field SC + yield + wilting score, Rainfed vs Irrigated, Aug 9 & Aug 30 2024
cleaned_data/   — Tidy CSVs ready for analysis (GC_combined_clean.csv is the primary analysis file)
scripts/        — R analysis scripts (run in order: GC_transpiration_analysis → GC_critical_analysis → GC_model_fitting)
figures/        — All PNG (150 dpi) and TIFF (300 dpi) publication figures
stats/          — Model comparison and statistics CSVs
references/     — Key papers (Chamarthi 2022, Fletcher 2007, King 2009, etc.)
reports/        — RIG proposal, interim report, final report

## Key data file
cleaned_data/GC_combined_clean.csv — 512 rows (GC1 n=180, GC2 n=192, GC3 n=140), Day 2 recommended for analysis

Columns: Experiment, Genotype, Quadrant, Expected_WUE, Plasticity, WUE_class,
         Replication, Pot, Day, Treatment, VPD_target, Temperature, RH, VPD,
         RelTranspiration, Transpiration_ml, TranspLeafArea, DryMatter_total,
         LeafArea_squarecm, StomatalConductance (GC1 only, NA for GC2/GC3)

## Genotype classification (Chamarthi et al. 2022 WUE quadrants)
Q1 — Plastic, High WUE:    PI471938 (GC1), PI398245 (GC2), PI398617 (GC3)  [expect LT-like response]
Q2 — Plastic, Low WUE:     PI507408 (GC1), PI567488B (GC2), PI070229 (GC3)
Q3 — Nonplastic, Low WUE:  PI548431 (GC1), PI377574 (GC2), PI423802 (GC3)
Q4 — Nonplastic, High WUE: PI603457A (GC1), PI424397 (GC2), PI603465C (GC3) [expect linear response]

## VPD treatments
GC1: T1=1.0, T2=1.5, T3=1.9, T4=2.4, T5=2.8, T6=3.3 kPa (6 levels)
GC2: T1=0.95, T2=1.52, T3=2.09, T4=2.66, T5=3.23, T6=3.80 kPa target (6 levels; T6 actual ~3.36–3.61 kPa)
GC3: T1=1.0, T2=1.5, T3=2.0, T4=2.7, T5=3.3 kPa (5 levels)

## Day effect — IMPORTANT
Day 1 shows systematic acclimation artifact (CV 0.32–0.39 vs 0.20–0.28 on Day 2 in GC1/GC3;
GC3 Day 1 mean RelTransp 23% lower than Day 2). GC2 showed a much smaller day effect
(Day 1 CV=0.229 vs Day 2 CV=0.219; mean RelTransp nearly identical), suggesting better
plant acclimation in that run. Always use Day == 2 for modelling for consistency.

## Running the analysis pipeline (Python, cloud-compatible)
```bash
python3 experimental_data/scripts/GC_figures.py       # generates all 13 PNG figures
python3 experimental_data/scripts/GC_build_report.py  # builds self-contained HTML report
```

GC2 data extraction (already done; GC2_clean.csv committed):
```bash
python3 experimental_data/scripts/GC2_process_raw.py   # re-extract from Data analysis 2.xlsx if needed
```

Original R scripts (require R + lme4 + segmented; not available in cloud):
```r
setwd("/home/user/LT-trait")
source("experimental_data/scripts/GC_transpiration_analysis.R")
source("experimental_data/scripts/GC_critical_analysis.R")
source("experimental_data/scripts/GC_model_fitting.R")
```

## Model fitting (GC_model_fitting.R)
Three mixed effects models per genotype (lme4, REML=FALSE, replication as random intercept):
- Linear: RelTransp ~ VPD + (1|Replication)
- Bilinear: breakpoint estimated via segmented package, profile likelihood
- Linear-Plateau: breakpoint by AIC grid search over VPD range
Selection by AIC. Results in stats/GC_model_comparison.csv.
Note: GC2 formal mixed-effects model fitting pending (requires R).

## Key results (Day 2)
GC1/GC3 mixed effects (R):
  Q1 genotypes (PI471938, PI398617): Linear-Plateau best fit, BP ~2.4–2.9 kPa — LT confirmed
  Q4 genotypes (PI603457A, PI603465C): Linear best fit, R²>0.91 — no VPD threshold
  ΔAIC < 2 for 6/8 genotypes → more replication needed to discriminate model forms

GC2 preliminary (Python linear regression, TranspLeafArea vs VPD):
  All genotypes peak at T4-T5 (VPD 2.67–3.12 kPa) then decline at T6 (3.60 kPa)
  PI398245 Q1: slope=5.19 mg m⁻² s⁻¹ kPa⁻¹, R²=0.12 — plateau-like
  PI567488B Q2: slope=8.42, R²=0.62 — stronger linear component
  PI377574 Q3: slope=5.27, R²=0.18
  PI424397 Q4: slope=6.02, R²=0.43
  No outliers detected (>2.5 SD criterion)

## Upcoming
- Run formal mixed-effects model fitting for GC2 in R (lme4/segmented)
- Extend VPD to 3.5–4.0 kPa in next run
- Increase replication to n ≥ 6 per treatment
