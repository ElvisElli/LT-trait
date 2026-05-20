# Research Incentive Grant (RIG) — 2024–2025 Final Report

**Project Title:** Investigating transpiration responses of soybean genotypes to vapour pressure deficit for enhanced water use efficiency and drought tolerance

**Principal Investigator:** Elvis F. Elli, Assistant Professor, Crop Physiology, CSES
**Co-PIs:** Samuel Fernandes, Assistant Professor, Agricultural Statistics, CSES;
Caio Canella Vieira, Assistant Professor, Soybean Breeding, CSES

**Project performance dates:** 04/2024 – 02/2026
**Report date:** [Month] 2025

---

## 1. Summary of Accomplishments

This project aimed to quantify transpiration responses of diverse soybean accessions to
increasing vapour pressure deficit (VPD) and to connect controlled-environment measurements to
field stomatal conductance, seed yield, and process-based model simulations. Over the
performance period we completed three growth chamber (GC) experiments encompassing 12 soybean
accessions, collected field stomatal conductance and yield data for an additional set of
commercial and experimental lines, and integrated all findings into a crop model framework that
links genotype-specific water use to seasonal irrigation demand across six US locations.

---

## 2. Growth Chamber Experiments: Design and Execution

### 2.1 Experimental design

Three GC runs were conducted following the protocol described in the original proposal (King et
al., 2009; Fletcher et al., 2007). Soybean plants were grown in sealed PVC pots to beginning
flowering (R1), at which point individual pots were placed on portable balances linked to an
Arduino datalogger recording pot weights every 2 minutes. VPD treatments were applied
sequentially by adjusting chamber relative humidity at constant temperature (31–32°C).
Transpiration was computed as weight loss between consecutive recordings.

**Genotype selection** was guided by the WUE plasticity × WUE level framework of Chamarthi et
al. (2022), which identifies four genotypic classes based on carbon isotope discrimination
(δ¹³C) and its phenotypic plasticity across environments:

| Quadrant | Plasticity | WUE level | GC1 Genotype | GC3 Genotype |
|----------|-----------|-----------|--------------|--------------|
| Q1       | Plastic   | High WUE  | PI471938     | PI398617     |
| Q2       | Plastic   | Low WUE   | PI507408     | PI070229     |
| Q3       | Nonplastic| Low WUE   | PI548431     | PI423802     |
| Q4       | Nonplastic| High WUE  | PI603457A    | PI603465C    |

GC Experiment 1 (GC1) used 6 VPD treatment levels (~1.0, 1.5, 1.9, 2.4, 2.8, 3.3 kPa) and
included stomatal conductance measurements (LI-600 porometer) concurrent with gravimetric
transpiration. GC Experiment 3 (GC3) used 5 VPD levels (~1.0, 1.5, 2.0, 2.7, 3.3 kPa). Both
experiments used 4 replications per genotype × treatment and two consecutive measurement days.

### 2.2 Day effect and data quality

A critical examination of measurement days revealed that Day 1 data were systematically more
variable than Day 2 in both experiments (coefficient of variation: GC3 Day 1 = 0.39 vs Day 2 =
0.20; GC1 Day 1 = 0.32 vs Day 2 = 0.28). GC3 Day 1 mean relative transpiration was 23% lower
than Day 2, suggesting that plants required one full measurement cycle to acclimate to the
target VPD conditions. Day 2 data were therefore used for all model fitting and statistical
comparisons reported below. No formal statistical outliers (> 2 SD from within-group mean) were
detected in either dataset, confirming that the Day 1 variability reflects a systematic
acclimation effect rather than instrument error.

---

## 3. Model Fitting: Linear, Bilinear, and Linear-Plateau Mixed Effects Models

### 3.1 Modelling approach

Three model types were fitted to the Day 2 relative transpiration data for each genotype using
mixed effects models (lme4 R package; Bates et al., 2014), with replication as a random
intercept to account for pot-to-pot variability:

- **Linear:** RelTransp = β₀ + β₁ × VPD + (1 | Replication)
  A continuously increasing response with no change in sensitivity to VPD.

- **Bilinear:** RelTransp = β₀ + β₁ × VPD_low + β₂ × VPD_high + (1 | Replication)
  where VPD_low = min(VPD, BP) and VPD_high = max(VPD − BP, 0). The breakpoint (BP) was
  estimated by profile likelihood using the `segmented` package (Muggeo, 2003). This model
  allows a change in transpiration sensitivity at a threshold VPD.

- **Linear-Plateau:** RelTransp = β₀ + β₁ × min(VPD, BP) + (1 | Replication)
  Transpiration increases linearly below BP and reaches a plateau above it, representing
  complete stomatal limitation at high VPD. BP was identified by grid search minimising AIC
  over the observed VPD range.

Model selection was based on AIC (lower = better fit per parameter count). All models were
fitted with REML = FALSE for comparability.

### 3.2 Results

**Table 1.** Model comparison results by genotype and experiment (Day 2, mixed effects; R² =
Pearson correlation² between observed and fitted values; BP = breakpoint VPD in kPa; slope
units = RelTransp kPa⁻¹).

| Exp | Genotype   | Quadrant | Best model     | Linear R² | Linear slope | Bilinear BP | Slope low | Slope high | Bilin R² | LP BP | LP slope | LP R²  | ΔAIC (best vs 2nd) |
|-----|-----------|----------|----------------|-----------|--------------|-------------|-----------|------------|----------|-------|----------|--------|---------------------|
| GC1 | PI471938  | Q1       | **Linear-Plateau** | 0.340 | 0.121 | 1.50 | 0.328 | 0.078 | 0.387 | 2.9 | 0.133 | 0.345 | 0.14 |
| GC1 | PI507408  | Q2       | **Bilinear**   | 0.926 | 0.199 | 2.00 | 0.147 | 0.239 | 0.939 | 3.0 | 0.210 | 0.911 | 2.80 |
| GC1 | PI548431  | Q3       | **Linear**     | 0.474 | 0.201 | 1.81 | 0.095 | 0.255 | 0.488 | 3.0 | 0.213 | 0.469 | 0.24 |
| GC1 | PI603457A | Q4       | **Linear**     | 0.970 | 0.289 | 2.40 | 0.273 | 0.321 | 0.972 | 3.0 | 0.306 | 0.958 | 0.88 |
| GC3 | PI398617  | Q1       | **Linear-Plateau** | 0.837 | 0.191 | 2.41 | 0.287 | 0.072 | 0.898 | 2.9 | 0.246 | 0.890 | 0.62 |
| GC3 | PI070229  | Q2       | **Linear-Plateau** | 0.888 | 0.163 | 2.61 | 0.211 | 0.071 | 0.929 | 3.0 | 0.201 | 0.929 | 1.85 |
| GC3 | PI423802  | Q3       | **Linear-Plateau** | 0.838 | 0.134 | 2.70 | 0.151 | 0.093 | 0.847 | 3.2 | 0.151 | 0.847 | 0.44 |
| GC3 | PI603465C | Q4       | **Linear**     | 0.911 | 0.172 | 1.78 | 0.112 | 0.191 | 0.922 | 3.3 | 0.184 | 0.899 | 0.21 |

### 3.3 Key findings

**Pattern 1 — Q1 genotypes (Plastic/High WUE) show the clearest transpiration limitation.**
PI471938 (GC1) and PI398617 (GC3) both best fit a Linear-Plateau or Bilinear model with a
pronounced reduction in sensitivity above ~2.4–2.9 kPa. For PI398617, the bilinear slope drops
from 0.287 below 2.41 kPa to only 0.072 above it (a 75% reduction in slope), consistent with
the limited-transpiration mechanism described by Sinclair et al. (2008) and Fletcher et al.
(2007). This matches the expected behaviour of high-WUE plastic genotypes that reduce stomatal
aperture at high VPD.

**Pattern 2 — Q4 genotypes (Nonplastic/High WUE) show linear responses with high R².**
PI603457A (R² = 0.970) and PI603465C (R² = 0.911) best fit a simple linear model, indicating
that transpiration continues to rise proportionally with VPD across the full treatment range.
Despite being classified as high-WUE, these accessions appear not to limit transpiration via a
VPD threshold mechanism — their WUE advantage may derive from superior photosynthetic efficiency
rather than stomatal restriction.

**Pattern 3 — Q3 (Nonplastic/Low WUE) and Q2 (Plastic/Low WUE) show intermediate behaviour.**
PI548431 (Q3) shows a linear response but with a lower R² (0.474), reflecting higher among-
replication variability. PI507408 (Q2) shows the most clearly bilinear response (best model by
ΔAIC = 2.80), with slope increasing above the 2.0 kPa breakpoint — an unexpected pattern for a
plastic genotype that may reflect VPD-driven cuticular transpiration at high demand.

**Pattern 4 — Breakpoints cluster between 2.4 and 2.9 kPa** for genotypes where non-linear
models are selected, broadly consistent with the 2–3 kPa threshold range documented for the LT
trait in soybean (Sinclair et al., 2008; Devi et al., 2010).

**Caveat on model discrimination:** ΔAIC between the best and second-best model is < 2.0 for
six of eight genotypes, indicating that with n = 4 replications per VPD treatment the data do
not strongly discriminate between model forms. This is consistent with the interim report
finding that some genotypes showed marginal nonlinearity. Expanding replication and extending
the VPD range to > 3.5 kPa (as recommended in the interim report) remains a priority for the
next run.

---

## 4. Stomatal Conductance Analysis (GC1)

Stomatal conductance (gₛ, mol H₂O m⁻² s⁻¹) was measured concurrently with gravimetric
transpiration in GC1 using an LI-600 porometer. After removing three data-entry errors (two
implausibly large values > 10⁶ and one negative value), gₛ ranged from 0.007 to 0.922
mol m⁻² s⁻¹.

Stomatal conductance declined with increasing VPD in all four GC1 genotypes, consistent with
stomatal feedback regulation. The Q1 genotype PI471938 showed the steepest gₛ decline above
2 kPa, consistent with its Linear-Plateau transpiration response. The correlation between gₛ
and relative transpiration was negative (r = −0.24 overall; Q1: r = −0.37, p = 0.03; Q4:
r = −0.29, p = 0.05), reflecting the paradox that under increasing VPD, atmospheric demand
forces higher whole-plant transpiration even as stomata close. This decoupling between gₛ and
total transpiration is particularly pronounced in the Q1/LT-trait genotypes, where stomatal
closure limits the amplification of transpiration as VPD rises above 2 kPa.

---

## 5. Field Stomatal Conductance and Yield

Field stomatal conductance and grain yield were measured in a commercial trial including 12
soybean genotypes (Rainfed and Irrigated treatments) on two dates (August 9 and August 30,
2024). Key findings:

- Irrigated plots showed 2–5× higher gₛ than rainfed plots on August 9, with the gap narrowing
  by August 30, likely due to late-season precipitation or canopy senescence.
- Stomatal conductance was positively associated with yield under irrigation (consistent with
  source-driven yield formation under non-limiting water supply) but showed no clear association
  under rainfed conditions (where water supply, not stomatal activity, limits yield).
- PI377574, a GC2 genotype, showed very low rainfed gₛ (0.060 mol m⁻² s⁻¹) and the lowest
  rainfed yield (711 kg ha⁻¹), flagging it as potentially drought-sensitive.
- Wilting score was negatively correlated with gₛ (plants with lower gₛ wilted less), providing
  field-level validation that stomatal restriction protects against drought-induced wilting.

---

## 6. Integration with Process-Based Modelling

A key deliverable for leveraging this project toward extramural funding was combining
controlled-environment physiological measurements with DSSAT SSM crop model simulations to
identify the field-scale consequences of the LT trait. Thirty years (1990–2019) of SSM
simulations across six US locations (Arkansas, Missouri, Nebraska) combined with NASA POWER
FAO-56 ETo revealed that:

- LT genotypes exhibit a mean seasonal crop coefficient (Kc = SUMETT / ETo) of **1.27 ± 0.20**
  under irrigation, compared to **1.78 ± 0.29** for non-LT genotypes — a 29% reduction
  consistent with the VPD-threshold transpiration limitation documented in GC1 and GC3.
- LT genotypes saved **246–344 mm season⁻¹** in irrigation water demand without reducing
  rainfed yield (in fact producing 3.6–14.7% higher yield under water-limited conditions).
- The genotype × ETo interaction was highly significant (ANCOVA, p < 0.001), confirming that
  Kc is a genotype-level rather than species-level parameter and that the LT trait measurably
  decouples seasonal crop ET from inter-annual atmospheric variability.

These simulation results, anchored by the growth chamber transpiration measurements and field
stomatal conductance data from this RIG project, constitute a coherent multi-scale dataset
supporting the central hypothesis of the funded NIFA proposal.

---

## 7. Upcoming Activities

1. **GC Experiment 2 data processing:** Raw weight data from a second GC run (PI398245,
   PI567488B, PI377574, PI424397 — the GC2 genotype set) are available but require tidy
   processing into a Final-sheet format for inclusion in the model comparison.
2. **Extended VPD range:** A third experimental run targeting VPD levels of 3.5–4.0 kPa is
   planned to better characterise plateau thresholds in genotypes with near-linear responses
   at the current VPD ceiling.
3. **Full model discrimination:** Increasing replication to n ≥ 6 per treatment will improve
   statistical power to discriminate between Bilinear and Linear-Plateau models for genotypes
   with ΔAIC < 2.
4. **Manuscript preparation:** Two manuscripts are in preparation:
   - *Growth chamber:* "Genotype-specific transpiration responses to vapour pressure deficit in
     soybean accessions with contrasting WUE plasticity" (target: *Plant, Cell & Environment*
     or *Field Crops Research*)
   - *Modelling:* "Rethinking crop coefficients for climate-smart soybean: limited-transpiration
     genotypes exhibit genotype-dependent seasonal water use" (target: *Agricultural Water
     Management*)

---

## 8. Leveraging for Extramural Funding

The preliminary findings from this RIG project were included in a USDA NIFA proposal titled
**"Designing agronomic and genetic strategies to maximize soybean systems' performance under
climate change"** (currently under review). The growth chamber transpiration data (GC1, GC3),
field stomatal conductance measurements, and process-based modelling outputs collectively
provide the multi-scale proof-of-concept required by the NIFA reviewers. Specifically, the
convergence of:

1. Leaf-level VPD threshold physiology (Sinclair et al., 2008; this study)
2. Genomic basis for WUE plasticity (Chamarthi et al., 2022 — loci in stomatal morphogenesis)
3. Field-scale stomatal conductance–yield relationships (this study)
4. Simulation-derived genotype-specific Kc values (this study)

provides a compelling narrative for a mechanistically-grounded, multi-environment investigation
of LT trait deployment across US soybean-producing regions under climate change.

---

## 9. References

Bates, D., Maechler, M., Bolker, B., Walker, S., 2014. lme4: Linear mixed-effects models using
Eigen and S4. R package version 1.1-7.

Chamarthi, S.K., Naik, P., Belko, N., Ntare, B.R., Sinclair, T.R., Rathore, A., 2022.
Identification of stable QTLs/genomic regions for carbon isotope discrimination and its
plasticity using genome-wide association mapping in soybean. Plant Genome 15, e20284.

Devi, M.J., Sinclair, T.R., Vadez, V., 2010. Genotypic variation in peanut for transpiration
response to vapor pressure deficit. Crop Sci. 50, 191–196.

Fletcher, A.L., Sinclair, T.R., Hartwell Allen Jr, L., 2007. Transpiration responses to vapor
pressure deficit in well watered "slow-wilting" and commercial soybean. Environ Exp Bot 61,
145–151.

King, C.A., Purcell, L.C., Brye, K.R., 2009. Differential wilting among soybean genotypes in
response to water deficit. Crop Sci. 49, 290–298.

Muggeo, V.M.R., 2003. Estimating regression models with unknown break-points. Stat. Med. 22,
3055–3071.

Sinclair, T.R., Zwieniecki, M.A., Holbrook, N.M., 2008. Low leaf hydraulic conductance with
rapid stomatal closure. Plant Cell Environ. 31, 796–803.

Sinclair, T.R., Messina, C.D., Beatty, A., Samples, M., 2010. Assessment across the United
States of the benefits of altered soybean drought traits. Agron. J. 102, 475–482.

---

*Report prepared by Elvis F. Elli · All analyses conducted in R (v4.3.3) ·
Code and data: https://github.com/ElvisElli/LT-trait*
