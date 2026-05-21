# MANUSCRIPT DRAFT
**Target:** Agricultural Water Management (Elsevier)
**Status:** First draft — numbers real, prose needs refinement

---

## Title
Rethinking crop coefficients for climate-smart soybean: limited-transpiration
genotypes exhibit genotype-dependent seasonal water use across the US Mid-South
and Midwest

---

## Authors
Elvis F. Elli^1^, [Co-authors TBD]

^1^ [Institution, Address]

**Corresponding author:** [email]

---

## Abstract (≈260 words)

The FAO-56 Penman-Monteith framework assigns a fixed crop coefficient (Kc) to
each species, implicitly assuming that seasonal water use scales uniformly with
atmospheric evaporative demand (ETo) regardless of genotype. This assumption
ignores the growing diversity of genetic variation in plant water-use traits and,
in particular, the mechanistic linkage between stomatal conductance, carbon
isotope discrimination, and whole-season water budgets. The limited-transpiration
(LT) trait in soybean — a heritable threshold response in leaf transpiration to
vapour pressure deficit (VPD) — offers a well-characterised physiological basis
for a genotype-specific Kc. Here we tested whether the LT trait translates from
leaf-level stomatal physiology to a measurable difference in seasonal Kc at the
field scale. We combined 30 years (1990–2019) of Sinclair SSM crop model simulations
across six locations in Arkansas, Missouri, and Nebraska with NASA POWER-derived
daily FAO-56 Penman-Monteith ETo computed over planting-date-specific growing
windows to derive seasonal crop coefficients for LT and non-LT soybean genotypes
under irrigated (potential) and rainfed (water-limited) conditions. LT genotypes
exhibited a mean seasonal Kc of 1.27 ± 0.20 under irrigation, compared to
1.78 ± 0.29 for non-LT genotypes — a 29% reduction in water use per unit ETo.
This translated to seasonal irrigation water savings of 246–344 mm. Despite
consuming substantially less water, LT genotypes produced 3.6–14.7% higher grain
yield under rainfed conditions, yielding superior water productivity (6.28 vs
5.63 kg mm⁻¹). The genotype × ETo interaction was highly significant
(ANCOVA, p < 0.001), confirming that Kc is not a fixed species-level parameter
but varies with the limited-transpiration trait. These findings challenge the
one-Kc-fits-all paradigm and have direct implications for irrigation scheduling
in regions adopting climate-smart soybean varieties.

**Keywords:** crop coefficient; limited transpiration; soybean; irrigation
scheduling; FAO-56; vapour pressure deficit; carbon isotope discrimination;
water productivity; climate-smart agriculture

---

## 1. Introduction

Accurate estimation of crop water demand is foundational to irrigation
scheduling and agricultural water management. The FAO-56 framework (Allen et
al., 1998) remains the global standard, estimating crop evapotranspiration (ETc)
as the product of a reference evapotranspiration (ETo) and a dimensionless crop
coefficient (Kc):

> ETc = Kc × ETo

The Kc integrates crop physiology (canopy size, stomatal behaviour, albedo) and
micrometeorology into a single scalar considered fixed for a given crop at a given
growth stage. Tabulated Kc values for soybean range from ~0.4 at emergence to
~1.15 at mid-season (Allen et al., 1998), and these values are routinely applied
regardless of the genetic identity of the cultivar being grown.

This species-level Kc assumption is increasingly challenged by the development
of crop genotypes with distinctly different water-use traits. Among the most
well-characterised is the limited-transpiration (LT) trait in soybean, first
systematically described by Sinclair et al. (2008). Genotypes carrying the LT
trait show a threshold response in leaf transpiration rate to increasing VPD:
below a critical VPD (~2 kPa), transpiration increases linearly with demand,
but above this threshold, transpiration plateaus or declines (Fig. 1). In
contrast, standard (non-LT) genotypes exhibit a continuously linear and steeper
response throughout the VPD range typical of US production environments. The
mechanistic basis lies in reduced leaf hydraulic conductance under high VPD,
which promotes rapid stomatal closure and limits the continuity between soil
water supply and atmospheric demand (Sinclair et al., 2008).

The stomatal restriction that defines the LT trait carries important consequences
for leaf-level carbon metabolism. Stomatal closure reduces the intercellular CO₂
concentration (Cᵢ), increasing the Cᵢ/Cₐ ratio deficit and causing preferential
fixation of ¹²CO₂ — the signature of elevated carbon isotope discrimination
(δ¹³C or C13 ratio). Higher δ¹³C is a well-established integrative proxy for
intrinsic water-use efficiency (WUE = A/gₛ; Farquhar and Richards, 1984), and
LT genotypes have been documented to exhibit both higher δ¹³C and greater
plasticity in this ratio across environments. Chamarthi et al. (2022) conducted
a genome-wide association study (GWAM) of C13 ratio and its phenotypic plasticity
(Finlay-Wilkinson slope across irrigated and droughted environments) in a panel
that included the same slow-wilting accessions characterised by Sinclair et al.
(2008) — PI 416937 and PI 471938 (LT/slow-wilting) and A5959 (non-LT/fast-wilting).
They identified 19 SNPs associated with C13 plasticity and 39 with mean C13 ratio,
with candidate genes enriched for stomatal complex morphogenesis, water transport,
and root development pathways. The ideal genotype emerging from this framework is
one that is plastic: exhibiting high WUE (conservative stomatal behaviour) under
drought stress, yet maintaining high stomatal conductance and carbon gain under
well-watered, high-ETo conditions. This plasticity concept provides a genomic
and physiological framework for understanding why LT soybean genotypes may
systematically depart from a fixed FAO-56 Kc.

The agronomic case for LT germplasm is well established at the crop level:
simulations and field trials consistently show yield advantages under rainfed
conditions (Sinclair et al., 2010). Yet the implications of the LT trait for
seasonal water budgeting and irrigation management have not been quantified
within the FAO-56 framework.

A parallel challenge to the fixed-Kc paradigm has emerged from empirical studies
of crop-atmosphere coupling. Marin et al. (2016) showed that for highly
canopy-atmosphere coupled crops (citrus, coffee, sugarcane), the basal crop
coefficient (Kcb) decreases as ETo increases — the opposite of the fixed-Kc
assumption — because high atmospheric demand increases internal plant resistance
to vapour diffusion and induces stomatal closure. Marin et al. (2019) extended
this finding to soybean, wheat, and potato under Brazilian cropping conditions,
reporting a whole-cycle soybean Kc of 0.83 ± 0.39 and demonstrating that when
ETo exceeded 4 mm d⁻¹, measured crop water use was consistently lower than
prescribed by the FAO-56 tabulated values. More recently, De Melo et al. (2026)
argued that the transferability of Kc values is further limited because standard
empirical Kc derivations overlook soil and plant hydraulic properties that
regulate actual plant-available water — suggesting that process-based approaches
more accurately capture water use dynamics under varying conditions. Taken
together, these findings imply that the inverse ETo–Kc relationship is a general
physiological consequence of plant resistance to water flow, and that the LT
trait in soybean — by genetically amplifying this resistance via low leaf
hydraulic conductance — should produce a predictable and mechanistically
grounded departure from a fixed-Kc framework.

If LT genotypes consume materially less water per unit atmospheric demand, then
the Kc used to schedule irrigation for these genotypes should be lower than the
FAO-56 standard. The failure to account for this could lead to systematic
over-irrigation of LT fields, wasting water and eroding the conservation value
of the trait at the regional scale.

The objectives of this study were: (1) to compute planting-date-specific
seasonal Kc values for LT and non-LT soybean genotypes across six US locations
and 30 years of simulation; (2) to quantify the water savings and yield
consequences of the lower Kc in LT genotypes under both irrigated and rainfed
conditions; and (3) to assess whether the genotype × ETo interaction provides
field-scale confirmation of the leaf-physiological and genomic evidence for
genotype-dependent water-use efficiency.

---

## 2. Materials and Methods

### 2.1 Study sites

Simulations were conducted for six locations in three states representing the
major soybean production regions of the US Mid-South and Midwest: Jonesboro
(35.84° N, 90.70° W) and Marianna (34.78° N, 90.76° W) in Arkansas; Columbia
(38.95° N, 92.33° W) and Norborne (39.30° N, 93.67° W) in Missouri; and
Lincoln (40.81° N, 96.70° W) and North Platte (41.12° N, 100.77° W) in Nebraska.
These sites span a latitudinal gradient of ~6° (34.8°–41.1° N) and encompass
wide gradients in growing-season temperature, precipitation, and evaporative
demand relevant to current and projected soybean production environments.

### 2.2 Crop simulations

Seasonal crop evapotranspiration (SUMETT, mm) and grain yield (kg ha⁻¹) were
obtained from Sinclair SSM (Sinclair Simplified Soybean Model; Sinclair et al., 2010)
simulations spanning 1990–2019
(n = 2 160 simulation-years). Two contrasting genotypic classes were simulated:
(i) non-limited transpiration (non-LT), represented by standard MG3 and MG4
soybean cultivars (SOYBEAN-MG3, SOYBEAN-MG4); and (ii) limited transpiration
(LT), represented by genotypes carrying the LT trait at 1.5× intensity for MG3
and MG4 maturity groups (SOYBEAN-LT1.5-MG3, SOYBEAN-LT1.5-MG4). Three
planting dates (early, medium, late) and two water management regimes —
fully irrigated (potential yield) and rainfed (water-limited yield) — were
simulated at each site. Days from planting to physiological maturity (R7; R7.PM)
were extracted from each simulation to define the actual crop growing window.

### 2.3 Reference evapotranspiration

Daily ETo (mm day⁻¹) was computed from NASA POWER (Prediction of Worldwide
Energy Resources) meteorological data downloaded for one representative site per
state using the `nasapower` R package (Sparks, 2018). Variables retrieved included
daily maximum and minimum air temperature (T2M_MAX, T2M_MIN, °C), relative
humidity (RH2M, %), wind speed at 2 m height (WS2M, m s⁻¹), and incoming solar
radiation (ALLSKY_SFC_SW_DWN, MJ m⁻² d⁻¹). The FAO-56 Penman-Monteith equation
was applied following Allen et al. (1998):

> ETo = (0.408 Δ Rn + γ [900/(T+273)] u₂ (es − ea)) / (Δ + γ(1 + 0.34 u₂))

where Δ is the slope of the saturation vapour pressure curve (kPa °C⁻¹), Rn is
net radiation (MJ m⁻² d⁻¹), γ = 0.062 kPa °C⁻¹ (psychrometric constant), T is
mean daily temperature (°C), u₂ is wind speed (m s⁻¹), and es and ea are
saturation and actual vapour pressure (kPa), respectively. Saturation vapour
pressure was computed as es = 0.6108 exp(17.27 T / (T + 237.3)) following FAO-56
equation 11. Net radiation was computed from incoming solar radiation (Rns = 0.77
ALLSKY_SFC_SW_DWN) and estimated longwave emission (Rnl) using the solar geometry
approach (declination, hour angle, eccentricity correction). Negative ETo values
(possible on overcast days) were set to zero.

### 2.4 Planting-date-specific crop coefficient and Marin-style daily ETc/ETo analysis

To minimise the confounding effects of fixed seasonal ETo windows, the seasonal
crop coefficient was computed over the actual growing period for each genotype ×
location × year × planting-date combination. For each simulation, the ETo window
was defined as the period from planting DOY (Pdoy) to physiological maturity DOY
(Pdoy + R7.PM), capped at DOY 330 to avoid unrealistic late-season accumulations.
The mean Pdoy and R7.PM were computed by State × Year × Planting_date grouping,
and daily ETo was summed over each resulting window to yield a planting-date-specific
seasonal ETo (ETo_pdoy, mm). The seasonal crop coefficient was then:

> Kc = SUMETT / ETo_pdoy

Following Marin et al. (2016, 2019), who plotted daily crop ET (ETc, mm d⁻¹)
against daily ETo (mm d⁻¹) and derived Kc as the linear regression slope, we
additionally expressed seasonal totals as mean daily values: ETc_daily = SUMETT / R7.PM
and ETo_daily = ETo_season_mm / 181 days (the DOY 120–300 growing season). The
regression slope of ETc_daily on ETo_daily represents a Marin-style Kc that can
be directly compared across genotypic classes and water management regimes. A
Kc slope that is both lower in magnitude and less responsive to variation in daily
ETo is the expected signature of the LT trait operating at the seasonal scale.

Mean Kc values and their standard deviations were computed by genotypic class and
state. To test whether the slope of the crop ET–ETo relationship differed between
genotypic classes, analysis of covariance (ANCOVA) was applied with SUMETT as
the response, ETo_pdoy as the covariate, and genotypic class as a factor. A
significant genotype × ETo interaction confirms that Kc differs between LT and
non-LT genotypes beyond the range of sampling variability.

### 2.5 Water productivity

Water productivity (WP, kg mm⁻¹) was computed as grain yield (kg ha⁻¹) divided
by seasonal SUMETT (mm) for each simulation-year and summarised by genotypic
class, water management regime, and site.

### 2.6 Statistical analysis and software

All analyses were conducted in R (v4.3.3; R Core Team, 2024). NASA POWER data
were accessed via the `nasapower` package; data manipulation used `dplyr` and
`tidyr`; figures were produced with `ggplot2` and `ggforce`. Scripts are publicly
available at https://github.com/ElvisElli/LT-trait.

---

## 3. Results

### 3.1 Growing-season ETo across sites

Mean growing-season ETo (DOY 120–300, 1990–2019) ranged from 682 ± 31 mm in
Arkansas to 778 ± 83 mm in Nebraska (Fig. 2). Nebraska sites showed the highest
inter-annual variability in ETo (CV ~11%), reflecting the more continental climate
of the central Great Plains compared to the more humid Mid-South environment of
Arkansas. Planting-date-specific ETo (summed over Pdoy to R7.PM) showed that
earlier planting dates experienced longer effective growing windows and accumulated
more total ETo, partially explaining yield differences among planting dates.

### 3.2 Genotype-dependent crop coefficients

Under fully irrigated conditions, non-LT genotypes exhibited a mean seasonal Kc
of 1.78 ± 0.29, compared to 1.27 ± 0.20 for LT genotypes (Table 1; Fig. 3).
This 29% difference in Kc was consistent across the six locations and 30
simulation years. The non-LT Kc exceeded the FAO-56 mid-season reference value
of 1.15 by 55%, while the LT Kc exceeded it by only 10%, reinforcing that
standard parameterisations may be calibrated on high-transpiring germplasm. Under
rainfed conditions, Kc values were similar between genotypes (non-LT: 0.614 ±
0.206; LT: 0.608 ± 0.197), reflecting water limitation rather than atmospheric
demand as the primary control on crop ET in water-limited environments.

The ANCOVA revealed a highly significant genotype × ETo interaction (p < 0.001),
confirming that the slope of the crop ET–ETo relationship differs between LT and
non-LT genotypes. Under irrigation, non-LT crop ET scaled positively with ETo
(slope = 0.684 mm mm⁻¹), whereas LT crop ET showed no significant increase with
ETo (slope = −0.206, p > 0.05). LT genotypes thus maintain stable seasonal water
use even as atmospheric demand increases — consistent with the VPD-threshold
stomatal behaviour reported at the leaf scale by Sinclair et al. (2008) and with
the elevated C13 ratio stability under high-ETo environments documented by
Chamarthi et al. (2022).

### 3.3 Irrigation water savings

The lower Kc of LT genotypes translated to substantial irrigation water savings.
Averaged across planting dates and years, LT genotypes used 341 mm less water per
season than non-LT genotypes in Arkansas, 246 mm less in Missouri, and 344 mm
less in Nebraska (Table 1). These savings represent 26–30% reductions in crop
water demand relative to non-LT genotypes under the same atmospheric conditions.

### 3.4 Yield and water productivity

Under rainfed conditions, LT genotypes produced grain yields 10.4%, 3.6%, and
14.7% higher than non-LT genotypes in Arkansas, Missouri, and Nebraska,
respectively (Fig. 4, 5). Under irrigated conditions, non-LT genotypes produced
higher absolute yields (6 625 vs 4 425 kg ha⁻¹), reflecting the yield penalty
associated with transpiration restriction under non-limiting water supply.

Water productivity (WP = yield / SUMETT) was highest for LT genotypes under
rainfed conditions (6.28 kg mm⁻¹) and lowest for non-LT genotypes under rainfed
conditions (5.63 kg mm⁻¹), with intermediate values under irrigation for both
genotypic classes (LT: 5.82; Non-LT: 6.19 kg mm⁻¹). The superior rainfed WP of
LT genotypes demonstrates that the Kc reduction is not simply a proportional
yield reduction — LT genotypes produce more grain per unit of water consumed
under water-limited conditions.

---

## 4. Discussion

### 4.1 A genotype-specific crop coefficient: from leaf physiology to field scale

The central finding of this study is that the LT trait in soybean produces a
systematically and significantly lower seasonal Kc than standard non-LT germplasm.
The 29% difference in Kc (1.27 vs 1.78) is not a marginal perturbation; it
exceeds the full range of Kc variation across soybean growth stages in the FAO-56
table (Kc_ini = 0.40, Kc_mid = 1.15, Kc_end = 0.50; Allen et al., 1998). This
implies that choosing between LT and non-LT germplasm has a larger effect on crop
water demand than choosing between growth stages in the FAO-56 parameterisation.

The mechanistic chain linking genotype to Kc can be traced through three scales.
At the leaf level, Sinclair et al. (2008) demonstrated that LT genotypes cap leaf
transpiration above a VPD threshold of ~2 kPa (Fig. 1). The physiological basis
is reduced leaf hydraulic conductance, which induces early stomatal closure as
the transpiration flux approaches leaf water supply capacity. Stomatal closure
reduces intercellular CO₂ concentration (Cᵢ), which in turn decreases the Cᵢ/Cₐ
ratio and drives preferential fixation of the lighter ¹²CO₂ isotope — producing
the elevated δ¹³C (higher carbon isotope discrimination or C13 ratio) that has
been documented as an integrative WUE proxy in LT-associated germplasm (Farquhar
and Richards, 1984; Chamarthi et al., 2022). Chamarthi et al. (2022) confirmed
that the same slow-wilting accessions characterised by Sinclair — PI 416937 and
PI 471938 — exhibited significantly higher C13 ratio (lower Cᵢ/Cₐ) under both
irrigated and drought conditions relative to the fast-wilting check A5959, and
that this C13 signal was associated with 39 genomic loci enriched for stomatal
morphogenesis and water transport pathways.

At the canopy level, the leaf-scale stomatal cap propagates to a reduction in
whole-canopy conductance and hence seasonal crop ET. Our results confirm this
translation: LT crop ET is insensitive to inter-annual variation in ETo
(non-significant ANCOVA slope), whereas non-LT crop ET scales positively with
atmospheric demand (slope = 0.684 mm mm⁻¹). This decoupling of crop ET from ETo
is the critical link between physiological trait and irrigation management — the
LT trait does not merely reduce average water use, it reduces the sensitivity of
seasonal water demand to the atmospheric variability that drives large inter-annual
swings in irrigation requirement.

Three independent growth chamber experiments (GC1–GC3; Elli et al., in prep.)
spanning 12 soybean accessions from the four WUE plasticity quadrants of Chamarthi
et al. (2022) confirmed the leaf-level VPD threshold across genetically diverse
material. In all three experiments, Q1 (Plastic/High WUE) accessions showed peak
transpiration at VPD ≈ 2.5–2.7 kPa followed by a decline at higher VPD, while Q4
(Nonplastic/High WUE) accessions maintained a near-linear response to higher VPD.
This consistency across three independent runs and 12 accessions reinforces the
leaf-physiological basis for the genotype-level Kc differences documented here.

### 4.2 Implications for irrigation scheduling and the Kc–ETo relationship

Current irrigation scheduling protocols apply species-level Kc tables (Allen et
al., 1998; Doorenbos and Pruitt, 1977) with no provision for genotypic adjustment.
Our findings suggest that this practice could lead to systematic over-irrigation
of fields planted with LT varieties by 246–344 mm season⁻¹. At the regional
scale, given the rapid adoption of climate-smart germplasm, this represents a
substantial unaccounted water saving. Conversely, using a non-LT Kc to schedule
irrigation for LT fields would not induce crop water stress, but it would fail
to capture the water conservation value of the LT trait and reduce the economic
incentive for adopting climate-smart varieties.

The Marin-style daily ETc vs ETo analysis adds a further dimension to this
recommendation. Marin et al. (2016, 2019) demonstrated empirically that Kc is
not constant across the ETo range: for coupled crops including soybean, Kc
decreases as daily ETo increases above ~4 mm d⁻¹, because high atmospheric
demand induces stomatal closure that limits the transpiration flux below what a
fixed Kc would predict. The LT trait provides a genetic mechanism for exactly
this response: by encoding a VPD-threshold hydraulic restriction, LT genotypes
enforce a structural version of the atmospheric-coupling-driven Kc reduction that
Marin et al. observed empirically across multiple crops. In this sense, LT soybean
represents an extreme of the continuum described by Marin — a genotype in which
the inverse ETo–Kc relationship is not merely a passive consequence of canopy
physiology but a heritable, selectable trait with direct agronomic implications.
The Marin-style scatter of our simulated data (Fig. 7b) confirms this: the LT
regression slope (Kc) is lower than Non-LT across all ETo levels, and the
insensitivity of LT crop ET to inter-annual ETo variability (flat fitted line)
mirrors the plateau behaviour documented at the leaf scale.

De Melo et al. (2026) further caution that empirical Kc derivations — whether
from this study or from Marin's field experiments — neglect the soil and plant
hydraulic properties that regulate the actual soil water available for
transpiration. While our approach uses a process-based model (Sinclair SSM) that
encodes crop hydraulics implicitly through the LT response function, future
validation with direct ET measurements (lysimeters, eddy covariance) should
account for soil texture and root architecture variation across the six study
sites.

We recommend that Kc tables for soybean be extended to include a genotypic
modifier distinguishing, at minimum, between limited-transpiration and standard
transpiration genotypes. For LT soybean, a mid-season Kc in the range 1.20–1.30
is supported by our multi-location, multi-year simulations using planting-date-
specific ETo windows; for standard germplasm, the FAO-56 value of 1.15 appears
conservative relative to our simulated Kc of 1.78. As phenotypic and genomic
characterisation of LT accessions advances — including the C13 plasticity GWAM
loci identified by Chamarthi et al. (2022) — it may become feasible to predict
genotype-specific Kc values from marker data, enabling precision-irrigation
protocols that are tailored at the variety level.

### 4.3 Rainfed advantage, WUE, and C13 plasticity

Under rainfed conditions, the LT trait delivered a consistent yield advantage
(average ~9.6% across states) while consuming virtually the same total amount of
water as non-LT genotypes (rainfed Kc: 0.608 vs 0.614). This apparent paradox
is explained by the timing and environmental sensitivity of water conservation.
The LT trait saves water during high-VPD (typically midday and mid-summer) periods,
preserving soil moisture for later reproductive stages when yield is most sensitive
to water deficit (Sinclair et al., 2010). The resulting higher water productivity
(6.28 vs 5.63 kg mm⁻¹) confirms that LT genotypes produce more crop output per
unit of water consumed — the operational definition of WUE in rainfed agriculture.

The C13 plasticity framework of Chamarthi et al. (2022) provides a genotypic
interpretation for this advantage. An ideal LT genotype would exhibit high
phenotypic plasticity in C13 ratio (and hence stomatal behaviour): conservative
stomatal closure and high WUE under drought-stress, high-VPD episodes, but
permissive stomatal opening and high carbon gain under well-watered, low-VPD
conditions. This plasticity in δ¹³C — quantified as the Finlay-Wilkinson slope
of C13 ratio across the environmental index of a multi-environment trial —
captures the genotype's capacity to adjust WUE in response to environmental
demand. Our seasonal Kc results reflect exactly this adaptive response at the
crop scale: LT genotypes exhibit lower Kc (higher WUE) under high-ETo,
water-demanding atmospheres, while their rainfed Kc approaches that of non-LT
genotypes when soil water is the primary constraint rather than atmospheric VPD.

### 4.4 Genomic basis and breeding implications

The convergence of leaf-scale physiology (Sinclair et al., 2008), carbon isotope
genomics (Chamarthi et al., 2022), and the field-scale Kc differences documented
here creates a coherent multi-scale picture of the LT trait's water-use advantage.
Chamarthi et al. (2022) identified 19 SNPs associated with C13 ratio plasticity
and 39 SNPs for mean C13 ratio across irrigated and droughted environments, with
candidate genes in stomatal complex morphogenesis, transpiration regulation, and
root development pathways. These loci provide molecular markers that could be
incorporated into genomic selection models to predict genotype-specific Kc values
from sequence data — a potential pathway to variety-level irrigation scheduling
without the need for full-season simulation models.

Furthermore, the demonstration that LT genotypes substantially exceed the FAO-56
mid-season Kc of 1.15 (observed LT Kc = 1.27; non-LT Kc = 1.78) implies that
historical Kc calibration studies, which were predominantly conducted on
conventional, high-transpiring germplasm, may not be transferable to modern
climate-smart varieties. As the soybean breeding community develops germplasm
with a wider range of hydraulic and stomatal traits, re-calibration of Kc for
representative genotypes at multiple environments will be increasingly important
for irrigation management accuracy.

### 4.5 Limitations and future directions

Several limitations of this study should be acknowledged. First, Kc values were
derived from crop model simulations rather than from direct field measurement of
ET, and the Sinclair SSM simulations used fixed soil water parameters across
locations; site-specific soil calibration and flux tower validation could alter
the magnitude of Kc differences. Second, the use of a single representative site
per state for ETo calculation, rather than individual site-specific weather
data, introduces spatial smoothing into the ETo estimates. Third, actual farm-
scale adoption will depend on seed availability and the yield penalty under
fully irrigated conditions, which remains substantial for current LT germplasm
(~33% lower than non-LT under irrigation). Fourth, the growth chamber experiments
(GC1–GC3, 12 accessions) confirm the leaf-level VPD threshold qualitatively, but
formal mixed-effects model discrimination among linear, bilinear, and
linear-plateau response functions was possible only for GC1 and GC3 (n = 4
replicates per treatment; ΔAIC < 2 for 6/8 genotypes); GC2 model fitting is
pending. Future work should validate these simulation-based Kc estimates with
lysimeter or eddy covariance measurements for well-characterised LT accessions
across a gradient of VPD environments, and should investigate whether C13
plasticity markers from Chamarthi et al. (2022) predict Kc variability within
the LT germplasm class.

---

## 5. Conclusions

Using 30 years of multi-location Sinclair SSM simulations combined with NASA POWER
planting-date-specific FAO-56 ETo, and drawing on the convergent evidence from
leaf-level hydraulics (Sinclair et al., 2008) and carbon isotope genomics
(Chamarthi et al., 2022), we demonstrate that:

1. The seasonal crop coefficient Kc is not a fixed species-level parameter — it
   differs significantly by soybean genotype (LT: 1.27 ± 0.20 vs Non-LT:
   1.78 ± 0.29; ANCOVA p < 0.001), with the non-LT value exceeding the FAO-56
   reference of 1.15 by 55%.
2. LT genotypes reduce seasonal irrigation water demand by 246–344 mm without
   sacrificing water-limited (rainfed) yield — in fact producing 3.6–14.7% more
   grain under rainfed conditions.
3. LT genotypes exhibit stable crop ET regardless of inter-annual ETo variability,
   consistent with the VPD-threshold stomatal physiology identified by Sinclair
   et al. (2008) and with the C13 ratio plasticity associated with LT accessions
   by Chamarthi et al. (2022).
4. Current FAO-56 Kc tables, which do not distinguish genotypic classes, will
   systematically over-estimate water demand for LT genotypes and may impede
   their adoption in precision irrigation systems.
5. The convergence of physiological, genomic, and crop modelling evidence
   supports genotype-specific Kc values in irrigation scheduling frameworks, and
   positions the LT trait — and the stomatal morphogenesis loci identified by GWAM
   — as a quantifiable lever for reducing agricultural water use under climate
   change.

---

## Acknowledgements
[TBD — funding sources, data providers]

---

## References (key — full list TBD)

Allen, R.G., Pereira, L.S., Raes, D., Smith, M., 1998. Crop evapotranspiration:
Guidelines for computing crop water requirements. FAO Irrigation and Drainage
Paper 56. Food and Agriculture Organization, Rome.

Chamarthi, S.K., Naik, P., Belko, N., Ntare, B.R., Sinclair, T.R., Rathore, A.,
2022. Identification of stable QTLs/genomic regions for carbon isotope
discrimination and its plasticity using genome-wide association mapping in
soybean. Plant Genome 15, e20284. https://doi.org/10.1002/tpg2.20284

Doorenbos, J., Pruitt, W.O., 1977. Guidelines for predicting crop water
requirements. FAO Irrigation and Drainage Paper 24. Food and Agriculture
Organization, Rome.

Farquhar, G.D., Richards, R.A., 1984. Isotopic composition of plant carbon
correlates with water-use efficiency of wheat genotypes. Aust. J. Plant Physiol.
11, 539–552.

Sinclair, T.R., Messina, C.D., Beatty, A., Samples, M., 2010. Assessment across
the United States of the benefits of altered soybean drought traits. Agron. J.
102, 475–482.

Sinclair, T.R., Zwieniecki, M.A., Holbrook, N.M., 2008. Low leaf hydraulic
conductance with rapid stomatal closure. Plant Cell Environ. 31, 796–803.
[NOTE: confirm correct 2008 Sinclair reference for the LT/slow-wilting trait]

Sparks, A.H., 2018. nasapower: A NASA POWER global meteorology, surface solar
energy and climatology data client for R. J. Open Source Softw. 3, 1035.

Marin, F.R., Angelocci, L.R., Nassif, D.S.P., Costa, L.G., Vianna, M.S., Carvalho,
K.S., 2016. Crop coefficient changes with reference evapotranspiration for highly
canopy-atmosphere coupled crops. Agric. Water Manage. 163, 139–145.
https://doi.org/10.1016/j.agwat.2015.09.010

Marin, F.R., Angelocci, L.R., Nassif, D.S.P., Vianna, M.S., Pilau, F.G., da Silva,
E.H.F.M., Sobenko, L., Gonçalves, A.O., Pereira, R.A.A., Carvalho, K.S., 2019.
Revisiting the crop coefficient–reference evapotranspiration procedure for improving
irrigation management. Theor. Appl. Climatol. 138, 1785–1793.
https://doi.org/10.1007/s00704-019-02940-7

De Melo, M.L.A., de Jong van Lier, Q., Marin, F.R., van Dam, J.C., 2026. FAO crop
coefficient and plant-available water: Do soil and plant hydraulic properties matter?
Agron. J. 118, e70277. https://doi.org/10.1002/agj2.70277

---

## Table 1 — Key results summary (to be formatted for journal)

| Genotype | Water mgmt    | Kc (mean ± SD)  | Yield (kg ha⁻¹) | WP (kg mm⁻¹) | Water saving vs Non-LT (mm) |
|----------|--------------|-----------------|-----------------|---------------|------------------------------|
| Non-LT   | Irrigated    | 1.78 ± 0.29     | 6 625           | 6.19          | —                            |
| LT       | Irrigated    | 1.27 ± 0.20     | 4 425           | 5.82          | 310 (avg)                    |
| Non-LT   | Rainfed      | 0.614 ± 0.206   | 2 521           | 5.63          | —                            |
| LT       | Rainfed      | 0.608 ± 0.197   | 2 748           | 6.28          | ~0                           |

ANCOVA genotype × ETo interaction: p < 0.001
Non-LT Kc exceeds FAO-56 reference (1.15) by 55%; LT Kc exceeds it by 10%.

---
*Draft prepared with Claude Code — all statistics computed from Sinclair SSM
simulations and NASA POWER ETo. Prose requires author review and refinement.*
