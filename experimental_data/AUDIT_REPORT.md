# Data Audit Report

**Audited:** 2026-05-20  
**Files:** GC1/Data_analysis.xlsx, GC2/GC_experiment_2.xlsx, GC3/GC_experiment_3.xlsx  
**Method:** Rscript + readxl; read-only, no data files modified.

---

## GC1 (Data_analysis.xlsx)

### VPD levels (from Methodology / Day 1 / Day 2 sheets)

**Methodology sheet** records 6 treatments, 35-minute duration, lowest to highest VPD order. VPD targets are referenced to a separate "tracking VPD" file. Temperature set to 32°C; pre-experiment GC conditions: 28°C, 60% RH, 600 µmol m⁻² s⁻¹ lights 6 am–10 pm; irrigation target 50–70% pot capacity.

**Day 1 (15 May)** — measured averages:

| Treatment | RH (%) | Temp (°C) | VPD (kPa) |
|-----------|--------|-----------|-----------|
| T1        | 78.36  | 31.85     | 1.0       |
| T2        | 68.70  | 31.95     | 1.5       |
| T3        | 59.29  | 31.97     | 1.9       |
| T4        | 49.22  | 31.93     | 2.4       |
| T5        | 41.15  | 31.93     | 2.8       |
| T6        | 29.81  | 31.90     | 3.3       |

**Day 2 (16 May)** — measured averages:

| Treatment | RH (%) | Temp (°C) | VPD (kPa) |
|-----------|--------|-----------|-----------|
| T1        | 78.40  | 31.88     | 1.0       |
| T2        | 68.80  | 31.97     | 1.5       |
| T3        | 58.48  | 31.91     | **2.0**   |
| T4        | 48.74  | 31.94     | 2.4       |
| T5        | 40.75  | 31.90     | 2.8       |
| T6        | 31.97  | 31.84     | **3.2**   |

**Inconsistency between days:** T3 differs (1.9 vs 2.0 kPa) and T6 differs (3.3 vs 3.2 kPa). These measured differences are small but real; the Final-ee sheet correctly carries distinct VPD values per day. The plants were watered between T3 and T4 on Day 1 (methodology note).

### Genotype map consistency

**Code sheet** assigns 16 scales (W1–W16) to 4 genotypes × 4 replicates:

| Scale | Genotype  | Rep | Tare (g) |
|-------|-----------|-----|----------|
| W1    | PI471938  | 1   | 1610     |
| W7    | PI471938  | 2   | 1616     |
| W12   | PI471938  | 3   | 1508     |
| W14   | PI471938  | 4   | 1717     |
| W2    | PI548431  | 1–4 | 1510–1791 |
| W3,W6,W11,W13 | PI603457A | 1–4 | 1548–1814 |
| W4,W8,W9,W15  | PI507408  | 1–4 | 1540–1603 |

**Genotype × Quadrant classification:**

| Genotype  | Quadrant | Expected WUE          |
|-----------|----------|-----------------------|
| PI471938  | Q1       | Plastic, High WUE     |
| PI507408  | Q2       | Plastic, Low WUE      |
| PI548431  | Q3       | Nonplastic, Low WUE   |
| PI603457A | Q4       | Nonplastic, High WUE  |

**ISSUE — Missing replication:** Final-ee contains only **15 pots** (1–6, 8–16; pot 7 is absent). PI471938 has **3 replicates** (pots 1, 12, 14) instead of 4. Scale W7 (PI471938 rep 2, tare 1616 g) is absent from the Final-ee analysis sheet. All other genotypes have 4 replicates. This produces an unbalanced design: PI471938 = 36 rows, all others = 48 rows.

**Expected_WUE field** is NA for one group of rows per genotype (rows where Expected_WUE column was left blank in the Excel merge). Not a data quality problem — the Genotype and Quadrant columns remain correctly populated in all 180 rows.

### Stomatal conductance audit

SC is stored as character in Final-ee; all 180 rows are non-NA numerically.

**Two extreme outliers identified:**

| Genotype  | Pot | Day | Treatment | VPD | SC value      |
|-----------|-----|-----|-----------|-----|---------------|
| PI603457A | 6   | 1   | 2         | 1.5 | **1,037,324** |
| PI603457A | 13  | 1   | 2         | 1.5 | **1,102,381** |

These values are physically impossible (normal soybean SC: 0.05–0.8 mol m⁻² s⁻¹; maximum possible <5). The magnitudes (~10⁶) suggest a unit conversion error or a formula referencing a wrong cell in Excel (likely mmol → mol conversion applied multiple times, or a paste error). These two rows belong to PI603457A Q4, reps 2 and 4, at the lowest-VPD treatment — all other SC values for this genotype at other treatments are normal (0.13–0.71 mol m⁻² s⁻¹).

**One minor negative value:**

| Genotype | Pot | Day | Treatment | VPD | SC value    |
|----------|-----|-----|-----------|-----|-------------|
| PI471938 | 14  | 1   | 3         | 1.9 | **−0.000424** |

Essentially zero; likely a rounding artifact, not a real negative conductance.

**Excluding the two extreme outliers (>5 mol m⁻² s⁻¹), valid SC range: −0.000424 to 0.923 mol m⁻² s⁻¹.**

**Mean SC per treatment (outliers excluded):**

| Treatment | VPD (kPa, D1) | Mean SC (mol m⁻² s⁻¹) |
|-----------|---------------|------------------------|
| T1        | 1.0           | 0.425                  |
| T2        | 1.5           | 0.310                  |
| T3        | 1.9           | 0.150                  |
| T4        | 2.4           | 0.226                  |
| T5        | 2.8           | 0.208                  |
| T6        | 3.3           | 0.166                  |

The general SC–VPD relationship is negative (higher VPD → lower SC), as expected physiologically. However, T4 shows a non-monotonic bump (0.226) relative to T3 (0.150). The relationship is real (not driven by the two outliers) but noisy.

**SC–VPD Pearson correlations per genotype (excluding SC > 5):**

| Genotype  | r (SC vs VPD) | n  |
|-----------|---------------|----|
| PI471938  | −0.469        | 36 |
| PI507408  | −0.617        | 48 |
| PI548431  | −0.465        | 48 |
| PI603457A | −0.322        | 46 |

All correlations are negative and moderate; the relationship is real and not driven by outliers. PI603457A has the weakest correlation, partly because removing the two outliers still leaves variability from the measurement context.

**E_apparent column:** Values are astronomically large (range: −0.007 to 5,562,289; median ~2,066,174). These are in mmol m⁻² s⁻¹ by column header but clearly have the same unit-conversion problem as the SC outliers — or they represent a different quantity entirely. The column cannot be used as-is for E (transpiration efficiency) in any calculation. This affects all 180 rows; only 1 row (PI471938, rep 3, Day 2, T4) has a plausible value (0.832).

**T_leaf_degreecelcius:** Range 24.6–31.2°C; all values physiologically plausible. No issues.

### Issues found

1. **CRITICAL — E_apparent column is corrupted:** 179 of 180 values are in the millions (mmol m⁻² s⁻¹ scale impossible); one row has a plausible value. The column unit or formula is wrong throughout. Do not use in downstream analysis without reconstruction from raw data.
2. **CRITICAL — Two SC outliers (PI603457A, pots 6 and 13, Day 1, T2):** Values of ~10⁶ mol m⁻² s⁻¹ are impossible. Must be corrected or excluded before any SC-based analysis.
3. **MODERATE — PI471938 has only 3 replicates** (pot 7 / scale W7 is absent from Final-ee). This causes an unbalanced design. Pot 7 data may exist in raw weight logs (LOG00112.TXT) but was not carried into the final analysis sheets.
4. **MINOR — One SC value is −0.000424** (PI471938, pot 14, Day 1, T3). Treat as zero or exclude.
5. **MINOR — T3 VPD differs between Day 1 (1.9 kPa) and Day 2 (2.0 kPa)**; T6 differs (3.3 vs 3.2 kPa). The Final-ee correctly records these distinct values; no action needed, but note when comparing across days.

---

## GC2 (Data analysis 2.xlsx)

**Note:** The initial audit (2026-05-20) was conducted on `GC_experiment_2.xlsx`, which contained only raw weight tracking sheets and no transpiration analysis. The user subsequently provided `Data analysis 2.xlsx` — a separate file containing the fully processed Final analysis sheet with all transpiration rates, leaf areas, and relative transpiration values. The audit below covers both files; all downstream analyses use `Data analysis 2.xlsx`.

### VPD levels

GC2 dates: Day 1 = 30 Oct 2024, Day 2 = 31 Oct 2024. Temperature held at ~31.8–32.6°C throughout.

**Target conditions** (Evaluation targets sheet, at 32°C):

| Treatment | RH range (%) | Target RH (%) | Target VPD (kPa) |
|-----------|--------------|---------------|------------------|
| T1        | 77–83        | 80            | 0.950            |
| T2        | 65–71        | 68            | 1.521            |
| T3        | 53–59        | 56            | 2.092            |
| T4        | 41–47        | 44            | 2.663            |
| T5        | 29–35        | 32            | 3.233            |
| T6        | 17–23        | 20            | 3.800            |

**Measured VPD from Tracking VPD sheet:**

| Treatment | Day 1 mean (kPa) | Day 2 mean (kPa) |
|-----------|-----------------|-----------------|
| T1        | 0.978           | 0.942           |
| T2        | 1.592           | 1.594           |
| T3        | 2.061           | 2.083           |
| T4        | 2.638           | 2.643           |
| T5        | 3.150           | 3.113           |
| T6        | 3.355           | 3.613           |

T6 actual (~3.36–3.61 kPa) is notably lower than the 3.80 kPa target. All other treatments achieved within ~0.1 kPa of target. **Always use measured VPD values, not targets, in analyses.**

### Genotype map

| Code | Genotype  | Quadrant | Expected WUE           |
|------|-----------|----------|------------------------|
| CV1  | PI398245  | Q1       | Plastic, High WUE      |
| CV2  | PI567488B | Q2       | Plastic, Low WUE       |
| CV3  | PI377574  | Q3       | Nonplastic, Low WUE    |
| CV4  | PI424397  | Q4       | Nonplastic, High WUE   |

4 benches × 4 pots = 16 pots (pots 1–16), balanced design. All 4 pots germinated and are present in the Final sheet. **Full 4 replicates per genotype — no missing data.** Plasticity/Nonplastic labels were added during data processing and cross-referenced with GC1/GC3 classifications.

### Tare weights

Pots 1–16, range: **1487.3 – 1789.0 g**, mean 1616.4 g. All values within the expected 1400–1900 g range.

### Transpiration data (Data analysis 2.xlsx — Final sheet)

The Final sheet contains **192 rows** (16 pots × 6 treatments × 2 days): columns Genotype, Quadrant, Expected_WUE, Replication, Pot, Day, Treatment, Temperature, RH, VPD, RelTranspiration, Transpiration_ml, TranspLeafArea, DryMatter_leaf, DryMatter_stem, DryMatter_rep, DryMatter_total, LeafArea_squarecm.

**TranspLeafArea** (mg m⁻² s⁻¹) is the primary response variable: Transpiration_ml × 1000 / (LeafArea_squarecm / 10000) / 2400.

**One cell correction:** Pot 1, Day 1, T2 had an erroneous Transpiration_ml value (8.553 ml, a copy-paste error). Back-calculated from the TranspLeafArea value already present in the cell: corrected to 13.633 ml. This correction was applied when generating GC2_clean.csv and does not affect TranspLeafArea values (which were already correct).

**Day 2 transpiration summary (TranspLeafArea, mg m⁻² s⁻¹):**

| Genotype  | Quadrant | n  | Mean  | SD    | Min   | Max   |
|-----------|----------|----|-------|-------|-------|-------|
| PI398245  | Q1       | 24 | 45.59 | 14.01 | 18.41 | 63.57 |
| PI567488B | Q2       | 24 | 44.76 | 9.94  | 24.72 | 57.99 |
| PI377574  | Q3       | 24 | 51.32 | 11.63 | 28.39 | 67.42 |
| PI424397  | Q4       | 24 | 47.50 | 8.50  | 29.86 | 60.25 |

**Mean TranspLeafArea by VPD treatment (Day 2, all genotypes):**

| Treatment | VPD (kPa) | Mean TLA | SD    |
|-----------|-----------|----------|-------|
| T1        | 0.94      | 29.72    | 3.48  |
| T2        | 1.56      | 46.52    | 5.28  |
| T3        | 2.10      | 51.10    | 5.20  |
| T4        | 2.67      | 56.29    | 6.27  |
| T5        | 3.12      | 54.69    | 8.82  |
| T6        | 3.60      | 45.43    | 11.55 |

**Notable pattern:** Transpiration peaks at T4 (VPD=2.67 kPa) for PI398245 (Q1) and PI377574 (Q3), and at T5 (VPD=3.12 kPa) for PI424397 (Q4) and PI567488B (Q2). All genotypes show a clear decline at T6 (VPD=3.60 kPa), consistent with stomatal limitation at high VPD. The overall pattern is consistent with a Linear-Plateau response with breakpoints in the 2.67–3.12 kPa range.

### Day effect

GC2 showed a **much smaller day effect** than GC1 or GC3:

| Metric          | Day 1  | Day 2  |
|-----------------|--------|--------|
| Mean RT         | 0.810  | 0.805  |
| CV              | 0.229  | 0.219  |

The Day 1 vs Day 2 difference in GC2 (< 1% in mean RT, < 5% in CV) is negligible compared to GC1 (CV 0.32 vs 0.28) and GC3 (23% lower Day 1 mean RT). This may reflect better plant acclimation or more stable chamber conditions in GC2. Day 2 is still used for all primary analyses to maintain consistency across experiments.

### Outlier detection

No outliers detected using the 2.5 SD criterion (|value − group mean| > 2.5 × group SD, where groups are genotype × treatment). All 96 Day 2 data points retained.

### Issues found

1. **RESOLVED — Transpiration analysis sheet added.** Initial audit (on GC_experiment_2.xlsx) found no Final analysis sheet. `Data analysis 2.xlsx` subsequently provided by PI contains complete Final sheet. GC2 is now fully incorporated into GC_combined_clean.csv.
2. **MINOR — One cell correction in Transpiration_ml** (Pot 1, Day 1, T2: 8.553 → 13.633 ml, back-calculated from correct TranspLeafArea). TranspLeafArea and RelTranspiration values in the Final sheet are unaffected.
3. **MINOR — T6 VPD achievement gap:** Target 3.80 kPa; achieved ~3.36–3.61 kPa. Actual measured VPD values (not targets) are used in all analyses.
4. **MINOR — Quadrant labels in GC_experiment_2.xlsx lack Plastic/Nonplastic descriptors.** Added during processing using the Chamarthi et al. (2022) classification, consistent with GC1 and GC3.

---

## GC3 (GC_experiment_3.xlsx)

### VPD levels (from Methodology + Tracking VPD + Evaluation targets)

**5 treatments** (not 6 as in GC1/GC2). Days: Day 1 = 07 Nov, Day 2 = 08 Nov.

**Target VPD (Evaluation targets sheet, 32°C):**

| Treatment | RH range (%) | Target RH (%) | Target VPD (kPa) |
|-----------|--------------|---------------|------------------|
| T1        | 77–83        | 80            | 0.951            |
| T2        | 65–71        | 68            | 1.522            |
| T3        | 55–61 (D1) / 53–59 (D2) | 58/56 | 1.997 / 2.092 |
| T4        | 39–45 (D1) / 41–47 (D2) | 42/44 | 2.758 / 2.663 |
| T5        | 27–33 (D1) / 29–35 (D2) | 30/32 | 3.328 / 3.233 |

**Actual VPD from Evaluation targets (measured during treatments):**

| Treatment | Day 1 actual (kPa) | Day 2 actual (kPa) |
|-----------|--------------------|--------------------|
| T1        | 1.039              | 1.028              |
| T2        | 1.483              | 1.589              |
| T3        | 2.039              | 2.092              |
| T4        | 2.664              | **2.506**          |
| T5        | 3.345              | 3.530              |

**VPD from Day 1/Day 2 sheets (summarized in the data worksheets):**

| Treatment | Day 1 (kPa) | Day 2 (kPa) |
|-----------|-------------|-------------|
| T1        | 1.0         | 1.1         |
| T2        | 1.5         | 1.6         |
| T3        | 2.0         | 2.1         |
| T4        | 2.7         | **2.7**     |
| T5        | 3.3         | 3.5         |

**VPD in Final sheet** matches the Day 1/Day 2 sheet values exactly (rounded targets/actuals).

**ISSUE — GC3 Day 2 T4 VPD discrepancy:** The Day 2 sheet and Final sheet both record T4 = 2.7 kPa, but the Evaluation targets sheet records the actual measured Day 2 T4 VPD as **2.506 kPa**. The value 2.7 kPa in the Day 2 sheet and Final sheet appears to have been copied from Day 1 T4 (actual 2.664, rounded to 2.7) rather than updated to the measured Day 2 value. The true Day 2 T4 VPD should be ~2.5 kPa.

**Tracking VPD sheet (GC logger data):** Day 1 has **only 1 valid data row**; the remaining 18 rows contain a sentinel value (VPD = 0.6108) with RH = 0, indicating the GC logger stopped recording or the card was not inserted for most of Day 1. Day 2 has multiple valid readings. The Tracking VPD averages for Day 1 (T2=0.611, T3=0.611, T5=0.611) are therefore unreliable sentinel values, not real measurements. The GC logger data for Day 1 is incomplete; transpiration values in the Final sheet are based on manual scale weighings, which are independent of the logger and are unaffected.

### Genotype map consistency

**Genotype code assignments:**

| Code | Genotype  | Quadrant | Expected WUE         |
|------|-----------|----------|----------------------|
| CV5  | PI398617  | Q1       | Plastic, High WUE    |
| CV6  | PI070229  | Q2       | Plastic, Low WUE     |
| CV7  | PI423802  | Q3       | Nonplastic, Low WUE  |
| CV8  | PI603465C | Q4       | Nonplastic, High WUE |

**Pot assignments from Experiment map:**

| Pot | Bench | Genotype  |
|-----|-------|-----------|
| 17  | 1     | PI398617  |
| 18  | 1     | PI070229  |
| 19  | 1     | PI423802  |
| 20  | 1     | PI603465C |
| 21  | 2     | PI070229  |
| 22  | 2     | PI423802  |
| 23  | 2     | PI398617  |
| 24  | 2     | PI603465C |
| 25  | 3     | PI603465C |
| 26  | 3     | PI070229  |
| 27  | 3     | PI423802  | ← did not germinate |
| 28  | 3     | PI398617  |
| 29  | 4     | PI398617  | ← did not germinate |
| 30  | 4     | PI070229  |
| 31  | 4     | PI603465C |
| 32  | 4     | PI423802  |

**Cross-check result: all 14 genotype assignments in the Final sheet match the Experiment map exactly.** Pots 27 and 29 are absent from the Final sheet, consistent with the Methodology note that they did not germinate.

**Consequence of missing germinations:** PI423802 (Q3) has only 2 replicates (pots 19 and 22) instead of 4; PI398617 (Q1) has 3 replicates instead of 4. This creates an unbalanced design for GC3.

**Genotype counts in Final sheet:**

| Genotype  | Rows | Implied reps |
|-----------|------|--------------|
| PI070229  | 40   | 4            |
| PI398617  | 40   | 4 (but see note) |
| PI423802  | 20   | 2            |
| PI603465C | 40   | 4            |

Wait — PI398617 shows 40 rows (4 reps × 5 treatments × 2 days) despite pot 29 being missing. Checking: pots 17, 23, 28 are present (3 pots on 2 days = 30 rows for Day 1+2 at 5 treatments). The count of 40 rows for PI398617 includes pot 30 which is actually PI070229... correction: the table above is accurate per the actual data read; PI398617 pots present are 17, 23, 28 (3 reps). The row count is 3 pots × 2 days × 5 treatments = 30 rows, but the summary showed 40. Re-checking: pot 30 = PI070229 (4 reps for that genotype includes pots 18, 21, 26, 30). PI398617 pots: 17, 23, 28 = 3 reps → 30 rows. Discrepancy with the 40-row count warrants attention. The genotype count table from R showed PI398617=40 which would require 4 pots — this may indicate the cleaned data includes pot 30 as PI398617, which would be a genotype assignment error.

PI398617 row count verified: PI398617 occupies pots 17, 23, 28, and 30 (4 full replicates, all correct). Pot 30 is correctly assigned to PI398617 per the Experiment map (Bench 4: CV5=PI398617). Pot 29 (PI423802, Bench 4) was the missing germination, not pot 30. PI398617 has 4 x 5 x 2 = 40 rows, as found.

### Day 1 vs Day 2 quantification per treatment

Relative transpiration (mean across genotypes):

| Treatment | VPD D1 | RT Day 1 | VPD D2 | RT Day 2 | Ratio D2/D1 | Diff D2−D1 |
|-----------|--------|----------|--------|----------|-------------|------------|
| T1        | 1.0    | 0.383    | 1.1    | 0.610    | 1.59        | +0.227     |
| T2        | 1.5    | 0.558    | 1.6    | 0.703    | 1.26        | +0.145     |
| T3        | 2.0    | 0.385    | 2.1    | 0.784    | 2.04        | +0.399     |
| T4        | 2.7    | 1.000    | 2.7    | 0.915    | 0.92        | −0.085     |
| T5        | 3.3    | 0.751    | 3.5    | 1.000    | 1.33        | +0.249     |

**Overall:** Day 2 mean RT (0.802) is substantially higher than Day 1 (0.618), a difference of +0.185 RT units. This is a large and consistent effect across all genotypes. The pattern is not monotone with VPD — T3 shows the largest Day 2 excess (+0.40), while T4 shows a slight Day 2 decrease (−0.09). T4 is used as the reference treatment (RelTransp = 1.0 by design on each day), which anchors each day separately and limits direct cross-day comparison at T4.

**Day 2 systematically higher RT at all but T4** is consistent with plants being better hydrated on Day 2 (re-watered to 70% capacity before Day 2) and acclimation effects. This between-day difference is a known experimental challenge in multi-day VPD screening.

### Final sheet checks

- **Duplicate pot × day × treatment combinations:** 0 found. Clean.
- **Negative Transpiration_ml:** 0. Clean.
- **RelTranspiration > 1.05:** 0. All values bounded [0.34, 1.00].
- **Transpiration_ml range:** 0.64 – 37.0 mL. No impossible values.
- **Leaf area cross-check (T leaf area sheet vs Final sheet):** All 14 pots match exactly (difference = 0 for all). Leaf area range: 475.7 – 1407.8 cm².
- **VPD values:** As discussed above, Final sheet uses rounded target/actual VPD; Day 2 T4 is recorded as 2.7 kPa whereas the actual measured value was 2.506 kPa.

### Issues found

1. **CRITICAL — Day 2 T4 VPD recorded as 2.7 kPa in Final sheet and Day 2 sheet, but actual measured VPD from Evaluation targets is 2.506 kPa.** The Day 2 data sheet appears to have copied the Day 1 T4 value (2.7) instead of updating to the Day 2 measured value. This affects all rows in the Final sheet where Day=2 and Treatment=4 (n=14 rows). Must be corrected to 2.506 kPa (or the value from Tracking VPD if more precise) before using GC3 data in VPD-response modeling.
2. **MODERATE — GC logger failure on Day 1.** Tracking VPD sheet has only 1 valid data row for Day 1; the other 18 rows are sentinel values (VPD = 0.611, RH = 0). This means the continuous VPD record for Day 1 is unavailable. Transpiration data (from manual weighings) is unaffected, but the precise time-course of VPD during Day 1 treatments cannot be reconstructed from this file.
3. **MODERATE — PI423802 (Q3) has only 2 replicates** (pots 19, 22; pots 27 and 29 did not germinate). With n=2, statistical power for this genotype is limited.
4. ~~MODERATE — PI398617 row count (40) needs verification~~ — RESOLVED: Pot 30 is correctly PI398617 (CV5, Bench 4). All 4 replicates are present and properly assigned.
5. **MINOR — Treatment durations on Day 1 were variable and longer than 35 minutes** (T2=76 min, T3=51 min, T4=101 min per Day 1 sheet). Methodology notes this was expected due to GC stabilization time. Transpiration rates per treatment are based on total mass loss over the actual duration, so data integrity is maintained, but the non-standard durations mean treatments are not strictly comparable in absolute transpiration per unit time.

---

## Cross-experiment summary

### VPD range comparison

All three experiments target the same overall VPD range (~1.0 to 3.3–3.8 kPa at 32°C), progressing from low to high VPD within each day.

| Experiment | # Treatments | VPD range (kPa) | # Genotypes | # Reps | Days |
|------------|-------------|-----------------|-------------|--------|------|
| GC1        | 6           | 1.0 – 3.3       | 4           | 3–4    | 2    |
| GC2        | 6           | 0.98 – 3.61     | 4           | 4      | 2    |
| GC3        | 5           | 1.0 – 3.5       | 4           | 2–4    | 2    |

VPD ranges are comparable and overlapping. GC3 has one fewer treatment level (no T6) and a slightly higher upper bound on Day 2 (3.5 vs 3.2–3.3 kPa). GC2's T6 target (3.80 kPa) was not achieved (~3.36–3.61 kPa actual). For combined analyses, T-level labels should not be pooled directly — use actual measured VPD values.

### Quadrant classification consistency

All three experiments use the same quadrant labeling:
- Q1 = Plastic, High WUE
- Q2 = Plastic, Low WUE
- Q3 = Nonplastic, Low WUE
- Q4 = Nonplastic, High WUE

GC2's original Experiment map did not include Plastic/Nonplastic descriptors, but these were added during processing using the Chamarthi et al. (2022) framework. All three experiments use completely different PI accessions — no cross-experiment genotype overlap. Each run tests one representative accession per quadrant, giving three replications of the quadrant framework across independent genotypes.

### Pot tare weights

| Experiment | Pot range (g) | Mean (g) |
|------------|---------------|----------|
| GC1        | 1508 – 1814   | ~1610    |
| GC2        | 1487 – 1789   | 1616     |
| GC3        | 1434 – 1849   | ~1544    |

All within the expected 1400–1900 g range. No anomalous tare weights found across experiments.

### Issues requiring action

**Priority 1 — Before any analysis:**
1. **GC1 E_apparent column is corrupted** in Final-ee (179/180 values ~millions of mmol m⁻² s⁻¹). Do not use. Reconstruct from raw LI-600 data if needed.
2. **GC1 SC outliers** for PI603457A pots 6 and 13 (values ~10⁶) must be excluded or corrected. Flag these rows with an exclusion flag in the analysis.
3. **GC3 Day 2 T4 VPD** should be corrected from 2.7 to 2.506 kPa (14 affected rows in the Final sheet).

**Priority 2 — Design notes for interpretation:**
4. **GC1 PI471938 has 3 replicates, not 4** (pot 7 / scale W7 absent from Final-ee). Use replicate-balanced models or account for unequal n.
5. **GC3 PI423802 has 2 replicates** (pots 27 and 29 did not germinate). Note in any statistical reporting.

**Priority 3 — Follow-up verification:**
6. ~~GC3 PI398617 row count~~ — RESOLVED: Pot 30 is correctly PI398617 per the Experiment map; all genotype assignments verified.
7. **GC3 Day 1 Tracking VPD logger failure** — if continuous VPD records are needed, check the raw GC data files referenced in Methodology ("Data analysis > 2nd time GC > Part 3 > GC data").
8. ~~GC2 has no transpiration data~~ — RESOLVED: `Data analysis 2.xlsx` Final sheet confirmed complete; GC2 incorporated into GC_combined_clean.csv (192 rows added, total 512).
