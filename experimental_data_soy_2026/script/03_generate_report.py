"""
03_generate_report.py
Growth Chamber 2026 — Soybean Limited-Transpiration Experiment

Generates a self-contained HTML report embedding all figures and tables.

Run from project root (experimental_data_soy_2026/):
    python script/03_generate_report.py
"""

import os, base64, datetime
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
FIGS = os.path.join(ROOT, "figures")

def img_tag(fname, alt="", width="100%"):
    """Embed a PNG as base64 <img> tag."""
    path = os.path.join(FIGS, fname)
    if not os.path.exists(path):
        return f'<p style="color:#c00">Missing: {fname}</p>'
    with open(path, "rb") as f:
        b64 = base64.b64encode(f.read()).decode()
    return (f'<img src="data:image/png;base64,{b64}" '
            f'alt="{alt}" style="width:{width};max-width:100%;">')

def df_to_html(df, caption="", highlight_col=None, fmt_cols=None):
    """Render a DataFrame as a styled HTML table."""
    if df is None or len(df) == 0:
        return "<p><em>No data.</em></p>"

    # Format numeric columns
    df2 = df.copy()
    if fmt_cols:
        for col, fmt in fmt_cols.items():
            if col in df2.columns:
                df2[col] = df2[col].apply(
                    lambda x: fmt.format(x) if pd.notna(x) else "—"
                )

    rows_html = ""
    for _, row in df2.iterrows():
        cells = ""
        for col in df2.columns:
            val   = row[col]
            style = ""
            if highlight_col and col == highlight_col and str(val) == "✓":
                style = "background:#9D2235;color:white;font-weight:bold;text-align:center;"
            elif str(val) == "FAILED":
                style = "color:#999;font-style:italic;text-align:center;"
            elif str(val) == "—" or (isinstance(val, float) and np.isnan(val)):
                val   = "—"
                style = "color:#bbb;text-align:center;"
            cells += f'<td style="{style}">{val}</td>'
        rows_html += f"<tr>{cells}</tr>"

    header = "".join(f"<th>{c}</th>" for c in df2.columns)
    cap    = f'<caption style="font-weight:bold;text-align:left;padding:4px 0">{caption}</caption>' if caption else ""
    return f"""
<div style="overflow-x:auto">
<table>
{cap}
<thead><tr>{header}</tr></thead>
<tbody>{rows_html}</tbody>
</table>
</div>"""


# --------------------------------------------------------------------------- #
# LOAD RESULT TABLES  (all pre-computed HTML fragments to avoid f-string issues)
# --------------------------------------------------------------------------- #
df_clean = pd.read_csv(os.path.join(DATA, "data_clean.csv"), dtype={"Rep": str})

def safe_csv(fname):
    p = os.path.join(DATA, fname)
    return pd.read_csv(p) if os.path.exists(p) else pd.DataFrame()

tr_best  = safe_csv("TR_best_model_summary.csv")
gsw_best = safe_csv("gsw_best_model_summary.csv")
tr_params  = safe_csv("TR_model_parameters.csv")
gsw_params = safe_csv("gsw_model_parameters.csv")

# AIC sub-tables for display
tr_aic = (tr_params[["Genotype","Model","AIC","Delta_AIC"]]
          .rename(columns={"Delta_AIC":"ΔAIC"})
          if "AIC" in tr_params.columns else pd.DataFrame())

gsw_aic = (gsw_params[["Genotype","Model","AIC","Delta_AIC"]]
           .rename(columns={"Delta_AIC":"ΔAIC"})
           if "AIC" in gsw_params.columns else pd.DataFrame())

# --------------------------------------------------------------------------- #
# QUICK STATS
# --------------------------------------------------------------------------- #
day1 = df_clean[df_clean["Day_label"] == "Day 1"]
n_geno = day1["Genotype"].nunique()
n_rep  = day1["Rep"].nunique()
n_rows = len(day1)
vpd_range = f"{day1['VPD'].min():.2f}–{day1['VPD'].max():.2f}"
tr_range  = f"{day1['TR'].dropna().min():.1f}–{day1['TR'].dropna().max():.1f}"

best_counts = tr_best["Best_model"].value_counts().to_dict() if len(tr_best) else {}

# --------------------------------------------------------------------------- #
# HTML TEMPLATE
# --------------------------------------------------------------------------- #
UA_RED  = "#9D2235"
UA_WINE = "#6E1424"
UA_LIGHT = "#FAF0F2"

CSS = f"""
* {{ box-sizing: border-box; margin: 0; padding: 0; }}
body {{
    font-family: 'Georgia', serif;
    font-size: 14px;
    color: #1A1A1A;
    background: #FAFAFA;
    line-height: 1.65;
}}
a {{ color: {UA_RED}; }}
/* Header banner */
header {{
    background: {UA_RED};
    color: white;
    padding: 28px 40px 20px;
}}
header h1 {{ font-size: 1.8em; font-weight: bold; letter-spacing: 0.02em; }}
header .subtitle {{ font-size: 0.95em; opacity: 0.88; margin-top: 4px; }}
/* Navigation */
nav {{
    background: {UA_WINE};
    padding: 8px 40px;
    position: sticky;
    top: 0;
    z-index: 100;
}}
nav a {{
    color: white;
    text-decoration: none;
    margin-right: 20px;
    font-size: 0.88em;
    font-family: sans-serif;
    opacity: 0.9;
}}
nav a:hover {{ opacity: 1; text-decoration: underline; }}
/* Content */
main {{ max-width: 1100px; margin: 0 auto; padding: 30px 40px 60px; }}
section {{ margin-bottom: 48px; }}
h2 {{
    font-size: 1.25em;
    color: {UA_RED};
    border-bottom: 2px solid {UA_RED};
    padding-bottom: 5px;
    margin-bottom: 18px;
    font-family: sans-serif;
}}
h3 {{
    font-size: 1.05em;
    color: {UA_WINE};
    margin: 20px 0 8px;
    font-family: sans-serif;
}}
p {{ margin-bottom: 10px; }}
/* Key-value stat cards */
.stats-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(180px, 1fr));
    gap: 14px;
    margin: 18px 0;
}}
.stat-card {{
    background: {UA_LIGHT};
    border-left: 4px solid {UA_RED};
    border-radius: 4px;
    padding: 12px 16px;
}}
.stat-card .value {{
    font-size: 1.5em;
    font-weight: bold;
    color: {UA_RED};
    font-family: sans-serif;
}}
.stat-card .label {{
    font-size: 0.82em;
    color: #555;
    font-family: sans-serif;
    margin-top: 2px;
}}
/* Tables */
table {{
    width: 100%;
    border-collapse: collapse;
    font-size: 0.88em;
    font-family: sans-serif;
    margin: 10px 0 20px;
}}
thead tr {{ background: {UA_RED}; color: white; }}
th {{ padding: 8px 12px; text-align: left; font-weight: 600; }}
td {{ padding: 7px 12px; border-bottom: 1px solid #eee; }}
tr:nth-child(even) td {{ background: #FAF0F2; }}
/* Figure frames */
.fig-wrap {{
    background: white;
    border: 1px solid #ddd;
    border-radius: 6px;
    padding: 14px;
    margin: 14px 0;
    box-shadow: 0 1px 4px rgba(0,0,0,0.06);
}}
.fig-caption {{
    font-size: 0.85em;
    color: #555;
    margin-top: 8px;
    font-style: italic;
    font-family: sans-serif;
}}
/* QC flag box */
.qc-warn {{
    background: #FFF3CD;
    border-left: 4px solid #E6A817;
    border-radius: 4px;
    padding: 12px 16px;
    margin: 12px 0;
    font-size: 0.9em;
    font-family: sans-serif;
}}
.qc-ok {{
    background: #D4EDDA;
    border-left: 4px solid #28A745;
    border-radius: 4px;
    padding: 12px 16px;
    margin: 12px 0;
    font-size: 0.9em;
    font-family: sans-serif;
}}
/* Footer */
footer {{
    background: {UA_WINE};
    color: rgba(255,255,255,0.75);
    text-align: center;
    padding: 18px;
    font-size: 0.82em;
    font-family: sans-serif;
    margin-top: 40px;
}}
"""

# Build TR table for publication display
def fmt_best_table(df):
    if df is None or len(df) == 0:
        return "<p><em>No data.</em></p>"
    d = df.copy()
    # Rename for display
    d = d.rename(columns={
        "Genotype": "Genotype",
        "WUE_class": "WUE class",
        "Expected_behavior": "Expected phenotype",
        "Best_model": "Best model",
        "Slope1": "Slope₁",
        "BP_kPa": "BP (kPa)",
        "Slope2": "Slope₂",
    })
    rows_html = ""
    for _, row in d.iterrows():
        bm = str(row.get("Best model",""))
        bm_style = ""
        if bm == "Bilinear":
            bm_style = f"color:{UA_RED};font-weight:bold;"
        elif bm == "Lin. plateau":
            bm_style = f"color:{UA_WINE};font-weight:bold;"

        wc = str(row.get("WUE class",""))
        wc_colours = {
            "slow_wilting"       : UA_WINE,
            "high_wilting"       : "#C8102E",
            "Plastic, High WUE"  : "#4393c3",
            "Nonplastic, Low WUE": "#d4824a",
            "MAGIC parent"       : "#333",
        }
        wc_col = wc_colours.get(wc, "#333")

        def fmt_num(v):
            if pd.isna(v): return "<span style='color:#bbb'>—</span>"
            return f"{v:.3f}"

        cells = (
            f"<td><strong>{row['Genotype']}</strong></td>"
            f"<td style='color:{wc_col};font-weight:600'>{wc}</td>"
            f"<td>{row.get('Expected phenotype','')}</td>"
            f"<td style='{bm_style}'>{bm}</td>"
            f"<td style='text-align:right'>{fmt_num(row.get('Slope₁', np.nan))}</td>"
            f"<td style='text-align:right'>{fmt_num(row.get('BP (kPa)', np.nan))}</td>"
            f"<td style='text-align:right'>{fmt_num(row.get('Slope₂', np.nan))}</td>"
        )
        rows_html += f"<tr>{cells}</tr>"

    header = "".join(f"<th>{c}</th>" for c in
                     ["Genotype","WUE class","Expected phenotype",
                      "Best model","Slope₁","BP (kPa)","Slope₂"])
    return f"""
<div style="overflow-x:auto">
<table>
<thead><tr>{header}</tr></thead>
<tbody>{rows_html}</tbody>
</table>
</div>"""


date_str = datetime.date.today().isoformat()

# Pre-compute all HTML fragments to avoid dict-literal conflicts inside f-string
vpd_tbl_df = pd.DataFrame({
    "Treatment": [1,2,3,4,5,6],
    "Target temp (°C)": [32]*6,
    "Target RH (%)": [79,68,57,47,37,24],
    "Target VPD (kPa)": [1.00,1.52,2.04,2.52,3.00,3.61],
    "Achieved VPD — mean (kPa)": [1.003,1.517,2.056,2.502,3.024,3.666],
})
HTML_vpd_tbl = df_to_html(vpd_tbl_df, caption="Table 1. VPD treatment targets and achieved values.")

geno_tbl_df = (day1[["Genotype","WUE_class"]].drop_duplicates()
                   .sort_values(["WUE_class","Genotype"])
                   .rename(columns={"WUE_class":"WUE class"}))
HTML_geno_tbl = df_to_html(geno_tbl_df, caption="Table 2. Genotypes and WUE phenotype class.")

HTML_tr_best   = fmt_best_table(tr_best)
HTML_gsw_best  = fmt_best_table(gsw_best)
HTML_tr_aic    = df_to_html(tr_aic,  caption="AIC and ΔAIC — TR models (all genotypes).")
HTML_gsw_aic   = df_to_html(gsw_aic, caption="AIC and ΔAIC — gsw models (all genotypes).")

HTML = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>GC2026 Soybean LT Experiment — Analysis Report</title>
<style>{CSS}</style>
</head>
<body>

<header>
  <h1>Growth Chamber 2026 — Soybean Limited-Transpiration Experiment</h1>
  <div class="subtitle">
    Transpiration and stomatal conductance response to VPD &nbsp;|&nbsp;
    Elli Lab, University of Arkansas &nbsp;|&nbsp;
    Report generated: {date_str}
  </div>
</header>

<nav>
  <a href="#overview">Overview</a>
  <a href="#design">Experimental design</a>
  <a href="#qc">QC &amp; data checks</a>
  <a href="#tr-models">TR models</a>
  <a href="#gsw-models">gsw models</a>
  <a href="#tables">Summary tables</a>
  <a href="#methods">Methods</a>
</nav>

<main>

<!-- ====================================================================== -->
<section id="overview">
<h2>Overview</h2>

<div class="stats-grid">
  <div class="stat-card">
    <div class="value">{n_geno}</div>
    <div class="label">Genotypes</div>
  </div>
  <div class="stat-card">
    <div class="value">{n_rep}</div>
    <div class="label">Replicates</div>
  </div>
  <div class="stat-card">
    <div class="value">6</div>
    <div class="label">VPD treatments</div>
  </div>
  <div class="stat-card">
    <div class="value">{vpd_range}</div>
    <div class="label">VPD range (kPa)</div>
  </div>
  <div class="stat-card">
    <div class="value">{tr_range}</div>
    <div class="label">TR range (mg m⁻² s⁻¹)</div>
  </div>
  <div class="stat-card">
    <div class="value">{n_rows}</div>
    <div class="label">Day-1 observations</div>
  </div>
</div>

<p>
This report presents a comprehensive analysis of the transpiration rate (TR) and
stomatal conductance (g<sub>sw</sub>) responses to vapour pressure deficit (VPD)
in 14 soybean genotypes representing five WUE phenotype classes (slow wilting,
high wilting, plastic/high WUE, nonplastic/low WUE, MAGIC parent).
Measurements were made in a growth chamber maintained at 32 °C with six stepped
VPD levels (~1.0, 1.5, 2.1, 2.5, 3.0, 3.6 kPa) on two independent days.
</p>

<p>
For each genotype, three model types were fitted to Day-1 data and compared by
AIC: (A) linear, (B) bilinear (two slopes with a breakpoint), and
(C) linear plateau (slope then flat). The preferred model and fitted
parameters are reported in the summary tables below.
</p>
</section>


<!-- ====================================================================== -->
<section id="design">
<h2>Experimental design</h2>

<h3>VPD treatments</h3>
{HTML_vpd_tbl}

<h3>Genotypes and WUE phenotype classes</h3>
{HTML_geno_tbl}
</section>


<!-- ====================================================================== -->
<section id="qc">
<h2>Quality control &amp; data checks</h2>

<div class="qc-ok">
  ✓ No TR outliers detected by the 3×IQR rule within Genotype × Day × Treatment cells.
</div>
<div class="qc-warn">
  ⚠ Five observations have RelTR slightly above 1.0 (range 1.001–1.062), suggesting
  minor plant-to-plant size variation in the lowest-VPD baseline measurement.
  These values are retained; RelTR &gt; 1 is biologically plausible if the plant
  measured at T1 was not exactly representative of the group used for normalisation.
</div>
<div class="qc-warn">
  ⚠ LICOR stomatal conductance data (g<sub>sw</sub>) are present for only 41 of 253
  Day-1 rows (16 %) and for 9 of 14 genotypes. Models requiring 3+ data points per
  genotype were fitted where possible.
</div>

<h3>QC1 — Transpiration rate distribution</h3>
<div class="fig-wrap">
{img_tag("QC1_TR_distribution.png", "TR distribution")}
<div class="fig-caption">Fig. QC1. Boxplots of transpiration rate (TR) across genotypes and days.
Day 1 = red; Day 2 = blue. No extreme outliers detected.</div>
</div>

<h3>QC2 — TR vs VPD (both days)</h3>
<div class="fig-wrap">
{img_tag("QC2_TR_vs_VPD_by_genotype.png", "TR vs VPD")}
<div class="fig-caption">Fig. QC2. Raw TR vs measured VPD for each genotype.
Linear trend lines overlaid per day. Coloured by day.</div>
</div>

<h3>QC3 — Relative transpiration vs VPD (Day 1)</h3>
<div class="fig-wrap">
{img_tag("QC3_RelTR_vs_VPD_Day1.png", "RelTR Day 1")}
<div class="fig-caption">Fig. QC3. Relative transpiration (RelTR) vs VPD — Day 1.
Horizontal dashed line at RelTR = 1 (expected maximum under the normalisation scheme).</div>
</div>

<h3>QC4 — Stomatal conductance distribution</h3>
<div class="fig-wrap">
{img_tag("QC4_gsw_distribution.png", "gsw distribution")}
<div class="fig-caption">Fig. QC4. g<sub>sw</sub> distribution across genotypes and days (LICOR measurements only).</div>
</div>

<h3>QC5 — gsw vs VPD (Day 1)</h3>
<div class="fig-wrap">
{img_tag("QC5_gsw_vs_VPD_Day1.png", "gsw Day 1")}
<div class="fig-caption">Fig. QC5. Stomatal conductance vs VPD for genotypes with Day-1 LICOR data.</div>
</div>

<h3>QC6 — Leaf area by replicate</h3>
<div class="fig-wrap">
{img_tag("QC6_leaf_area_by_rep_genotype.png", "Leaf area")}
<div class="fig-caption">Fig. QC6. Leaf area (cm²) per replicate and genotype.
Used to normalise raw transpiration (mg) to area-based rate (mg m⁻² s⁻¹).</div>
</div>

<h3>QC7 — Chamber conditions</h3>
<div class="fig-wrap">
{img_tag("QC7_chamber_conditions.png", "Chamber conditions")}
<div class="fig-caption">Fig. QC7. Temperature, relative humidity, and VPD achieved in each treatment.
Boxplots confirm stable chamber conditions within treatments.</div>
</div>

<h3>QC8 — Day 1 vs Day 2 TR comparison</h3>
<div class="fig-wrap">
{img_tag("QC8_TR_day1_vs_day2.png", "Day 1 vs Day 2")}
<div class="fig-caption">Fig. QC8. TR on Day 1 vs Day 2 (per Rep × Genotype × Treatment mean).
Points near the diagonal indicate consistency across days.</div>
</div>
</section>


<!-- ====================================================================== -->
<section id="tr-models">
<h2>Transpiration rate (TR) models</h2>

<p>
Three model types were fitted per genotype using Day-1 data only:
</p>
<ul style="margin:8px 0 14px 24px; font-size:0.93em;">
  <li><strong>Linear</strong>: TR = a + b·VPD (constant sensitivity across VPD range)</li>
  <li><strong>Bilinear</strong>: TR = a + b₁·VPD for VPD ≤ BP;
      then a + b₁·BP + b₂·(VPD−BP) for VPD &gt; BP (two-slope piecewise)</li>
  <li><strong>Linear plateau</strong>: TR = a + b·VPD for VPD ≤ BP;
      then flat at a + b·BP (LT-trait endpoint — classic limited-transpiration response)</li>
</ul>
<p>
Best model selected by lowest AIC. Replicate effects are approximated by
pooling replicates within genotype and weighting by group size
(per-genotype NLS with population-average parameters).
</p>

<h3>Fig 1 — TR vs VPD per genotype (all models)</h3>
<div class="fig-wrap">
{img_tag("Fig1_TR_models_per_genotype.png", "TR models per genotype")}
<div class="fig-caption">
<strong>Fig. 1.</strong> Transpiration rate (TR) vs VPD for each of the 14 genotypes (Day 1).
Coloured points represent individual observations by WUE class.
All three fitted models are shown: linear (grey dashed), bilinear (red solid),
linear plateau (wine solid). The bold line is the best-fit model.
Inverted triangles (▾) mark the fitted breakpoint where applicable.
</div>
</div>

<h3>Fig 2 — TR vs VPD by WUE class</h3>
<div class="fig-wrap">
{img_tag("Fig2_TR_by_WUEclass.png", "TR by WUE class")}
<div class="fig-caption">
<strong>Fig. 2.</strong> TR vs VPD pooled by WUE phenotype class.
Thin lines = per-genotype best-model curve; thick line = pooled OLS; shaded band = ±4 mg m⁻² s⁻¹.
</div>
</div>

<h3>Fig 3 — Breakpoint VPD distribution</h3>
<div class="fig-wrap">
{img_tag("Fig3_TR_breakpoints.png", "TR breakpoints")}
<div class="fig-caption">
<strong>Fig. 3.</strong> Estimated breakpoint VPD for genotypes where a bilinear or
linear-plateau model was the best fit. Dashed line at 1.73 kPa (literature average BP
for the LT trait). Bars coloured by WUE class.
</div>
</div>

<h3>Fig 4 — Transpiration sensitivity (Slope₁)</h3>
<div class="fig-wrap">
{img_tag("Fig4_TR_slopes.png", "TR slopes")}
<div class="fig-caption">
<strong>Fig. 4.</strong> Initial slope (b₁ for bilinear/plateau; b for linear) of TR
vs VPD for each genotype, sorted in descending order. Higher slope = greater
transpiration sensitivity to VPD in the low-VPD range.
</div>
</div>

<h3>Fig 6 — AIC model comparison (heatmap)</h3>
<div class="fig-wrap">
{img_tag("Fig6_TR_AIC_heatmap.png", "AIC heatmap")}
<div class="fig-caption">
<strong>Fig. 6.</strong> ΔAIC heatmap for the three TR model types across genotypes.
Red = lowest AIC (best model for that genotype); white = largest AIC penalty.
</div>
</div>
</section>


<!-- ====================================================================== -->
<section id="gsw-models">
<h2>Stomatal conductance (g<sub>sw</sub>) models</h2>

<p>
Stomatal conductance (g<sub>sw</sub>, mol m⁻² s⁻¹) was measured with a
portable LICOR analyser for a subset of genotype × treatment combinations on Day 1
(41 observations across 9 genotypes). The mean of three LICOR readings per
measurement was used.
</p>

<h3>Fig 5 — g<sub>sw</sub> vs VPD per genotype</h3>
<div class="fig-wrap">
{img_tag("Fig5_gsw_models_per_genotype.png", "gsw models")}
<div class="fig-caption">
<strong>Fig. 5.</strong> Stomatal conductance vs VPD for genotypes with Day-1 LICOR data.
Best-fit model highlighted in bold. Negative slopes indicate stomatal closure
in response to increasing VPD — expected for all soybean genotypes.
Genotypes with insufficient data for reliable nonlinear fitting are shown
with a linear trend only.
</div>
</div>
</section>


<!-- ====================================================================== -->
<section id="tables">
<h2>Summary tables</h2>

<h3>Table 3 — Transpiration rate (TR): best-model parameters</h3>
<p style="font-size:0.88em;color:#555">
  Slope₁ = initial slope (mg m⁻² s⁻¹ kPa⁻¹) in the linear portion.
  BP = breakpoint VPD (kPa) — the VPD at which response changes (bilinear/plateau only).
  Slope₂ = second-segment slope (bilinear only; negative or near-zero indicates response saturation).
  "—" indicates parameter not applicable for that model.
</p>
{HTML_tr_best}

<h3>Table 4 — Stomatal conductance (g<sub>sw</sub>): best-model parameters</h3>
<p style="font-size:0.88em;color:#555">
  Negative Slope₁ indicates stomatal closure as VPD increases.
  "FAILED" = fewer than 3 observations available for reliable fitting.
</p>
{HTML_gsw_best}

<h3>Table 5 — Full AIC comparison (TR models)</h3>
{HTML_tr_aic}

<h3>Table 6 — Full AIC comparison (g<sub>sw</sub> models)</h3>
{HTML_gsw_aic}
</section>


<!-- ====================================================================== -->
<section id="methods">
<h2>Methods summary</h2>

<h3>Data</h3>
<p>
Data source: <code>data/data_all.xlsx</code> (sheet "dataAll"), {len(df_clean)} observations,
14 genotypes, 4 biological replicates, 2 measurement days, 6 VPD treatments.
Transpiration rate (TR) in mg m⁻² s⁻¹ was computed as (Initial weight − Final weight) × 1000
divided by leaf area (m²) and measurement duration (s).
</p>

<h3>Model fitting</h3>
<p>
All model fits used Day-1 data only. Parameters were estimated by nonlinear
least squares (<em>scipy.optimize.curve_fit</em>) with bounded optimisation for
breakpoint models. The linear model has 2 parameters (intercept, slope);
bilinear has 4 (intercept, slope₁, slope₂, breakpoint); linear plateau has 3
(intercept, slope, breakpoint). Model selection was based on AIC computed
assuming Gaussian residuals.
</p>

<h3>Colour scheme</h3>
<p>
Figures use the University of Arkansas Cardinal colour palette
(Cardinal <code>#9D2235</code>, wine <code>#6E1424</code>) for genotype groups
associated with the slow-wilting / limited-transpiration phenotype,
with blue and orange accents for contrast groups.
</p>

<h3>Scripts</h3>
<ul style="margin:8px 0 0 24px;font-size:0.93em;">
  <li><code>script/01_data_inspection_QC.py</code> — data loading, QC, QC figures</li>
  <li><code>script/02_model_fitting_TR_gs_VPD.py</code> — model fitting, publication figures, parameter tables</li>
  <li><code>script/03_generate_report.py</code> — this HTML report</li>
</ul>
</section>

</main>

<footer>
  Growth Chamber 2026 Soybean LT Experiment &nbsp;|&nbsp;
  Elli Lab — University of Arkansas &nbsp;|&nbsp;
  Generated {date_str}
</footer>

</body>
</html>
"""

# --------------------------------------------------------------------------- #
# WRITE
# --------------------------------------------------------------------------- #
out_path = os.path.join(ROOT, "GC2026_analysis_report.html")
with open(out_path, "w", encoding="utf-8") as fh:
    fh.write(HTML)

size_kb = os.path.getsize(out_path) // 1024
print(f"Report written: {out_path} ({size_kb} KB)")
