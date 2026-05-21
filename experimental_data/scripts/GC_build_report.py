"""
Build GC_analysis_report.html from figures + data.
Run from /home/user/LT-trait:
    python3 experimental_data/scripts/GC_build_report.py
"""

import base64, os, pandas as pd, numpy as np
from datetime import date

FIGDIR = "experimental_data/figures"
OUT    = "experimental_data/GC_analysis_report.html"

def b64(name):
    path = os.path.join(FIGDIR, name)
    if not os.path.exists(path):
        return None
    with open(path, "rb") as f:
        return base64.b64encode(f.read()).decode()

def img(name, caption="", width="100%"):
    data = b64(name)
    if data is None:
        return f'<p class="missing">[Figure not found: {name}]</p>'
    return (f'<figure>'
            f'<img src="data:image/png;base64,{data}" style="width:{width};max-width:1400px;" alt="{caption}">'
            f'<figcaption>{caption}</figcaption></figure>')

# ── Data summary ───────────────────────────────────────────────────────────────
df = pd.read_csv("experimental_data/cleaned_data/GC_combined_clean.csv")
df["Genotype"] = df["Genotype"].str.strip()
df["Day"]      = df["Day"].astype(str)

def flag_outliers(sub, col, n_sd=2.5):
    grp = sub.groupby(["Experiment","Genotype","Treatment"])[col]
    mu  = grp.transform("mean"); sd = grp.transform("std")
    return (sub[col] - mu).abs() > n_sd * sd.fillna(np.inf)

d2 = df[df["Day"] == "2"].copy()
d2["out_tla"] = flag_outliers(d2, "TranspLeafArea")
d2["out_rel"] = flag_outliers(d2, "RelTranspiration")
d2_clean = d2[~d2["out_tla"] & d2["TranspLeafArea"].notna()]

stats_tbl = d2_clean.groupby(["Experiment","Genotype","Quadrant"]).agg(
    n       =("TranspLeafArea","count"),
    TLA_mean=("TranspLeafArea","mean"),
    TLA_sd  =("TranspLeafArea","std"),
    RT_mean =("RelTranspiration","mean"),
    RT_sd   =("RelTranspiration","std"),
).round(2).reset_index()

# Day effect
day_means = df.groupby(["Experiment","Genotype","Day"])["RelTranspiration"].mean().round(3).unstack("Day")
day_means.columns = ["Day 1 mean RT","Day 2 mean RT"]
day_means["Δ (D2−D1)"] = (day_means["Day 2 mean RT"] - day_means["Day 1 mean RT"]).round(3)
day_means = day_means.reset_index()

# Out counts
out_n = d2["out_tla"].sum()

QUAD_LAB = {
    "Q1": "Q1 · Plastic / High WUE",
    "Q2": "Q2 · Plastic / Low WUE",
    "Q3": "Q3 · Nonplastic / Low WUE",
    "Q4": "Q4 · Nonplastic / High WUE",
}

def df_to_html(df, cls="table"):
    return df.to_html(index=False, classes=cls, border=0, na_rep="NA")

# ── Experiment metadata ────────────────────────────────────────────────────────
exp_meta = {
    "GC1": {"date":"Oct 2024","n_pots":16,"n_geno":4,"vpd":"1.0–3.3 kPa (6 levels)","days":2,"has_tla":True,"has_dm":True,"has_sc":True,
            "genotypes":["PI471938 (Q1)","PI507408 (Q2)","PI548431 (Q3)","PI603457A (Q4)"]},
    "GC2": {"date":"Oct 30–31 2024","n_pots":16,"n_geno":4,"vpd":"0.95–3.6 kPa (6 levels)","days":2,"has_tla":True,"has_dm":True,"has_sc":False,
            "genotypes":["PI398245 (Q1)","PI567488B (Q2)","PI377574 (Q3)","PI424397 (Q4)"]},
    "GC3": {"date":"Nov 2024","n_pots":16,"n_geno":4,"vpd":"1.0–3.3 kPa (5 levels)","days":2,"has_tla":True,"has_dm":True,"has_sc":False,
            "genotypes":["PI398617 (Q1)","PI070229 (Q2)","PI423802 (Q3)","PI603465C (Q4)"]},
}

# ── HTML ───────────────────────────────────────────────────────────────────────
css = """
* { box-sizing: border-box; }
body { font-family: "Segoe UI", sans-serif; max-width: 1600px; margin: 0 auto;
       padding: 20px 30px; color: #222; background: #fff; }
h1 { color: #1a5276; border-bottom: 3px solid #1a5276; padding-bottom: 8px; }
h2 { color: #1f618d; margin-top: 2em; border-bottom: 1px solid #aed6f1; }
h3 { color: #2874a6; margin-top: 1.5em; }
h4 { color: #555; margin-top: 1em; }
figure { margin: 1.2em 0; text-align: center; }
figcaption { color: #555; font-size: 0.88em; margin-top: 4px; font-style: italic; }
.table { border-collapse: collapse; width: 100%; font-size: 0.88em; margin: 1em 0; }
.table th { background: #1f618d; color: #fff; padding: 6px 10px; text-align: left; }
.table td { padding: 5px 10px; border-bottom: 1px solid #dde; }
.table tr:nth-child(even) td { background: #f4f8fb; }
.callout { background: #eaf4fb; border-left: 4px solid #1f618d; padding: 10px 16px;
           margin: 1em 0; border-radius: 3px; }
.warn    { background: #fef9e7; border-left: 4px solid #f39c12; padding: 10px 16px;
           margin: 1em 0; border-radius: 3px; }
.missing { color: #e74c3c; font-style: italic; }
code { background: #f4f4f4; padding: 1px 5px; border-radius: 3px; font-size: 0.87em; }
.toc { background: #f8f9fa; border: 1px solid #dde; padding: 14px 20px; margin: 1em 0;
       border-radius: 4px; display: inline-block; min-width: 280px; }
.toc a { color: #1a5276; text-decoration: none; }
.toc a:hover { text-decoration: underline; }
.toc li { margin: 4px 0; }
"""

# Build stats table for body
def stats_section(exp):
    s = stats_tbl[stats_tbl["Experiment"] == exp].copy()
    s = s.rename(columns={
        "Genotype":"Genotype","Quadrant":"Quadrant","n":"n (Day 2)",
        "TLA_mean":"TLA mean\n(mg m⁻² s⁻¹)","TLA_sd":"TLA SD",
        "RT_mean":"Rel.Transp mean","RT_sd":"Rel.Transp SD"
    })
    return s.drop(columns=["Experiment"]).to_html(index=False, classes="table", border=0)

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>GC Transpiration Analysis — GC1, GC2, GC3</title>
<style>{css}</style>
</head>
<body>

<h1>Growth Chamber (GC) Soybean Transpiration Analysis</h1>
<p><strong>PI:</strong> Elvis F. Elli · University of Arkansas &nbsp;|&nbsp;
   <strong>Project:</strong> RIG 2024–2026 &nbsp;|&nbsp;
   <strong>Report generated:</strong> {date.today().isoformat()}</p>
<p>Three growth chamber experiments testing VPD-response in soybean accessions classified
by WUE plasticity (Chamarthi et al. 2022). Primary response variable:
<strong>transpiration per unit leaf area</strong> (mg m⁻² s⁻¹), measured as pot
weight loss over a 40-min treatment window normalised by leaf area.
Relative transpiration (0–1) retained for cross-experiment comparison.</p>

<div class="toc">
  <strong>Contents</strong>
  <ol>
    <li><a href="#overview">Experiment overview</a></li>
    <li><a href="#methods">Methods &amp; QC</a></li>
    <li><a href="#gc1">GC1 results</a></li>
    <li><a href="#gc2">GC2 results</a></li>
    <li><a href="#gc3">GC3 results</a></li>
    <li><a href="#overall">Overall / combined analysis</a></li>
    <li><a href="#summary">Key findings</a></li>
    <li><a href="#data">Data notes</a></li>
  </ol>
</div>

<!-- ═════════════════════════════════════════ OVERVIEW ════════════════════════ -->
<h2 id="overview">1. Experiment Overview</h2>

<table class="table">
<thead><tr><th>Experiment</th><th>Date</th><th>Pots</th><th>Genotypes</th><th>VPD range</th>
<th>Days</th><th>Leaf area</th><th>Dry matter</th><th>Stomatal cond.</th></tr></thead>
<tbody>
{"".join(
  f'<tr><td><strong>{e}</strong></td>'
  f'<td>{m["date"]}</td>'
  f'<td>{m["n_pots"]} ({m["n_geno"]} geno × 4 rep)</td>'
  f'<td>{"<br>".join(m["genotypes"])}</td>'
  f'<td>{m["vpd"]}</td>'
  f'<td>{m["days"]}</td>'
  f'<td>{"✓" if m["has_tla"] else "✗"}</td>'
  f'<td>{"✓" if m["has_dm"] else "✗"}</td>'
  f'<td>{"✓" if m["has_sc"] else "✗"}</td></tr>'
  for e,m in exp_meta.items()
)}
</tbody>
</table>

<h3>Genotype classification (Chamarthi et al. 2022)</h3>
<table class="table">
<thead><tr><th>Quadrant</th><th>WUE class</th><th>GC1</th><th>GC2</th><th>GC3</th></tr></thead>
<tbody>
<tr><td>Q1</td><td>Plastic, High WUE</td><td>PI471938</td><td>PI398245</td><td>PI398617</td></tr>
<tr><td>Q2</td><td>Plastic, Low WUE</td><td>PI507408</td><td>PI567488B</td><td>PI070229</td></tr>
<tr><td>Q3</td><td>Nonplastic, Low WUE</td><td>PI548431</td><td>PI377574</td><td>PI423802</td></tr>
<tr><td>Q4</td><td>Nonplastic, High WUE</td><td>PI603457A</td><td>PI424397</td><td>PI603465C</td></tr>
</tbody>
</table>

<!-- ═══════════════════════════════════════════ METHODS ══════════════════════ -->
<h2 id="methods">2. Methods &amp; QC</h2>

<h3>Transpiration measurement</h3>
<p>Pots were weighed at the beginning and end of each 40-min VPD treatment window.
Weight loss (g ≈ ml) = water transpired. Three weighings were taken per treatment;
values are averaged before computing transpiration (as in the processed spreadsheets).</p>
<p>Transpiration per leaf area (mg m⁻² s⁻¹) = Transpiration (ml) × 1000 / leaf area (m²) / 2400 s.</p>
<p>Relative transpiration = Transpiration(T<sub>i</sub>) / max(Transpiration) per pot and day,
normalised so the treatment with maximum weight loss = 1.0.</p>

<h3>VPD treatments</h3>
<ul>
  <li>GC1 &amp; GC2: 6 levels, target RH = 80, 68, 56, 44, 32, 20% at 32°C</li>
  <li>GC3: 5 levels, target RH = 80, 68, 56, 44, 32% at 32°C</li>
  <li>Actual VPD measured from chamber probes; values in the cleaned CSV are measured averages per treatment.</li>
</ul>

<h3>Day effect and primary analysis day</h3>
<div class="callout">
  <strong>Day 2 is used for all primary analysis.</strong>
  Day 1 shows a systematic acclimation artifact: plants not yet at equilibrium after
  transitioning from growth conditions to the measurement protocol.
  GC1 and GC3 show Day 1 mean relative transpiration 10–30% below Day 2.
  GC2 shows a smaller Day effect (plants may have been closer to equilibrium).
  All figures and statistics below use Day 2 unless stated otherwise.
</div>

{img("GC1_Fig3_day_comparison.png", "GC1 — Day 1 vs Day 2 comparison. Dashed lines = Day 1 (acclimation); solid = Day 2 (stable).")}

<h3>Outlier detection</h3>
<p>Outliers defined as observations &gt; 2.5 SD from the genotype × treatment group mean
(within each experiment, Day 2 only). Flagged values excluded from summary statistics
and figures; raw data retained in <code>GC_combined_clean.csv</code>.</p>
<p><strong>Outliers flagged (Day 2, TranspLeafArea): {int(out_n)}</strong></p>

{df_to_html(day_means)}

<!-- ═══════════════════════════════════════════ GC1 ══════════════════════════ -->
<h2 id="gc1">3. GC1 Results</h2>
<p>GC1 tested four genotypes (one per WUE quadrant), with 4 replications across 6 VPD levels,
measured on 2 consecutive days. Stomatal conductance (LI-600 porometer) was recorded concurrently.</p>

{stats_section("GC1")}

{img("GC1_Fig1_actual_transpiration.png",
     "GC1 — Transpiration per leaf area vs VPD (Day 2). Points = individual pots; line = mean ± SE; dashed = Sinclair 2008 VPD threshold (2.0 kPa).")}

{img("GC1_Fig2_relative_transpiration.png",
     "GC1 — Relative transpiration vs VPD (Day 2). Max = 1.0 per pot.")}

<!-- ═══════════════════════════════════════════ GC2 ══════════════════════════ -->
<h2 id="gc2">4. GC2 Results</h2>
<p>GC2 tested a second set of four genotypes (one per quadrant), 4 replications,
6 VPD levels over 2 days (Oct 30–31, 2024). Leaf area, dry matter, and transpiration
per leaf area were measured.</p>

<div class="warn">
  <strong>Data note — GC2 anomaly:</strong> Pot 1, Day 1, Treatment 6 (highest VPD)
  recorded a net weight <em>gain</em> (−7.4 ml transpiration), indicating irrigation
  or measurement error during this window. This value is set to <code>NA</code> in
  <code>GC2_clean.csv</code>. All other 191 records are valid.
</div>

{stats_section("GC2")}

{img("GC2_Fig1_actual_transpiration.png",
     "GC2 — Transpiration per leaf area vs VPD (Day 2). Same axes and colour convention as GC1.")}

{img("GC2_Fig2_relative_transpiration.png",
     "GC2 — Relative transpiration vs VPD (Day 2).")}

{img("GC2_Fig3_day_comparison.png",
     "GC2 — Day 1 vs Day 2 comparison. GC2 shows smaller day effect than GC1/GC3.")}

<!-- ═══════════════════════════════════════════ GC3 ══════════════════════════ -->
<h2 id="gc3">5. GC3 Results</h2>
<p>GC3 tested a third set of four genotypes, 4 replications, 5 VPD levels over 2 days
(November 2024). No stomatal conductance data.</p>

{stats_section("GC3")}

{img("GC3_Fig1_actual_transpiration.png",
     "GC3 — Transpiration per leaf area vs VPD (Day 2).")}

{img("GC3_Fig2_relative_transpiration.png",
     "GC3 — Relative transpiration vs VPD (Day 2).")}

{img("GC3_Fig3_day_comparison.png",
     "GC3 — Day 1 vs Day 2 comparison. GC3 shows a strong Day 1 deficit (~20–30% lower).")}

<!-- ══════════════════════════════════════════ OVERALL ═══════════════════════ -->
<h2 id="overall">6. Overall / Combined Analysis</h2>
<p>All three experiments combined: 12 genotypes, 512 observations total (Day 2 n ≈ 256).
Quadrant colour coding: Q1 = orange, Q2 = teal, Q3 = purple, Q4 = pink.
Experiment line style: GC1 = solid, GC2 = dashed, GC3 = dotted.</p>

{img("Overall_Fig1_actual_transpiration.png",
     "Overall — Actual transpiration per leaf area vs VPD, faceted by WUE quadrant. Each line = one genotype from one experiment.")}

{img("Overall_Fig2_relative_transpiration.png",
     "Overall — Relative transpiration vs VPD, faceted by quadrant.")}

{img("Overall_Fig3_all_genotypes.png",
     "Overall — All genotypes on a single panel. Colour = quadrant, line style = experiment.")}

{img("Overall_Fig4_experiment_comparison.png",
     "Overall — Experiment-level comparison within each quadrant. Lines represent quadrant-average transpiration per experiment.")}

<!-- ══════════════════════════════════════════ SUMMARY ═══════════════════════ -->
<h2 id="summary">7. Key Findings</h2>

<div class="callout">
<h4>Response pattern</h4>
<ul>
  <li>Across all three experiments, transpiration (actual and relative) increases with VPD
      up to approximately 2–3 kPa, consistent with the Sinclair (2008) limited-transpiration
      threshold framework.</li>
  <li>Q1 genotypes (Plastic, High WUE — PI471938, PI398245, PI398617) tend toward
      a plateau or slight decline in transpiration at high VPD, consistent with a
      VPD-threshold response and the limited-transpiration (LT) trait.</li>
  <li>Q4 genotypes (Nonplastic, High WUE) show a more linear increase with VPD
      across the tested range.</li>
  <li>Q2 and Q3 genotypes show variable responses, with Q3 generally displaying
      lower overall transpiration per leaf area.</li>
</ul>
</div>

<div class="callout">
<h4>Experiment consistency</h4>
<ul>
  <li>GC2 genotypes show somewhat higher absolute transpiration per leaf area than
      GC1/GC3 genotypes at equivalent VPDs, likely due to larger plant size at time
      of measurement (later vegetative stage).</li>
  <li>Relative transpiration patterns are more consistent across experiments,
      supporting cross-experiment comparability.</li>
  <li>The Day 1 acclimation artifact is consistent across GC1 and GC3;
      GC2 shows a smaller effect, possibly because plants were measured at a more
      advanced and stable growth stage.</li>
</ul>
</div>

<h3>Descriptive statistics — Day 2, all experiments</h3>
{df_to_html(stats_tbl)}

<!-- ══════════════════════════════════════════ DATA NOTES ════════════════════ -->
<h2 id="data">8. Data Notes</h2>

<h3>File structure</h3>
<table class="table">
<thead><tr><th>File</th><th>Description</th></tr></thead>
<tbody>
<tr><td><code>cleaned_data/GC2_clean.csv</code></td>
    <td>GC2 primary analysis file extracted from <em>Data analysis 2.xlsx</em> Final sheet.
        192 rows; one NA (Pot 1, Day 1, T6 weight-gain anomaly).</td></tr>
<tr><td><code>cleaned_data/GC_combined_clean.csv</code></td>
    <td>All three experiments combined. 512 rows × 25 columns. Primary analysis file.</td></tr>
<tr><td><code>raw_data/GC2/GC_experiment_2.xlsx</code></td>
    <td>Raw weight tracking spreadsheet (Tracking weights sheet). Source for verification.</td></tr>
<tr><td><code>raw_data/Data/Data analysis 2.xlsx</code></td>
    <td>Processed GC2 data (transpiration per pot per treatment, leaf area, dry matter).
        Source for <em>GC2_clean.csv</em>.</td></tr>
</tbody>
</table>

<h3>GC2 data processing</h3>
<p>GC2 was processed from the <em>Final</em> sheet of <code>Data analysis 2.xlsx</code>,
which matches the format of GC3's Final sheet. One cell error was detected and corrected:
Pot 1, Day 1, Treatment 2 had <code>Transpiration_ml = 8.553</code> (incorrect, copy/paste
error from Pot 2 T1) — corrected to <code>13.633</code> by back-calculation from
<code>TranspLeafArea × LeafArea_squarecm</code>.</p>

<h3>Scripts</h3>
<ul>
  <li><code>scripts/GC2_process_raw.py</code> — alternative raw-weights extraction
      (from Tracking weights sheet; superseded by Final sheet extraction)</li>
  <li><code>scripts/GC_figures.py</code> — generates all PNG figures in this report</li>
  <li><code>scripts/GC_build_report.py</code> — generates this HTML report</li>
  <li><code>scripts/GC_comprehensive_analysis.R</code> — R-based analysis (GC1 + GC3;
      update needed to include GC2)</li>
</ul>

<hr>
<p style="color:#888;font-size:0.83em;">Generated {date.today().isoformat()} ·
LT-trait project · Elvis F. Elli, University of Arkansas</p>
</body>
</html>
"""

with open(OUT, "w", encoding="utf-8") as f:
    f.write(html)

print(f"Report written: {OUT}  ({os.path.getsize(OUT)//1024} KB)")
