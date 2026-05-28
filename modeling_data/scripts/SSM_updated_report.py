"""
SSM_updated_report.py — Detailed HTML report for updated 10-location SSM simulations.
Run from /home/user/LT-trait:
    python3 modeling_data/scripts/SSM_updated_report.py
"""

import os, base64, io, warnings
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats as scipy_stats

warnings.filterwarnings("ignore")

YEARLY = "modeling_data/data/updated-outputs/all_yearly.csv"
OUTDIR = "modeling_data/figures"
REPDIR = "modeling_data/reports"
os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(REPDIR, exist_ok=True)

# ── Colour conventions ────────────────────────────────────────────────────────
COL_GENO = {
    "Non-LT":  "#31a354",
    "LT1.5":   "#d95f02",
    "LT2":     "#7570b3",
    "LT2.5":   "#e7298a",
}
STATE_ORDER = ["Arkansas", "Missouri", "Nebraska"]
LOC_ORDER = [
    "Jonesboro","Keiser","Marianna","Rohwer",         # AR
    "Albany","Mount Vernon","Novelty",                 # MO
    "Eustis","Lincoln","North Platte",                 # NE
]
GENO_ORDER = ["Non-LT","LT1.5","LT2","LT2.5"]

plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 9,
    "figure.dpi": 150, "savefig.dpi": 150,
})

# ── Load data ─────────────────────────────────────────────────────────────────
print("Loading yearly data …")
yr = pd.read_csv(YEARLY)
yr["Genotype_class"] = pd.Categorical(yr["Genotype_class"], categories=GENO_ORDER, ordered=True)
yr["Location"] = yr["Location"].str.strip()

# Fix North Platte soft-hyphen in name
yr["Location"] = yr["Location"].str.replace("­", "", regex=False)

# Numeric cols
for c in ["SUMETT","Yield_kgha","WP","R7PM","MXLAI","HI","CRAIN","CIRGW"]:
    yr[c] = pd.to_numeric(yr[c], errors="coerce")

yr_irri = yr[yr["water_mgmt"] == "Irrigated"]
yr_rfed = yr[yr["water_mgmt"] == "Rainfed"]

print(f"  {len(yr)} rows | {yr['Location'].nunique()} locations | {yr['Crop'].nunique()} crops")

# ── Helper: fig → base64 ──────────────────────────────────────────────────────
def fig_to_b64():
    buf = io.BytesIO()
    plt.savefig(buf, format="png", bbox_inches="tight", facecolor="white", dpi=150)
    plt.close()
    return base64.b64encode(buf.getvalue()).decode()

# ── Figure 1: Mean irrigated yield by location and genotype ──────────────────
print("Fig 1: Irrigated yield by location …")
agg = (yr_irri.groupby(["Location","Genotype_class"])["Yield_kgha"]
       .agg(mean="mean", se=lambda x: x.std()/np.sqrt(len(x))).reset_index())
agg["Location"] = pd.Categorical(agg["Location"], categories=LOC_ORDER, ordered=True)
agg = agg.sort_values(["Location","Genotype_class"])

fig, ax = plt.subplots(figsize=(14, 5))
locs = [l for l in LOC_ORDER if l in agg["Location"].values]
x = np.arange(len(locs))
w = 0.18
for i, g in enumerate(GENO_ORDER):
    sub = agg[agg["Genotype_class"]==g]
    sub = sub.set_index("Location").reindex(locs)
    ax.bar(x + (i-1.5)*w, sub["mean"], w, yerr=sub["se"],
           color=COL_GENO[g], label=g, alpha=0.85, capsize=3, error_kw={"linewidth":0.7})
ax.set_xticks(x); ax.set_xticklabels(locs, rotation=35, ha="right")
ax.set_ylabel("Grain yield (kg ha⁻¹)")
ax.legend(title=None, frameon=True)
ax.axvline(3.5, color="grey", linewidth=0.7, linestyle="--")
ax.axvline(6.5, color="grey", linewidth=0.7, linestyle="--")
ax.text(1.75, ax.get_ylim()[1]*0.97, "Arkansas", ha="center", fontsize=9, color="grey")
ax.text(5.0,  ax.get_ylim()[1]*0.97, "Missouri",  ha="center", fontsize=9, color="grey")
ax.text(8.5,  ax.get_ylim()[1]*0.97, "Nebraska",  ha="center", fontsize=9, color="grey")
ax.grid(axis="y", alpha=0.3); ax.set_axisbelow(True)
fig_irri_yield = fig_to_b64()

# ── Figure 2: Mean rainfed yield by location and genotype ────────────────────
print("Fig 2: Rainfed yield by location …")
agg2 = (yr_rfed.groupby(["Location","Genotype_class"])["Yield_kgha"]
        .agg(mean="mean", se=lambda x: x.std()/np.sqrt(len(x))).reset_index())
agg2["Location"] = pd.Categorical(agg2["Location"], categories=LOC_ORDER, ordered=True)
agg2 = agg2.sort_values(["Location","Genotype_class"])

fig, ax = plt.subplots(figsize=(14, 5))
for i, g in enumerate(GENO_ORDER):
    sub = agg2[agg2["Genotype_class"]==g].set_index("Location").reindex(locs)
    ax.bar(x + (i-1.5)*w, sub["mean"], w, yerr=sub["se"],
           color=COL_GENO[g], label=g, alpha=0.85, capsize=3, error_kw={"linewidth":0.7})
ax.set_xticks(x); ax.set_xticklabels(locs, rotation=35, ha="right")
ax.set_ylabel("Grain yield (kg ha⁻¹)")
ax.legend(title=None, frameon=True)
ax.axvline(3.5, color="grey", linewidth=0.7, linestyle="--")
ax.axvline(6.5, color="grey", linewidth=0.7, linestyle="--")
ax.text(1.75, ax.get_ylim()[1]*0.97, "Arkansas", ha="center", fontsize=9, color="grey")
ax.text(5.0,  ax.get_ylim()[1]*0.97, "Missouri",  ha="center", fontsize=9, color="grey")
ax.text(8.5,  ax.get_ylim()[1]*0.97, "Nebraska",  ha="center", fontsize=9, color="grey")
ax.grid(axis="y", alpha=0.3); ax.set_axisbelow(True)
fig_rfed_yield = fig_to_b64()

# ── Figure 3: Irrigated SUMETT by location ───────────────────────────────────
print("Fig 3: Irrigated SUMETT …")
agg3 = (yr_irri.groupby(["Location","Genotype_class"])["SUMETT"]
        .agg(mean="mean", se=lambda x: x.std()/np.sqrt(len(x))).reset_index())
agg3["Location"] = pd.Categorical(agg3["Location"], categories=LOC_ORDER, ordered=True)
agg3 = agg3.sort_values(["Location","Genotype_class"])

fig, ax = plt.subplots(figsize=(14, 5))
for i, g in enumerate(GENO_ORDER):
    sub = agg3[agg3["Genotype_class"]==g].set_index("Location").reindex(locs)
    ax.bar(x + (i-1.5)*w, sub["mean"], w, yerr=sub["se"],
           color=COL_GENO[g], label=g, alpha=0.85, capsize=3, error_kw={"linewidth":0.7})
ax.set_xticks(x); ax.set_xticklabels(locs, rotation=35, ha="right")
ax.set_ylabel("Seasonal ETc (mm)")
ax.legend(title=None, frameon=True)
ax.axvline(3.5, color="grey", linewidth=0.7, linestyle="--")
ax.axvline(6.5, color="grey", linewidth=0.7, linestyle="--")
ax.text(1.75, ax.get_ylim()[1]*0.97, "Arkansas", ha="center", fontsize=9, color="grey")
ax.text(5.0,  ax.get_ylim()[1]*0.97, "Missouri",  ha="center", fontsize=9, color="grey")
ax.text(8.5,  ax.get_ylim()[1]*0.97, "Nebraska",  ha="center", fontsize=9, color="grey")
ax.grid(axis="y", alpha=0.3); ax.set_axisbelow(True)
fig_sumett = fig_to_b64()

# ── Figure 4: Water productivity by water_mgmt and genotype (violin) ─────────
print("Fig 4: Water productivity violin …")
fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
for ax, (mgmt, sub) in zip(axes, yr.groupby("water_mgmt")):
    vp = ax.violinplot(
        [sub[sub["Genotype_class"]==g]["WP"].dropna().values for g in GENO_ORDER],
        positions=range(len(GENO_ORDER)), showmedians=True, showextrema=False
    )
    for i, (body, g) in enumerate(zip(vp["bodies"], GENO_ORDER)):
        body.set_facecolor(COL_GENO[g]); body.set_alpha(0.5)
    vp["cmedians"].set_color("black"); vp["cmedians"].set_linewidth(1.5)
    ax.set_xticks(range(len(GENO_ORDER)))
    ax.set_xticklabels(GENO_ORDER, rotation=20, ha="right")
    ax.set_title(mgmt)
    ax.grid(axis="y", alpha=0.3); ax.set_axisbelow(True)
axes[0].set_ylabel("Water productivity (kg mm⁻¹)")
fig.tight_layout()
fig_wp = fig_to_b64()

# ── Figure 5: Yield advantage LT vs Non-LT (rainfed) ────────────────────────
print("Fig 5: LT yield advantage rainfed …")
pivot = (yr_rfed.groupby(["Location","Pyear","Genotype_class"])["Yield_kgha"]
         .mean().unstack("Genotype_class").reset_index())
for g in ["LT1.5","LT2","LT2.5"]:
    pivot[f"adv_{g}"] = (pivot[g] - pivot["Non-LT"]) / pivot["Non-LT"] * 100

pivot["Location"] = pd.Categorical(pivot["Location"], categories=LOC_ORDER, ordered=True)
adv_mean = (pivot.groupby("Location")[[f"adv_{g}" for g in ["LT1.5","LT2","LT2.5"]]]
            .mean().reindex(locs))

fig, ax = plt.subplots(figsize=(12, 5))
w2 = 0.22
for i, g in enumerate(["LT1.5","LT2","LT2.5"]):
    ax.bar(x + (i-1)*w2, adv_mean[f"adv_{g}"], w2,
           color=COL_GENO[g], label=g, alpha=0.85)
ax.axhline(0, color="black", linewidth=0.8)
ax.set_xticks(x); ax.set_xticklabels(locs, rotation=35, ha="right")
ax.set_ylabel("Rainfed yield advantage vs Non-LT (%)")
ax.legend(title="LT intensity", frameon=True)
ax.axvline(3.5, color="grey", linewidth=0.7, linestyle="--")
ax.axvline(6.5, color="grey", linewidth=0.7, linestyle="--")
ax.grid(axis="y", alpha=0.3); ax.set_axisbelow(True)
fig_adv = fig_to_b64()

# ── Figure 6: SUMETT savings LT vs Non-LT (irrigated) ───────────────────────
print("Fig 6: Irrigated SUMETT savings …")
pivot_et = (yr_irri.groupby(["Location","Pyear","Genotype_class"])["SUMETT"]
            .mean().unstack("Genotype_class").reset_index())
for g in ["LT1.5","LT2","LT2.5"]:
    pivot_et[f"sav_{g}"] = pivot_et["Non-LT"] - pivot_et[g]
pivot_et["Location"] = pd.Categorical(pivot_et["Location"], categories=LOC_ORDER, ordered=True)
sav_mean = (pivot_et.groupby("Location")[[f"sav_{g}" for g in ["LT1.5","LT2","LT2.5"]]]
            .mean().reindex(locs))

fig, ax = plt.subplots(figsize=(12, 5))
for i, g in enumerate(["LT1.5","LT2","LT2.5"]):
    ax.bar(x + (i-1)*w2, sav_mean[f"sav_{g}"], w2,
           color=COL_GENO[g], label=g, alpha=0.85)
ax.axhline(0, color="black", linewidth=0.8)
ax.set_xticks(x); ax.set_xticklabels(locs, rotation=35, ha="right")
ax.set_ylabel("Irrigation water saving vs Non-LT (mm season⁻¹)")
ax.legend(title="LT intensity", frameon=True)
ax.axvline(3.5, color="grey", linewidth=0.7, linestyle="--")
ax.axvline(6.5, color="grey", linewidth=0.7, linestyle="--")
ax.grid(axis="y", alpha=0.3); ax.set_axisbelow(True)
fig_savings = fig_to_b64()

# ── Figure 7: Season length (R7PM) by location and genotype (irrigated) ──────
print("Fig 7: Season length …")
agg7 = (yr_irri.groupby(["Location","Genotype_class"])["R7PM"]
        .agg(mean="mean", se=lambda x: x.std()/np.sqrt(len(x))).reset_index())
agg7["Location"] = pd.Categorical(agg7["Location"], categories=LOC_ORDER, ordered=True)
agg7 = agg7.sort_values(["Location","Genotype_class"])

fig, ax = plt.subplots(figsize=(14, 5))
for i, g in enumerate(GENO_ORDER):
    sub = agg7[agg7["Genotype_class"]==g].set_index("Location").reindex(locs)
    ax.bar(x + (i-1.5)*w, sub["mean"], w, yerr=sub["se"],
           color=COL_GENO[g], label=g, alpha=0.85, capsize=3, error_kw={"linewidth":0.7})
ax.set_xticks(x); ax.set_xticklabels(locs, rotation=35, ha="right")
ax.set_ylabel("Days to physiological maturity (R7.PM)")
ax.legend(title=None, frameon=True)
ax.axvline(3.5, color="grey", linewidth=0.7, linestyle="--")
ax.axvline(6.5, color="grey", linewidth=0.7, linestyle="--")
ax.grid(axis="y", alpha=0.3); ax.set_axisbelow(True)
fig_r7pm = fig_to_b64()

# ── Summary statistics table ──────────────────────────────────────────────────
print("Computing summary stats …")
def fmt(v, decimals=1):
    if pd.isna(v): return "—"
    return f"{v:.{decimals}f}"

summary_rows = []
for loc in locs:
    for mgmt in ["Irrigated","Rainfed"]:
        sub = yr[(yr["Location"]==loc) & (yr["water_mgmt"]==mgmt)]
        for g in GENO_ORDER:
            sg = sub[sub["Genotype_class"]==g]
            if len(sg) == 0: continue
            summary_rows.append({
                "Location": loc,
                "Water mgmt": mgmt,
                "Genotype": g,
                "n (years)": len(sg["Pyear"].unique()),
                "Yield mean (kg/ha)": f"{sg['Yield_kgha'].mean():.0f} ± {sg['Yield_kgha'].std():.0f}",
                "SUMETT mean (mm)":   f"{sg['SUMETT'].mean():.0f} ± {sg['SUMETT'].std():.0f}",
                "WP mean (kg/mm)":    f"{sg['WP'].mean():.2f} ± {sg['WP'].std():.2f}",
                "R7PM mean (days)":   f"{sg['R7PM'].mean():.0f}",
            })
summary_df = pd.DataFrame(summary_rows)

# Per-location LT advantage summary
adv_rows = []
for loc in locs:
    sub_i = yr_irri[yr_irri["Location"]==loc]
    sub_r = yr_rfed[yr_rfed["Location"]==loc]
    nlt_irri_y = sub_i[sub_i["Genotype_class"]=="Non-LT"]["Yield_kgha"].mean()
    nlt_irri_e = sub_i[sub_i["Genotype_class"]=="Non-LT"]["SUMETT"].mean()
    nlt_rfed_y = sub_r[sub_r["Genotype_class"]=="Non-LT"]["Yield_kgha"].mean()
    for g in ["LT1.5","LT2","LT2.5"]:
        lt_irri_y = sub_i[sub_i["Genotype_class"]==g]["Yield_kgha"].mean()
        lt_irri_e = sub_i[sub_i["Genotype_class"]==g]["SUMETT"].mean()
        lt_rfed_y = sub_r[sub_r["Genotype_class"]==g]["Yield_kgha"].mean()
        et_sav = nlt_irri_e - lt_irri_e
        rfed_adv = (lt_rfed_y - nlt_rfed_y) / nlt_rfed_y * 100 if nlt_rfed_y > 0 else np.nan
        irri_pen = (lt_irri_y - nlt_irri_y) / nlt_irri_y * 100 if nlt_irri_y > 0 else np.nan
        adv_rows.append({
            "Location": loc,
            "LT genotype": g,
            "ET saving irrigated (mm)": f"{et_sav:.0f}",
            "Rainfed yield advantage (%)": f"{rfed_adv:+.1f}",
            "Irrigated yield penalty (%)": f"{irri_pen:+.1f}",
        })
adv_df = pd.DataFrame(adv_rows)

# ── Build HTML ────────────────────────────────────────────────────────────────
print("Building HTML report …")

def df_to_html(df, idx=False):
    return df.to_html(index=idx, border=0, classes="tbl")

CSS = """
body{font-family:Arial,sans-serif;max-width:1400px;margin:auto;padding:20px;color:#222}
h1{color:#1a5276;border-bottom:3px solid #1a5276;padding-bottom:8px}
h2{color:#2874a6;border-bottom:1px solid #aed6f1;margin-top:40px}
h3{color:#1f618d}
.tbl{border-collapse:collapse;width:100%;font-size:13px;margin:12px 0}
.tbl th{background:#2874a6;color:#fff;padding:6px 10px;text-align:left}
.tbl td{padding:5px 10px;border-bottom:1px solid #eee}
.tbl tr:nth-child(even) td{background:#f5f8fa}
img{max-width:100%;border:1px solid #ddd;border-radius:4px;margin:10px 0}
.note{background:#fef9e7;border-left:4px solid #f39c12;padding:10px 14px;margin:12px 0;font-size:13px}
.grid2{display:grid;grid-template-columns:1fr 1fr;gap:16px}
.kv{display:inline-block;background:#eaf4fb;border-radius:4px;padding:4px 10px;margin:4px;font-size:13px}
"""

# Per-location section
loc_sections = ""
for loc in locs:
    state = yr[yr["Location"]==loc]["State"].iloc[0]
    sub_i = yr_irri[yr_irri["Location"]==loc]
    sub_r = yr_rfed[yr_rfed["Location"]==loc]
    years = sorted(yr[yr["Location"]==loc]["Pyear"].dropna().astype(int).unique())

    # stats
    rows = []
    for mgmt, sub in [("Irrigated", sub_i), ("Rainfed", sub_r)]:
        for g in GENO_ORDER:
            sg = sub[sub["Genotype_class"]==g]
            if len(sg) == 0: continue
            rows.append({
                "Management": mgmt, "Genotype": g,
                "n yrs": int(sg["Pyear"].nunique()),
                "Yield (kg/ha)": f"{sg['Yield_kgha'].mean():.0f} ± {sg['Yield_kgha'].std():.0f}",
                "ETc (mm)":      f"{sg['SUMETT'].mean():.0f} ± {sg['SUMETT'].std():.0f}",
                "WP (kg/mm)":    f"{sg['WP'].mean():.2f} ± {sg['WP'].std():.2f}",
                "R7.PM (d)":     f"{sg['R7PM'].mean():.0f}",
            })
    loc_df = pd.DataFrame(rows)

    loc_sections += f"""
<h2>{loc} ({state})</h2>
<p><span class='kv'>Years: {years[0]}–{years[-1]}</span>
<span class='kv'>N simulations: {len(sub_i.drop_duplicates(['Crop','Manag','Pyear']))} irrigated · {len(sub_r.drop_duplicates(['Crop','Manag','Pyear']))} rainfed</span></p>
{df_to_html(loc_df)}
"""

REPORT_HTML = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>SSM Updated Simulations — Report</title>
<style>{CSS}</style>
</head>
<body>
<h1>SSM Updated Simulations — Detailed Report</h1>
<p><b>10 locations · 4 genotypes (Non-LT, LT1.5, LT2, LT2.5) · 6 managements · 1990–2019 (30 years)</b></p>
<p>Generated: 2026-05-28 · Data: <code>modeling_data/data/updated-outputs/</code></p>

<div class="note">
<b>Genotype classification:</b>
Non-LT = SOYBEAN-MG3/MG4 (standard); LT1.5 / LT2 / LT2.5 = limited-transpiration trait at increasing
intensity (VPD threshold multiplier). Managements: RFD = rainfed, IRRI = irrigated; ELY/MID/LTE = early/mid/late planting.
Locations by state — Arkansas: Jonesboro, Keiser, Marianna, Rohwer; Missouri: Albany, Mount Vernon, Novelty;
Nebraska: Eustis, Lincoln, North Platte.
</div>

<h2>Overview: all locations</h2>

<h3>Irrigated grain yield by location and genotype</h3>
<img src="data:image/png;base64,{fig_irri_yield}" alt="Irrigated yield">

<h3>Rainfed grain yield by location and genotype</h3>
<img src="data:image/png;base64,{fig_rfed_yield}" alt="Rainfed yield">

<h3>Seasonal ETc (SUMETT) under irrigation by location and genotype</h3>
<img src="data:image/png;base64,{fig_sumett}" alt="SUMETT">

<h3>Rainfed yield advantage of LT genotypes vs Non-LT (%)</h3>
<img src="data:image/png;base64,{fig_adv}" alt="Yield advantage">

<h3>Irrigation water saving (SUMETT reduction) of LT vs Non-LT</h3>
<img src="data:image/png;base64,{fig_savings}" alt="ET savings">

<h3>Water productivity (yield / ETc) by water management and genotype</h3>
<img src="data:image/png;base64,{fig_wp}" alt="Water productivity">

<h3>Season length (days to R7.PM) under irrigation</h3>
<img src="data:image/png;base64,{fig_r7pm}" alt="Season length">

<h2>LT advantage summary across all locations</h2>
<p>Mean irrigated ET saving, rainfed yield advantage, and irrigated yield penalty relative to Non-LT.
All values are averages over planting dates and 30 simulation years.</p>
{df_to_html(adv_df)}

<h2>Per-location summaries</h2>
{loc_sections}

<h2>Overall summary statistics</h2>
{df_to_html(summary_df)}

<hr>
<p style="font-size:11px;color:#888">
Report generated by <code>SSM_updated_report.py</code> ·
Data: modeling_data/data/updated-outputs/ ·
Code: https://github.com/ElvisElli/LT-trait
</p>
</body>
</html>"""

out = os.path.join(REPDIR, "SSM_updated_report.html")
with open(out, "w", encoding="utf-8") as fh:
    fh.write(REPORT_HTML)
print(f"Report saved: {out}  ({os.path.getsize(out)//1024} KB)")
