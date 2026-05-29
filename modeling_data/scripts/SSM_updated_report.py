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
DAILY  = "modeling_data/data/updated-outputs/all_daily.csv"
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

# ══════════════════════════════════════════════════════════════════════════════
# DAILY ETc / GLAI / DDMP ANALYSIS
# Reference ET: NASA POWER Penman-Monteith (nasapower_ETo_pdoy.csv / seasonal)
# SSM internal PET column is NOT used for reference ET here.
# ══════════════════════════════════════════════════════════════════════════════
print("Loading daily data for ETc / GLAI / DDMP analysis (large file, ~30 s) …")
DAILY_COLS = ["Location","Genotype_class","water_mgmt","DAP","TR","SEVP","ETc","GLAI","DDMP"]
chunks = []
for chunk in pd.read_csv(DAILY, usecols=DAILY_COLS, chunksize=200_000):
    chunk = chunk[chunk["DAP"] > 0]
    for c in ["ETc","GLAI","DDMP"]:
        chunk[c] = pd.to_numeric(chunk[c], errors="coerce")
    chunk["Location"] = chunk["Location"].str.strip().str.replace("­", "", regex=False)
    chunk = chunk[chunk["ETc"].notna()]
    chunk = chunk[chunk["ETc"] >= 0]
    chunks.append(chunk)
dy = pd.concat(chunks, ignore_index=True)
dy["Genotype_class"] = pd.Categorical(dy["Genotype_class"], categories=GENO_ORDER, ordered=True)
state_map = {l: "Arkansas" for l in ["Jonesboro","Keiser","Marianna","Rohwer"]}
state_map.update({l: "Missouri"  for l in ["Albany","Mount Vernon","Novelty"]})
state_map.update({l: "Nebraska"  for l in ["Eustis","Lincoln","North Platte"]})
dy["State"] = dy["Location"].map(state_map)
print(f"  {len(dy):,} valid growing-season daily rows loaded")

# 5-day DAP bins covering full cycle (up to 245 DAP)
MAX_DAP = int(dy["DAP"].max()) + 6
BIN_EDGES  = list(range(0, MAX_DAP, 5))
BIN_LABELS = [f"{i+1}–{i+5}" for i in range(0, MAX_DAP - 5, 5)]
dy["DAP_bin"] = pd.cut(dy["DAP"], bins=BIN_EDGES, labels=BIN_LABELS)
dy["DAP_mid"] = dy["DAP_bin"].astype(str).str.split("–").str[0].astype(float) + 2.5

locs_d = [l for l in LOC_ORDER if l in dy["Location"].values]

# ── Fig D1: Mean daily ETc — bar chart by location × genotype × management ───
print("Fig D1: Mean daily ETc by location, management, genotype …")
daily_mean = (dy.groupby(["Location","Genotype_class","water_mgmt"])["ETc"]
                .mean().reset_index(name="mean_daily_ETc"))
daily_se   = (dy.groupby(["Location","Genotype_class","water_mgmt"])["ETc"]
                .sem().reset_index(name="se_daily_ETc"))
daily_agg  = daily_mean.merge(daily_se, on=["Location","Genotype_class","water_mgmt"])
daily_agg["Location"] = pd.Categorical(daily_agg["Location"], categories=LOC_ORDER, ordered=True)
daily_agg = daily_agg.sort_values(["Location","Genotype_class"])

fig, axes = plt.subplots(2, 1, figsize=(14, 9), sharey=False)
for ax, mgmt in zip(axes, ["Irrigated","Rainfed"]):
    sub = daily_agg[daily_agg["water_mgmt"] == mgmt]
    x = np.arange(len(locs_d)); w = 0.18
    for i, g in enumerate(GENO_ORDER):
        sg = sub[sub["Genotype_class"] == g].set_index("Location").reindex(locs_d)
        ax.bar(x + (i-1.5)*w, sg["mean_daily_ETc"], w, yerr=sg["se_daily_ETc"],
               color=COL_GENO[g], label=g, alpha=0.85, capsize=3, error_kw={"linewidth":0.7})
    ax.set_xticks(x); ax.set_xticklabels(locs_d, rotation=35, ha="right")
    ax.set_ylabel("Mean daily ETc (mm day⁻¹)")
    ax.set_title(mgmt)
    ax.legend(title=None, frameon=True)
    ax.axvline(3.5, color="grey", linewidth=0.7, linestyle="--")
    ax.axvline(6.5, color="grey", linewidth=0.7, linestyle="--")
    ax.text(1.75, ax.get_ylim()[1]*0.97, "Arkansas", ha="center", fontsize=9, color="grey")
    ax.text(5.0,  ax.get_ylim()[1]*0.97, "Missouri",  ha="center", fontsize=9, color="grey")
    ax.text(8.5,  ax.get_ylim()[1]*0.97, "Nebraska",  ha="center", fontsize=9, color="grey")
    ax.grid(axis="y", alpha=0.3); ax.set_axisbelow(True)
plt.suptitle("Mean daily crop ET (ETc = TR + Soil evaporation) by location, management, and genotype",
             fontsize=12, y=1.01)
plt.tight_layout()
fig_d1 = fig_to_b64()

# ── Generic progression plotter (shared for ETc, GLAI, DDMP) ─────────────────
def make_progression_fig(var, ylabel, title_suffix):
    dap_state = (dy.groupby(["State","Genotype_class","water_mgmt","DAP_bin","DAP_mid"])[var]
                   .mean().reset_index(name=f"mean_{var}"))
    fig, axes = plt.subplots(2, 3, figsize=(18, 9), sharey=False, sharex=False)
    for col, state in enumerate(STATE_ORDER):
        for row, mgmt in enumerate(["Irrigated","Rainfed"]):
            ax = axes[row][col]
            sub = dap_state[(dap_state["State"] == state) & (dap_state["water_mgmt"] == mgmt)]
            for g in GENO_ORDER:
                sg = sub[sub["Genotype_class"] == g].sort_values("DAP_mid")
                sg = sg[sg[f"mean_{var}"].notna()]
                if len(sg) == 0: continue
                ax.plot(sg["DAP_mid"], sg[f"mean_{var}"],
                        color=COL_GENO[g], linewidth=1.8, label=g, marker="o", markersize=2.5)
            ax.set_title(f"{state} — {mgmt}", fontsize=10)
            ax.set_ylim(bottom=0)
            ax.grid(alpha=0.25); ax.set_axisbelow(True)
            # reference lines at R1~50, R5~90, PM~end
            for xv, lbl in [(50,"V→R1"),(90,"R1→R5")]:
                ax.axvline(xv, color="grey", linestyle=":", linewidth=0.8, alpha=0.5)
                ax.text(xv+1, ax.get_ylim()[1]*0.95, lbl, fontsize=7, color="grey")
            if col == 0: ax.set_ylabel(ylabel)
            if row == 1: ax.set_xlabel("Days after planting (DAP)")
    axes[0][0].legend(fontsize=8, frameon=True, loc="upper left")
    plt.suptitle(f"{title_suffix}\n(5-day DAP bins, averaged across 30 years and planting dates)",
                 y=1.01, fontsize=12)
    plt.tight_layout()
    return fig_to_b64()

# ── Fig D2: ETc progression ───────────────────────────────────────────────────
print("Fig D2: Seasonal ETc progression (5-day bins, full cycle) …")
fig_d2 = make_progression_fig("ETc", "Mean daily ETc (mm day⁻¹)",
                               "Seasonal progression of daily ETc by state and management")

# ── Fig D3: GLAI progression ──────────────────────────────────────────────────
print("Fig D3: GLAI seasonal progression …")
fig_d3 = make_progression_fig("GLAI", "Mean green LAI (m² m⁻²)",
                               "Seasonal progression of green leaf area index (GLAI) by state and management")

# ── Fig D4: DDMP progression ──────────────────────────────────────────────────
print("Fig D4: DDMP seasonal progression …")
fig_d4 = make_progression_fig("DDMP", "Mean daily dry matter production (g m⁻²)",
                               "Seasonal progression of daily dry matter production (DDMP) by state and management")

# ── Fig D5: Daily ETc reduction vs Non-LT ────────────────────────────────────
print("Fig D5: Daily ETc reduction LT vs Non-LT …")
daily_pivot = (daily_agg.pivot_table(index=["Location","water_mgmt"],
                                      columns="Genotype_class",
                                      values="mean_daily_ETc").reset_index())
for g in ["LT1.5","LT2","LT2.5"]:
    daily_pivot[f"red_{g}"] = daily_pivot["Non-LT"] - daily_pivot[g]
daily_pivot["Location"] = pd.Categorical(daily_pivot["Location"],
                                          categories=LOC_ORDER, ordered=True)

fig, axes = plt.subplots(1, 2, figsize=(15, 5), sharey=True)
for ax, mgmt in zip(axes, ["Irrigated","Rainfed"]):
    sub = daily_pivot[daily_pivot["water_mgmt"] == mgmt].sort_values("Location")
    locs_s = [l for l in LOC_ORDER if l in sub["Location"].values]
    x2 = np.arange(len(locs_s)); w3 = 0.22
    for i, g in enumerate(["LT1.5","LT2","LT2.5"]):
        vals = sub.set_index("Location").reindex(locs_s)[f"red_{g}"]
        ax.bar(x2 + (i-1)*w3, vals, w3, color=COL_GENO[g], label=g, alpha=0.85)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(x2); ax.set_xticklabels(locs_s, rotation=35, ha="right")
    ax.set_title(mgmt)
    ax.legend(title="LT intensity", fontsize=8, frameon=True)
    ax.axvline(3.5, color="grey", linewidth=0.7, linestyle="--")
    ax.axvline(6.5, color="grey", linewidth=0.7, linestyle="--")
    ax.grid(axis="y", alpha=0.3); ax.set_axisbelow(True)
axes[0].set_ylabel("Daily ETc reduction vs Non-LT (mm day⁻¹)")
plt.suptitle("Daily ETc reduction of LT genotypes relative to Non-LT\n"
             "(positive = water saved per day)", fontsize=12, y=1.01)
plt.tight_layout()
fig_d5 = fig_to_b64()

# ── Daily ETc summary table ───────────────────────────────────────────────────
print("Computing daily ETc summary table …")
daily_table_rows = []
for loc in locs_d:
    for mgmt in ["Irrigated","Rainfed"]:
        sub = dy[(dy["Location"] == loc) & (dy["water_mgmt"] == mgmt)]
        nlt  = sub[sub["Genotype_class"] == "Non-LT"]["ETc"].mean()
        for g in GENO_ORDER:
            sg = sub[sub["Genotype_class"] == g]["ETc"]
            if len(sg) == 0: continue
            m = sg.mean(); sd = sg.std()
            red = nlt - m if g != "Non-LT" else 0.0
            daily_table_rows.append({
                "Location": loc, "Management": mgmt, "Genotype": g,
                "n days": len(sg),
                "Mean daily ETc (mm)": f"{m:.3f}",
                "SD": f"{sd:.3f}",
                "Reduction vs Non-LT (mm day⁻¹)": f"{red:+.3f}" if g != "Non-LT" else "—",
            })
daily_table_df = pd.DataFrame(daily_table_rows)


# ── Build HTML ────────────────────────────────────────────────────────────────
print("Building HTML report …")

def df_to_html(df, idx=False):
    return df.to_html(index=idx, border=0, classes="tbl")

daily_table_html = df_to_html(daily_table_df)

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
<p>Generated: 2026-05-29 · Data: <code>modeling_data/data/updated-outputs/</code></p>

<h2>Experiment Design</h2>
<p>
The Sinclair Simplified Soybean Model (SSM) was used to simulate soybean growth and water use across
<b>10 locations in three US states</b> (Arkansas, Missouri, and Nebraska) for the 30-year period
1990–2019. Each simulation represents one combination of location, genotype, water management, and
planting date, run for each year independently.
</p>

<h3>Locations</h3>
<p>Locations were selected to span the major US soybean production belt, from the humid Mid-South to
the semi-arid Great Plains:</p>
<table class="tbl" style="width:auto">
<tr><th>State</th><th>Locations</th><th>Climate</th></tr>
<tr><td>Arkansas</td><td>Jonesboro, Keiser, Marianna, Rohwer</td><td>Humid subtropical; high VPD summers; drought risk during reproductive stages</td></tr>
<tr><td>Missouri</td><td>Albany, Mount Vernon, Novelty</td><td>Humid continental; moderate summer VPD; variable rainfall</td></tr>
<tr><td>Nebraska</td><td>Eustis, Lincoln, North Platte</td><td>Semi-arid to sub-humid; high summer VPD; irrigation-dependent production</td></tr>
</table>

<h3>Genotypes and the LT Trait</h3>
<p>
Four genotype classes were simulated using maturity groups matched to each region:
<b>MG 4</b> for Arkansas and Missouri, and <b>MG 3</b> for Nebraska — chosen to approximate
locally adapted commercial varieties and keep phenology comparable within each region:
</p>
<table class="tbl" style="width:auto">
<tr><th>Genotype class</th><th>SSM crop name</th><th>MG (AR &amp; MO / NE)</th><th>LT trait description</th><th>VPD threshold</th></tr>
<tr><td><b>Non-LT</b></td><td>SOYBEAN-MG4 / SOYBEAN-MG3</td><td>MG4 (AR, MO) / MG3 (NE)</td><td>Standard cultivar — transpiration increases linearly with VPD; no stomatal limit</td><td>None</td></tr>
<tr><td><b>LT1.5</b></td><td>SOYBEAN-LT1.5-MG4/3</td><td>MG4 (AR, MO) / MG3 (NE)</td><td>Mild limited transpiration — stomata restrict TR above ~1.5 kPa VPD</td><td>~1.5 kPa</td></tr>
<tr><td><b>LT2</b></td><td>SOYBEAN-LT2-MG4/3</td><td>MG4 (AR, MO) / MG3 (NE)</td><td>Moderate LT — TR limited above ~2.0 kPa (canonical PI 416937 threshold)</td><td>~2.0 kPa</td></tr>
<tr><td><b>LT2.5</b></td><td>SOYBEAN-LT2.5-MG4/3</td><td>MG4 (AR, MO) / MG3 (NE)</td><td>Strong LT — TR limited above ~2.5 kPa; most conservative water use</td><td>~2.5 kPa</td></tr>
</table>
<p>
The LT trait is encoded in the SSM by a VPD-threshold response: when daily VPD exceeds the
genotype-specific threshold, transpiration is capped at the rate reached at that VPD level,
conserving soil water. Above-threshold days are common in all three states during July–August,
making the LT trait agronomically relevant across the simulation domain.
</p>

<h3>Water Management</h3>
<p>
Each genotype was simulated under two water management regimes:
<ul>
<li><b>Rainfed (RFD)</b> — no supplemental water; crop relies entirely on rainfall and stored soil water.
    This represents the majority of US soybean production.</li>
<li><b>Irrigated (IRRI)</b> — water is supplied on demand to keep the fraction of transpirable soil
    water (FTSW) above stress thresholds. This isolates the atmospheric demand effect (VPD) from
    soil water limitation, providing a cleaner signal of the LT trait's transpiration response.</li>
</ul>
</p>

<h3>Planting Dates</h3>
<p>
Three planting dates were simulated per location to capture agronomic flexibility:
<b>Early (ELY)</b>, <b>Mid (MID)</b>, and <b>Late (LTE)</b>. Results reported here are averaged
across planting dates unless stated otherwise.
</p>

<h3>Output Variables</h3>
<table class="tbl" style="width:auto">
<tr><th>Variable</th><th>Description</th><th>Unit</th></tr>
<tr><td>SUMETT</td><td>Seasonal total crop evapotranspiration (ETc)</td><td>mm season⁻¹</td></tr>
<tr><td>Yield_kgha</td><td>Grain yield (WGRN × 10)</td><td>kg ha⁻¹</td></tr>
<tr><td>WP</td><td>Water productivity (Yield / SUMETT)</td><td>kg mm⁻¹</td></tr>
<tr><td>R7.PM</td><td>Days from planting to physiological maturity</td><td>days</td></tr>
<tr><td>ETc (daily)</td><td>Daily crop ET = TR (transpiration) + SEVP (soil evaporation)</td><td>mm day⁻¹</td></tr>
<tr><td>GLAI</td><td>Green leaf area index</td><td>m² m⁻²</td></tr>
<tr><td>DDMP</td><td>Daily dry matter production</td><td>g m⁻² day⁻¹</td></tr>
</table>
<p>Note: reference ET (ETo) used for Kc calculations is the Penman-Monteith ETo from NASA POWER,
not the SSM internal PET output. Seasonal Kc results are reported in the yearly analysis section.</p>

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

<h2>Average Daily Crop ET (ETc) Analysis</h2>
<p>
Daily ETc = transpiration (TR) + soil evaporation (SEVP), using growing-season days only (DAP &gt; 0).
All values are based on {len(dy):,} valid daily rows from the 10-location SSM simulation (1990–2019).
</p>

<h3>D1 — Mean daily ETc by location, management, and genotype</h3>
<p>Bars show the grand mean ± SE of all growing-season daily ETc values pooled across years and planting dates
for each location × genotype × management combination.</p>
<img src="data:image/png;base64,{fig_d1}" alt="Mean daily ETc by location">

<h3>D2 — Seasonal progression of daily ETc (5-day DAP bins)</h3>
<p>ETc averaged within consecutive 5-day bins from planting to physiological maturity (full cycle),
split by state and management. Dotted reference lines mark approximate vegetative→reproductive
transition (~DAP 50) and early reproductive→grain fill transition (~DAP 90).
LT genotypes diverge from Non-LT most strongly during mid-season when VPD peaks.
Season length varies by region: up to ~200 DAP in Arkansas, ~240 DAP in Missouri and Nebraska.</p>
<img src="data:image/png;base64,{fig_d2}" alt="Seasonal ETc progression">

<h3>D3 — Seasonal progression of Green Leaf Area Index (GLAI)</h3>
<p>Mean GLAI (m² m⁻²) within 5-day DAP bins by state, management, and genotype.
GLAI trajectories reflect canopy development and senescence; differences among LT classes indicate
the indirect effect of water conservation on leaf area maintenance during reproductive stages.</p>
<img src="data:image/png;base64,{fig_d3}" alt="Seasonal GLAI progression">

<h3>D4 — Seasonal progression of daily dry matter production (DDMP)</h3>
<p>Mean daily dry matter production (g m⁻² day⁻¹) within 5-day DAP bins by state, management, and genotype.
DDMP integrates radiation interception and radiation use efficiency; comparing LT and Non-LT reveals
whether water conservation translates into maintained or improved carbon assimilation.</p>
<img src="data:image/png;base64,{fig_d4}" alt="Seasonal DDMP progression">

<h3>D5 — Daily ETc reduction of LT genotypes relative to Non-LT</h3>
<p>Positive values indicate water saved per day relative to Non-LT.
LT2 and LT2.5 show the largest daily savings under irrigation in all three states.</p>
<img src="data:image/png;base64,{fig_d5}" alt="Daily ETc reduction">

<h3>D6 — Mean daily ETc summary table</h3>
<p>Complete numerical summary by location, management, and genotype.
Reduction column reports the absolute daily water saving relative to Non-LT (mm day⁻¹).</p>
{daily_table_html}

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
