"""
01_data_inspection_QC.py
Growth Chamber 2026 — Soybean Limited-Transpiration Experiment

Purpose : Load data_all.xlsx, inspect structure, flag QC issues,
          produce 8 diagnostic figures, save data_clean.csv.

Run from project root (experimental_data_soy_2026/):
    python script/01_data_inspection_QC.py
"""

import warnings
warnings.filterwarnings("ignore")

import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# --------------------------------------------------------------------------- #
# PATHS
# --------------------------------------------------------------------------- #
ROOT   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA   = os.path.join(ROOT, "data")
FIGS   = os.path.join(ROOT, "figures")
os.makedirs(FIGS, exist_ok=True)

# --------------------------------------------------------------------------- #
# STYLE — University of Arkansas theme
# --------------------------------------------------------------------------- #
UA_RED    = "#9D2235"
UA_DARK   = "#1A1A1A"
UA_GRAY   = "#888888"
UA_LIGHT  = "#F5F5F5"

plt.rcParams.update({
    "font.family"       : "DejaVu Sans",
    "font.size"         : 9,
    "axes.spines.top"   : False,
    "axes.spines.right" : False,
    "axes.linewidth"    : 0.8,
    "xtick.major.width" : 0.8,
    "ytick.major.width" : 0.8,
    "figure.facecolor"  : "white",
    "axes.facecolor"    : "white",
    "savefig.dpi"       : 300,
    "savefig.bbox"      : "tight",
    "savefig.facecolor" : "white",
})

# WUE-class colour palette
WUE_COLS = {
    "slow_wilting"        : "#9D2235",
    "high_wilting"        : "#E8002D",
    "Plastic, High WUE"   : "#4393c3",
    "Nonplastic, Low WUE" : "#f4a582",
    "MAGIC parent"        : "#555555",
}
DAY_COLS = {"Day 1": "#9D2235", "Day 2": "#4393c3"}

# --------------------------------------------------------------------------- #
# 1. LOAD DATA
# --------------------------------------------------------------------------- #
print("Loading data …")
raw = pd.read_excel(os.path.join(DATA, "data_all.xlsx"), sheet_name="dataAll")
print(f"  Raw shape: {raw.shape}")

# Drop trailing unnamed column(s)
raw = raw.loc[:, ~raw.columns.str.startswith("Unnamed")]

# --------------------------------------------------------------------------- #
# 2. RENAME / RECODE
# --------------------------------------------------------------------------- #
rename = {
    "RepNo"              : "Rep",
    "day"                : "Day",
    "TR_mg_sqm_seg"      : "TR",
    "RelativeTransp"     : "RelTR",
    "avg_vpd"            : "VPD",
    "avg_temperature"    : "Temp_C",
    "avg_rh"             : "RH_pct",
    "LeafArea_squarecm"  : "LeafArea",
    "treatment_time_mins": "duration_min",
}
df = raw.rename(columns=rename).copy()

# Average stomatal conductance
df["gsw"] = df[["gsw1","gsw2","gsw3"]].mean(axis=1)
# Flag where ALL gsw readings are NaN (not measured)
all_gsw_nan = df[["gsw1","gsw2","gsw3"]].isna().all(axis=1)
df.loc[all_gsw_nan, "gsw"] = np.nan

# Encode factors
df["Rep"]      = df["Rep"].astype(str)
df["Day_label"]= df["Day"].map({1: "Day 1", 2: "Day 2"})
df["Genotype"] = df["Genotype"].astype(str)
WUE_ORDER = ["slow_wilting","high_wilting","Plastic, High WUE",
             "Nonplastic, Low WUE","MAGIC parent"]
df["WUE_class"] = pd.Categorical(df["WUE_class"], categories=WUE_ORDER, ordered=True)
df["Treatment"] = df["Treatment"].astype(int)

print(f"  Columns after renaming: {list(df.columns)}")

# --------------------------------------------------------------------------- #
# 3. MISSING VALUE REPORT
# --------------------------------------------------------------------------- #
print("\n--- Missing values ---")
miss = df.isnull().sum()
print(miss[miss > 0].to_string())

print("\n--- gsw availability by Day ---")
gsw_avail = (df.groupby("Day_label")
               .apply(lambda x: pd.Series({
                   "n_total" : len(x),
                   "n_gsw"   : x["gsw"].notna().sum(),
                   "pct_gsw" : f'{100*x["gsw"].notna().mean():.1f}%'
               }))
               .reset_index())
print(gsw_avail.to_string(index=False))

# --------------------------------------------------------------------------- #
# 4. SAMPLE DESIGN SUMMARY
# --------------------------------------------------------------------------- #
print("\n--- Observations per Genotype × Day ---")
pivot = (df.groupby(["Genotype","Day_label"])
           .size()
           .unstack(fill_value=0)
           .reset_index())
print(pivot.to_string(index=False))

print("\n--- VPD achieved per Treatment ---")
vpd_sum = (df.groupby("Treatment")["VPD"]
             .agg(n="count", mean=np.mean, min=np.min, max=np.max)
             .round(3)
             .reset_index())
print(vpd_sum.to_string(index=False))

# --------------------------------------------------------------------------- #
# 5. OUTLIER FLAGS
# --------------------------------------------------------------------------- #
print("\n--- RelTR > 1 ---")
print(df.loc[df["RelTR"] > 1, ["Rep","Day_label","Genotype","Treatment","VPD","TR","RelTR"]]
        .to_string(index=False))

# IQR flag (3× rule within Genotype × Day × Treatment)
def compute_outlier(series):
    q1, q3 = series.quantile(0.25), series.quantile(0.75)
    iqr = q3 - q1
    return (series < q1 - 3*iqr) | (series > q3 + 3*iqr)

df["TR_outlier"] = (df.groupby(["Genotype","Day_label","Treatment"])["TR"]
                      .transform(compute_outlier))

print("\n--- TR outliers (3×IQR) ---")
out_rows = df[df["TR_outlier"] == True][["Rep","Day_label","Genotype","Treatment","VPD","TR","RelTR"]]
if len(out_rows):
    print(out_rows.to_string(index=False))
else:
    print("  None detected")

# --------------------------------------------------------------------------- #
# 6. SAVE CLEAN DATA
# --------------------------------------------------------------------------- #
df.to_csv(os.path.join(DATA, "data_clean.csv"), index=False)
print(f"\nSaved data/data_clean.csv ({len(df)} rows)")

# --------------------------------------------------------------------------- #
# 7. FIGURES
# --------------------------------------------------------------------------- #

genotypes = sorted(df["Genotype"].unique())
n_geno    = len(genotypes)
NCOLS     = 5
NROWS     = int(np.ceil(n_geno / NCOLS))


def strip_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# --- QC1: TR distribution by Genotype × Day ----------------------------------
fig, ax = plt.subplots(figsize=(13, 4))
positions = np.arange(n_geno)
w = 0.32
for di, (dlabel, dc) in enumerate(DAY_COLS.items()):
    sub = df[df["Day_label"] == dlabel]
    data_by_geno = [sub[sub["Genotype"] == g]["TR"].dropna().values for g in genotypes]
    bp = ax.boxplot(data_by_geno, positions=positions + (di - 0.5)*w, widths=w*0.85,
                    patch_artist=True, sym=".",
                    boxprops=dict(facecolor=dc, alpha=0.7, linewidth=0.8),
                    medianprops=dict(color="white", linewidth=1.5),
                    whiskerprops=dict(linewidth=0.8),
                    capprops=dict(linewidth=0.8),
                    flierprops=dict(marker=".", markersize=3, color=dc, alpha=0.5))

ax.set_xticks(positions)
ax.set_xticklabels(genotypes, rotation=40, ha="right", fontsize=8)
ax.set_ylabel("TR (mg m⁻² s⁻¹)")
ax.set_title("QC1  Transpiration rate distribution by genotype and day", fontweight="bold")
legend_handles = [Patch(facecolor=c, label=d, alpha=0.7) for d, c in DAY_COLS.items()]
ax.legend(handles=legend_handles, frameon=False, fontsize=8)
strip_axes(ax)
fig.tight_layout()
fig.savefig(os.path.join(FIGS, "QC1_TR_distribution.png"))
plt.close(fig)
print("Saved QC1")


# --- QC2: TR vs VPD by genotype, both days -----------------------------------
fig, axes = plt.subplots(NROWS, NCOLS, figsize=(14, NROWS*2.8), sharey=False)
axes_flat = axes.flatten()
for i, g in enumerate(genotypes):
    ax = axes_flat[i]
    for dlabel, dc in DAY_COLS.items():
        sub = df[(df["Genotype"] == g) & (df["Day_label"] == dlabel)]
        ax.scatter(sub["VPD"], sub["TR"], c=dc, s=20, alpha=0.7, label=dlabel)
        if len(sub) > 3:
            z = np.polyfit(sub["VPD"].dropna(), sub["TR"].dropna(), 1)
            xr = np.linspace(sub["VPD"].min(), sub["VPD"].max(), 80)
            ax.plot(xr, np.polyval(z, xr), c=dc, lw=1.2, alpha=0.8)
    ax.set_title(g, fontsize=7, fontweight="bold")
    ax.set_xlabel("VPD (kPa)", fontsize=7)
    ax.set_ylabel("TR", fontsize=7)
    strip_axes(ax)

for j in range(i+1, len(axes_flat)):
    axes_flat[j].set_visible(False)

handles = [Line2D([0],[0], color=c, lw=2, label=d) for d, c in DAY_COLS.items()]
fig.legend(handles=handles, loc="lower right", frameon=False, fontsize=8)
fig.suptitle("QC2  TR vs VPD by genotype (both days)", fontweight="bold", y=1.01)
fig.tight_layout()
fig.savefig(os.path.join(FIGS, "QC2_TR_vs_VPD_by_genotype.png"))
plt.close(fig)
print("Saved QC2")


# --- QC3: RelTR vs VPD, Day 1 ------------------------------------------------
day1 = df[df["Day_label"] == "Day 1"].copy()
fig, axes = plt.subplots(NROWS, NCOLS, figsize=(14, NROWS*2.8), sharey=True)
axes_flat = axes.flatten()
for i, g in enumerate(genotypes):
    ax = axes_flat[i]
    sub = day1[day1["Genotype"] == g]
    wc  = sub["WUE_class"].iloc[0] if len(sub) else "MAGIC parent"
    col = WUE_COLS.get(str(wc), UA_GRAY)
    ax.axhline(1.0, lw=0.8, ls="--", color=UA_GRAY)
    ax.scatter(sub["VPD"], sub["RelTR"], c=col, s=22, alpha=0.8, zorder=3)
    ax.set_title(g, fontsize=7, fontweight="bold")
    ax.set_xlabel("VPD (kPa)", fontsize=7)
    ax.set_ylabel("RelTR", fontsize=7)
    strip_axes(ax)

for j in range(i+1, len(axes_flat)):
    axes_flat[j].set_visible(False)

legend_handles = [Patch(facecolor=c, label=wc) for wc, c in WUE_COLS.items()]
fig.legend(handles=legend_handles, loc="lower right", frameon=False, fontsize=7, ncol=2)
fig.suptitle("QC3  Relative transpiration vs VPD — Day 1", fontweight="bold", y=1.01)
fig.tight_layout()
fig.savefig(os.path.join(FIGS, "QC3_RelTR_vs_VPD_Day1.png"))
plt.close(fig)
print("Saved QC3")


# --- QC4: gsw distribution ---------------------------------------------------
gsw_df = df[df["gsw"].notna()].copy()
fig, ax = plt.subplots(figsize=(10, 4))
geno_gsw = sorted(gsw_df["Genotype"].unique())
positions = np.arange(len(geno_gsw))
w = 0.32
for di, (dlabel, dc) in enumerate(DAY_COLS.items()):
    sub = gsw_df[gsw_df["Day_label"] == dlabel]
    data_by_geno = [sub[sub["Genotype"] == g]["gsw"].dropna().values for g in geno_gsw]
    ax.boxplot(data_by_geno, positions=positions + (di - 0.5)*w, widths=w*0.85,
               patch_artist=True, sym=".",
               boxprops=dict(facecolor=dc, alpha=0.7, linewidth=0.8),
               medianprops=dict(color="white", linewidth=1.5),
               whiskerprops=dict(linewidth=0.8),
               capprops=dict(linewidth=0.8),
               flierprops=dict(marker=".", markersize=3, color=dc, alpha=0.5))
ax.set_xticks(positions)
ax.set_xticklabels(geno_gsw, rotation=40, ha="right", fontsize=8)
ax.set_ylabel("g$_{sw}$ (mol m⁻² s⁻¹)")
ax.set_title("QC4  Stomatal conductance distribution", fontweight="bold")
legend_handles = [Patch(facecolor=c, label=d, alpha=0.7) for d, c in DAY_COLS.items()]
ax.legend(handles=legend_handles, frameon=False, fontsize=8)
strip_axes(ax)
fig.tight_layout()
fig.savefig(os.path.join(FIGS, "QC4_gsw_distribution.png"))
plt.close(fig)
print("Saved QC4")


# --- QC5: gsw vs VPD, Day 1 --------------------------------------------------
day1_gsw = day1[day1["gsw"].notna()].copy()
geno_gsw_d1 = sorted(day1_gsw["Genotype"].unique())
nc = min(5, len(geno_gsw_d1))
nr = int(np.ceil(len(geno_gsw_d1) / nc))
fig, axes = plt.subplots(nr, nc, figsize=(nc*2.6, nr*2.8), sharey=True)
axes_flat = np.array(axes).flatten()
for i, g in enumerate(geno_gsw_d1):
    ax = axes_flat[i]
    sub = day1_gsw[day1_gsw["Genotype"] == g]
    wc  = sub["WUE_class"].iloc[0]
    col = WUE_COLS.get(str(wc), UA_GRAY)
    ax.scatter(sub["VPD"], sub["gsw"], c=col, s=35, alpha=0.85, zorder=3)
    if len(sub) > 2:
        try:
            z = np.polyfit(sub["VPD"].values, sub["gsw"].values, 1)
            xr = np.linspace(sub["VPD"].min(), sub["VPD"].max(), 80)
            ax.plot(xr, np.polyval(z, xr), c=col, lw=1.2)
        except Exception:
            pass
    ax.set_title(g, fontsize=7, fontweight="bold")
    ax.set_xlabel("VPD (kPa)", fontsize=7)
    ax.set_ylabel("g$_{sw}$", fontsize=7)
    strip_axes(ax)
for j in range(i+1, len(axes_flat)):
    axes_flat[j].set_visible(False)
legend_handles = [Patch(facecolor=c, label=wc) for wc, c in WUE_COLS.items()]
fig.legend(handles=legend_handles, loc="lower right", frameon=False, fontsize=7, ncol=2)
fig.suptitle("QC5  Stomatal conductance vs VPD — Day 1 (LICOR)", fontweight="bold", y=1.01)
fig.tight_layout()
fig.savefig(os.path.join(FIGS, "QC5_gsw_vs_VPD_Day1.png"))
plt.close(fig)
print("Saved QC5")


# --- QC6: Leaf area by Rep × Genotype ----------------------------------------
fig, axes = plt.subplots(NROWS, NCOLS, figsize=(14, NROWS*2.8))
axes_flat = axes.flatten()
rep_cols = {"1":"#9D2235","2":"#4393c3","3":"#2ca25f","4":"#f4a582"}
for i, g in enumerate(genotypes):
    ax = axes_flat[i]
    sub = df[df["Genotype"] == g]
    for r in sorted(sub["Rep"].unique()):
        vals = sub[sub["Rep"] == r]["LeafArea"].dropna()
        ax.scatter([float(r)]*len(vals), vals, c=rep_cols.get(r, UA_GRAY),
                   s=20, alpha=0.7)
    ax.set_xticks([1,2,3,4])
    ax.set_xticklabels(["R1","R2","R3","R4"], fontsize=7)
    ax.set_title(g, fontsize=7, fontweight="bold")
    ax.set_ylabel("Leaf area (cm²)", fontsize=7)
    strip_axes(ax)
for j in range(i+1, len(axes_flat)):
    axes_flat[j].set_visible(False)
fig.suptitle("QC6  Leaf area by replicate and genotype", fontweight="bold", y=1.01)
fig.tight_layout()
fig.savefig(os.path.join(FIGS, "QC6_leaf_area_by_rep_genotype.png"))
plt.close(fig)
print("Saved QC6")


# --- QC7: Chamber conditions -------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(11, 4))
for ax, (col, lbl) in zip(axes, [("Temp_C","Temperature (°C)"),
                                   ("RH_pct","Relative humidity (%)"),
                                   ("VPD","VPD (kPa)")]):
    treats = sorted(df["Treatment"].unique())
    data_by_treat = [df[df["Treatment"]==t][col].dropna().values for t in treats]
    ax.boxplot(data_by_treat, patch_artist=True,
               boxprops=dict(facecolor=UA_RED, alpha=0.6, linewidth=0.8),
               medianprops=dict(color="white", linewidth=1.5),
               whiskerprops=dict(linewidth=0.8),
               capprops=dict(linewidth=0.8),
               flierprops=dict(marker=".", markersize=3, color=UA_RED))
    ax.set_xticklabels([f"T{t}" for t in treats])
    ax.set_xlabel("Treatment")
    ax.set_ylabel(lbl)
    strip_axes(ax)
fig.suptitle("QC7  Chamber conditions by treatment", fontweight="bold")
fig.tight_layout()
fig.savefig(os.path.join(FIGS, "QC7_chamber_conditions.png"))
plt.close(fig)
print("Saved QC7")


# --- QC8: Day 1 vs Day 2 TR --------------------------------------------------
pivot_days = (df.groupby(["Rep","Genotype","Treatment"])["TR"]
                .mean()
                .reset_index())
pivot_days2 = (df.groupby(["Rep","Genotype","Treatment","Day_label"])["TR"]
                 .mean()
                 .unstack("Day_label")
                 .reset_index())

fig, ax = plt.subplots(figsize=(6, 5))
for wc, col in WUE_COLS.items():
    genos_wc = df[df["WUE_class"] == wc]["Genotype"].unique()
    sub = pivot_days2[pivot_days2["Genotype"].isin(genos_wc)]
    if "Day 1" in sub and "Day 2" in sub:
        ax.scatter(sub["Day 1"], sub["Day 2"], c=col, s=20, alpha=0.6, label=wc)

lim_max = max(df["TR"].max(), df["TR"].max()) * 1.05
ax.plot([0, lim_max], [0, lim_max], "k--", lw=0.8, alpha=0.5)
ax.set_xlabel("TR Day 1 (mg m⁻² s⁻¹)")
ax.set_ylabel("TR Day 2 (mg m⁻² s⁻¹)")
ax.set_title("QC8  TR Day 1 vs Day 2 (diagonal = no day effect)", fontweight="bold")
ax.legend(frameon=False, fontsize=7, loc="upper left")
strip_axes(ax)
fig.tight_layout()
fig.savefig(os.path.join(FIGS, "QC8_TR_day1_vs_day2.png"))
plt.close(fig)
print("Saved QC8")

print("\nQC complete — all figures in figures/")
print(f"Clean data: data/data_clean.csv ({len(df)} rows)\n")
