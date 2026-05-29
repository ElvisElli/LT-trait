"""
Fig_marin_daily_ETc_ETo.py — Marin-style daily ETc vs ETo scatter plots.
Uses daily SSM ETc (TR + SEVP) on y-axis and NASA POWER PM ETo on x-axis.
PM ETo computed by compute_pm_eto_daily.py (FAO-56 Penman-Monteith, state-level sites).
Produces plots per location, per water management, and overall.
Run from /home/user/LT-trait:
    python3 modeling_data/scripts/Fig_marin_daily_ETc_ETo.py
"""

import os, warnings
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats as scipy_stats

warnings.filterwarnings("ignore")

DAILY  = "modeling_data/data/updated-outputs/all_daily.csv"
YEARLY = "modeling_data/data/updated-outputs/all_yearly.csv"
FIGDIR = "modeling_data/figures"
os.makedirs(FIGDIR, exist_ok=True)

COL_GENO = {
    "Non-LT": "#31a354",
    "LT1.5":  "#d95f02",
    "LT2":    "#7570b3",
    "LT2.5":  "#e7298a",
}
GENO_ORDER = ["Non-LT", "LT1.5", "LT2", "LT2.5"]
LOC_ORDER = [
    "Jonesboro","Keiser","Marianna","Rohwer",
    "Albany","Mount Vernon","Novelty",
    "Eustis","Lincoln","North Platte",
]

plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 10, "axes.labelsize": 10,
    "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 8,
    "figure.dpi": 150, "savefig.dpi": 150,
})

def save(name):
    path = os.path.join(FIGDIR, name + ".png")
    plt.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  saved {name}.png")

# ── Load and filter daily data ─────────────────────────────────────────────
print("Loading daily data (this may take a moment) …")
dy = pd.read_csv(DAILY, usecols=[
    "Location","Manag","Crop","Pyear","DOY","DAP",
    "TR","SEVP","Genotype_class","water_mgmt","State","Planting"
])

dy["Location"] = dy["Location"].str.strip().str.replace("­","",regex=False)

# ETc = transpiration + soil evaporation
dy["ETc"] = pd.to_numeric(dy["TR"],   errors="coerce") + \
            pd.to_numeric(dy["SEVP"], errors="coerce")

# Keep only growing-season days with valid ETc
dy = dy[(dy["DAP"] > 0) & dy["ETc"].notna() & (dy["ETc"] >= 0)]

# ── Join NASA POWER PM ETo (replaces SSM internal PET as reference ET) ────
ETO_DAILY = "modeling_data/data/nasapower_ETo_daily_py.csv"
eto = pd.read_csv(ETO_DAILY)
state_abbr = {"AR": "Arkansas", "MO": "Missouri", "NE": "Nebraska"}
eto["State"] = eto["State"].map(state_abbr)
eto = eto.rename(columns={"Year": "Pyear"})
dy = dy.merge(eto[["State","Pyear","DOY","ETo_mm"]], on=["State","Pyear","DOY"], how="left")
dy = dy[dy["ETo_mm"].notna() & (dy["ETo_mm"] >= 0)]

print(f"  {len(dy):,} valid daily rows (PM ETo joined)")

def kc_slope(x, y):
    """OLS slope forced through origin."""
    if len(x) < 10:
        return np.nan, np.nan, np.nan
    slope, intercept, r, p, se = scipy_stats.linregress(x, y)
    return slope, r**2, len(x)

def scatter_panel(ax, sub, color, label, alpha=0.03, s=1):
    ax.scatter(sub["ETo_mm"], sub["ETc"], c=color, alpha=alpha, s=s, rasterized=True)
    if len(sub) > 10:
        sl, r2, n = kc_slope(sub["ETo_mm"].values, sub["ETc"].values)
        xr = np.linspace(0, sub["ETo_mm"].quantile(0.99), 100)
        slope, intercept, *_ = scipy_stats.linregress(sub["ETo_mm"].values, sub["ETc"].values)
        yr = slope * xr + intercept
        ax.plot(xr, yr, color=color, linewidth=1.5,
                label=f"{label}  Kc={slope:.2f}  R²={r2:.2f}")
    ax.set_xlabel("PM ETo — NASA POWER (mm day⁻¹)")
    ax.set_ylabel("Crop ET (ETc, mm day⁻¹)")
    ax.set_xlim(left=0); ax.set_ylim(bottom=0)

# ─────────────────────────────────────────────────────────────────────────────
# 1. OVERALL: all locations combined, irrigated only
# ─────────────────────────────────────────────────────────────────────────────
print("Fig 1: Overall irrigated ETc vs ETo by genotype …")
fig, axes = plt.subplots(1, 4, figsize=(16, 4), sharey=True, sharex=True)
irri = dy[dy["water_mgmt"] == "Irrigated"]
for ax, g in zip(axes, GENO_ORDER):
    sub = irri[irri["Genotype_class"] == g]
    scatter_panel(ax, sub, COL_GENO[g], g)
    ax.set_title(g)
    ax.legend(fontsize=7, loc="upper left")
    ax.plot([0, ax.get_xlim()[1]], [0, ax.get_xlim()[1]], "k--", linewidth=0.5, alpha=0.4)
axes[1].set_ylabel("")
axes[2].set_ylabel("")
axes[3].set_ylabel("")
plt.suptitle("Daily ETc vs ETo — all locations, irrigated (Marin-style)", y=1.01)
plt.tight_layout()
save("Fig_marin_overall_irrigated")

# ─────────────────────────────────────────────────────────────────────────────
# 2. OVERALL: rainfed
# ─────────────────────────────────────────────────────────────────────────────
print("Fig 2: Overall rainfed ETc vs ETo by genotype …")
fig, axes = plt.subplots(1, 4, figsize=(16, 4), sharey=True, sharex=True)
rfed = dy[dy["water_mgmt"] == "Rainfed"]
for ax, g in zip(axes, GENO_ORDER):
    sub = rfed[rfed["Genotype_class"] == g]
    scatter_panel(ax, sub, COL_GENO[g], g)
    ax.set_title(g)
    ax.legend(fontsize=7, loc="upper left")
    ax.plot([0, ax.get_xlim()[1]], [0, ax.get_xlim()[1]], "k--", linewidth=0.5, alpha=0.4)
plt.suptitle("Daily ETc vs ETo — all locations, rainfed (Marin-style)", y=1.01)
plt.tight_layout()
save("Fig_marin_overall_rainfed")

# ─────────────────────────────────────────────────────────────────────────────
# 3. OVERALL: overlay all genotypes, irrigated — single panel per management
# ─────────────────────────────────────────────────────────────────────────────
print("Fig 3: Overall overlay by genotype …")
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for ax, (mgmt, sub_m) in zip(axes, [("Irrigated", irri), ("Rainfed", rfed)]):
    for g in GENO_ORDER:
        sub = sub_m[sub_m["Genotype_class"] == g]
        if len(sub) < 10: continue
        slope, intercept, r, p, se = scipy_stats.linregress(sub["ETo_mm"].values, sub["ETc"].values)
        xr = np.linspace(0, sub["ETo_mm"].quantile(0.99), 200)
        ax.scatter(sub["ETo_mm"], sub["ETc"], c=COL_GENO[g], alpha=0.015, s=1, rasterized=True)
        ax.plot(xr, slope*xr + intercept, color=COL_GENO[g], linewidth=2,
                label=f"{g}  Kc={slope:.2f} (R²={r**2:.2f})")
    ax.set_xlabel("PM ETo — NASA POWER (mm day⁻¹)")
    ax.set_ylabel("Crop ET (ETc, mm day⁻¹)")
    ax.set_title(mgmt)
    ax.set_xlim(left=0); ax.set_ylim(bottom=0)
    ax.legend(fontsize=8)
    ax.plot([0,8],[0,8], "k--", linewidth=0.5, alpha=0.3, label="1:1")
plt.tight_layout()
save("Fig_marin_overlay_geno")

# ─────────────────────────────────────────────────────────────────────────────
# 4. PER LOCATION: irrigated — faceted 2×5 grid
# ─────────────────────────────────────────────────────────────────────────────
print("Fig 4: Per-location irrigated ETc vs ETo …")
locs = [l for l in LOC_ORDER if l in dy["Location"].values]
fig, axes = plt.subplots(2, 5, figsize=(20, 8), sharex=True, sharey=True)
axes = axes.flatten()
for ax, loc in zip(axes, locs):
    sub_loc = irri[irri["Location"] == loc]
    for g in GENO_ORDER:
        sub = sub_loc[sub_loc["Genotype_class"] == g]
        if len(sub) < 10: continue
        slope, intercept, r, *_ = scipy_stats.linregress(sub["ETo_mm"].values, sub["ETc"].values)
        xr = np.linspace(0, sub["ETo_mm"].quantile(0.99), 100)
        ax.scatter(sub["ETo_mm"], sub["ETc"], c=COL_GENO[g], alpha=0.02, s=0.8, rasterized=True)
        ax.plot(xr, slope*xr + intercept, color=COL_GENO[g], linewidth=1.5,
                label=f"{g}: Kc={slope:.2f}")
    ax.set_title(loc, fontsize=9)
    ax.set_xlim(left=0); ax.set_ylim(bottom=0)
    ax.legend(fontsize=6, loc="upper left")
    ax.plot([0,8],[0,8],"k--",linewidth=0.4,alpha=0.3)
for ax in axes[len(locs):]: ax.set_visible(False)
for ax in axes[5:]: ax.set_xlabel("ETo (mm day⁻¹)")
for ax in axes[::5]: ax.set_ylabel("ETc (mm day⁻¹)")
plt.suptitle("Daily ETc vs ETo per location — irrigated", y=1.01)
plt.tight_layout()
save("Fig_marin_per_location_irrigated")

# ─────────────────────────────────────────────────────────────────────────────
# 5. PER LOCATION: rainfed
# ─────────────────────────────────────────────────────────────────────────────
print("Fig 5: Per-location rainfed ETc vs ETo …")
fig, axes = plt.subplots(2, 5, figsize=(20, 8), sharex=True, sharey=True)
axes = axes.flatten()
for ax, loc in zip(axes, locs):
    sub_loc = rfed[rfed["Location"] == loc]
    for g in GENO_ORDER:
        sub = sub_loc[sub_loc["Genotype_class"] == g]
        if len(sub) < 10: continue
        slope, intercept, r, *_ = scipy_stats.linregress(sub["ETo_mm"].values, sub["ETc"].values)
        xr = np.linspace(0, sub["ETo_mm"].quantile(0.99), 100)
        ax.scatter(sub["ETo_mm"], sub["ETc"], c=COL_GENO[g], alpha=0.02, s=0.8, rasterized=True)
        ax.plot(xr, slope*xr + intercept, color=COL_GENO[g], linewidth=1.5,
                label=f"{g}: Kc={slope:.2f}")
    ax.set_title(loc, fontsize=9)
    ax.set_xlim(left=0); ax.set_ylim(bottom=0)
    ax.legend(fontsize=6, loc="upper left")
    ax.plot([0,8],[0,8],"k--",linewidth=0.4,alpha=0.3)
for ax in axes[len(locs):]: ax.set_visible(False)
for ax in axes[5:]: ax.set_xlabel("ETo (mm day⁻¹)")
for ax in axes[::5]: ax.set_ylabel("ETc (mm day⁻¹)")
plt.suptitle("Daily ETc vs ETo per location — rainfed", y=1.01)
plt.tight_layout()
save("Fig_marin_per_location_rainfed")

# ─────────────────────────────────────────────────────────────────────────────
# 6. PER LOCATION AND GENOTYPE: Kc heatmap
# ─────────────────────────────────────────────────────────────────────────────
print("Fig 6: Kc heatmap by location × genotype × management …")
kc_rows = []
for loc in locs:
    for mgmt in ["Irrigated","Rainfed"]:
        sub_m = dy[(dy["Location"]==loc) & (dy["water_mgmt"]==mgmt)]
        for g in GENO_ORDER:
            sub = sub_m[sub_m["Genotype_class"]==g]
            if len(sub) < 10: continue
            sl, intercept, r, p, se = scipy_stats.linregress(sub["ETo_mm"].values, sub["ETc"].values)
            kc_rows.append({"Location":loc,"Genotype":g,"Management":mgmt,
                            "Kc":sl,"R2":r**2,"n":len(sub)})
kc_df = pd.DataFrame(kc_rows)

fig, axes = plt.subplots(1, 2, figsize=(14, 5))
for ax, mgmt in zip(axes, ["Irrigated","Rainfed"]):
    sub = kc_df[kc_df["Management"]==mgmt].pivot(index="Location", columns="Genotype", values="Kc")
    sub = sub.reindex(locs)[GENO_ORDER]
    im = ax.imshow(sub.values, cmap="RdYlGn_r", aspect="auto", vmin=0.4, vmax=1.4)
    ax.set_xticks(range(len(GENO_ORDER))); ax.set_xticklabels(GENO_ORDER, rotation=20)
    ax.set_yticks(range(len(sub))); ax.set_yticklabels(sub.index, fontsize=9)
    ax.set_title(mgmt)
    for i in range(len(sub)):
        for j in range(len(GENO_ORDER)):
            v = sub.values[i,j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=8,
                        color="white" if v > 1.1 or v < 0.6 else "black")
    plt.colorbar(im, ax=ax, label="Kc (slope ETc/ETo)")
plt.suptitle("Daily Kc (ETc/ETo slope) by location × genotype")
plt.tight_layout()
save("Fig_marin_kc_heatmap")

# ─────────────────────────────────────────────────────────────────────────────
# 7. STATE-LEVEL comparison: overlay regressions
# ─────────────────────────────────────────────────────────────────────────────
print("Fig 7: State-level overlay …")
fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharex=True, sharey=True)
states = ["Arkansas","Missouri","Nebraska"]
for ax, state in zip(axes, states):
    sub_s = irri[irri["State"]==state]
    for g in GENO_ORDER:
        sub = sub_s[sub_s["Genotype_class"]==g]
        if len(sub) < 10: continue
        slope, intercept, r, *_ = scipy_stats.linregress(sub["ETo_mm"].values, sub["ETc"].values)
        xr = np.linspace(0, sub["ETo_mm"].quantile(0.99), 200)
        ax.scatter(sub["ETo_mm"], sub["ETc"], c=COL_GENO[g], alpha=0.015, s=1, rasterized=True)
        ax.plot(xr, slope*xr+intercept, color=COL_GENO[g], linewidth=2,
                label=f"{g}: Kc={slope:.2f}")
    ax.set_title(state); ax.set_xlabel("ETo (mm day⁻¹)")
    ax.set_xlim(left=0); ax.set_ylim(bottom=0)
    ax.legend(fontsize=8)
    ax.plot([0,8],[0,8],"k--",linewidth=0.5,alpha=0.3)
axes[0].set_ylabel("ETc (mm day⁻¹)")
plt.suptitle("Daily ETc vs ETo by state — irrigated", y=1.01)
plt.tight_layout()
save("Fig_marin_by_state_irrigated")

# Save Kc summary table
kc_df.to_csv("modeling_data/data/updated-outputs/kc_daily_summary.csv", index=False)
print("Kc summary saved: modeling_data/data/updated-outputs/kc_daily_summary.csv")
print("\nDone.")
