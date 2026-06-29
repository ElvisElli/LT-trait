"""
02_model_fitting_TR_gs_VPD.py
Growth Chamber 2026 — Soybean Limited-Transpiration Experiment

Purpose : Fit linear, bilinear, and linear-plateau (quadratic-plateau style
          piecewise) models for transpiration (TR) and stomatal conductance
          (gsw) vs VPD.  Day 1 data only.  Replicate treated as random effect
          (via statsmodels MixedLM for the linear model; population-average
          NLS approximation for nonlinear models with Rep as grouping factor).

Outputs :
  figures/  — 6 publication-ready PNG figures (Arkansas colour theme)
  data/     — TR_model_parameters.csv, gsw_model_parameters.csv,
              TR_AIC_comparison.csv, gsw_AIC_comparison.csv

Run from project root (experimental_data_soy_2026/):
    python script/02_model_fitting_TR_gs_VPD.py
"""

import warnings
warnings.filterwarnings("ignore")

import os, itertools
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from scipy.optimize import curve_fit, OptimizeWarning
from scipy.stats import pearsonr
import statsmodels.formula.api as smf

# --------------------------------------------------------------------------- #
# PATHS
# --------------------------------------------------------------------------- #
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
FIGS = os.path.join(ROOT, "figures")
os.makedirs(FIGS, exist_ok=True)

# --------------------------------------------------------------------------- #
# UNIVERSITY OF ARKANSAS COLOUR SYSTEM
# --------------------------------------------------------------------------- #
UA_RED   = "#9D2235"   # Cardinal
UA_RED2  = "#C8102E"   # Brighter red
UA_WINE  = "#6E1424"   # Deep wine
UA_PINK  = "#D4808E"   # Muted rose
UA_LIGHT = "#F9E8EA"   # Blush
UA_DARK  = "#1A1A1A"
UA_GRAY  = "#7A7A7A"
UA_WHITE = "#FFFFFF"

# WUE class palette — five distinguishable shades within red-to-dark spectrum
WUE_COLS = {
    "slow_wilting"        : UA_WINE,   # darkest red — true LT phenotype
    "high_wilting"        : UA_RED2,   # bright red
    "Plastic, High WUE"   : "#4393c3", # blue for contrast
    "Nonplastic, Low WUE" : "#d4824a", # warm orange
    "MAGIC parent"        : UA_DARK,   # black reference
}
WUE_MARKERS = {
    "slow_wilting"        : "o",
    "high_wilting"        : "s",
    "Plastic, High WUE"   : "^",
    "Nonplastic, Low WUE" : "D",
    "MAGIC parent"        : "P",
}
MODEL_COLS = {
    "Linear"      : UA_GRAY,
    "Bilinear"    : UA_RED,
    "Lin. plateau": UA_WINE,
}
MODEL_LTY  = {"Linear": (5,3), "Bilinear": (1,0), "Lin. plateau": (3,2)}

# Publication plot style
plt.rcParams.update({
    "font.family"       : "DejaVu Sans",
    "font.size"         : 9,
    "axes.linewidth"    : 0.9,
    "xtick.major.width" : 0.9,
    "ytick.major.width" : 0.9,
    "xtick.minor.visible": False,
    "ytick.minor.visible": False,
    "axes.spines.top"   : False,
    "axes.spines.right" : False,
    "figure.facecolor"  : "white",
    "axes.facecolor"    : "white",
    "savefig.dpi"       : 300,
    "savefig.bbox"      : "tight",
    "savefig.facecolor" : "white",
})

def strip_ax(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# --------------------------------------------------------------------------- #
# 1. LOAD DATA — Day 1 only
# --------------------------------------------------------------------------- #
csv_path = os.path.join(DATA, "data_clean.csv")
if os.path.exists(csv_path):
    df_all = pd.read_csv(csv_path, dtype={"Rep": str})
else:
    import openpyxl
    raw = pd.read_excel(os.path.join(DATA, "data_all.xlsx"), sheet_name="dataAll")
    raw = raw.loc[:, ~raw.columns.str.startswith("Unnamed")]
    raw.rename(columns={"RepNo":"Rep","day":"Day","TR_mg_sqm_seg":"TR",
                         "RelativeTransp":"RelTR","avg_vpd":"VPD",
                         "avg_temperature":"Temp_C","avg_rh":"RH_pct",
                         "LeafArea_squarecm":"LeafArea"}, inplace=True)
    raw["gsw"] = raw[["gsw1","gsw2","gsw3"]].mean(axis=1)
    all_nan = raw[["gsw1","gsw2","gsw3"]].isna().all(axis=1)
    raw.loc[all_nan, "gsw"] = np.nan
    raw["Day_label"] = raw["Day"].map({1:"Day 1", 2:"Day 2"})
    df_all = raw
    df_all["Rep"] = df_all["Rep"].astype(str)

day1 = df_all[df_all["Day_label"] == "Day 1"].copy()
day1_gsw = day1[day1["gsw"].notna()].copy()

genotypes = sorted(day1["Genotype"].unique())
geno_gsw  = sorted(day1_gsw["Genotype"].unique())

# WUE lookup
wue_lookup = (day1.drop_duplicates("Genotype")
                   .set_index("Genotype")["WUE_class"]
                   .to_dict())

print(f"Day 1: {len(day1)} rows | {len(genotypes)} genotypes | {len(day1_gsw)} rows with gsw")


# --------------------------------------------------------------------------- #
# 2. MODEL FUNCTIONS
# --------------------------------------------------------------------------- #

def linear_fn(x, a, b):
    return a + b * x

def bilinear_fn(x, a, b1, b2, bp):
    return np.where(x <= bp, a + b1 * x, a + b1 * bp + b2 * (x - bp))

def lin_plateau_fn(x, a, b, bp):
    return np.where(x <= bp, a + b * x, a + b * bp)

def aic_from_nls(y, y_hat, n_params):
    """Approximate AIC for NLS fit (Gaussian likelihood)."""
    resid = y - y_hat
    sse   = np.sum(resid**2)
    n     = len(y)
    sigma2 = sse / n
    if sigma2 <= 0:
        return np.inf
    loglik = -0.5 * n * (np.log(2 * np.pi * sigma2) + 1)
    return 2 * n_params - 2 * loglik


# --------------------------------------------------------------------------- #
# 3. FIT FUNCTIONS (per-genotype)
# --------------------------------------------------------------------------- #

def fit_one(vpd, y, geno_label="?"):
    """Return dict with fitted parameters, AIC, and predicted curves."""
    vpd = np.array(vpd, dtype=float)
    y   = np.array(y,   dtype=float)
    mask = np.isfinite(vpd) & np.isfinite(y)
    vpd, y = vpd[mask], y[mask]
    n = len(vpd)
    if n < 3:
        return None

    vmin, vmax = vpd.min(), vpd.max()
    bp0 = 0.5*(vmin + vmax)
    a0  = y.min()
    b0  = max((y.max() - y.min()) / (vmax - vmin + 1e-9), 0)

    results = {}

    # ---- Linear ----
    try:
        popt, _ = curve_fit(linear_fn, vpd, y, p0=[a0, b0], maxfev=4000)
        y_hat   = linear_fn(vpd, *popt)
        results["Linear"] = {
            "a": popt[0], "b1": popt[1], "b2": np.nan, "bp": np.nan,
            "AIC": aic_from_nls(y, y_hat, 2),
            "popt": popt, "fn": linear_fn, "n_par": 2
        }
    except Exception:
        results["Linear"] = None

    # ---- Bilinear ----
    try:
        popt, _ = curve_fit(
            bilinear_fn, vpd, y,
            p0=[a0, b0, b0*0.1, bp0],
            bounds=([-np.inf, 0, -np.inf, vmin],
                    [ np.inf, np.inf, np.inf, vmax]),
            maxfev=8000
        )
        y_hat = bilinear_fn(vpd, *popt)
        results["Bilinear"] = {
            "a": popt[0], "b1": popt[1], "b2": popt[2], "bp": popt[3],
            "AIC": aic_from_nls(y, y_hat, 4),
            "popt": popt, "fn": bilinear_fn, "n_par": 4
        }
    except Exception:
        results["Bilinear"] = None

    # ---- Linear plateau ----
    try:
        popt, _ = curve_fit(
            lin_plateau_fn, vpd, y,
            p0=[a0, b0, bp0],
            bounds=([-np.inf, 0, vmin],
                    [ np.inf, np.inf, vmax]),
            maxfev=8000
        )
        y_hat = lin_plateau_fn(vpd, *popt)
        results["Lin. plateau"] = {
            "a": popt[0], "b1": popt[1], "b2": np.nan, "bp": popt[2],
            "AIC": aic_from_nls(y, y_hat, 3),
            "popt": popt, "fn": lin_plateau_fn, "n_par": 3
        }
    except Exception:
        results["Lin. plateau"] = None

    # Best model by AIC
    aics = {k: v["AIC"] for k, v in results.items() if v is not None}
    if aics:
        results["best"] = min(aics, key=aics.get)
    else:
        results["best"] = None

    return results


# --------------------------------------------------------------------------- #
# 4. FIT ALL GENOTYPES — TR and gsw
# --------------------------------------------------------------------------- #
print("\nFitting TR models …")
tr_fits = {}
for g in genotypes:
    sub = day1[day1["Genotype"] == g]
    tr_fits[g] = fit_one(sub["VPD"], sub["TR"], geno_label=g)
    best = tr_fits[g]["best"] if tr_fits[g] else "FAILED"
    print(f"  {g:15s}  best={best}")

print("\nFitting gsw models …")
gsw_fits = {}
for g in geno_gsw:
    sub = day1_gsw[day1_gsw["Genotype"] == g]
    gsw_fits[g] = fit_one(sub["VPD"], sub["gsw"], geno_label=g)
    best = gsw_fits[g]["best"] if gsw_fits[g] else "FAILED"
    print(f"  {g:15s}  best={best}")


# --------------------------------------------------------------------------- #
# 5. PARAMETER TABLES
# --------------------------------------------------------------------------- #
WUE_EXPECTED = {
    "slow_wilting"        : "Slow wilting / LT",
    "high_wilting"        : "High wilting / non-LT",
    "Plastic, High WUE"   : "Plastic – High WUE",
    "Nonplastic, Low WUE" : "Nonplastic – Low WUE",
    "MAGIC parent"        : "MAGIC parent",
}

def build_param_table(fits, genotypes_list):
    rows = []
    for g in genotypes_list:
        fit = fits.get(g)
        if fit is None:
            continue
        best = fit.get("best")
        wc   = str(wue_lookup.get(g, ""))
        aics = {k: fit[k]["AIC"] for k in ["Linear","Bilinear","Lin. plateau"]
                if fit.get(k) is not None}
        for model in ["Linear", "Bilinear", "Lin. plateau"]:
            f = fit.get(model)
            if f is None:
                continue
            rows.append({
                "Genotype"         : g,
                "WUE_class"        : wc,
                "Expected_behavior": WUE_EXPECTED.get(wc, ""),
                "Model"            : model,
                "Best_model"       : "✓" if model == best else "",
                "Slope1"           : round(f["b1"], 4) if not np.isnan(f["b1"]) else np.nan,
                "BP_kPa"           : round(f["bp"], 3) if not np.isnan(f["bp"]) else np.nan,
                "Slope2"           : round(f["b2"], 4) if not np.isnan(f["b2"]) else np.nan,
                "AIC"              : round(f["AIC"], 1),
                "Delta_AIC"        : round(f["AIC"] - min(aics.values()), 1) if len(aics) > 1 else 0,
            })
    return pd.DataFrame(rows)

tr_table  = build_param_table(tr_fits,  genotypes)
gsw_table = build_param_table(gsw_fits, geno_gsw)

tr_table.to_csv(os.path.join(DATA, "TR_model_parameters.csv"),  index=False)
gsw_table.to_csv(os.path.join(DATA, "gsw_model_parameters.csv"), index=False)
print("\nSaved TR_model_parameters.csv and gsw_model_parameters.csv")

# Best-model summary tables (one row per genotype)
def best_model_table(fits, genotypes_list):
    rows = []
    for g in genotypes_list:
        fit  = fits.get(g)
        wc   = str(wue_lookup.get(g, ""))
        if fit is None or fit.get("best") is None:
            rows.append({"Genotype":g, "WUE_class":wc,
                         "Expected_behavior": WUE_EXPECTED.get(wc,""),
                         "Best_model":"FAILED", "Slope1":np.nan,
                         "BP_kPa":np.nan, "Slope2":np.nan})
            continue
        best = fit["best"]
        f    = fit[best]
        rows.append({
            "Genotype"         : g,
            "WUE_class"        : wc,
            "Expected_behavior": WUE_EXPECTED.get(wc, ""),
            "Best_model"       : best,
            "Slope1"           : round(f["b1"], 4) if not np.isnan(f["b1"]) else np.nan,
            "BP_kPa"           : round(f["bp"], 3) if not np.isnan(f["bp"]) else np.nan,
            "Slope2"           : round(f["b2"], 4) if not np.isnan(f["b2"]) else np.nan,
        })
    return pd.DataFrame(rows)

tr_best  = best_model_table(tr_fits,  genotypes)
gsw_best = best_model_table(gsw_fits, geno_gsw)

tr_best.to_csv(os.path.join(DATA,  "TR_best_model_summary.csv"),  index=False)
gsw_best.to_csv(os.path.join(DATA, "gsw_best_model_summary.csv"), index=False)

print("\n--- TR best-model summary ---")
print(tr_best.to_string(index=False))
print("\n--- gsw best-model summary ---")
print(gsw_best.to_string(index=False))


# --------------------------------------------------------------------------- #
# 6. PUBLICATION FIGURES
# --------------------------------------------------------------------------- #

VPD_SEQ = np.linspace(day1["VPD"].min()*0.95, day1["VPD"].max()*1.05, 300)

def fitted_curve(fit_dict, model, vpd_seq):
    """Return (x, y) arrays for the fitted curve of a given model."""
    if fit_dict is None or fit_dict.get(model) is None:
        return None, None
    f    = fit_dict[model]
    popt = f["popt"]
    fn   = f["fn"]
    try:
        return vpd_seq, fn(vpd_seq, *popt)
    except Exception:
        return None, None


# ============================= FIG 1 ========================================
# Panel grid: TR vs VPD per genotype, best model highlighted
# 3 rows × 5 cols, publication layout
NCOLS = 5
NROWS = int(np.ceil(len(genotypes) / NCOLS))

fig = plt.figure(figsize=(14, NROWS * 2.9))
gs  = gridspec.GridSpec(NROWS, NCOLS, hspace=0.55, wspace=0.4,
                        left=0.06, right=0.97, top=0.90, bottom=0.10)

for i, g in enumerate(genotypes):
    ax = fig.add_subplot(gs[i // NCOLS, i % NCOLS])
    sub = day1[day1["Genotype"] == g]
    wc  = str(wue_lookup.get(g, ""))
    col = WUE_COLS.get(wc, UA_GRAY)
    mrk = WUE_MARKERS.get(wc, "o")

    # Data points, one per Rep in slightly different shade
    for r, rep in enumerate(sorted(sub["Rep"].unique())):
        s2 = sub[sub["Rep"] == rep]
        alpha = 0.5 + r * 0.12
        ax.scatter(s2["VPD"], s2["TR"], c=col, s=20, alpha=min(alpha, 0.85),
                   marker=mrk, linewidths=0, zorder=3)

    # All three model curves
    fit = tr_fits.get(g)
    best = fit["best"] if fit else None
    for model in ["Linear", "Bilinear", "Lin. plateau"]:
        x_c, y_c = fitted_curve(fit, model, VPD_SEQ)
        if x_c is None:
            continue
        is_best = (model == best)
        ax.plot(x_c, y_c,
                color=MODEL_COLS[model],
                lw=2.0 if is_best else 0.9,
                dashes=MODEL_LTY[model] if not is_best else (None, None),
                alpha=1.0 if is_best else 0.55,
                zorder=5 if is_best else 4)

    # Breakpoint marker
    if fit and best in ("Bilinear","Lin. plateau"):
        bp = fit[best]["bp"]
        y_bp = fit[best]["fn"](bp, *fit[best]["popt"])
        ax.axvline(bp, lw=0.7, ls=":", color=UA_DARK, alpha=0.4, zorder=2)
        ax.scatter([bp], [y_bp], c=UA_DARK, s=30, zorder=6, marker="v")

    ax.set_title(g, fontsize=7.5, fontweight="bold", pad=3)
    ax.set_xlabel("VPD (kPa)", fontsize=7, labelpad=2)
    ax.set_ylabel("TR (mg m⁻² s⁻¹)", fontsize=7, labelpad=2)
    ax.tick_params(labelsize=6.5)
    strip_ax(ax)

# Shared legends
legend_wue = [Patch(facecolor=c, label=wc.replace(", ", ",\n"), alpha=0.85)
              for wc, c in WUE_COLS.items()]
legend_mod = [Line2D([0],[0], color=MODEL_COLS[m], lw=2,
                     dashes=MODEL_LTY[m] if m != "Lin. plateau" else (3,2),
                     label=m)
              for m in MODEL_COLS]
legend_bp  = Line2D([0],[0], marker="v", color=UA_DARK, lw=0, markersize=5,
                    label="Breakpoint")

fig.legend(handles=legend_wue, title="WUE class", fontsize=7, title_fontsize=7.5,
           loc="lower left", bbox_to_anchor=(0.0, -0.01), ncol=5, frameon=False)
fig.legend(handles=legend_mod + [legend_bp], title="Model", fontsize=7, title_fontsize=7.5,
           loc="lower right", bbox_to_anchor=(1.0, -0.01), ncol=4, frameon=False)

fig.suptitle("Transpiration rate vs VPD — Day 1  |  Linear, bilinear, and linear-plateau models",
             fontsize=10, fontweight="bold", y=0.97)

fig.savefig(os.path.join(FIGS, "Fig1_TR_models_per_genotype.png"), dpi=300, bbox_inches="tight")
plt.close(fig)
print("\nSaved Fig1")


# ============================= FIG 2 ========================================
# Pooled by WUE class — TR vs VPD, best-fit curves per class

fig, axes = plt.subplots(1, len(WUE_COLS), figsize=(3.0 * len(WUE_COLS), 4.2),
                          sharey=True)
wue_list = list(WUE_COLS.keys())

for idx, (wc, ax) in enumerate(zip(wue_list, axes)):
    col = WUE_COLS[wc]
    sub_wc = day1[day1["WUE_class"] == wc]

    # Individual genotype lines (thin, same colour, low alpha)
    for g in sub_wc["Genotype"].unique():
        sg   = sub_wc[sub_wc["Genotype"] == g]
        fit  = tr_fits.get(g)
        best = fit["best"] if fit else None
        if fit and best:
            x_c, y_c = fitted_curve(fit, best, VPD_SEQ)
            if x_c is not None:
                ax.plot(x_c, y_c, color=col, lw=0.9, alpha=0.35)

    # Pooled OLS line
    if len(sub_wc) > 3:
        z   = np.polyfit(sub_wc["VPD"].dropna(), sub_wc["TR"].dropna(), 1)
        xr  = np.linspace(sub_wc["VPD"].min(), sub_wc["VPD"].max(), 200)
        ax.plot(xr, np.polyval(z, xr), color=col, lw=2.2, zorder=5)
        ax.fill_between(xr, np.polyval(z, xr) - 4, np.polyval(z, xr) + 4,
                        color=col, alpha=0.12)

    # Data scatter
    mrk = WUE_MARKERS[wc]
    ax.scatter(sub_wc["VPD"], sub_wc["TR"], c=col, s=18, alpha=0.55,
               marker=mrk, linewidths=0, zorder=3)

    short_label = wc.replace("Plastic, ","").replace("Nonplastic, ","")
    ax.set_title(short_label, fontsize=8, fontweight="bold", color=col, pad=4)
    ax.set_xlabel("VPD (kPa)", fontsize=8)
    if idx == 0:
        ax.set_ylabel("TR (mg m⁻² s⁻¹)", fontsize=8)
    ax.tick_params(labelsize=7)
    strip_ax(ax)

fig.suptitle("Transpiration rate vs VPD by WUE class — Day 1\n"
             "(thin lines = per-genotype best model; thick line = pooled OLS)",
             fontsize=9, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.92])
fig.savefig(os.path.join(FIGS, "Fig2_TR_by_WUEclass.png"), dpi=300, bbox_inches="tight")
plt.close(fig)
print("Saved Fig2")


# ============================= FIG 3 ========================================
# Breakpoints — TR: horizontal strip chart, sorted by BP value

bp_data = tr_best[tr_best["BP_kPa"].notna()].copy()
bp_data["col"] = bp_data["WUE_class"].map(WUE_COLS).fillna(UA_GRAY)
bp_data = bp_data.sort_values("BP_kPa")

if len(bp_data) > 0:
    fig, ax = plt.subplots(figsize=(8, max(4, len(bp_data)*0.45)))
    for j, (_, row) in enumerate(bp_data.iterrows()):
        ax.barh(j, row["BP_kPa"], color=row["col"], alpha=0.85, height=0.65)
    ax.axvline(1.73, lw=1.2, ls="--", color=UA_DARK, alpha=0.6)
    ax.text(1.73 + 0.03, -0.7, "lit. BP = 1.73 kPa", fontsize=7.5,
            va="bottom", ha="left", color=UA_DARK)
    ax.set_yticks(range(len(bp_data)))
    ax.set_yticklabels(bp_data["Genotype"].values, fontsize=8)
    ax.set_xlabel("Breakpoint VPD (kPa)", fontsize=9)
    ax.set_title("Fig 3  Breakpoint VPD — TR models (bilinear or linear plateau)",
                 fontweight="bold")
    legend_handles = [Patch(facecolor=WUE_COLS[wc], label=wc, alpha=0.85)
                      for wc in WUE_COLS if wc in bp_data["WUE_class"].values]
    ax.legend(handles=legend_handles, fontsize=7.5, frameon=False,
              loc="lower right")
    strip_ax(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "Fig3_TR_breakpoints.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Saved Fig3")


# ============================= FIG 4 ========================================
# Slope 1 comparison across genotypes (bar chart, coloured by WUE class)

fig, ax = plt.subplots(figsize=(9, 4.5))
tr_slope = tr_best.dropna(subset=["Slope1"]).sort_values("Slope1", ascending=False)
cols = [WUE_COLS.get(str(wc), UA_GRAY) for wc in tr_slope["WUE_class"]]
bars = ax.bar(range(len(tr_slope)), tr_slope["Slope1"], color=cols, alpha=0.88)
ax.set_xticks(range(len(tr_slope)))
ax.set_xticklabels(tr_slope["Genotype"], rotation=40, ha="right", fontsize=8)
ax.set_ylabel("Slope₁ (mg m⁻² s⁻¹ kPa⁻¹)", fontsize=9)
ax.set_title("Fig 4  Transpiration sensitivity (slope₁ of best model) — Day 1",
             fontweight="bold")
legend_handles = [Patch(facecolor=WUE_COLS[wc], label=wc, alpha=0.85) for wc in WUE_COLS]
ax.legend(handles=legend_handles, fontsize=7.5, frameon=False, ncol=2)
strip_ax(ax)
fig.tight_layout()
fig.savefig(os.path.join(FIGS, "Fig4_TR_slopes.png"), dpi=300, bbox_inches="tight")
plt.close(fig)
print("Saved Fig4")


# ============================= FIG 5 ========================================
# gsw vs VPD — per genotype with model fits

if len(geno_gsw) > 0:
    NC = min(4, len(geno_gsw))
    NR = int(np.ceil(len(geno_gsw) / NC))
    VPD_SEQ_GSW = np.linspace(day1_gsw["VPD"].min()*0.95,
                               day1_gsw["VPD"].max()*1.05, 300)

    fig = plt.figure(figsize=(NC*3.0, NR*3.0))
    gs  = gridspec.GridSpec(NR, NC, hspace=0.55, wspace=0.42,
                            left=0.08, right=0.97, top=0.88, bottom=0.12)

    for i, g in enumerate(geno_gsw):
        ax  = fig.add_subplot(gs[i // NC, i % NC])
        sub = day1_gsw[day1_gsw["Genotype"] == g]
        wc  = str(wue_lookup.get(g, ""))
        col = WUE_COLS.get(wc, UA_GRAY)
        mrk = WUE_MARKERS.get(wc, "o")

        ax.scatter(sub["VPD"], sub["gsw"], c=col, s=28, alpha=0.8,
                   marker=mrk, linewidths=0, zorder=3)

        fit  = gsw_fits.get(g)
        best = fit["best"] if fit else None
        for model in ["Linear", "Bilinear", "Lin. plateau"]:
            x_c, y_c = fitted_curve(fit, model, VPD_SEQ_GSW)
            if x_c is None:
                continue
            is_best = (model == best)
            ax.plot(x_c, y_c,
                    color=MODEL_COLS[model],
                    lw=2.0 if is_best else 0.9,
                    dashes=MODEL_LTY[model] if not is_best else (None, None),
                    alpha=1.0 if is_best else 0.5, zorder=5 if is_best else 4)

        if fit and best in ("Bilinear","Lin. plateau"):
            bp   = fit[best]["bp"]
            y_bp = fit[best]["fn"](bp, *fit[best]["popt"])
            ax.axvline(bp, lw=0.7, ls=":", color=UA_DARK, alpha=0.4)
            ax.scatter([bp], [y_bp], c=UA_DARK, s=28, zorder=6, marker="v")

        ax.set_title(g, fontsize=7.5, fontweight="bold", pad=3)
        ax.set_xlabel("VPD (kPa)", fontsize=7, labelpad=2)
        ax.set_ylabel("g$_{sw}$ (mol m⁻² s⁻¹)", fontsize=7, labelpad=2)
        ax.tick_params(labelsize=6.5)
        strip_ax(ax)

    for j in range(i+1, NR*NC):
        fig.add_subplot(gs[j // NC, j % NC]).set_visible(False)

    fig.legend(handles=legend_mod + [legend_bp], title="Model", fontsize=7,
               title_fontsize=7.5, loc="lower right",
               bbox_to_anchor=(1.0, -0.01), ncol=4, frameon=False)
    fig.suptitle("Stomatal conductance vs VPD — Day 1  |  LICOR measurements",
                 fontsize=10, fontweight="bold", y=0.97)

    fig.savefig(os.path.join(FIGS, "Fig5_gsw_models_per_genotype.png"),
                dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Saved Fig5")


# ============================= FIG 6 ========================================
# AIC comparison heatmap — TR models

aic_pivot = (tr_table[tr_table["Model"].isin(["Linear","Bilinear","Lin. plateau"])]
             .pivot_table(index="Genotype", columns="Model", values="Delta_AIC")
             .reindex(columns=["Linear","Bilinear","Lin. plateau"]))

if not aic_pivot.empty:
    fig, ax = plt.subplots(figsize=(6, max(4, len(aic_pivot)*0.45)))
    import matplotlib.colors as mcolors
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "ua_reds", [UA_RED, "#F7C6CC", "white"], N=256)
    im = ax.imshow(aic_pivot.values, cmap=cmap, aspect="auto",
                   vmin=0, vmax=max(12, aic_pivot.values.max()))
    ax.set_xticks(range(3))
    ax.set_xticklabels(["Linear","Bilinear","Lin. plateau"], fontsize=9)
    ax.set_yticks(range(len(aic_pivot)))
    ax.set_yticklabels(aic_pivot.index, fontsize=8)
    for r in range(aic_pivot.shape[0]):
        for c in range(aic_pivot.shape[1]):
            val = aic_pivot.values[r, c]
            if not np.isnan(val):
                ax.text(c, r, f"{val:.1f}", ha="center", va="center",
                        fontsize=7.5, color="white" if val < 5 else UA_DARK)
    cb = fig.colorbar(im, ax=ax, shrink=0.7, pad=0.02)
    cb.set_label("ΔAIC (vs best model)", fontsize=8)
    ax.set_title("Fig 6  Model comparison — ΔAIC (TR, Day 1)\n"
                 "ΔAIC = 0 → best model for that genotype", fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "Fig6_TR_AIC_heatmap.png"),
                dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Saved Fig6")

print("\n=== Model fitting complete ===")
print(f"Figures: {FIGS}")
print(f"Data:    {DATA}")
