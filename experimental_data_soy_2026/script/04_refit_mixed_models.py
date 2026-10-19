#!/usr/bin/env python3
"""
04_refit_mixed_models.py
Growth Chamber 2026 - Soybean Limited-Transpiration Experiment

Python re-implementation of the R pipeline (02_model_fitting_TR_gs_VPD.R):
  * Day 1 data only, per genotype
  * Linear, bilinear and linear-plateau TR ~ VPD models (and gsw ~ VPD)
  * Replicate as a random intercept. Because the random effect is additive,
    the marginal likelihood is exactly Gaussian with a block-diagonal
    covariance (sigma^2 I + sigma_u^2 J per Rep), so all three models are
    fitted by maximum (marginal) likelihood with scipy and compared by AIC.
    This is the same likelihood nlme::lme uses for the linear model and the
    one nlme::nlme approximates for the nonlinear ones.

Inputs : data/data_clean.csv (made by 01_data_inspection_QC.py)
Outputs: data/TR_best_model_summary.csv, data/gsw_best_model_summary.csv,
         data/TR_model_parameters.csv, data/TR_quadrant_classification.csv,
         data/predicted_TR_gsw_at_target_VPD.csv
         figures/paper_figures/Figure1_TR_response_VPD.png (14 panels)
         figures/paper_figures/Figure2_TR_breakpoints.png
         figures/paper_figures/Figure3_intercept_slope_quadrants.png
         figures/paper_figures/Figure5_gsw_predicted.png

Run from experimental_data_soy_2026/:
    python3 script/04_refit_mixed_models.py
Optional: --out-suffix _5rep  (appends to every output file name)
"""
import os, sys, argparse, warnings
import numpy as np
import pandas as pd
from scipy.optimize import minimize, curve_fit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

warnings.filterwarnings("ignore")
ap = argparse.ArgumentParser()
ap.add_argument("--out-suffix", default="")
ap.add_argument("--input", default="data/data_clean.csv")
ap.add_argument("--min-side", type=int, default=3, help="min observations on each side of a breakpoint")
args = ap.parse_args()
SUF = args.out_suffix

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
DATA = "data"
FIGS = "figures/paper_figures"
os.makedirs(FIGS, exist_ok=True)

GREEN = "#31a354"
WUE_ORDER = ["slow_wilting", "high_wilting", "Plastic, High WUE", "Nonplastic, Low WUE", "MAGIC parent"]
WUE_LABELS = {"slow_wilting": "Slow wilting", "high_wilting": "High wilting",
              "Plastic, High WUE": "High WUE plasticity", "Nonplastic, Low WUE": "Low WUE plasticity",
              "MAGIC parent": "MAGIC parent"}
STRIP_COLS = {"slow_wilting": "#6E1424", "high_wilting": "#C8102E", "Plastic, High WUE": "#4393c3",
              "Nonplastic, Low WUE": "#d4824a", "MAGIC parent": "#555555"}
WUE_EXPECTED = {"slow_wilting": "Slow wilting / LT", "high_wilting": "High wilting / non-LT",
                "Plastic, High WUE": "Plastic - High WUE", "Nonplastic, Low WUE": "Nonplastic - Low WUE",
                "MAGIC parent": "MAGIC parent"}

# ----------------------------------------------------------------------------- data
df = pd.read_csv(args.input, dtype={"Rep": str})
day1 = df[df["Day_label"] == "Day 1"].copy()
# gsw sanity filter: leaf-level conductance above 3 mol m-2 s-1 is physically impossible
# (one Rep 5 entry for A5959 reads 29.3, an entry/unit error) -> excluded from the gsw fits
bad = day1["gsw"] > 3
if bad.any():
    print("Excluded impossible gsw values:\n", day1.loc[bad, ["Rep", "Genotype", "Treatment", "VPD", "gsw"]].to_string(index=False))
day1_gsw = day1[day1["gsw"].notna() & ~bad].copy()
wue = day1.drop_duplicates("Genotype").set_index("Genotype")["WUE_class"].to_dict()
genotypes = sorted(day1["Genotype"].unique())
geno_gsw = sorted(day1_gsw["Genotype"].unique())
print(f"Day 1: {len(day1)} rows, {len(genotypes)} genotypes, {len(day1_gsw)} gsw rows ({len(geno_gsw)} genotypes)")

# ----------------------------------------------------------------------------- models
def f_linear(x, a, b): return a + b * x
def f_bilinear(x, a, b1, b2, bp): return np.where(x <= bp, a + b1 * x, a + b1 * bp + b2 * (x - bp))
def f_plateau(x, a, b, bp): return np.where(x <= bp, a + b * x, a + b * bp)

MODELS = {"Linear": (f_linear, 2), "Bilinear": (f_bilinear, 4), "Lin. plateau": (f_plateau, 3)}


def marginal_nll(params, fn, k, x, y, groups):
    """-log L of y = fn(x; theta) + u_g + e, u_g ~ N(0, su2), e ~ N(0, s2)."""
    theta, ls2, lsu2 = params[:k], params[k], params[k + 1]
    s2, su2 = np.exp(ls2), np.exp(lsu2)
    r = y - fn(x, *theta)
    nll = 0.0
    for g in np.unique(groups):
        m = groups == g
        rg = r[m]; n = m.sum()
        # Sigma = s2 I + su2 J ; use Sherman-Morrison / determinant lemma
        # logdet = (n-1) log s2 + log(s2 + n su2)
        # quad = (rg'rg)/s2 - su2 (sum rg)^2 / (s2 (s2 + n su2))
        logdet = (n - 1) * np.log(s2) + np.log(s2 + n * su2)
        quad = rg @ rg / s2 - su2 * rg.sum() ** 2 / (s2 * (s2 + n * su2))
        nll += 0.5 * (n * np.log(2 * np.pi) + logdet + quad)
    return nll


def fit_mixed(fn, k, x, y, groups, theta0, bounds_theta):
    s2_0 = np.var(y - fn(x, *theta0)) + 1e-6
    p0 = np.r_[theta0, np.log(s2_0), np.log(s2_0 * 0.2)]
    bnds = list(bounds_theta) + [(np.log(1e-6), None), (np.log(1e-8), None)]
    best = None
    for start_scale in (1.0, 0.5, 2.0):
        p = p0.copy(); p[k + 1] += np.log(start_scale)
        try:
            res = minimize(marginal_nll, p, args=(fn, k, x, y, groups), method="L-BFGS-B", bounds=bnds,
                           options={"maxiter": 2000})
            if best is None or res.fun < best.fun:
                best = res
        except Exception:
            pass
    if best is None or not np.isfinite(best.fun):
        return None
    theta = best.x[:k]
    return {"theta": theta, "logL": -best.fun, "AIC": 2 * (k + 2) + 2 * best.fun,
            "s2": np.exp(best.x[k]), "su2": np.exp(best.x[k + 1])}


def fit_one(sub, col):
    d = sub[["VPD", col, "Rep"]].dropna()
    if len(d) < 4:
        return None
    x = d["VPD"].values.astype(float); y = d[col].values.astype(float); g = d["Rep"].values
    vmin, vmax = x.min(), x.max(); bp0 = 0.5 * (vmin + vmax)
    sl0 = np.polyfit(x, y, 1)[0]; a0 = y.mean() - sl0 * x.mean()
    out = {}
    # Linear
    out["Linear"] = fit_mixed(f_linear, 2, x, y, g, [a0, sl0], [(None, None), (None, None)])
    # Bilinear: NLS start then mixed
    try:
        popt, _ = curve_fit(f_bilinear, x, y, p0=[a0, max(sl0, 0.1), sl0 * 0.1, bp0],
                            bounds=([-np.inf, 0, -np.inf, vmin], [np.inf, np.inf, np.inf, vmax]), maxfev=8000)
        out["Bilinear"] = fit_mixed(f_bilinear, 4, x, y, g, popt,
                                    [(0, None), (0, None), (None, None), (vmin, vmax)])
    except Exception:
        out["Bilinear"] = None
    # Linear plateau
    try:
        popt, _ = curve_fit(f_plateau, x, y, p0=[a0, max(sl0, 0.1), bp0],
                            bounds=([-np.inf, 0, vmin], [np.inf, np.inf, vmax]), maxfev=8000)
        out["Lin. plateau"] = fit_mixed(f_plateau, 3, x, y, g, popt, [(0, None), (0, None), (vmin, vmax)])
    except Exception:
        out["Lin. plateau"] = None
    # A breakpoint is only admissible when it is supported by data on both sides
    # (>= 2 observations below and above). Nonlinear intercepts are bounded at 0
    # (transpiration cannot be negative at zero VPD); both rules mirror what the
    # R nls/nlme fits produced on the 4-rep data.
    for m in ("Bilinear", "Lin. plateau"):
        r = out.get(m)
        if r is not None:
            bp = r["theta"][-1]
            if (x < bp).sum() < args.min_side or (x > bp).sum() < args.min_side:
                out[m] = None
    aics = {m: r["AIC"] for m, r in out.items() if r is not None}
    out["best"] = min(aics, key=aics.get) if aics else None
    out["x"] = x; out["y"] = y
    return out


def params(r, model):
    th = r["theta"]
    if model == "Linear":
        return dict(a=th[0], b1=th[1], b2=np.nan, bp=np.nan)
    if model == "Bilinear":
        return dict(a=th[0], b1=th[1], b2=th[2], bp=th[3])
    return dict(a=th[0], b1=th[1], b2=np.nan, bp=th[2])


def predict(fit, model, xs):
    fn = MODELS[model][0]
    return fn(np.asarray(xs, float), *fit[model]["theta"])


print("\nFitting TR models (mixed, ML) ...")
tr = {g: fit_one(day1[day1.Genotype == g], "TR") for g in genotypes}
print("Fitting gsw models ...")
gs = {g: fit_one(day1_gsw[day1_gsw.Genotype == g], "gsw") for g in geno_gsw}

# ----------------------------------------------------------------------------- tables
def best_table(fits):
    rows = []
    for g, r in fits.items():
        if r is None or r["best"] is None:
            continue
        p = params(r[r["best"]], r["best"])
        rows.append(dict(Genotype=g, WUE_class=wue[g], Expected_behavior=WUE_EXPECTED.get(wue[g], ""),
                         Best_model=r["best"], Intercept=round(p["a"], 3), Slope1=round(p["b1"], 4),
                         BP_kPa=round(p["bp"], 3) if np.isfinite(p["bp"]) else np.nan,
                         Slope2=round(p["b2"], 4) if np.isfinite(p["b2"]) else np.nan,
                         AIC_best=round(r[r["best"]]["AIC"], 1), n=len(r["y"])))
    return pd.DataFrame(rows)


def param_table(fits):
    rows = []
    for g, r in fits.items():
        if r is None:
            continue
        aics = {m: r[m]["AIC"] for m in MODELS if r.get(m) is not None}
        for m in MODELS:
            if r.get(m) is None:
                continue
            p = params(r[m], m)
            rows.append(dict(Genotype=g, WUE_class=wue[g], Model=m, Best_model="yes" if m == r["best"] else "",
                             a=round(p["a"], 3), Slope1=round(p["b1"], 4),
                             Slope2=round(p["b2"], 4) if np.isfinite(p["b2"]) else np.nan,
                             BP_kPa=round(p["bp"], 3) if np.isfinite(p["bp"]) else np.nan,
                             AIC=round(r[m]["AIC"], 1), Delta_AIC=round(r[m]["AIC"] - min(aics.values()), 1)))
    return pd.DataFrame(rows)


tr_best = best_table(tr); gs_best = best_table(gs)
tr_best.to_csv(f"{DATA}/TR_best_model_summary{SUF}.csv", index=False)
gs_best.to_csv(f"{DATA}/gsw_best_model_summary{SUF}.csv", index=False)
param_table(tr).to_csv(f"{DATA}/TR_model_parameters{SUF}.csv", index=False)
param_table(gs).to_csv(f"{DATA}/gsw_model_parameters{SUF}.csv", index=False)

x_mid = tr_best["Intercept"].median(); s_mid = tr_best["Slope1"].median()
quad = tr_best.copy()
quad["Quadrant"] = np.where(quad.Intercept >= x_mid, "High intercept", "Low intercept") + " / " + \
                   np.where(quad.Slope1 >= s_mid, "High sensitivity", "Low sensitivity")
quad[["Genotype", "WUE_class", "Best_model", "Intercept", "Slope1", "Quadrant"]].to_csv(
    f"{DATA}/TR_quadrant_classification{SUF}.csv", index=False)

pred_rows = []
for g, r in tr.items():
    if r and r["best"]:
        for v in (1.5, 3.5):
            pred_rows.append(dict(Genotype=g, WUE_class=wue[g], VPD_target=v, variable="TR",
                                  value=float(predict(r, r["best"], [v])[0])))
for g, r in gs.items():
    if r and r["best"]:
        for v in (1.5, 3.5):
            pred_rows.append(dict(Genotype=g, WUE_class=wue[g], VPD_target=v, variable="gsw",
                                  value=float(predict(r, r["best"], [v])[0])))
pred = pd.DataFrame(pred_rows)
pred.to_csv(f"{DATA}/predicted_TR_gsw_at_target_VPD{SUF}.csv", index=False)

print("\n--- TR best-model summary ---")
print(tr_best.to_string(index=False))
print("\n--- gsw best-model summary ---")
print(gs_best.to_string(index=False))
print(f"\nQuadrant medians: intercept {x_mid:.2f}, slope {s_mid:.3f}")

# ----------------------------------------------------------------------------- figures
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "axes.linewidth": 1.0,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": "white", "savefig.facecolor": "white", "savefig.dpi": 300})


def panel_label(g):
    lab = WUE_LABELS.get(wue[g], wue[g])
    if g == "PI398198":  # treated as low plasticity in the paper figures
        lab = "Low WUE plasticity"
    return f"{g}\n({lab})"


# Figure 1: TR vs VPD, all genotypes, 3 x 5
order = ["PI398492", "PI398617", "PI424231", "PI442012A", "PI603418D",
         "PI377574", "PI398272", "PI404191", "PI603171", "PI398198",
         "A5959", "AG43XF2", "PI416937", "PI471938"]
order = [g for g in order if g in tr] + [g for g in genotypes if g not in order]
ncol = 5; nrow = int(np.ceil(len(order) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(15, 3.3 * nrow), sharex=True, sharey=True)
xs = np.linspace(0.9, 3.9, 200)
for ax, g in zip(axes.flat, order):
    r = tr[g]
    ax.scatter(r["x"], r["y"], s=40, facecolor=GREEN, edgecolor="black", linewidth=0.8, zorder=3)
    if r["best"]:
        ax.plot(xs, predict(r, r["best"], xs), color=GREEN, lw=2.2)
    ax.set_title(panel_label(g), fontsize=11, pad=4)
    ax.set_ylim(0, 125); ax.set_xlim(0.8, 4.0); ax.set_xticks([1, 2, 3, 4])
    for s in ("top", "right"):
        ax.spines[s].set_visible(True)
for ax in axes.flat[len(order):]:
    ax.axis("off")
fig.supxlabel("Vapor Pressure Deficit (kPa)", fontsize=13)
fig.supylabel("Transpiration (mg m$^{-2}$ s$^{-1}$)", fontsize=13)
fig.tight_layout(rect=(0.02, 0.02, 1, 1))
fig.savefig(f"{FIGS}/Figure1_TR_response_VPD{SUF}.png", bbox_inches="tight"); plt.close(fig)

# Figure 2: breakpoints
bp = tr_best.dropna(subset=["BP_kPa"]).sort_values("BP_kPa")
fig, ax = plt.subplots(figsize=(9, 6))
ax.bar(bp.Genotype, bp.BP_kPa, color=GREEN, width=0.65)
m = bp.BP_kPa.mean()
ax.axhline(m, ls="--", color="grey", lw=2)
ax.text(len(bp) - 0.5, m + 0.05, f"mean = {m:.2f} kPa", ha="right", color="grey", style="italic")
ax.set_ylabel("Breakpoint VPD (kPa)", fontsize=14); ax.set_ylim(0, max(3.4, bp.BP_kPa.max() + 0.2))
ax.tick_params(axis="x", rotation=35)
for s in ("top", "right"):
    ax.spines[s].set_visible(True)
fig.tight_layout(); fig.savefig(f"{FIGS}/Figure2_TR_breakpoints{SUF}.png", bbox_inches="tight"); plt.close(fig)

# Figure 3: intercept vs slope quadrants
markers = {"Linear": "o", "Bilinear": "s", "Lin. plateau": "^"}
# manual label nudges (points) for crowded labels
LABEL_OFFSET = {"A5959": (0, -20), "AG43XF2": (-38, 10), "PI424231": (0, -20), "PI398492": (-42, 12), "PI398198": (0, -20), "PI404191": (-10, 14)}
fig, ax = plt.subplots(figsize=(10, 7))
xr = (tr_best.Intercept.min() - 3, tr_best.Intercept.max() + 3)
yr = (tr_best.Slope1.min() - 2, tr_best.Slope1.max() + 2)
ax.add_patch(Rectangle((xr[0], s_mid), x_mid - xr[0], yr[1] - s_mid, color="#f3f3f3", zorder=0))
ax.add_patch(Rectangle((x_mid, s_mid), xr[1] - x_mid, yr[1] - s_mid, color="#f3f3f3", zorder=0))
ax.axvline(x_mid, ls="--", color="grey", lw=1.8); ax.axhline(s_mid, ls="--", color="grey", lw=1.8)
for _, r in tr_best.iterrows():
    ax.scatter(r.Intercept, r.Slope1, marker=markers[r.Best_model], s=160, facecolor=GREEN,
               edgecolor="black", linewidth=1, zorder=3)
    dx, dy = LABEL_OFFSET.get(r.Genotype, (0, 12))
    ax.annotate(r.Genotype, (r.Intercept, r.Slope1), xytext=(dx, dy), textcoords="offset points",
                ha="center", fontsize=9, bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="black", lw=0.5))
ax.text(x_mid - 0.5, yr[1] - 0.6, "Low water use\nunder high VPD", ha="right", va="top", fontsize=10)
ax.text(xr[1] - 0.5, yr[1] - 0.6, "High water use\nunder high VPD", ha="right", va="top", fontsize=10)
ax.text(xr[0] + 0.5, yr[0] + 0.4, "Low water use\nunder low VPD", ha="left", va="bottom", fontsize=10)
ax.text(xr[1] - 0.5, yr[0] + 0.4, "High water use\nunder low VPD", ha="right", va="bottom", fontsize=10)
ax.set_xlim(*xr); ax.set_ylim(*yr)
ax.set_xlabel("Intercept (mg m$^{-2}$ s$^{-1}$)", fontsize=14)
ax.set_ylabel("Slope (mg m$^{-2}$ s$^{-1}$ kPa$^{-1}$)", fontsize=14)
for s in ("top", "right"):
    ax.spines[s].set_visible(True)
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([], [], marker=mk, ls="", mfc=GREEN, mec="black", ms=11, label=lab)
                   for lab, mk in markers.items()], loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=3,
          frameon=True, edgecolor="black")
fig.tight_layout(); fig.savefig(f"{FIGS}/Figure3_intercept_slope_quadrants{SUF}.png", bbox_inches="tight"); plt.close(fig)

# Figure 5 (gsw part): predicted gsw at 1.5 and 3.5 kPa
gp = pred[pred.variable == "gsw"].pivot(index="Genotype", columns="VPD_target", values="value")
gorder = [g for g in order if g in gp.index]
ncol = 5; nrow = int(np.ceil(len(gorder) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(15, 2.6 * nrow + 0.6), sharey=True)
for ax, g in zip(axes.flat, gorder):
    v15, v35 = gp.loc[g, 1.5], gp.loc[g, 3.5]
    ax.bar([0, 1], [max(v15, 0), max(v35, 0)], color=["#4393c3", "#6E1424"], width=0.6)
    lab = {"Plastic, High WUE": "Plastic / High WUE", "Nonplastic, Low WUE": "Nonplastic / Low WUE",
           "slow_wilting": "Slow wilting", "high_wilting": "High wilting", "MAGIC parent": "MAGIC parent"}[wue[g]]
    ax.set_title(f"{g}\n[{lab}]", fontsize=9, color="white", fontweight="bold",
                 bbox=dict(boxstyle="square,pad=0.35", fc=STRIP_COLS[wue[g]], ec="none"), pad=6)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["1.5 kPa", "3.5 kPa"]); ax.set_ylim(0, 1.0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(True)
for ax in axes.flat[len(gorder):]:
    ax.axis("off")
fig.supylabel("g$_{sw}$ (mol m$^{-2}$ s$^{-1}$)", fontsize=13)
fig.tight_layout(rect=(0.02, 0, 1, 1))
fig.savefig(f"{FIGS}/Figure5_gsw_predicted{SUF}.png", bbox_inches="tight"); plt.close(fig)

# gsw decline summary
gp["decline_%"] = 100 * (gp[3.5] - gp[1.5]) / gp[1.5]
print("\n--- gsw change 1.5 -> 3.5 kPa (%) ---")
print(gp.round(3).to_string())
tp = pred[pred.variable == "TR"].pivot(index="Genotype", columns="VPD_target", values="value")
print("\n--- TR change 1.5 -> 3.5 kPa (%) ---")
print((100 * (tp[3.5] - tp[1.5]) / tp[1.5]).round(1).to_string())
print("\nDone.")
