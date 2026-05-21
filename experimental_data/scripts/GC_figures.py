"""
GC_figures.py — Generate all growth-chamber figures for GC1, GC2, GC3 and overall.
Run from /home/user/LT-trait:
    python3 experimental_data/scripts/GC_figures.py
"""

import os, warnings
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats

warnings.filterwarnings("ignore")

# ── Paths ──────────────────────────────────────────────────────────────────────
FIGDIR  = "experimental_data/figures"
os.makedirs(FIGDIR, exist_ok=True)

df = pd.read_csv("experimental_data/cleaned_data/GC_combined_clean.csv")
df["Genotype"] = df["Genotype"].str.strip()
df["Day"] = df["Day"].astype(str)

# ── Colour conventions ─────────────────────────────────────────────────────────
COL_QUAD = {"Q1": "#d95f02", "Q2": "#1b9e77", "Q3": "#7570b3", "Q4": "#e7298a"}
COL_EXP  = {"GC1": "#1f78b4", "GC2": "#ff7f00", "GC3": "#e31a1c"}
COL_DAY  = {"1": "#aec7e8", "2": "#1f77b4"}
VPD_BREAK = 2.0  # Sinclair 2008 reference breakpoint (kPa)

# Quadrant labels
QUAD_LAB = {
    "Q1": "Q1 · Plastic / High WUE",
    "Q2": "Q2 · Plastic / Low WUE",
    "Q3": "Q3 · Nonplastic / Low WUE",
    "Q4": "Q4 · Nonplastic / High WUE",
}

# Genotype → quadrant lookup
GENO_QUAD = df.drop_duplicates("Genotype").set_index("Genotype")["Quadrant"].to_dict()

plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 9,
    "figure.dpi": 150, "savefig.dpi": 150,
})

# ── Helpers ───────────────────────────────────────────────────────────────────
def save(name):
    path = os.path.join(FIGDIR, name + ".png")
    plt.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  saved {name}.png")


def outlier_mask(sub, col, n_sd=2.5):
    """True where col is a valid, non-outlier value."""
    if col not in sub.columns:
        return pd.Series(True, index=sub.index)
    grp = sub.groupby(["Genotype", "Day", "Treatment"])[col]
    mu  = grp.transform("mean")
    sd  = grp.transform("std")
    ok  = sub[col].notna() & ((sub[col] - mu).abs() <= n_sd * sd.fillna(np.inf))
    return ok


def plot_exp(ax, sub_day2, geno, col, color, vpd_col="VPD"):
    """Scatter + mean±SE for one genotype on one axes."""
    pts = sub_day2[sub_day2["Genotype"] == geno].dropna(subset=[col])
    if pts.empty:
        return
    mn = pts.groupby("VPD")[col].agg(["mean", "sem"]).reset_index()
    ax.axvline(VPD_BREAK, color="grey", lw=0.8, ls="--", zorder=0)
    ax.scatter(pts["VPD"], pts[col], color=color, alpha=0.30, s=18, zorder=2)
    ax.plot(mn["VPD"], mn["mean"], color=color, lw=1.6, zorder=3)
    ax.fill_between(mn["VPD"],
                    mn["mean"] - mn["sem"], mn["mean"] + mn["sem"],
                    color=color, alpha=0.20, zorder=1)
    ax.errorbar(mn["VPD"], mn["mean"], mn["sem"],
                fmt="o", color="white", mec=color, mew=1.2, ms=6, ecolor=color,
                elinewidth=0.9, zorder=4)


def exp_title(exp, genotypes, quadrants):
    parts = [f"{g} ({QUAD_LAB.get(q,'?')})" for g, q in zip(genotypes, quadrants)]
    return f"{exp} — " + " · ".join(parts)


# ═══════════════════════════════════════════════════════════════════════════════
# Per-experiment figures
# ═══════════════════════════════════════════════════════════════════════════════
for exp in ["GC1", "GC2", "GC3"]:
    sub = df[df["Experiment"] == exp].copy()
    sub["Genotype"] = sub["Genotype"].str.strip()

    # Day 2, non-outliers
    ok_tla = outlier_mask(sub, "TranspLeafArea")
    ok_rel = outlier_mask(sub, "RelTranspiration")
    sub["outlier_tla"] = ~ok_tla
    sub["outlier_rel"] = ~ok_rel
    d2       = sub[sub["Day"] == "2"].copy()
    d2_clean = d2[ok_tla[d2.index]].copy()

    genotypes = sorted(d2["Genotype"].dropna().unique())
    quadrants = [GENO_QUAD.get(g, "?") for g in genotypes]
    n         = len(genotypes)

    # ── Fig 1: Actual transpiration per leaf area (Day 2) ────────────────────
    has_tla = d2_clean["TranspLeafArea"].notna().any()
    if has_tla:
        fig, axes = plt.subplots(1, n, figsize=(5 * n, 4.5), sharey=True)
        if n == 1:
            axes = [axes]
        for ax, g in zip(axes, genotypes):
            q = GENO_QUAD.get(g, "Q1")
            plot_exp(ax, d2_clean, g, "TranspLeafArea", COL_QUAD[q])
            ax.set_title(f"{g}\n{QUAD_LAB.get(q,'')}", fontsize=9)
            ax.set_xlabel("VPD (kPa)")
        axes[0].set_ylabel("Transpiration (mg m⁻² s⁻¹)")
        fig.suptitle(f"{exp} — Transpiration per leaf area vs VPD  (Day 2, mean ± SE)",
                     fontsize=11, fontweight="bold")
        plt.tight_layout()
        save(f"{exp}_Fig1_actual_transpiration")
    else:
        print(f"  [skip] {exp} has no TranspLeafArea data")

    # ── Fig 2: Relative transpiration (Day 2) ───────────────────────────────
    d2_rel = d2[ok_rel[d2.index]].copy()
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 4.5), sharey=True)
    if n == 1:
        axes = [axes]
    for ax, g in zip(axes, genotypes):
        q = GENO_QUAD.get(g, "Q1")
        plot_exp(ax, d2_rel, g, "RelTranspiration", COL_QUAD[q])
        ax.axhline(1.0, color="grey", lw=0.6, ls=":", zorder=0)
        ax.set_ylim(0, 1.15)
        ax.set_title(f"{g}\n{QUAD_LAB.get(q,'')}", fontsize=9)
        ax.set_xlabel("VPD (kPa)")
    axes[0].set_ylabel("Relative transpiration (0–1)")
    fig.suptitle(f"{exp} — Relative transpiration vs VPD  (Day 2, mean ± SE)",
                 fontsize=11, fontweight="bold")
    plt.tight_layout()
    save(f"{exp}_Fig2_relative_transpiration")

    # ── Fig 3: Day 1 vs Day 2 diagnostic ────────────────────────────────────
    col3 = "TranspLeafArea" if has_tla else "RelTranspiration"
    ylabel3 = "Transpiration (mg m⁻² s⁻¹)" if has_tla else "Relative transpiration (0–1)"
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 4.5), sharey=True)
    if n == 1:
        axes = [axes]
    for ax, g in zip(axes, genotypes):
        q = GENO_QUAD.get(g, "Q1")
        for day, ls, lw in [("1", "--", 1.0), ("2", "-", 1.6)]:
            pts = sub[sub["Genotype"] == g][sub["Day"] == day].dropna(subset=[col3])
            if pts.empty:
                continue
            mn = pts.groupby("VPD")[col3].agg(["mean", "sem"]).reset_index()
            c = COL_QUAD[q]
            alpha = 0.5 if day == "1" else 0.9
            ax.scatter(pts["VPD"], pts[col3], color=c, alpha=0.15 if day == "1" else 0.30,
                       s=15, zorder=2)
            ax.plot(mn["VPD"], mn["mean"], color=c, lw=lw, ls=ls, alpha=alpha, zorder=3,
                    label=f"Day {day}")
        # Flag outliers (Day 2)
        out = d2[d2["Genotype"] == g][sub["Day"] == "2"] if "TranspLeafArea" in d2.columns else pd.DataFrame()
        try:
            out = d2[(d2["Genotype"] == g) & d2["outlier_tla"]].dropna(subset=[col3])
            if not out.empty:
                ax.scatter(out["VPD"], out[col3], marker="x", color="red", s=40, lw=1.5,
                           zorder=5, label="Outlier")
        except Exception:
            pass
        ax.axvline(VPD_BREAK, color="grey", lw=0.8, ls="--", zorder=0)
        ax.set_title(f"{g}\n{QUAD_LAB.get(q,'')}", fontsize=9)
        ax.set_xlabel("VPD (kPa)")
        ax.legend(fontsize=8)
    axes[0].set_ylabel(ylabel3)
    fig.suptitle(f"{exp} — Day 1 vs Day 2 comparison  (Day 1 dashed, Day 2 solid)",
                 fontsize=11, fontweight="bold")
    plt.tight_layout()
    save(f"{exp}_Fig3_day_comparison")

print()
print("Per-experiment figures done.")

# ═══════════════════════════════════════════════════════════════════════════════
# Overall / combined figures
# ═══════════════════════════════════════════════════════════════════════════════

# Day 2, all experiments, non-outliers
ok_tla_all = outlier_mask(df, "TranspLeafArea")
ok_rel_all  = outlier_mask(df, "RelTranspiration")
df["outlier_tla"] = ~ok_tla_all
d2_all       = df[df["Day"] == "2"].copy()
d2_all_clean = d2_all[ok_tla_all[d2_all.index]].copy()

all_genos = df.groupby(["Experiment", "Genotype", "Quadrant"]).size().reset_index()
all_genos = all_genos.sort_values(["Quadrant", "Experiment"])

# ── Overall Fig 1: Actual transpiration per leaf area — facets by quadrant ──
fig, axes = plt.subplots(1, 4, figsize=(22, 5), sharey=True)
for ax, (quad, qlab) in zip(axes, QUAD_LAB.items()):
    rows = all_genos[all_genos["Quadrant"] == quad]
    for _, row in rows.iterrows():
        g   = row["Genotype"]
        exp = row["Experiment"]
        pts = d2_all_clean[(d2_all_clean["Genotype"] == g) &
                            (d2_all_clean["Experiment"] == exp)].dropna(subset=["TranspLeafArea"])
        if pts.empty:
            continue
        mn = pts.groupby("VPD")["TranspLeafArea"].agg(["mean", "sem"]).reset_index()
        c  = COL_EXP[exp]
        ax.scatter(pts["VPD"], pts["TranspLeafArea"], color=c, alpha=0.22, s=14, zorder=2)
        ax.plot(mn["VPD"], mn["mean"], color=c, lw=1.5, zorder=3, label=f"{exp}: {g}")
        ax.fill_between(mn["VPD"], mn["mean"] - mn["sem"], mn["mean"] + mn["sem"],
                        color=c, alpha=0.15, zorder=1)
        ax.errorbar(mn["VPD"], mn["mean"], mn["sem"], fmt="o",
                    color="white", mec=c, mew=1.2, ms=5.5, ecolor=c,
                    elinewidth=0.8, zorder=4)
    ax.axvline(VPD_BREAK, color="grey", lw=0.8, ls="--", zorder=0)
    ax.set_title(qlab, fontsize=9.5, fontweight="bold")
    ax.set_xlabel("VPD (kPa)")
    ax.legend(fontsize=7.5, loc="upper left")
axes[0].set_ylabel("Transpiration (mg m⁻² s⁻¹)")
fig.suptitle("Overall — Actual transpiration per leaf area vs VPD  (Day 2, GC1 + GC2 + GC3)",
             fontsize=12, fontweight="bold")
plt.tight_layout()
save("Overall_Fig1_actual_transpiration")

# ── Overall Fig 2: Relative transpiration — facets by quadrant ───────────────
d2_all_rel = d2_all[ok_rel_all[d2_all.index]].copy()
fig, axes = plt.subplots(1, 4, figsize=(22, 5), sharey=True)
for ax, (quad, qlab) in zip(axes, QUAD_LAB.items()):
    rows = all_genos[all_genos["Quadrant"] == quad]
    for _, row in rows.iterrows():
        g   = row["Genotype"]
        exp = row["Experiment"]
        pts = d2_all_rel[(d2_all_rel["Genotype"] == g) &
                          (d2_all_rel["Experiment"] == exp)].dropna(subset=["RelTranspiration"])
        if pts.empty:
            continue
        mn = pts.groupby("VPD")["RelTranspiration"].agg(["mean", "sem"]).reset_index()
        c  = COL_EXP[exp]
        ax.scatter(pts["VPD"], pts["RelTranspiration"], color=c, alpha=0.22, s=14, zorder=2)
        ax.plot(mn["VPD"], mn["mean"], color=c, lw=1.5, zorder=3, label=f"{exp}: {g}")
        ax.fill_between(mn["VPD"], mn["mean"] - mn["sem"], mn["mean"] + mn["sem"],
                        color=c, alpha=0.15, zorder=1)
        ax.errorbar(mn["VPD"], mn["mean"], mn["sem"], fmt="o",
                    color="white", mec=c, mew=1.2, ms=5.5, ecolor=c,
                    elinewidth=0.8, zorder=4)
    ax.axvline(VPD_BREAK, color="grey", lw=0.8, ls="--", zorder=0)
    ax.axhline(1.0, color="grey", lw=0.6, ls=":", zorder=0)
    ax.set_ylim(0, 1.2)
    ax.set_title(qlab, fontsize=9.5, fontweight="bold")
    ax.set_xlabel("VPD (kPa)")
    ax.legend(fontsize=7.5, loc="upper left")
axes[0].set_ylabel("Relative transpiration (0–1)")
fig.suptitle("Overall — Relative transpiration vs VPD  (Day 2, GC1 + GC2 + GC3)",
             fontsize=12, fontweight="bold")
plt.tight_layout()
save("Overall_Fig2_relative_transpiration")

# ── Overall Fig 3: All genotypes in one panel — coloured by quadrant ─────────
fig, ax = plt.subplots(figsize=(12, 6))
for _, row in all_genos.iterrows():
    g   = row["Genotype"]
    exp = row["Experiment"]
    q   = row["Quadrant"]
    pts = d2_all_clean[(d2_all_clean["Genotype"] == g) &
                        (d2_all_clean["Experiment"] == exp)].dropna(subset=["TranspLeafArea"])
    if pts.empty:
        continue
    mn = pts.groupby("VPD")["TranspLeafArea"].agg(["mean", "sem"]).reset_index()
    c  = COL_QUAD[q]
    ls = {"GC1": "-", "GC2": "--", "GC3": ":"}.get(exp, "-")
    ax.plot(mn["VPD"], mn["mean"], color=c, lw=1.6, ls=ls, label=f"{exp}: {g} ({q})")
    ax.errorbar(mn["VPD"], mn["mean"], mn["sem"], fmt="o",
                color="white", mec=c, mew=1.1, ms=5, ecolor=c, elinewidth=0.8, zorder=4)
ax.axvline(VPD_BREAK, color="grey", lw=0.9, ls="--", zorder=0)
ax.annotate("Sinclair 2008\nVPD threshold", xy=(VPD_BREAK + 0.05, ax.get_ylim()[1] * 0.95),
            fontsize=8, color="grey")
# Quadrant legend patches
quad_patches = [mpatches.Patch(color=v, label=QUAD_LAB[k]) for k, v in COL_QUAD.items()]
exp_lines = [plt.Line2D([0], [0], color="grey", lw=1.5, ls=ls, label=exp)
             for exp, ls in [("GC1", "-"), ("GC2", "--"), ("GC3", ":")]]
leg1 = ax.legend(handles=quad_patches, loc="upper left", title="Quadrant", fontsize=8.5)
ax.add_artist(leg1)
ax.legend(handles=exp_lines, loc="center left", title="Experiment", fontsize=8.5)
ax.set_xlabel("VPD (kPa)", fontsize=11)
ax.set_ylabel("Transpiration (mg m⁻² s⁻¹)", fontsize=11)
ax.set_title("All genotypes — Actual transpiration per leaf area vs VPD  (Day 2)",
             fontsize=12, fontweight="bold")
plt.tight_layout()
save("Overall_Fig3_all_genotypes")

# ── Overall Fig 4: GC1 vs GC2 vs GC3 mean by quadrant ──────────────────────
fig, axes = plt.subplots(1, 4, figsize=(22, 5), sharey=True)
for ax, (quad, qlab) in zip(axes, QUAD_LAB.items()):
    for exp in ["GC1", "GC2", "GC3"]:
        pts = d2_all_clean[(d2_all_clean["Quadrant"] == quad) &
                            (d2_all_clean["Experiment"] == exp)].dropna(subset=["TranspLeafArea"])
        if pts.empty:
            continue
        mn = pts.groupby("VPD")["TranspLeafArea"].agg(["mean", "sem"]).reset_index()
        c  = COL_EXP[exp]
        ax.plot(mn["VPD"], mn["mean"], color=c, lw=1.8, label=exp)
        ax.fill_between(mn["VPD"], mn["mean"] - mn["sem"], mn["mean"] + mn["sem"],
                        color=c, alpha=0.18)
        ax.errorbar(mn["VPD"], mn["mean"], mn["sem"], fmt="o",
                    color="white", mec=c, mew=1.2, ms=6, ecolor=c, elinewidth=0.9)
    ax.axvline(VPD_BREAK, color="grey", lw=0.8, ls="--")
    ax.set_title(qlab, fontsize=9.5, fontweight="bold")
    ax.set_xlabel("VPD (kPa)")
    ax.legend(fontsize=8.5)
axes[0].set_ylabel("Transpiration (mg m⁻² s⁻¹)")
fig.suptitle("Overall — Experiment comparison by quadrant  (Day 2, mean ± SE)",
             fontsize=12, fontweight="bold")
plt.tight_layout()
save("Overall_Fig4_experiment_comparison")

print()
print("Overall figures done.")
print(f"\nAll figures written to {FIGDIR}/")
