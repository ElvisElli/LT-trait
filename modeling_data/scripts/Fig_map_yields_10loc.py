#!/usr/bin/env python3
"""
Fig_map_yields_10loc.py
Water-limited yield vs potential yield at the 10 updated SSM locations.

Panel A: pie per site (green = water-limited/rainfed yield as a share of the
         irrigated potential, yellow = yield gap by water deficit), one map for
         the non-LT genotype and one for the LT1.5 genotype.
Panel B: mean attainable yield + yield gap across all locations.

Input : modeling_data/data/updated-outputs/all_yearly.csv
Output: modeling_data/figures/Fig_map_yields_10loc.png (A + B)
        modeling_data/figures/Fig_map_yields_10loc_A.png, _B.png
        modeling_data/data/yield_gap_10loc_summary.csv
Base map: a US states GeoJSON (path via --geojson; default downloads none).

Run from the repo root:
    python3 modeling_data/scripts/Fig_map_yields_10loc.py --geojson us-states.json
"""
import os, json, argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Wedge, Rectangle

ap = argparse.ArgumentParser()
ap.add_argument("--geojson", required=True, help="US states GeoJSON with a 'name' property")
ap.add_argument("--yearly", default="modeling_data/data/updated-outputs/all_yearly.csv")
ap.add_argument("--outdir", default="modeling_data/figures")
args = ap.parse_args()

GREEN, YELLOW = "#31a354", "#ffffb2"
SITES = pd.DataFrame({
    "Location": ["Jonesboro", "Keiser", "Marianna", "Rohwer", "Albany", "Mount Vernon", "Novelty",
                 "Eustis", "Lincoln", "North Platte"],
    "lat": [35.84, 35.67, 34.77, 33.80, 40.25, 37.10, 40.01, 40.66, 40.81, 41.12],
    "lon": [-90.70, -90.09, -90.76, -91.27, -94.33, -93.82, -92.21, -100.03, -96.70, -100.77],
})

# ---------------------------------------------------------------- data
y = pd.read_csv(args.yearly)
y["Location"] = y["Location"].str.replace("­", "", regex=False).str.strip()
y["Yield_kgha"] = pd.to_numeric(y["Yield_kgha"], errors="coerce")
y = y[y["Genotype_class"].isin(["Non-LT", "LT1.5"])].copy()
y["Genotype"] = y["Genotype_class"].map({"Non-LT": "Non-limited transpiration trait (Non-LT)",
                                         "LT1.5": "Limited transpiration trait (LT)"})
site = (y.groupby(["Location", "Genotype", "water_mgmt"])["Yield_kgha"].mean()
        .unstack("water_mgmt").reset_index())
site["WL_share"] = site["Rainfed"] / site["Irrigated"]
site["Gap"] = site["Irrigated"] - site["Rainfed"]
site = site.merge(SITES, on="Location")
site.to_csv("modeling_data/data/yield_gap_10loc_summary.csv", index=False)

allloc = y.groupby(["Genotype", "water_mgmt", "Location"])["Yield_kgha"].mean().unstack("water_mgmt")
summ = allloc.groupby("Genotype").agg(Rainfed=("Rainfed", "mean"), Irrigated=("Irrigated", "mean"),
                                      Rainfed_se=("Rainfed", lambda s: s.std() / np.sqrt(len(s))),
                                      Irrigated_se=("Irrigated", lambda s: s.std() / np.sqrt(len(s))))
summ["Gap_%"] = 100 * (summ.Irrigated - summ.Rainfed) / summ.Irrigated
print(summ.round(1))
nl, lt = summ.loc["Non-limited transpiration trait (Non-LT)"], summ.loc["Limited transpiration trait (LT)"]
print(f"LT water-limited yield vs Non-LT: {100*(lt.Rainfed/nl.Rainfed-1):+.1f}%")
print(f"Non-LT potential yield vs LT:     {100*(nl.Irrigated/lt.Irrigated-1):+.1f}%")

# ---------------------------------------------------------------- base map
gj = json.load(open(args.geojson))
FOCUS = {"Arkansas", "Missouri", "Nebraska"}


def polys(feature):
    geom = feature["geometry"]
    if geom["type"] == "Polygon":
        return [geom["coordinates"][0]]
    return [p[0] for p in geom["coordinates"]]


def draw_base(ax):
    for f in gj["features"]:
        name = f["properties"].get("name")
        focus = name in FOCUS
        for ring in polys(f):
            ax.add_patch(Polygon(ring, closed=True, facecolor="#9e9e9e" if focus else "#e8e8e8",
                                 edgecolor="black" if focus else "#c8c8c8", linewidth=1.6 if focus else 0.6,
                                 zorder=2 if focus else 1))
    ax.set_xlim(-104.5, -87.5); ax.set_ylim(32.5, 43.5)
    ax.set_aspect(1.25); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_linewidth(1.2)
    ax.set_facecolor("#f4f4f4")


def draw_pies(ax, d):
    rmax = 0.8
    for _, r in d.iterrows():
        rad = rmax * np.sqrt(r.Irrigated / site.Irrigated.max())
        ang = 90 - 360 * r.WL_share
        ax.add_patch(Wedge((r.lon, r.lat), rad, ang, 90, facecolor=GREEN, edgecolor="black", lw=1.0, zorder=5))
        ax.add_patch(Wedge((r.lon, r.lat), rad, 90, ang + 360, facecolor=YELLOW, edgecolor="black", lw=1.0, zorder=5))


def panel_A(fig_or_axes):
    axes = fig_or_axes
    for ax, geno in zip(axes, ["Non-limited transpiration trait (Non-LT)", "Limited transpiration trait (LT)"]):
        draw_base(ax); draw_pies(ax, site[site.Genotype == geno])
        ax.set_title(geno, fontsize=13, pad=6, backgroundcolor="white",
                     bbox=dict(boxstyle="square,pad=0.4", fc="white", ec="black", lw=1.2))


# ---------------------------------------------------------------- Panel A figure
figA, axes = plt.subplots(1, 2, figsize=(12.5, 6.0))
panel_A(axes)
figA.legend(handles=[Rectangle((0, 0), 1, 1, fc=GREEN, ec="black"), Rectangle((0, 0), 1, 1, fc=YELLOW, ec="black")],
            labels=["Water-limited yield", "Yield gap by water deficit"], loc="upper center", ncol=2,
            frameon=False, fontsize=12, bbox_to_anchor=(0.5, 1.0))
figA.text(0.01, 0.97, "A", fontsize=16, fontweight="bold")
figA.tight_layout(rect=(0, 0, 1, 0.94))
figA.savefig(f"{args.outdir}/Fig_map_yields_10loc_A.png", dpi=300, bbox_inches="tight"); plt.close(figA)

# ---------------------------------------------------------------- Panel B figure
figB, ax = plt.subplots(figsize=(4.2, 6.0))
labels = ["Non-LT", "LT"]
rows = [nl, lt]
for i, (lab, r) in enumerate(zip(labels, rows)):
    ax.bar(i, r.Rainfed, color=GREEN, edgecolor="black", width=0.6, zorder=3)
    ax.bar(i, r.Irrigated - r.Rainfed, bottom=r.Rainfed, color=YELLOW, edgecolor="black", width=0.6, zorder=3)
    ax.errorbar(i, r.Rainfed, yerr=r.Rainfed_se, color="black", capsize=0, lw=1.5, zorder=4)
    ax.errorbar(i, r.Irrigated, yerr=r.Irrigated_se, color="black", capsize=0, lw=1.5, zorder=4)
ax.set_xticks([0, 1]); ax.set_xticklabels(labels, fontsize=13)
ax.set_ylabel("Yield (kg ha$^{-1}$)", fontsize=13)
ax.set_ylim(0, 8000); ax.set_yticks([0, 2000, 4000, 6000]); ax.set_yticklabels(["0", "2,000", "4,000", "6,000"])
ax.set_title("All locations", fontsize=13, pad=6, bbox=dict(boxstyle="square,pad=0.4", fc="white", ec="black", lw=1.2))
ax.grid(axis="y", color="#dddddd", lw=0.8, zorder=0)
for s in ax.spines.values():
    s.set_linewidth(1.2)
figB.text(0.01, 0.97, "B", fontsize=16, fontweight="bold")
figB.tight_layout(rect=(0, 0, 1, 0.95))
figB.savefig(f"{args.outdir}/Fig_map_yields_10loc_B.png", dpi=300, bbox_inches="tight"); plt.close(figB)

# ---------------------------------------------------------------- combined
fig = plt.figure(figsize=(17, 6.0))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.42], wspace=0.08)
axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])]
panel_A(axes)
axb = fig.add_subplot(gs[0, 2])
for i, (lab, r) in enumerate(zip(labels, rows)):
    axb.bar(i, r.Rainfed, color=GREEN, edgecolor="black", width=0.6, zorder=3)
    axb.bar(i, r.Irrigated - r.Rainfed, bottom=r.Rainfed, color=YELLOW, edgecolor="black", width=0.6, zorder=3)
    axb.errorbar(i, r.Rainfed, yerr=r.Rainfed_se, color="black", lw=1.5, zorder=4)
    axb.errorbar(i, r.Irrigated, yerr=r.Irrigated_se, color="black", lw=1.5, zorder=4)
axb.set_xticks([0, 1]); axb.set_xticklabels(labels, fontsize=13); axb.set_ylabel("Yield (kg ha$^{-1}$)", fontsize=13)
axb.set_ylim(0, 8000); axb.set_yticks([0, 2000, 4000, 6000]); axb.set_yticklabels(["0", "2,000", "4,000", "6,000"])
axb.set_title("All locations", fontsize=13, pad=6, bbox=dict(boxstyle="square,pad=0.4", fc="white", ec="black", lw=1.2))
axb.grid(axis="y", color="#dddddd", lw=0.8, zorder=0)
fig.legend(handles=[Rectangle((0, 0), 1, 1, fc=GREEN, ec="black"), Rectangle((0, 0), 1, 1, fc=YELLOW, ec="black")],
           labels=["Water-limited yield", "Yield gap by water deficit"], loc="upper center", ncol=2,
           frameon=False, fontsize=12, bbox_to_anchor=(0.4, 1.0))
fig.text(0.01, 0.95, "A", fontsize=16, fontweight="bold"); fig.text(0.72, 0.95, "B", fontsize=16, fontweight="bold")
fig.savefig(f"{args.outdir}/Fig_map_yields_10loc.png", dpi=300, bbox_inches="tight"); plt.close(fig)
print("Saved figures to", args.outdir)
