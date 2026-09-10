#!/usr/bin/env python3
"""
Generate Figure 1 (9-panel and 14-panel) with new Rep 5 data marked in red.

Uses the new data_all_9_10_2026.xlsx with model fitting for VPD response curves.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rcParams
from scipy.optimize import curve_fit
from scipy.stats import linregress
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# Settings
# ============================================================================
DATA_FILE = "/home/user/LT-trait/experimental_data_soy_2026/data/data_all_9_10_2026.xlsx"
OUTPUT_DIR = "."
DPI = 150

# Color system
COLOR_OLD_REP = "#31a354"  # Green for reps 1-4
COLOR_NEW_REP = "#e41a1c"  # Red for rep 5
COLOR_CURVE = "#4393c3"    # Blue for fitted curves

# 9 genotypes for first figure (custom order)
GENOTYPES_9 = ["PI398492", "PI398617", "PI424231", "PI603418D",
               "PI377574", "PI603171", "PI398198", "A5959", "PI471938"]

# ============================================================================
# Model functions
# ============================================================================

def linear_fn(x, a, b):
    return a + b * x

def bilinear_fn(x, a, b1, b2, bp):
    """Bilinear with breakpoint"""
    return np.where(x <= bp,
                   a + b1 * x,
                   a + b1 * bp + b2 * (x - bp))

def lin_plateau_fn(x, a, b, bp):
    """Linear-plateau: rises then flattens"""
    return np.where(x <= bp,
                   a + b * x,
                   a + b * bp)

# ============================================================================
# Model fitting
# ============================================================================

def fit_linear(x, y):
    """Fit linear model"""
    try:
        # Use robust linear regression via numpy polyfit
        popt, _ = curve_fit(linear_fn, x, y, p0=[np.mean(y), 1.0])
        return {'params': {'a': popt[0], 'b': popt[1]}, 'success': True}
    except:
        return {'success': False}

def fit_bilinear(x, y):
    """Fit bilinear model with starting values"""
    try:
        x_range = x.max() - x.min()
        bp_init = x.min() + x_range * 0.5
        b_init = (y.max() - y.min()) / x_range
        popt, _ = curve_fit(bilinear_fn, x, y,
                           p0=[y.min(), b_init, b_init * 0.1, bp_init],
                           maxfev=5000, ftol=1e-5)
        return {'params': {'a': popt[0], 'b1': popt[1], 'b2': popt[2], 'bp': popt[3]}, 'success': True}
    except:
        return {'success': False}

def fit_lin_plateau(x, y):
    """Fit linear-plateau model"""
    try:
        x_range = x.max() - x.min()
        bp_init = x.min() + x_range * 0.5
        b_init = (y.max() - y.min()) / x_range
        popt, _ = curve_fit(lin_plateau_fn, x, y,
                           p0=[y.min(), b_init, bp_init],
                           maxfev=5000, ftol=1e-5)
        return {'params': {'a': popt[0], 'b': popt[1], 'bp': popt[2]}, 'success': True}
    except:
        return {'success': False}

def calc_aic(y_obs, y_pred, k):
    """Calculate AIC: 2k + n*ln(RSS/n)"""
    n = len(y_obs)
    rss = np.sum((y_obs - y_pred) ** 2)
    if rss <= 0:
        return np.inf
    return 2 * k + n * np.log(rss / n)

def fit_all_models(x, y):
    """Fit all three models and return best via AIC"""
    results = {}

    # Linear
    lin_fit = fit_linear(x, y)
    if lin_fit['success']:
        y_pred = linear_fn(x, lin_fit['params']['a'], lin_fit['params']['b'])
        results['Linear'] = {
            'params': lin_fit['params'],
            'aic': calc_aic(y, y_pred, 2),
            'y_pred': y_pred
        }

    # Bilinear
    bil_fit = fit_bilinear(x, y)
    if bil_fit['success']:
        p = bil_fit['params']
        y_pred = bilinear_fn(x, p['a'], p['b1'], p['b2'], p['bp'])
        results['Bilinear'] = {
            'params': bil_fit['params'],
            'aic': calc_aic(y, y_pred, 4),
            'y_pred': y_pred
        }

    # Linear-plateau
    lp_fit = fit_lin_plateau(x, y)
    if lp_fit['success']:
        p = lp_fit['params']
        y_pred = lin_plateau_fn(x, p['a'], p['b'], p['bp'])
        results['Lin. plateau'] = {
            'params': lp_fit['params'],
            'aic': calc_aic(y, y_pred, 3),
            'y_pred': y_pred
        }

    if not results:
        return None

    # Select best model by AIC
    best_model = min(results.keys(), key=lambda k: results[k]['aic'])
    return {
        'best_model': best_model,
        'aic': {k: results[k]['aic'] for k in results},
        'params': results[best_model]['params'],
        'y_pred': results[best_model]['y_pred'],
        'all_models': results
    }

def predict_curve(params, model_type, vpd_seq):
    """Generate prediction curve"""
    if model_type == 'Linear':
        return linear_fn(vpd_seq, params['a'], params['b'])
    elif model_type == 'Bilinear':
        return bilinear_fn(vpd_seq, params['a'], params['b1'], params['b2'], params['bp'])
    elif model_type == 'Lin. plateau':
        return lin_plateau_fn(vpd_seq, params['a'], params['b'], params['bp'])
    return np.zeros_like(vpd_seq)

# ============================================================================
# Load and prepare data
# ============================================================================

print("Loading data...")
df = pd.read_excel(DATA_FILE, sheet_name='dataAll')
print(f"  Loaded {len(df)} rows")

# Rename columns
df = df.rename(columns={
    'RepNo': 'Rep',
    'day': 'Day',
    'TR_mg_sqm_seg': 'TR',
    'avg_vpd': 'VPD'
})

df['IsNewRep'] = (df['Rep'] == 5)
day1 = df[df['Day'] == 1].copy()

# Identify genotypes with >2 reps
genotype_counts = day1.groupby('Genotype')['Rep'].nunique()
genotypes_9panel = [g for g in GENOTYPES_9 if genotype_counts.get(g, 0) > 2]
genotypes_all = sorted([g for g in genotype_counts[genotype_counts > 2].index])

print(f"  Day 1 data: {len(day1)} observations")
print(f"  Genotypes (>2 reps): {len(genotypes_all)}")
print(f"  9-panel genotypes available: {len(genotypes_9panel)}/{len(GENOTYPES_9)}")
print(f"  New Rep 5 points: {day1['IsNewRep'].sum()}")

# Prepare VPD sequence for curves
vpd_seq = np.linspace(day1['VPD'].min() * 0.95, day1['VPD'].max() * 1.05, 300)

# ============================================================================
# Fit models for each genotype
# ============================================================================

print("\nFitting models...")
model_results = {}
aic_table = []

for genotype in genotypes_all:
    geno_data = day1[day1['Genotype'] == genotype]
    x = geno_data['VPD'].values
    y = geno_data['TR'].values

    if len(x) < 4:
        continue

    # Remove NaN
    mask = ~(np.isnan(x) | np.isnan(y))
    x = x[mask]
    y = y[mask]

    fit = fit_all_models(x, y)
    if fit:
        model_results[genotype] = fit
        aic_table.append({
            'Genotype': genotype,
            'Linear': fit['aic'].get('Linear', np.nan),
            'Bilinear': fit['aic'].get('Bilinear', np.nan),
            'Lin.plateau': fit['aic'].get('Lin. plateau', np.nan),
            'Best': fit['best_model']
        })
        print(f"  {genotype:15s}: {fit['best_model']:15s} (AIC: {fit['aic'][fit['best_model']]:.1f})")

# ============================================================================
# Generate Figure 1 (9-panel)
# ============================================================================

print("\nGenerating Figure 1 (9-panel)...")

fig, axes = plt.subplots(3, 3, figsize=(12, 12), dpi=DPI)
axes = axes.flatten()

for idx, genotype in enumerate(genotypes_9panel):
    if genotype not in model_results:
        print(f"  Skipping {genotype} (no model)")
        continue

    ax = axes[idx]
    geno_data = day1[day1['Genotype'] == genotype]
    fit = model_results[genotype]

    # Plot old reps (green) and new rep (red)
    old_data = geno_data[~geno_data['IsNewRep']]
    new_data = geno_data[geno_data['IsNewRep']]

    ax.scatter(old_data['VPD'], old_data['TR'], s=100, marker='o',
               facecolor=COLOR_OLD_REP, edgecolor='black', linewidth=0.8,
               alpha=0.9, label='Reps 1-4', zorder=3)

    if len(new_data) > 0:
        ax.scatter(new_data['VPD'], new_data['TR'], s=100, marker='o',
                   facecolor=COLOR_NEW_REP, edgecolor='black', linewidth=0.8,
                   alpha=0.9, label='Rep 5 (NEW)', zorder=4)

    # Plot fitted curve
    y_curve = predict_curve(fit['params'], fit['best_model'], vpd_seq)
    y_curve = np.maximum(y_curve, 0)  # Clamp to 0
    ax.plot(vpd_seq, y_curve, color=COLOR_CURVE, linewidth=1.2, zorder=2)

    ax.set_xlim(-0.1, day1['VPD'].max() * 1.05)
    ax.set_ylim(0, 120)
    ax.set_xlabel('Vapor Pressure Deficit (kPa)', fontsize=11)
    ax.set_ylabel('Transpiration (mg m⁻² s⁻¹)', fontsize=11)
    ax.set_title(f'{genotype}\n[{fit["best_model"]}]', fontsize=10, fontweight='bold')
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# Hide unused subplots
for idx in range(len(genotypes_9panel), len(axes)):
    axes[idx].set_visible(False)

plt.tight_layout()
fig.savefig(f'{OUTPUT_DIR}/Figure1_9panel_FITTED_NEW_REP_RED.png', dpi=DPI, bbox_inches='tight')
print(f"  Saved: Figure1_9panel_FITTED_NEW_REP_RED.png")

plt.close()

# ============================================================================
# Generate Figure 1B (14-panel - all genotypes)
# ============================================================================

print("\nGenerating Figure 1B (14-panel - all genotypes)...")

n_genos = len(genotypes_all)
n_cols = 7
n_rows = (n_genos + n_cols - 1) // n_cols

fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 6), dpi=DPI)
if n_rows == 1:
    axes = axes.reshape(1, -1)
axes = axes.flatten()

for idx, genotype in enumerate(genotypes_all):
    if genotype not in model_results:
        print(f"  Skipping {genotype} (no model)")
        continue

    ax = axes[idx]
    geno_data = day1[day1['Genotype'] == genotype]
    fit = model_results[genotype]

    # Plot old reps (green) and new rep (red)
    old_data = geno_data[~geno_data['IsNewRep']]
    new_data = geno_data[geno_data['IsNewRep']]

    ax.scatter(old_data['VPD'], old_data['TR'], s=80, marker='o',
               facecolor=COLOR_OLD_REP, edgecolor='black', linewidth=0.7,
               alpha=0.9, label='Reps 1-4', zorder=3)

    if len(new_data) > 0:
        ax.scatter(new_data['VPD'], new_data['TR'], s=80, marker='o',
                   facecolor=COLOR_NEW_REP, edgecolor='black', linewidth=0.7,
                   alpha=0.9, label='Rep 5 (NEW)', zorder=4)

    # Plot fitted curve
    y_curve = predict_curve(fit['params'], fit['best_model'], vpd_seq)
    y_curve = np.maximum(y_curve, 0)
    ax.plot(vpd_seq, y_curve, color=COLOR_CURVE, linewidth=1.0, zorder=2)

    ax.set_xlim(-0.1, day1['VPD'].max() * 1.05)
    ax.set_ylim(0, 120)
    ax.set_xlabel('Vapor Pressure Deficit (kPa)', fontsize=9)
    ax.set_ylabel('Transpiration (mg m⁻² s⁻¹)', fontsize=9)
    ax.set_title(f'{genotype}\n[{fit["best_model"]}]', fontsize=8, fontweight='bold')
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
    ax.tick_params(labelsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# Hide unused subplots
for idx in range(len(genotypes_all), len(axes)):
    axes[idx].set_visible(False)

plt.tight_layout()
fig.savefig(f'{OUTPUT_DIR}/Figure1B_14genotypes_FITTED_NEW_REP_RED.png', dpi=DPI, bbox_inches='tight')
print(f"  Saved: Figure1B_14genotypes_FITTED_NEW_REP_RED.png")

plt.close()

# ============================================================================
# Summary
# ============================================================================

print("\n" + "="*60)
print("ANALYSIS SUMMARY")
print("="*60)

aic_df = pd.DataFrame(aic_table)
print("\nModel AIC comparison:")
print(aic_df.to_string(index=False))

print(f"\nGenotypes by best model:")
for model in ['Linear', 'Bilinear', 'Lin. plateau']:
    genos = aic_df[aic_df['Best'] == model]['Genotype'].tolist()
    print(f"  {model:15s}: {len(genos)} genotypes")

print(f"\nData changes with Rep 5:")
old_total = len(day1[~day1['IsNewRep']])
new_total = len(day1[day1['IsNewRep']])
print(f"  Previous observations (Reps 1-4): {old_total}")
print(f"  New observations (Rep 5): {new_total}")
print(f"  Total increase: {new_total} observations (+{100*new_total/(old_total):.1f}%)")

print("\n✓ Figures generated successfully!")
