#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_distributions.py

Changes in this version
-----------------------

Notes
-----
• Adjusted p-values are NOT drawn on plots.
• Input columns expected (per TSV):

[
    "seqname", "chr", "repeat_length", "PAR", "CEN",
    "global_PID_maskN", "global_P_maskN", "global_Q_maskN", "global_kimura2p_maskN",
    "global_PID_countN", "global_P_countN", "global_Q_countN", "global_kimura2p_countN",
    "takki_PID_maskN", "takki_P_maskN", "takki_Q_maskN", "takki_kimura2p_maskN",
    "takki_PID_countN", "takki_P_countN", "takki_Q_countN", "takki_kimura2p_countN",
    "chr_PID_maskN", "chr_P_maskN", "chr_Q_maskN", "chr_kimura2p_maskN",
    "chr_PID_countN", "chr_P_countN", "chr_Q_countN", "chr_kimura2p_countN",
]

Output
------
Plots under $wd/plots/:
  plots/repeat_length_*.{pdf,svg,png}
  plots/{global,takki,chr}/*.{pdf,svg,png}

Stats under $wd/plots/:
  plots/Tgut716A_stats.tsv
  plots/Tgut191A_stats.tsv
"""

import os
import math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # headless backend
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from scipy.stats import gaussian_kde, mannwhitneyu

# ───────────────────────────────────────────── Config / paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
wd = os.path.normpath(os.path.join(SCRIPT_DIR, ".."))
in_191 = os.path.join(wd, "Tgut191A_pid_kimura.tsv")
in_716 = os.path.join(wd, "Tgut716A_pid_kimura.tsv")

plot_root = os.path.join(wd, "plots")
os.makedirs(plot_root, exist_ok=True)
for sub in ("global", "takki", "chr"):
    os.makedirs(os.path.join(plot_root, sub), exist_ok=True)

families = ["Tgut716A", "Tgut191A"]  # plotting order
palette = {"Tgut191A": "#009E73", "Tgut716A": "#CC79A7"}

# Family-specific length cutoffs (ignore repeats shorter than these values)
# Set to >0 to filter; 0 disables filtering.
Tgut716A_cutoff = 0
Tgut191A_cutoff = 0

# ---- toggles
RUN_RAINCLOUDS = True
# Produce both grouped (4 panels) and boolean-highlight (2 panels)
PLOT_GROUPED = True         # boolean_grouped = False case
PLOT_BOOLEAN_HIGHLIGHT = True  # boolean_grouped = True case

# ── Tunables (you can tweak these to control geometry/widths) ────────────────
VIOLIN_WIDTH_GROUPED = 0.35     # half-violin width for 4-group panels
VIOLIN_WIDTH_BOOL    = 0.32     # half-violin width for 2-group panels  ← thinner
SCATTER_JITTER_SD_GROUPED = 0.05
SCATTER_JITTER_SD_BOOL    = 0.045  # tighter spread for 2-group panels   ← thinner
SCATTER_X_OFFSET_GROUPED  = 0.30
SCATTER_X_OFFSET_BOOL     = 0.30  # draw points closer to center        ← thinner
SCATTER_SIZE_GROUPED      = 18
SCATTER_SIZE_BOOL         = 15    # slightly smaller points in 2-group panels
FIGSIZE_GROUPED           = (4, 6)
FIGSIZE_BOOL              = (4, 6)  # narrower canvas for 2-group panels

# >>> NEW: keep vectors for violins/boxes, rasterize only the heavy scatter layers
RASTERIZE_SCATTER = True
RASTERIZE_MIN_POINTS = 2000  # rasterize a group's scatter only if it has ≥ this many points

# # Variables to plot
# comparisons = ["global", "takki", "chr"]
# methods = ["maskN", "countN"]

# pid_vars = [f"{c}_PID_{m}" for c in comparisons for m in methods]          # 3×2 = 6
# k2p_vars = [f"{c}_kimura2p_{m}" for c in comparisons for m in methods]     # 3×2 = 6
# all_vars = ["repeat_length"] + pid_vars + k2p_vars

# ───────────────────────────────────────────── Fast plotting
# Variables to plot
comparisons = ["chr"]
methods = ["countN"]
k2p_vars = [f"{c}_kimura2p_{m}" for c in comparisons for m in methods]     # 2×1 = 2
all_vars = k2p_vars


# ───────────────────────────────────────────── Helpers

def lighten_color(hex_color, amount=0.55):
    """Pastel version by mixing with white."""
    rgb = np.array(mcolors.to_rgb(hex_color))
    white = np.array([1.0, 1.0, 1.0])
    mixed = rgb * (1 - amount) + white * amount
    return mcolors.to_hex(mixed, keep_alpha=False)

def bh_adjust(pvals):
    """
    Benjamini–Hochberg FDR adjustment.
    pvals: list of floats (assumed finite here). Returns list of adj p in original order.
    """
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = [pvals[i] for i in order]
    adj = [np.nan] * n
    running_min = 1.0
    m = n
    for k in reversed(range(m)):
        rank = k + 1  # 1-based
        pv = ranked[k]
        q = pv * m / rank
        running_min = min(running_min, q)
        adj_val = min(1.0, running_min)
        adj[order[k]] = adj_val
    return adj

def savefig(fig, outpath_stem):
    for ext in ("pdf", "svg", "png"):
        fig.savefig(f"{outpath_stem}.{ext}", dpi=600, bbox_inches="tight")

def _format_values(vals, precision=6):
    return "[" + ", ".join(f"{float(x):.{precision}g}" for x in vals) + "]"

def _median(arr):
    return float(np.median(arr)) if len(arr) else float("nan")

def _iqr(arr):
    if len(arr) == 0:
        return float("nan")
    q1, q3 = np.percentile(arr, [25, 75])
    return float(q3 - q1)

def _cliffs_delta_from_u(u, n, m):
    # rank-biserial correlation == Cliff's delta
    return float(2.0 * u / (n * m) - 1.0) if (n > 0 and m > 0) else float("nan")

def _ylabel_for(var: str) -> str:
    """Minimal, de-redundant Y-axis labels per variable."""
    if var.endswith("_kimura2p_maskN") or var.endswith("_kimura2p_countN"):
        return "Kimura 2p"
    if var.endswith("_PID_maskN") or var.endswith("_PID_countN"):
        return "% Identity"
    if var == "repeat_length":
        return "Repeat length (bp)"
    return var

# ───────────────────────────────────────────── Data

df191 = pd.read_csv(in_191, sep="\t")
df191 = df191[df191["repeat_length"] >= Tgut191A_cutoff].copy()
df191["family"] = "Tgut191A"

df716 = pd.read_csv(in_716, sep="\t")
df716 = df716[df716["repeat_length"] >= Tgut716A_cutoff].copy()
df716["family"] = "Tgut716A"

df = pd.concat([df716, df191], ignore_index=True)

# bool_flag used only for boolean_highlight plots:
#   • Tgut716A: CEN
#   • Tgut191A: PAR
df["bool_flag"] = np.where(df["family"].eq("Tgut716A"), df["CEN"], df["PAR"]).astype(bool)

# ───────────────────────────────────────────── Stats (compute & save)

def compute_pair_pvals(var):
    """Return raw p-values for three comparisons on a variable."""
    # 4.1 Tgut716A vs Tgut191A
    a = df.loc[df["family"].eq("Tgut716A"), var].dropna().values
    b = df.loc[df["family"].eq("Tgut191A"), var].dropna().values
    p1 = mannwhitneyu(a, b, alternative="two-sided").pvalue if (len(a) and len(b)) else np.nan

    # 4.2 Tgut716A CEN vs non-CEN
    a1 = df.loc[(df["family"].eq("Tgut716A")) & (df["CEN"].astype(bool)), var].dropna().values
    a0 = df.loc[(df["family"].eq("Tgut716A")) & (~df["CEN"].astype(bool)), var].dropna().values
    p2 = mannwhitneyu(a1, a0, alternative="two-sided").pvalue if (len(a1) and len(a0)) else np.nan

    # 4.3 Tgut191A PAR vs non-PAR
    b1 = df.loc[(df["family"].eq("Tgut191A")) & (df["PAR"].astype(bool)), var].dropna().values
    b0 = df.loc[(df["family"].eq("Tgut191A")) & (~df["PAR"].astype(bool)), var].dropna().values
    p3 = mannwhitneyu(b1, b0, alternative="two-sided").pvalue if (len(b1) and len(b0)) else np.nan

    return p1, p2, p3

def _mean(arr):
    return float(np.mean(arr)) if len(arr) else float("nan")

# Gather raw p-values across variables for BH per comparison
raw_p_716_vs_191 = []
raw_p_716_cen_vs_non = []
raw_p_191_par_vs_non = []
for v in all_vars:
    p1, p2, p3 = compute_pair_pvals(v)
    raw_p_716_vs_191.append(p1)
    raw_p_716_cen_vs_non.append(p2)
    raw_p_191_par_vs_non.append(p3)

# BH-adjusted per comparison set
adj_716_vs_191 = bh_adjust(raw_p_716_vs_191)
adj_716_cen_vs_non = bh_adjust(raw_p_716_cen_vs_non)
adj_191_par_vs_non = bh_adjust(raw_p_191_par_vs_non)

# Assemble means/medians/IQR and effect sizes per variable
rows_716 = []
rows_191 = []
for i, v in enumerate(all_vars):
    # 716A sets
    vals_716       = df.loc[df["family"].eq("Tgut716A"), v].dropna().values
    vals_716_cen   = df.loc[(df["family"].eq("Tgut716A")) & (df["CEN"].astype(bool)), v].dropna().values
    vals_716_non   = df.loc[(df["family"].eq("Tgut716A")) & (~df["CEN"].astype(bool)), v].dropna().values
    # 191A sets
    vals_191       = df.loc[df["family"].eq("Tgut191A"), v].dropna().values
    vals_191_par   = df.loc[(df["family"].eq("Tgut191A")) & (df["PAR"].astype(bool)), v].dropna().values
    vals_191_non   = df.loc[(df["family"].eq("Tgut191A")) & (~df["PAR"].astype(bool)), v].dropna().values

    # Effect sizes via U-statistic
    if len(vals_716) and len(vals_191):
        u_716_191 = mannwhitneyu(vals_716, vals_191, alternative="two-sided").statistic
        delta_716_191 = _cliffs_delta_from_u(u_716_191, len(vals_716), len(vals_191))
    else:
        delta_716_191 = float("nan")

    if len(vals_716_cen) and len(vals_716_non):
        u_cen_non = mannwhitneyu(vals_716_cen, vals_716_non, alternative="two-sided").statistic
        delta_cen_non = _cliffs_delta_from_u(u_cen_non, len(vals_716_cen), len(vals_716_non))
    else:
        delta_cen_non = float("nan")

    if len(vals_191_par) and len(vals_191_non):
        u_par_non = mannwhitneyu(vals_191_par, vals_191_non, alternative="two-sided").statistic
        delta_par_non = _cliffs_delta_from_u(u_par_non, len(vals_191_par), len(vals_191_non))
    else:
        delta_par_non = float("nan")

    rows_716.append({
        "variable": v,
        "mean_Tgut716A": _mean(vals_716),
        "median_Tgut716A": _median(vals_716),
        "iqr_Tgut716A": _iqr(vals_716),
        "mean_Tgut716A_CEN": _mean(vals_716_cen),
        "median_Tgut716A_CEN": _median(vals_716_cen),
        "iqr_Tgut716A_CEN": _iqr(vals_716_cen),
        "mean_Tgut716A_nonCEN": _mean(vals_716_non),
        "median_Tgut716A_nonCEN": _median(vals_716_non),
        "iqr_Tgut716A_nonCEN": _iqr(vals_716_non),
        "BH_p_Tgut716A_vs_Tgut191A": adj_716_vs_191[i],
        "CliffsDelta_Tgut716A_vs_Tgut191A": delta_716_191,
        "BH_p_Tgut716A_CEN_vs_nonCEN": adj_716_cen_vs_non[i],
        "CliffsDelta_Tgut716A_CEN_vs_nonCEN": delta_cen_non,
    })
    rows_191.append({
        "variable": v,
        "mean_Tgut191A": _mean(vals_191),
        "median_Tgut191A": _median(vals_191),
        "iqr_Tgut191A": _iqr(vals_191),
        "mean_Tgut191A_PAR": _mean(vals_191_par),
        "median_Tgut191A_PAR": _median(vals_191_par),
        "iqr_Tgut191A_PAR": _iqr(vals_191_par),
        "mean_Tgut191A_nonPAR": _mean(vals_191_non),
        "median_Tgut191A_nonPAR": _median(vals_191_non),
        "iqr_Tgut191A_nonPAR": _iqr(vals_191_non),
        "BH_p_Tgut716A_vs_Tgut191A": adj_716_vs_191[i],
        "CliffsDelta_Tgut716A_vs_Tgut191A": delta_716_191,
        "BH_p_Tgut191A_PAR_vs_nonPAR": adj_191_par_vs_non[i],
        "CliffsDelta_Tgut191A_PAR_vs_nonPAR": delta_par_non,
    })

# Write TSVs
pd.DataFrame(rows_716).to_csv(os.path.join(plot_root, "Tgut716A_stats.tsv"), sep="\t", index=False)
pd.DataFrame(rows_191).to_csv(os.path.join(plot_root, "Tgut191A_stats.tsv"), sep="\t", index=False)

# ───────────────────────────────────────────── Plotting primitives

plt.style.use("default")
rng = np.random.default_rng(7)

def _draw_half_violin(ax, center_x, vals, color, width_violin=0.35):
    """Right half-violin with KDE."""
    if len(vals) > 1 and np.max(vals) > np.min(vals):
        kde = gaussian_kde(vals, bw_method=0.3)
        y = np.linspace(vals.min(), vals.max(), 300)
        dens = kde(y)
        if dens.max() > 0:
            scale = width_violin / dens.max()
            x_right = center_x + dens * scale
            ax.fill_betweenx(y, center_x, x_right, alpha=0.6, linewidth=0, color=color, zorder=1)
            ax.plot(x_right, y, linewidth=1.0, color=color, zorder=2)

def _draw_box(ax, center_x, vals, color, box_width=0.12):
    sns.boxplot(
        x=np.full_like(vals, center_x, dtype=float),
        y=vals,
        width=box_width,
        showcaps=True,
        boxprops={'facecolor': color, 'edgecolor': 'black', 'linewidth': 1.0},
        whiskerprops={'color': 'black'},
        capprops={'color': 'black'},
        medianprops={'color': 'black', 'linewidth': 1.0},
        showfliers=False,
        orient='v',
        ax=ax,
        zorder=3
    )

def raincloud_plot(ax, data, value_col, palette_map, boolean_grouped=False):
    """
    Draw rainclouds for one variable.

    boolean_grouped=False -> 4 groups: 716A CEN, 716A non-CEN, 191A PAR, 191A non-PAR (all hollow points)
    boolean_grouped=True  -> 2 groups: 716A, 191A (bool_flag points filled)
    """
    # thinner geometry for boolean panels
    jitter_sd = SCATTER_JITTER_SD_BOOL if boolean_grouped else SCATTER_JITTER_SD_GROUPED
    sizes = SCATTER_SIZE_BOOL if boolean_grouped else SCATTER_SIZE_GROUPED
    x_offset = SCATTER_X_OFFSET_BOOL if boolean_grouped else SCATTER_X_OFFSET_GROUPED
    violin_w = VIOLIN_WIDTH_BOOL if boolean_grouped else VIOLIN_WIDTH_GROUPED

    if not boolean_grouped:
        # Four groups
        groups = [
            ("Tgut716A", "CEN", True,  "Tgut716A CEN"),
            ("Tgut716A", "CEN", False, "Tgut716A non-CEN"),
            ("Tgut191A", "PAR", True,  "Tgut191A PAR"),
            ("Tgut191A", "PAR", False, "Tgut191A non-PAR"),
        ]
        pos = np.arange(len(groups), dtype=float)
        xticklabels = []

        for i, (fam, flag_col, flag_val, label) in enumerate(groups):
            sub = data[(data["family"].eq(fam)) & (data[flag_col].astype(bool) == flag_val)]
            vals = sub[value_col].dropna().values
            if len(vals) == 0:
                xticklabels.append(label.replace(" ", "\n"))
                continue
            deep = palette_map[fam]
            _draw_half_violin(ax, pos[i], vals, deep, width_violin=violin_w)
            _draw_box(ax, pos[i], vals, deep)
            # scatter: all hollow
            x_jit = rng.normal(loc=pos[i] - x_offset, scale=jitter_sd, size=len(vals))
            ras = RASTERIZE_SCATTER and (len(vals) >= RASTERIZE_MIN_POINTS)
            ax.scatter(
                x_jit, vals, s=sizes, linewidths=0.6, facecolors='none', edgecolors=deep,
                zorder=4, rasterized=ras
            )
            xticklabels.append(label.replace(" ", "\n"))

        ax.set_xticks(pos, xticklabels, fontsize=9)  # flat, two-line labels
        ax.set_xlim(-0.6, len(groups) - 0.4)

    else:
        # Two groups with boolean highlight
        fam_order = ["Tgut716A", "Tgut191A"]
        pos = np.arange(len(fam_order), dtype=float)
        for i, fam in enumerate(fam_order):
            sub = data[data["family"].eq(fam)]
            vals = sub[value_col].dropna().values
            if len(vals) == 0:
                continue
            deep = palette_map[fam]
            pastel = lighten_color(deep, 0.55)

            _draw_half_violin(ax, pos[i], vals, deep, width_violin=violin_w)
            _draw_box(ax, pos[i], vals, deep)

            # align boolean mask to non-NaN of value_col
            mask_idx = sub.index[sub[value_col].notna()]
            mask_true = sub.loc[mask_idx, "bool_flag"].astype(bool).values
            x_jit = rng.normal(loc=pos[i] - x_offset, scale=jitter_sd, size=len(vals))
            ras = RASTERIZE_SCATTER and (len(vals) >= RASTERIZE_MIN_POINTS)

            # hollow for FALSE, filled pastel for TRUE
            ax.scatter(
                x_jit[~mask_true], vals[~mask_true],
                s=sizes, linewidths=0.6, facecolors='none', edgecolors=deep, zorder=4, rasterized=ras
            )
            ax.scatter(
                x_jit[mask_true], vals[mask_true],
                s=sizes, linewidths=0.6, facecolors=pastel, edgecolors=deep, zorder=5, alpha=0.85, rasterized=ras
            )

        ax.set_xticks(pos, fam_order, fontsize=10)
        ax.set_xlim(-0.6, len(fam_order) - 0.4)

    # Aesthetics
    ax.set_ylabel(_ylabel_for(value_col))
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

# ───────────────────────────────────────────── Rainclouds

def plot_raincloud_variable(var, boolean_grouped=False):
    fig_size = FIGSIZE_BOOL if boolean_grouped else FIGSIZE_GROUPED
    fig, ax = plt.subplots(figsize=fig_size)
    fig.subplots_adjust(bottom=0.16)  # give room for two-line x-labels
    # choose out dir
    if var == "repeat_length":
        out_dir = plot_root
    else:
        comp = var.split("_", 1)[0]  # global / takki / chr
        out_dir = os.path.join(plot_root, comp)

    raincloud_plot(ax, df, var, palette, boolean_grouped=boolean_grouped)

    suffix = "bool" if boolean_grouped else "grouped"
    stem = os.path.join(out_dir, f"{var}_raincloud_{suffix}")
    savefig(fig, stem)
    plt.close(fig)

if RUN_RAINCLOUDS:
    for v in all_vars:
        if PLOT_GROUPED:
            plot_raincloud_variable(v, boolean_grouped=False)   # 4 groups: CEN/non & PAR/non
        if PLOT_BOOLEAN_HIGHLIGHT:
            plot_raincloud_variable(v, boolean_grouped=True)    # 2 groups with boolean highlight

print("All plots and stats written to:", plot_root)
