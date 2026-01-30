#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress, t, gaussian_kde, mannwhitneyu
from itertools import combinations
from matplotlib.ticker import FormatStrFormatter
from matplotlib.lines import Line2D

def significance_label(p):
    if p <= 0.001:
        return '***'
    elif p <= 0.01:
        return '**'
    elif p <= 0.05:
        return '*'
    else:
        return 'n.s.'

def fdr_bh(pvals):
    """Benjamini–Hochberg FDR adjustment for a list of p-values."""
    pvals = np.asarray(pvals, dtype=float)
    m = len(pvals)
    order = np.argsort(pvals)
    ranked = np.empty(m, dtype=float)
    ranked[order] = np.arange(1, m+1)
    adj = pvals * m / ranked

    # enforce monotonicity
    adj_sorted = np.minimum.accumulate(adj[order][::-1])[::-1]
    adj_corrected = np.empty_like(adj_sorted)
    adj_corrected[order] = np.clip(adj_sorted, 0, 1)
    return adj_corrected

# Load data: paths, telomeres BED file and chromosome classification TSV
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

plot_dir = os.path.join(script_dir, "plots")
os.makedirs(plot_dir, exist_ok=True)

bed_fp = os.path.join(script_dir, "bTaeGut7_T2T", "bTaeGut7.fa.gz_terminal_telomeres.bed")
bed_cols = ['chr', 'start', 'end', 'length', 'label', 'fwdCounts', 'revCounts', 'canCounts', 'nonCanCounts', 'chrSize']
df_bed = pd.read_csv(bed_fp, sep="\t", header=None, names=bed_cols)

class_fp = os.path.join(script_dir, "chr_classification.tsv")
df_class = pd.read_csv(class_fp, sep="\t", header=0)
df_class = df_class.rename(columns={'Name': 'chr'})

# Merge the two dataframes
df = pd.merge(df_bed, df_class[['chr', 'Chr_type', 'Classification_size', 'Dot', 'TCHEST']], on='chr', how='left')

# Normalize TCHEST labels just in case of trailing punctuation/case
if 'TCHEST' in df.columns:
    df['TCHEST'] = df['TCHEST'].astype(str).str.strip().str.replace(r'\.$', '', regex=True).str.capitalize()

# derive haplotype & sizes
df['haplotype'] = df['chr'].str.rsplit('_', n=1).str[-1]
df['lengthKb']  = df['length'] / 1e3
df['chrMb']     = df['chrSize'] / 1e6

# outlier bounds (used for plotting and stats)
q1, q3   = df['lengthKb'].quantile([0.25, 0.75])
iqr      = q3 - q1
lb, ub   = q1 - 1.5 * iqr, q3 + 1.5 * iqr
df_clean = df[(df['lengthKb'] >= lb) & (df['lengthKb'] <= ub)]
df_outliers = df[(df['lengthKb'] < lb) | (df['lengthKb'] > ub)]

# shape mapping
shape_map = {'mat': 'o', 'pat': 's'}

# variables and titles
vars_to_plot = ["Chr_type", "Classification_size", "Dot", "TCHEST"]
titles = {
    "Chr_type": "Chromosome type",
    "Classification_size": "Classification size",
    "Dot": "Dot chromosomes",
    "TCHEST": "TCHEST type"
}

# TCHEST palette
tchest_colors = {
    "Strict":  "#D55E00",  # dot microchromosome
    "Absent":  "#0072B2",  # macrochromosome
    "Relaxed": "#F0E442",  # non-dot microchromosome
}
tchest_order = ["Strict", "Relaxed", "Absent"]

# use minimal seaborn style (no grid)
sns.set_style("white")

# Dimensions:
SCATTER_W, RAIN_W, HEIGHT = 4.0, 2.25, 4.0

for var in vars_to_plot:
    cats_found = df[var].dropna().unique().tolist()

    # enforce order for TCHEST; otherwise keep sorted
    if var == "TCHEST":
        cats = [c for c in tchest_order if c in cats_found]
    else:
        cats = sorted(cats_found)

    n = len(cats)
    if n == 0:
        continue

    # figure
    fig, (ax_scatter, ax_rain) = plt.subplots(
        1, 2,
        figsize=(SCATTER_W + RAIN_W, HEIGHT),
        sharey=True,
        gridspec_kw={'width_ratios': [SCATTER_W, RAIN_W], 'wspace': 0.1}
    )


    # color palettes
    if var == "TCHEST":
        solid_colors = [tchest_colors[c] for c in cats]
        pastel_colors = solid_colors  # use same hues with alpha for fill
    else:
        pastel_colors = sns.color_palette("pastel", n_colors=n)
        solid_colors  = sns.color_palette("deep",   n_colors=n)

    # --- Scatter panel ---
    for i, val in enumerate(cats):
        color = solid_colors[i]
        for hap, marker in shape_map.items():
            sub_in = df_clean[(df_clean[var] == val) & (df_clean['haplotype'] == hap)]
            ax_scatter.scatter(
                sub_in['chrMb'], sub_in['lengthKb'],
                marker=marker, facecolor=color, edgecolor=color,
                alpha=0.7, s=40
            )
            sub_out = df_outliers[(df_outliers[var] == val) & (df_outliers['haplotype'] == hap)]
            ax_scatter.scatter(
                sub_out['chrMb'], sub_out['lengthKb'],
                marker=marker, facecolor='none', edgecolor=color, s=40
            )

    ax_scatter.set_xscale('log')
    ax_scatter.get_xaxis().set_major_formatter(FormatStrFormatter('%g'))
    ax_scatter.set_xlabel('Chromosome Size (Mbp)')
    ax_scatter.set_ylabel('Telomere Length (Kbp)')
    ax_scatter.grid(False)
    ax_scatter.tick_params(direction='in', bottom=True, right=True)
    ax_scatter.set_ylim(0, 35)
    ax_scatter.set_yticks([0, 10, 20, 30])

    # regression + CI (overall)
    x = df_clean['chrMb']
    y = df_clean['lengthKb']
    logx = np.log10(x)
    slope, intercept, r_value, p_value, std_err = linregress(logx, y)

    x_line    = np.logspace(np.log10(x.min()), np.log10(x.max()), 100)
    logx_line = np.log10(x_line)
    y_line    = intercept + slope * logx_line

    n_d       = len(logx)
    dfree     = n_d - 2
    t_stat    = t.ppf(0.975, dfree)
    mean_logx = np.mean(logx)
    resid     = y - (intercept + slope * logx)
    sigma     = np.sqrt(np.sum(resid**2) / dfree)
    ssx       = np.sum((logx - mean_logx)**2)
    se_fit    = sigma * np.sqrt(1/n_d + (logx_line - mean_logx)**2 / ssx)
    y_low     = y_line - t_stat * se_fit
    y_up      = y_line + t_stat * se_fit

    ax_scatter.plot(x_line, y_line, linewidth=1, color='grey', zorder=2)
    ax_scatter.fill_between(x_line, y_low, y_up, color='lightgrey', alpha=0.4, zorder=1)

    # Regression p is not adjusted; keep "p", but set fontsize to match legend
    ann_fs = 10
    p_label = "p"  # change to "p-adj" only if you adjust this regression p-value elsewhere
    ax_scatter.text(
        0.05, 0.95,
        f"$R^2$ = {r_value**2:.2f}\n{p_label} = {p_value:.2g}",
        transform=ax_scatter.transAxes,
        verticalalignment='top',
        fontsize=ann_fs
    )

    # Shape-only legend in dark grey (circle=mat, square=pat), larger font
    legend_handles = [
        Line2D([0], [0], marker='o', linestyle='None', color='dimgray', label='Maternal', markersize=6),
        Line2D([0], [0], marker='s', linestyle='None', color='dimgray', label='Paternal', markersize=6),
    ]
    ax_scatter.legend(handles=legend_handles, fontsize=ann_fs, frameon=True)

    # --- Raincloud panel: half-violin + box ---
    positions = np.arange(n)
    for i, cat in enumerate(cats):
        sub = df_clean[df_clean[var] == cat]
        clean_vals = sub['lengthKb']

        # half-violin
        if len(clean_vals) > 1:
            kde     = gaussian_kde(clean_vals, bw_method=0.5)
            y_vals  = np.linspace(clean_vals.min(), clean_vals.max(), 200)
            dens    = kde(y_vals)
            scale   = 0.4 / dens.max()
            offsets = dens * scale

            ax_rain.fill_betweenx(
                y_vals, positions[i], positions[i] + offsets,
                color=pastel_colors[i], alpha=0.6, linewidth=0, zorder=1
            )
            ax_rain.plot(
                positions[i] + offsets, y_vals,
                color=pastel_colors[i], linewidth=1, zorder=2
            )

        # boxplot
        sns.boxplot(
            y='lengthKb',
            data=sub,
            orient='v',
            width=0.2,
            showcaps=True,
            boxprops={'facecolor': solid_colors[i], 'edgecolor': 'black'},
            whiskerprops={'color': 'black'},
            capprops={'color': 'black'},
            medianprops={'color': 'black'},
            showfliers=False,
            positions=[positions[i]],
            ax=ax_rain,
            zorder=3
        )

    # --- Statistics: Mann–Whitney U with BH-FDR within each variable ---
    group_vals = [df[df[var] == cat]['lengthKb'] for cat in cats]  # Full data
    pairs = list(combinations(range(n), 2))
    raw_ps = []
    U_stats = []

    for (i1, i2) in pairs:
        v1 = group_vals[i1]
        v2 = group_vals[i2]
        U, p = mannwhitneyu(v1, v2, alternative='two-sided')
        U_stats.append(U)
        raw_ps.append(p)

    if len(raw_ps) > 0:
        if len(raw_ps) == 1:
            padjs = raw_ps[:]  # no correction needed for a single comparison
        else:
            padjs = fdr_bh(raw_ps)
        p_adj_map = {pairs[k]: padjs[k] for k in range(len(pairs))}
    else:
        p_adj_map = {}

    # significance bars use adjusted p-values
    y_min, y_max = ax_rain.get_ylim()
    margin      = (y_max - y_min) * 0.1
    bar_y_start = y_max - margin

    if len(pairs) == 1:
        (i1, i2) = pairs[0]
        p_adj = p_adj_map[(i1, i2)]
        label = significance_label(p_adj)
        x1, x2 = positions[i1], positions[i2]
        ax_rain.plot([x1, x2], [bar_y_start, bar_y_start], color='black', lw=1, zorder=5)
        ax_rain.text((x1 + x2) / 2, bar_y_start + margin*0.2, label,
                     ha='center', va='bottom', fontsize=10, zorder=5)
        new_y_max = bar_y_start + margin
    else:
        for k, (i1, i2) in enumerate(pairs):
            p_adj = p_adj_map[(i1, i2)]
            label = significance_label(p_adj)
            y_bar = bar_y_start + k * margin
            ax_rain.plot([i1, i2], [y_bar, y_bar], color='black', lw=1, zorder=5)
            ax_rain.text((i1 + i2) / 2, y_bar + margin*0.2, label,
                         ha='center', va='bottom', fontsize=10, zorder=5)
        new_y_max = bar_y_start + (len(pairs)-1) * margin + margin

    ax_rain.set_ylim(y_min, new_y_max)
    ax_rain.set_xticks(positions)
    ax_rain.set_xticklabels(cats, fontsize=10, rotation=45, ha='right')
    ax_rain.set_xlabel('')
    sns.despine(ax=ax_rain, left=False, top=False, right=False, bottom=False)
    ax_rain.yaxis.grid(True, linestyle='--', linewidth=0.5, alpha=0.5)
    ax_rain.xaxis.grid(False)

    # Title should belong to the raincloud (top-to-bottom: Title > bars > raincloud)
    ax_rain.set_title(titles[var], fontsize=12, pad=6)

    # save figure
    for ext, dpi in [('pdf', None), ('svg', None), ('png', 600)]:
        fp = os.path.join(plot_dir, f"chr_telo_sizes_{var}.{ext}")
        fig.savefig(fp, dpi=dpi, bbox_inches='tight')
    plt.close(fig)

    # --- Console reporting ---
    print(f"\n[{var}] group means (Kb):")
    for cat, vals in zip(cats, group_vals):
        print(f"  {cat}: mean={vals.mean():.2f} Kb, n={len(vals)}")

    if len(pairs) > 0:
        print(f"[{var}] pairwise Mann-Whitney U tests with BH-FDR within {var}:")
        for (U, p_raw, (i1, i2)) in zip(U_stats, raw_ps, pairs):
            p_adj = p_adj_map[(i1, i2)]
            print(f"  {cats[i1]} vs {cats[i2]}: U={U:.1f}, p={p_raw:.3e}, p-adj={p_adj:.3e}")

print("\nCombined scatter + raincloud plots saved to", plot_dir)
