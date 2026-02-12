#!/usr/bin/env python3
import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# Load data: paths and telomeres BED files
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

plot_dir = os.path.join(script_dir, "plots", "GGswu")
os.makedirs(plot_dir, exist_ok=True)

# Column names for BED files
bed_cols = ['chr', 'start', 'end', 'length', 'label', 'fwdCounts', 'revCounts', 'canCounts', 'nonCanCounts', 'chrSize']

# Load terminal telomeres
terminal_fp = os.path.join(script_dir, "GCA_024206055.2", "GCA_024206055.2_GGswu_genomic.fna_terminal_telomeres.bed")
df_terminal = pd.read_csv(terminal_fp, sep="\t", header=None, names=bed_cols)

# Load interstitial telomeres
interstitial_fp = os.path.join(script_dir, "GCA_024206055.2", "GCA_024206055.2_GGswu_genomic.fna_interstitial_telomeres.bed")
df_interstitial = pd.read_csv(interstitial_fp, sep="\t", header=None, names=bed_cols)

# Exclude mitochondrial sequence (CP115610.1)
MITO_ACC = "CP115610.1"
df_terminal = df_terminal[df_terminal['chr'] != MITO_ACC].copy()
df_interstitial = df_interstitial[df_interstitial['chr'] != MITO_ACC].copy()

# Derive common columns for both dataframes
for df in [df_terminal, df_interstitial]:
    df['lengthKb'] = df['length'] / 1e3
    df['canonicalProp'] = (df['canCounts'] * 6 * 100) / df['length']
    df['distanceToEnd'] = df[['start', 'end', 'chrSize']].apply(lambda row: min(row['start'], row['chrSize'] - row['end']), axis=1)

# Extract accession number for labeling (e.g., "CP100555.1" -> "555")
for df in [df_terminal, df_interstitial]:
    df['shortLabel'] = df['chr'].str.extract(r'CP\d*?(\d{3}\.\d)$')

# Additional columns for interstitial
df_interstitial['fwdProp'] = df_interstitial['fwdCounts'] * 100 / (df_interstitial['fwdCounts'] + df_interstitial['revCounts'])
df_interstitial['repeatCoverage'] = ((df_interstitial['fwdCounts'] + df_interstitial['revCounts']) * 6 * 100) / df_interstitial['length']

# Default scatter color
POINT_COLOR = '#0072B2'

# Label category palette
label_colors = {
    'p': '#56A3FF',   # light blue - CCCTAA dominant
    'q': '#FCCF14',   # yellow - TTAGGG dominant
    'u': '#417F32',   # green - fusion
}
label_order = ['p', 'q', 'u']

# Use minimal seaborn style
sns.set_style("white")

# Common plot parameters
POINT_SIZE = 25
POINT_ALPHA = 0.5
TICK_FONTSIZE = 8
LABEL_FONTSIZE = 9
LEGEND_FONTSIZE = 7
GRID_ALPHA = 0.5
GRID_LINEWIDTH = 0.5

# ITS filtering criteria
ITS_CANONICAL_THRESHOLD = 50  # Canonical proportion >= 50%
ITS_COVERAGE_THRESHOLD = 0   # Repeat coverage >= 0%

# Create filtered interstitial dataframe
df_its_filtered = df_interstitial[
    (df_interstitial['canonicalProp'] >= ITS_CANONICAL_THRESHOLD) &
    (df_interstitial['repeatCoverage'] >= ITS_COVERAGE_THRESHOLD)
].copy()

print(f"Total ITS: {len(df_interstitial)}")
print(f"Filtered ITS (canonical >= {ITS_CANONICAL_THRESHOLD}%, coverage >= {ITS_COVERAGE_THRESHOLD}%): {len(df_its_filtered)}")

def setup_axes(ax):
    """Common axis setup: ticks, spines, gridlines."""
    ax.tick_params(direction='in', which='both', labelsize=TICK_FONTSIZE)
    for spine in ax.spines.values():
        spine.set_visible(True)

def save_plot(fig, name):
    """Save figure in PDF, SVG, and PNG formats."""
    for ext, dpi in [('pdf', None), ('svg', None), ('png', 600)]:
        fp = os.path.join(plot_dir, f"{name}.{ext}")
        fig.savefig(fp, dpi=dpi, bbox_inches='tight')
    plt.close(fig)

# =============================================================================
# TERMINAL TELOMERE PLOTS
# =============================================================================

# 1.1 Terminal canonical proportion
# X: Telomere length (kbp), Y: Canonical proportion (%)
fig, ax = plt.subplots(figsize=(3.5, 3.5))

ax.scatter(df_terminal['lengthKb'], df_terminal['canonicalProp'],
           c=POINT_COLOR, alpha=POINT_ALPHA, s=POINT_SIZE)

# Add median lines
median_x = df_terminal['lengthKb'].median()
median_y = df_terminal['canonicalProp'].median()
ax.axvline(x=median_x, color='grey', linestyle='--', linewidth=0.5)
ax.axhline(y=median_y, color='grey', linestyle='--', linewidth=0.5)

ax.set_xlabel('Telomere length (kbp)', fontsize=LABEL_FONTSIZE)
ax.set_ylabel('Canonical proportion (%)', fontsize=LABEL_FONTSIZE)
setup_axes(ax)

save_plot(fig, 'terminal_canonical_proportion')
print("Saved: terminal_canonical_proportion")


# 1.2 Terminal chr positioning
# X: Distance to end (bp) in log scale, Y: Telomere length (kbp)
fig, ax = plt.subplots(figsize=(3.5, 3.5))

ax.scatter(df_terminal['distanceToEnd'], df_terminal['lengthKb'],
           c=POINT_COLOR, alpha=POINT_ALPHA, s=POINT_SIZE)

# Add 100bp dashed line
ax.axvline(x=100, color='grey', linestyle='--', linewidth=0.5)

# Count datapoints <= 100bp and > 100bp
n_le_100 = (df_terminal['distanceToEnd'] <= 100).sum()
n_gt_100 = (df_terminal['distanceToEnd'] > 100).sum()

# Add count labels at the bottom
ylim = ax.get_ylim()
ax.text(30, ylim[0] + 0.5, f'n={n_le_100}', ha='center', va='bottom', fontsize=LEGEND_FONTSIZE)
ax.text(300, ylim[0] + 0.5, f'n={n_gt_100}', ha='center', va='bottom', fontsize=LEGEND_FONTSIZE)

ax.set_xscale('log')
ax.set_xlabel('Distance to end (bp)', fontsize=LABEL_FONTSIZE)
ax.set_ylabel('Telomere length (kbp)', fontsize=LABEL_FONTSIZE)
setup_axes(ax)

save_plot(fig, 'terminal_chr_positioning')
print("Saved: terminal_chr_positioning")

# =============================================================================
# INTERSTITIAL TELOMERE PLOTS
# =============================================================================

# 2.1 Interstitial canonical proportion (filtered)
def plot_interstitial_canonical_proportion(df, suffix=''):
    """Plot interstitial canonical proportion."""
    fig, ax = plt.subplots(figsize=(3.5, 3.5))

    ax.scatter(df['lengthKb'], df['canonicalProp'],
               c=POINT_COLOR, alpha=POINT_ALPHA, s=POINT_SIZE)

    ax.axvline(x=1, color='grey', linestyle='--', linewidth=0.5)

    ax.set_xlabel('Telomere length (kbp)', fontsize=LABEL_FONTSIZE)
    ax.set_ylabel('Canonical proportion (%)', fontsize=LABEL_FONTSIZE)
    setup_axes(ax)

    name = f'interstitial_canonical_proportion{suffix}'
    save_plot(fig, name)
    print(f"Saved: {name}")

plot_interstitial_canonical_proportion(df_its_filtered, '_filtered')


# 2.2 Interstitial chr positioning by label (filtered)
def plot_interstitial_chr_positioning_label(df, suffix=''):
    """Plot interstitial chr positioning by label category (p/q/u)."""
    fig, ax = plt.subplots(figsize=(3.5, 3.5))

    for label_cat in label_order:
        sub = df[df['label'] == label_cat]
        if len(sub) > 0:
            ax.scatter(sub['distanceToEnd'], sub['lengthKb'],
                    c=label_colors[label_cat], label=label_cat, alpha=POINT_ALPHA, s=POINT_SIZE)

    ax.axhline(y=1, color='grey', linestyle='--', linewidth=0.5)

    ax.set_xscale('log')
    ax.set_xlabel('Distance to end (bp)', fontsize=LABEL_FONTSIZE)
    ax.set_ylabel('Telomere length (kbp)', fontsize=LABEL_FONTSIZE)
    ax.legend(title='Label type', fontsize=LEGEND_FONTSIZE, title_fontsize=LEGEND_FONTSIZE, frameon=True, loc='upper right')
    setup_axes(ax)

    name = f'interstitial_chr_positioning_label{suffix}'
    save_plot(fig, name)
    print(f"Saved: {name}")

plot_interstitial_chr_positioning_label(df_its_filtered, '_filtered')


# 2.3 Interstitial chr positioning by orientation (filtered)
def plot_interstitial_chr_positioning_orientation(df, suffix=''):
    """Plot interstitial chr positioning colored by forward proportion."""
    fig, ax = plt.subplots(figsize=(3.5, 3))

    scatter = ax.scatter(df['distanceToEnd'], df['lengthKb'],
                         c=df['fwdProp'], cmap='viridis',
                         vmin=0, vmax=100, alpha=POINT_ALPHA, s=POINT_SIZE)

    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Forward proportion (%)', fontsize=LABEL_FONTSIZE)
    cbar.ax.tick_params(labelsize=TICK_FONTSIZE)

    ax.axhline(y=1, color='grey', linestyle='--', linewidth=0.5)
    ax.set_xscale('log')
    ax.set_xlabel('Distance to end (bp)', fontsize=LABEL_FONTSIZE)
    ax.set_ylabel('Telomere length (kbp)', fontsize=LABEL_FONTSIZE)
    setup_axes(ax)

    name = f'interstitial_chr_positioning_orientation{suffix}'
    save_plot(fig, name)
    print(f"Saved: {name}")

plot_interstitial_chr_positioning_orientation(df_its_filtered, '_filtered')

# 2.4 Interstitial repeat coverage (filtered)
def plot_interstitial_repeat_coverage(df, suffix=''):
    """Plot interstitial repeat coverage with hollow points for >100%."""
    fig, ax = plt.subplots(figsize=(3.5, 3.5))

    # Split into normal (<= 100%) and outliers (> 100%)
    normal = df[df['repeatCoverage'] <= 100]
    outliers = df[df['repeatCoverage'] > 100]

    # Plot normal points (filled)
    if len(normal) > 0:
        ax.scatter(normal['lengthKb'], normal['repeatCoverage'],
                   c=POINT_COLOR, alpha=POINT_ALPHA, s=POINT_SIZE)

    # Plot outliers collapsed to 100% (hollow)
    if len(outliers) > 0:
        ax.scatter(outliers['lengthKb'], [100] * len(outliers),
                   facecolors='none', edgecolors=POINT_COLOR,
                   alpha=POINT_ALPHA, s=POINT_SIZE, label='>100%')
        ax.legend(fontsize=LEGEND_FONTSIZE, frameon=True, loc='lower right')

    ax.set_xlabel('Telomere length (kbp)', fontsize=LABEL_FONTSIZE)
    ax.set_ylabel('Repeat coverage (%)', fontsize=LABEL_FONTSIZE)
    setup_axes(ax)

    name = f'interstitial_repeat_coverage{suffix}'
    save_plot(fig, name)
    print(f"Saved: {name}")

plot_interstitial_repeat_coverage(df_its_filtered, '_filtered')

# 2.5 Interstitial repeat coverage vs canonical proportion (filtered)
def plot_interstitial_coverage_vs_canonical(df, suffix=''):
    """Plot interstitial repeat coverage vs canonical proportion."""
    fig, ax = plt.subplots(figsize=(3.5, 3.5))

    # Split into normal (<= 100%) and outliers (> 100%)
    normal = df[df['repeatCoverage'] <= 100]
    outliers = df[df['repeatCoverage'] > 100]

    # Plot normal points (filled)
    if len(normal) > 0:
        ax.scatter(normal['canonicalProp'], normal['repeatCoverage'],
                   c=POINT_COLOR, alpha=POINT_ALPHA, s=POINT_SIZE)

    # Plot outliers collapsed to 100% (hollow)
    if len(outliers) > 0:
        ax.scatter(outliers['canonicalProp'], [100] * len(outliers),
                   facecolors='none', edgecolors=POINT_COLOR,
                   alpha=POINT_ALPHA, s=POINT_SIZE, label='>100%')
        ax.legend(fontsize=LEGEND_FONTSIZE, frameon=True, loc='lower right')

    # Add median lines
    median_x = df['canonicalProp'].median()
    median_y = df['repeatCoverage'].median()
    ax.axvline(x=median_x, color='grey', linestyle='--', linewidth=0.5)
    ax.axhline(y=median_y, color='grey', linestyle='--', linewidth=0.5)

    ax.set_xlabel('Canonical proportion (%)', fontsize=LABEL_FONTSIZE)
    ax.set_ylabel('Repeat coverage (%)', fontsize=LABEL_FONTSIZE)
    setup_axes(ax)

    name = f'interstitial_coverage_vs_canonical{suffix}'
    save_plot(fig, name)
    print(f"Saved: {name}")

plot_interstitial_coverage_vs_canonical(df_its_filtered, '_filtered')

print(f"\nAll plots saved to {plot_dir}")
