import matplotlib.pyplot as plt
import pandas as pd
import argparse
import numpy as np
import seaborn as sns
from scipy.stats import gaussian_kde
from matplotlib.lines import Line2D
import os
import matplotlib.ticker as mticker

# Set up argument parser
parser = argparse.ArgumentParser(description='Scatter/raincloud plots of "Tgut716A" and "Tgut191A" repeats from two CSV files.')
parser.add_argument('mat_csv', type=str, help='Path to the CSV file for "Tgut716A" repeats')
parser.add_argument('pat_csv', type=str, help='Path to the CSV file for "Tgut191A" repeats')

# Parse the arguments
args = parser.parse_args()

# Read the "Tgut716A" and "Tgut191A" CSV files
mat_repeats = pd.read_csv(args.mat_csv, header=None, names=['Chromosome', 'Value'])
pat_repeats = pd.read_csv(args.pat_csv, header=None, names=['Chromosome', 'Value'])

# Sort both "Tgut716A" and "Tgut191A" by chromosome name to match corresponding repeats
mat_repeats = mat_repeats.sort_values(by='Chromosome')
pat_repeats = pat_repeats.sort_values(by='Chromosome')

# Add identifiers used for coloring/faceting
mat_repeats['Repeat'] = 'Tgut716A'
pat_repeats['Repeat'] = 'Tgut191A'
# Expect labels like chrN_mat / chrN_pat
mat_repeats['Haplotype'] = mat_repeats['Chromosome'].astype(str).str.extract(r'_(mat|pat)$', expand=False)
pat_repeats['Haplotype'] = pat_repeats['Chromosome'].astype(str).str.extract(r'_(mat|pat)$', expand=False)
df = pd.concat([mat_repeats, pat_repeats], ignore_index=True)

# Raincloud helpers (half-violin + thin box + jitter)
def _draw_box(ax, center_x, vals, color, box_width=0.12):
    sns.boxplot(
        x=np.full(len(vals), center_x, dtype=float),
        y=np.asarray(vals, dtype=float),
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

def _jitter_points(ax, center_x, vals, edge_color, rng, x_offset=0.25, jitter_sd=0.045, size=15):
    x_jit = rng.normal(loc=center_x - x_offset, scale=jitter_sd, size=len(vals))
    ax.scatter(x_jit, vals, s=size, linewidths=0.6, facecolors='none', edgecolors=edge_color, zorder=4)

def _legend_proxies(palette):
    return [
        Line2D([0], [0], marker='o', linestyle='None', markerfacecolor=palette['Tgut716A'],
               markeredgecolor=palette['Tgut716A'], markeredgewidth=2.0, label='Tgut716A'),
        Line2D([0], [0], marker='o', linestyle='None', markerfacecolor=palette['Tgut191A'],
               markeredgecolor=palette['Tgut191A'], markeredgewidth=2.0, label='Tgut191A'),
    ]

outdir = 'plots'  # Output directory
def _save_all(fig, base):
    os.makedirs(outdir, exist_ok=True)
    fig.savefig(f'{outdir}/{base}.png', dpi=600, bbox_inches='tight', pad_inches=0.02)
    fig.savefig(f'{outdir}/{base}.svg',           bbox_inches='tight', pad_inches=0.02)
    fig.savefig(f'{outdir}/{base}.pdf',           bbox_inches='tight', pad_inches=0.02)

# Palette and geometry shared across plots
palette = {"Tgut716A": "#CC79A7", "Tgut191A": "#009E73"}
rng = np.random.default_rng(7)
y_ticks = [1e3, 1e4, 1e5, 1e6, 1e7]  # Only show these Y grid/ticks

# Two adjacent subpanels (Maternal, Paternal) sharing Y; X-axis shows ONLY these two groups
fig2A, (axL, axR) = plt.subplots(1, 2, figsize=(4, 6), sharey=True, layout='constrained')

# Minimal, centralized params for readability
BOX_WIDTH = 0.12
X_OFFSET   = 0.25
JITTER_SD  = 0.045
POINT_SIZE = 15
X_POS      = [0.0, 1.0]  # positions for 716A, 191A inside each panel
X_LIM      = (-0.55, 1.35)

# --- Left panel: Maternal ---
for pos, rep in zip(X_POS, ['Tgut716A', 'Tgut191A']):
    vals = df.loc[df['Haplotype'].eq('mat') & df['Repeat'].eq(rep), 'Value'].dropna().values
    if len(vals) == 0:
        continue
    _draw_box(axL, pos, vals, palette[rep], box_width=BOX_WIDTH)
    _jitter_points(axL, pos, vals, palette[rep], rng, x_offset=X_OFFSET, jitter_sd=JITTER_SD, size=POINT_SIZE)

# --- Right panel: Paternal ---
for pos, rep in zip(X_POS, ['Tgut716A', 'Tgut191A']):
    vals = df.loc[df['Haplotype'].eq('pat') & df['Repeat'].eq(rep), 'Value'].dropna().values
    if len(vals) == 0:
        continue
    _draw_box(axR, pos, vals, palette[rep], box_width=BOX_WIDTH)
    _jitter_points(axR, pos, vals, palette[rep], rng, x_offset=X_OFFSET, jitter_sd=JITTER_SD, size=POINT_SIZE)

# Axes styling
for ax in (axL, axR):
    ax.set_xlim(*X_LIM)
    ax.set_xticks([])
    ax.set_yscale('log')
    ax.set_yticks(y_ticks)
    ax.yaxis.set_minor_locator(mticker.NullLocator())
    ax.set_ylabel('Total repeat content (bp)')
    ax.grid(True, axis='y', which='major')
    ax.axhline(y=10000, color='black', linestyle='--', linewidth=1.0)

# Move “Maternal”/“Paternal” labels to the bottom
axL.set_xlabel('Maternal', labelpad=6)
axR.set_xlabel('Paternal', labelpad=6)

# Legend: simple and opaque white frame (no background lines bleeding through)
leg = axR.legend(handles=_legend_proxies(palette), loc='upper right', frameon=True)
leg.get_frame().set_facecolor('white')
# leg.get_frame().set_edgecolor('black')
leg.get_frame().set_alpha(1.0)

_save_all(fig2A, 'presence_by_hap')

# Show the plots (optional)
plt.show()
