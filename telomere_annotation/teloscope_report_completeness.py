# teloscope_report_completeness.py
import matplotlib
matplotlib.use('Agg') 
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import MultipleLocator
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch

# ──────────────────────────────────────────────────────────────────────────────
# Function to annotate generic stacked horizontal bars with their widths (integers)
def annotate_bars(ax):
    for patch in ax.patches:
        w = patch.get_width()
        if w > 0:
            ax.text(
                patch.get_x() + w/2,
                patch.get_y() + patch.get_height()/2,
                f"{int(round(w))}",
                ha='center', va='center',
                fontsize=8, color='black'
            )

# Annotate stacked horizontal bars with custom (absolute) labels for each patch
def annotate_bars_with_values(ax, values):
    """
    values: list/iterable of numbers (absolute counts) in the same order as ax.patches
    """
    for patch, val in zip(ax.patches, values):
        w = patch.get_width()
        if w > 0:
            ax.text(
                patch.get_x() + w/2,
                patch.get_y() + patch.get_height()/2,
                f"{int(round(val))}",
                ha='center', va='center',
                fontsize=8, color='black'
            )

# ──────────────────────────────────────────────────────────────────────────────
# Paths
# Use compiled CSV from parse_telo_reports.py
reports_file = "teloscope_reports_compiled.csv"

# Output dir and formats
plot_dir = "./plots"
os.makedirs(plot_dir, exist_ok=True)

# ──────────────────────────────────────────────────────────────────────────────
# Global style
plt.rcParams.update({
    'font.family': 'Arial',
    'font.weight': 'regular',
    'font.size': 9,
    'axes.linewidth': 0.5,
    'axes.titlesize': 9,
    'axes.labelsize': 9,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'legend.frameon': False,
    'figure.dpi': 600,
    'savefig.dpi': 600
})

# ──────────────────────────────────────────────────────────────────────────────
# Load data
# CSV → rows per accession with required columns
df = pd.read_csv(reports_file).fillna(0)
if 'accession' in df.columns:
    df.set_index('accession', inplace=True)

# Optional manual list of accessions that INCLUDE mitochondrial DNA.
# If provided, we subtract 1 from: total_paths, zero_telomeres, no_telomeres
# for those accessions. If empty, auto-correct any row with no_telomeres>0.
mito_accessions = []  # e.g., ["GCA_048771995.1", "GCA_051427915.1"]

def _safe_sub_one(s: pd.Series) -> pd.Series:
    s = s.astype(int)
    return (s - 1).clip(lower=0)

if mito_accessions:
    mask = df.index.isin(mito_accessions)
    for col in ['total_paths', 'zero_telomeres', 'no_telomeres']:
        if col in df.columns:
            df.loc[mask, col] = _safe_sub_one(df.loc[mask, col])
else:
    # Auto rule: assembled CM → only mitochondrion lacks telomeres/gaps.
    # Subtract 1 where there is at least one "no_telomeres".
    if 'no_telomeres' in df.columns:
        mask = df['no_telomeres'].astype(int) > 0
        for col in ['total_paths', 'zero_telomeres', 'no_telomeres']:
            if col in df.columns:
                df.loc[mask, col] = _safe_sub_one(df.loc[mask, col])

# ──────────────────────────────────────────────────────────────────────────────
# Display order + year tags for x-axis
order = [
    "GCA_000151805.2",
    "GCA_008822115.3",
    "GCA_003957565.4",
    "GCA_051427915.1",
    "GCA_048771995.1",
]
years = {
    "GCA_000151805.2": "2013",
    "GCA_008822115.3": "2020",
    "GCA_003957565.4": "2021",
    "GCA_051427915.1": "2022",
    "GCA_048771995.1": "2025",
}
# Keep only those present and reindex in desired order
order = [acc for acc in order if acc in df.index]
df = df.reindex(order)

# Pretty x tick labels (accession + year). Remove xlabel text later.
xticklabels = [f"{acc}\n({years.get(acc,'')})" for acc in df.index]

# ──────────────────────────────────────────────────────────────────────────────
# Series per accession (ensure integer dtype where applicable)
paths                 = df['total_paths'].astype(int)
total_telomeres       = df['total_telomeres'].astype(int)
two_telomeres         = df['two_telomeres'].astype(int)
one_telomere          = df['one_telomere'].astype(int)
zero_telomeres        = df['zero_telomeres'].astype(int)

# Completeness classes (full palette)
t2t                   = df['t2t'].astype(int)
gapped_t2t            = df['gapped_t2t'].astype(int)
missassembled         = df['missassembled'].astype(int)
gapped_missassembled  = df['gapped_missassembled'].astype(int)
incomplete            = df['incomplete'].astype(int)
gapped_incomplete     = df['gapped_incomplete'].astype(int)
discordant            = df.get('discordant', 0).astype(int)
gapped_discordant     = df.get('gapped_discordant', 0).astype(int)
no_telomeres_cls      = df.get('no_telomeres', 0).astype(int)
gapped_no_telomeres   = df.get('gapped_no_telomeres', 0).astype(int)

# ──────────────────────────────────────────────────────────────────────────────
# Color map (A uses Found + Missing; C uses full palette; "Missing" = #F0F0F0)
colors = {
    'Found':                   '#4DCCBD',
    'Missing':                 '#F0F0F0',

    't2t':                     '#1A9641',
    'gapped_t2t':              '#9CCF60',
    'incomplete':              '#FFC754',
    'gapped_incomplete':       '#FFE885',
    'missassembled':           '#D6594C',
    'gapped_missassembled':    '#F58B6D',
    'discordant':              '#8278F4',
    'gapped_discordant':       '#B395EB',
    'no_telomeres':            '#C8C8C8',
    'gapped_no_telomeres':     '#F0F0F0',
}

# ──────────────────────────────────────────────────────────────────────────────
# Helper to place legend at bottom with 2 rows
def bottom_legend(ax, handles, title=None, ncol=None):
    if ncol is None:
        ncol = (len(handles) + 1) // 2  # 2 rows
    ax.legend(handles=handles, title=title, frameon=False,
              loc='upper center', bbox_to_anchor=(0.5, -0.2),
              ncol=ncol, borderaxespad=0.)

# ──────────────────────────────────────────────────────────────────────────────
# MODE 1 — Absolute counts vs expected from assembled chromosomes
#   • Plot A: Found telomeres vs Expected–Found (expected = 2×paths)
#   • Plot C: Absolute category counts (no "Missing" category)
abs_fig = plt.figure(figsize=(8, 3))
abs_fig.subplots_adjust(bottom=0.05)
abs_gs  = gridspec.GridSpec(1, 2, wspace=0.15, width_ratios=[1, 1])
axA_abs = plt.subplot(abs_gs[0, 0])
axC_abs = plt.subplot(abs_gs[0, 1])

axA_abs.xaxis.set_major_locator(MultipleLocator(20))
axC_abs.xaxis.set_major_locator(MultipleLocator(10))

baseline_tel = (2 * paths).astype(int)
missing_tel_abs = (baseline_tel - total_telomeres).clip(lower=0).astype(int)

dataA_abs = pd.DataFrame({
    'Found':   total_telomeres,
    'Missing': missing_tel_abs
}, index=df.index)

dataA_abs.plot.barh(stacked=True, ax=axA_abs,
                    color=[colors['Found'], colors['Missing']],
                    width=0.7, edgecolor='white', linewidth=0.5)
axA_abs.set_title('Assembly telomeres')
axA_abs.set_xlabel('Total telomeres')
axA_abs.set_yticklabels(xticklabels)
axA_abs.set_ylabel('')
axA_abs.invert_yaxis()
bottom_legend(axA_abs, [
    Patch(facecolor=colors['Found'],   label='Found'),
    Patch(facecolor=colors['Missing'], label='Missing'),
], ncol=2)
# annotate with absolute numbers per segment
abs_vals_A = list(total_telomeres.values) + list(missing_tel_abs.values)
annotate_bars_with_values(axA_abs, abs_vals_A)

# Panel C (absolute)
dataC_abs = pd.DataFrame({
    't2t':                     t2t,
    'gapped_t2t':              gapped_t2t,
    'incomplete':              incomplete,
    'gapped_incomplete':       gapped_incomplete,
    'missassembled':           missassembled,
    'gapped_missassembled':    gapped_missassembled,
    'discordant':              discordant,
    'gapped_discordant':       gapped_discordant,
    'no_telomeres':            no_telomeres_cls,
    'gapped_no_telomeres':     gapped_no_telomeres,
}, index=df.index)

dataC_abs.plot.barh(stacked=True, ax=axC_abs,
                    color=[colors[k] for k in dataC_abs.columns],
                    width=0.7, edgecolor='white', linewidth=0.5)
axC_abs.set_title('Chrs. by telomere completeness')
axC_abs.set_xlabel('Assembled chromosomes')
axC_abs.set_yticklabels([])  # shared y-axis labels with panel A
axC_abs.set_ylabel('')
axC_abs.invert_yaxis()
bottom_legend(axC_abs, [Patch(facecolor=colors[k], label=k.replace('_',' '))
                        for k in dataC_abs.columns if dataC_abs[k].sum() > 0], ncol=2)
annotate_bars(axC_abs)

# Save MODE 1
for ext, dpi in [('pdf', None), ('svg', None), ('png', 600)]:
    save_kwargs = {'format': ext, 'bbox_inches': 'tight'}
    if dpi is not None:
        save_kwargs['dpi'] = dpi
    abs_fig.savefig(os.path.join(plot_dir, f"telomere_completeness_absolute.{ext}"), **save_kwargs)

# ──────────────────────────────────────────────────────────────────────────────
# MODE 2 — Percent of expected from assembled chromosomes
#   • Plot A: Observed% (= found / (2*paths) × 100), remainder = Expected − observed
#   • Plot C: Each category / paths × 100
pct_fig = plt.figure(figsize=(8, 3))
pct_fig.subplots_adjust(bottom=0.05)
pct_gs  = gridspec.GridSpec(1, 2, wspace=0.15, width_ratios=[1, 1])
axA_pct = plt.subplot(pct_gs[0, 0])
axC_pct = plt.subplot(pct_gs[0, 1])
# MODE 2 (percent)
for ax in (axA_pct, axC_pct):
    ax.xaxis.set_major_locator(MultipleLocator(20))


# A (percent)
with np.errstate(divide='ignore', invalid='ignore'):
    baseline_tel_nonzero = baseline_tel.replace(0, np.nan)
    found_pct = (total_telomeres / baseline_tel_nonzero * 100).fillna(0)
found_pct = found_pct.astype(float)
missing_pct = (100 - found_pct).clip(lower=0)

dataA_pct = pd.DataFrame({
    'Found':   found_pct,
    'Missing': missing_pct
}, index=df.index)

dataA_pct.plot.barh(stacked=True, ax=axA_pct,
                    color=[colors['Found'], colors['Missing']],
                    width=0.7, edgecolor='white', linewidth=0.5)
axA_pct.set_title('Assembly telomeres')
axA_pct.set_xlabel('Total telomeres %')
axA_pct.set_yticklabels(xticklabels)
axA_pct.set_ylabel('')
axA_pct.invert_yaxis()
# show absolute values in labels
abs_vals_Apct = list(total_telomeres.values) + list(missing_tel_abs.values)
annotate_bars_with_values(axA_pct, abs_vals_Apct)
bottom_legend(axA_pct, [
    Patch(facecolor=colors['Found'],   label='Found'),
    Patch(facecolor=colors['Missing'], label='Missing'),
], ncol=2)

# C (percent)
with np.errstate(divide='ignore', invalid='ignore'):
    denom = paths.replace(0, np.nan).astype(float)
    dataC_pct = pd.DataFrame({
        't2t':                   (t2t / denom * 100),
        'gapped_t2t':            (gapped_t2t / denom * 100),
        'incomplete':            (incomplete / denom * 100),
        'gapped_incomplete':     (gapped_incomplete / denom * 100),
        'missassembled':         (missassembled / denom * 100),
        'gapped_missassembled':  (gapped_missassembled / denom * 100),
        'discordant':            (discordant / denom * 100),
        'gapped_discordant':     (gapped_discordant / denom * 100),
        'no_telomeres':          (no_telomeres_cls / denom * 100),
        'gapped_no_telomeres':   (gapped_no_telomeres / denom * 100),
    }, index=df.index).fillna(0.0)

dataC_pct.plot.barh(stacked=True, ax=axC_pct,
                    color=[colors[k] for k in dataC_pct.columns],
                    width=0.7, edgecolor='white', linewidth=0.5)
axC_pct.set_title('Chrs. by telomere completeness')
axC_pct.set_xlabel('Assembled chromosomes %')
axC_pct.set_yticklabels([])  # shared y-axis labels with panel A
axC_pct.set_ylabel('')
axC_pct.invert_yaxis()
# annotate with ABSOLUTE counts, not percents
abs_vals_Cpct = []
for col in ['t2t','gapped_t2t','incomplete','gapped_incomplete',
            'missassembled','gapped_missassembled','discordant','gapped_discordant',
            'no_telomeres','gapped_no_telomeres']:
    abs_vals_Cpct.extend(list(df[col].astype(int).values))
annotate_bars_with_values(axC_pct, abs_vals_Cpct)
bottom_legend(axC_pct, [Patch(facecolor=colors[k], label=k.replace('_',' '))
                        for k in dataC_pct.columns if dataC_pct[k].sum() > 0], ncol=2)

# Save MODE 2
for ext, dpi in [('pdf', None), ('svg', None), ('png', 600)]:
    save_kwargs = {'format': ext, 'bbox_inches': 'tight'}
    if dpi is not None:
        save_kwargs['dpi'] = dpi
    pct_fig.savefig(os.path.join(plot_dir, f"telomere_completeness_percent.{ext}"), **save_kwargs)

print("Saved figures to", os.path.abspath(plot_dir))
