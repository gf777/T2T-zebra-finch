#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generates supplementary figure:
  - Figure A: Chromosome vs. Telomere Length (kbp), sorted by mean
  - Figure B: Chromosome vs. Canonical repeats (%), sorted by mean

Input
-----
  - bTaeGut7_T2T/bTaeGut7.fa.gz_terminal_telomeres.bed
  - chr_classification.tsv

Output
------
Figures (PDF, SVG, PNG at 600 DPI) under {script_dir}/plots/
  - supfig16_1A_sorted_mean.pdf / .svg / .png
  - supfig16_1B_sorted_mean.pdf / .svg / .png
"""
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ----------------------------- Config & Style ------------------------------ #

# Science-style minimal aesthetics
mpl.rcParams.update({
    "figure.dpi": 300,
    "savefig.dpi": 600,
    "pdf.fonttype": 42,     # embed fonts properly
    "ps.fonttype": 42,
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial"],
    "font.weight": "normal",
    "axes.linewidth": 1.0,
    "axes.labelsize": 9,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "xtick.major.size": 3.2,
    "ytick.major.size": 3.2,
    "xtick.major.width": 0.9,
    "ytick.major.width": 0.9,
    "legend.frameon": True,
})

# Exact palette for TCHEST classes
TCHEST_COLORS = {
    "Strict":  "#D55E00",  # dot microchromosome
    "Absent":  "#0072B2",  # macrochromosome
    "Relaxed": "#F0E442",  # non-dot microchromosome
}

TCHEST_PRIORITY = ["Strict", "Relaxed", "Absent"]


# ------------------------------- Utilities -------------------------------- #

def pick_tchest_class(series: pd.Series) -> str:
    """Collapse per-row TCHEST into one class. Priority: Strict > Relaxed > Absent."""
    vals = set(series.dropna().astype(str))
    for c in TCHEST_PRIORITY:
        if c in vals:
            return c
    return "Absent"


# ------------------------------- Data Loading ------------------------------ #

def load_data(script_dir: str) -> pd.DataFrame:
    """Load terminal telomeres BED and chromosome classification, merge them."""
    
    # Column names for BED files
    bed_cols = ['chr', 'start', 'end', 'length', 'label', 'fwdCounts', 'revCounts', 
                'canCounts', 'nonCanCounts', 'chrSize']

    # Load terminal telomeres
    terminal_fp = os.path.join(script_dir, "bTaeGut7_T2T", "bTaeGut7.fa.gz_terminal_telomeres.bed")
    if not os.path.exists(terminal_fp):
        raise FileNotFoundError(f"Terminal telomeres file not found: {terminal_fp}")
    df_terminal = pd.read_csv(terminal_fp, sep="\t", header=None, names=bed_cols)
    print(f"Loaded {len(df_terminal)} terminal telomere records from {terminal_fp}")

    # Load chromosome classification
    class_fp = os.path.join(script_dir, "chr_classification.tsv")
    if not os.path.exists(class_fp):
        raise FileNotFoundError(f"Chromosome classification file not found: {class_fp}")
    df_class = pd.read_csv(class_fp, sep="\t", header=0)
    df_class = df_class.rename(columns={'Name': 'chr'})
    print(f"Loaded {len(df_class)} chromosome classification records from {class_fp}")

    # Merge classification data
    df_terminal = pd.merge(
        df_terminal, 
        df_class[['chr', 'Chr_type', 'Classification_size', 'Dot', 'TCHEST']], 
        on='chr', 
        how='left'
    )

    # Normalize TCHEST labels
    if 'TCHEST' in df_terminal.columns:
        df_terminal['TCHEST'] = (df_terminal['TCHEST']
                                 .astype(str)
                                 .str.strip()
                                 .str.replace(r'\.$', '', regex=True)
                                 .str.capitalize())

    # Derive common columns
    df_terminal['lengthKb'] = df_terminal['length'] / 1e3
    df_terminal['canonicalProp'] = (df_terminal['canCounts'] * 6 * 100) / df_terminal['length']
    df_terminal['distanceToEnd'] = df_terminal[['start', 'end', 'chrSize']].apply(
        lambda row: min(row['start'], row['chrSize'] - row['end']), axis=1
    )

    # Extract chromosome number and haplotype for labeling
    df_terminal['chrNum'] = df_terminal['chr'].str.extract(r'chr(\d+[A-Z]?|[A-Z]+)')
    df_terminal['haplotype'] = df_terminal['chr'].str.rsplit('_', n=1).str[-1]
    df_terminal['shortLabel'] = df_terminal['chrNum'] + df_terminal['haplotype'].str[0]

    # Check for missing data
    missing_tchest = df_terminal['TCHEST'].isna().sum()
    if missing_tchest > 0:
        print(f"Warning: {missing_tchest} records have missing TCHEST classification")
        df_terminal['TCHEST'] = df_terminal['TCHEST'].fillna('Absent')

    print(f"Final dataset: {len(df_terminal)} terminal telomere records")
    print(f"Unique chromosomes (chrNum): {df_terminal['chrNum'].nunique()}")
    
    return df_terminal


def prepare_stats(df: pd.DataFrame, value_col: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Prepare per-chromosome statistics and values dataframe.
    
    Args:
        df: Input dataframe with terminal telomeres
        value_col: Column name for values ('lengthKb' or 'canonicalProp')
    
    Returns:
        values_df: Long-form values per chromosome
        stats_df: Per-chromosome statistics
    """
    # Group by chrNum, haplotype, label (p/q arm)
    agg = (df.groupby(['chrNum', 'haplotype', 'label'], as_index=False)
             .agg({
                 value_col: 'max',
                 'TCHEST': lambda s: pick_tchest_class(s)
             }))
    agg = agg.rename(columns={value_col: 'value', 'TCHEST': 'tchest_any'})

    # Per-chromosome TCHEST class
    tchest_per_chr = (agg.groupby('chrNum')['tchest_any']
                        .apply(pick_tchest_class)
                        .rename('TCHEST_chr'))

    # Values DF (long form)
    values_df = agg[['chrNum', 'haplotype', 'label', 'value']].copy()
    values_df['end'] = values_df['haplotype'].str[0] + values_df['label'].str.upper()

    # Stats per chromosome
    stats = (values_df.groupby('chrNum')['value']
                      .agg(mean_val='mean', median_val='median', n_points='count')
                      .reset_index())
    stats = stats.merge(tchest_per_chr.reset_index(), on='chrNum', how='left')
    stats['TCHEST_chr'] = stats['TCHEST_chr'].fillna('Absent')

    return values_df, stats


def order_chromosomes_by_mean(stats_df: pd.DataFrame) -> list[str]:
    """Return chromosome order sorted by mean value."""
    srt = stats_df.sort_values(['mean_val', 'chrNum'], ascending=[True, True])
    return srt['chrNum'].tolist()


# ------------------------------- Plotting --------------------------------- #

def plot_chr_telomere_length(values_df: pd.DataFrame,
                              stats_df: pd.DataFrame,
                              order: list[str],
                              out_pdf: str,
                              out_svg: str,
                              out_png: str) -> None:
    """
    Scatter 4 values per chromosome with mean ticks. X tick labels colored by TCHEST.
    Legend in upper left.
    """
    tchest_map = dict(zip(stats_df['chrNum'], stats_df['TCHEST_chr']))
    chr_to_x = {ch: i for i, ch in enumerate(order, start=1)}
    n_chr = len(order)

    width = 0.16 * 41 + 1.2
    height = 3.0
    fig, ax = plt.subplots(figsize=(width, height), constrained_layout=True)

    # Light grey horizontal line at overall mean (kbp)
    overall_mean = values_df['value'].mean()
    ax.axhline(overall_mean, color="0.85", lw=1.0, zorder=1)
    print(f"1A: Overall mean telomere length = {overall_mean:.2f} kbp")

    slots = np.array([-0.15, -0.05, 0.05, 0.15])

    # Track which TCHEST classes are present for legend
    legend_handles = {}

    for ch in order:
        sub = values_df[values_df['chrNum'] == ch].copy()
        x0 = chr_to_x[ch]

        end_order = ['mP', 'mQ', 'pP', 'pQ']
        sub['end_key'] = pd.Categorical(sub['end'], categories=end_order, ordered=True)
        sub = sub.sort_values('end_key').reset_index(drop=True)

        sub['x'] = x0 + slots[:len(sub)]
        tchest_class = tchest_map.get(ch, 'Absent')
        col = TCHEST_COLORS.get(tchest_class, 'black')
        rgb = np.array(mpl.colors.to_rgb(col))
        pastel = 0.6 + 0.4 * rgb

        mp = sub['label'] == 'p'
        mq = sub['label'] == 'q'
        
        if mp.any():
            h = ax.scatter(sub.loc[mp, 'x'], sub.loc[mp, 'value'],
                          s=18, marker='s', linewidths=0.8, edgecolors=col,
                          facecolors=[pastel], zorder=3)
            if tchest_class not in legend_handles:
                legend_handles[tchest_class] = h
        if mq.any():
            h = ax.scatter(sub.loc[mq, 'x'], sub.loc[mq, 'value'],
                          s=18, marker='o', linewidths=0.8, edgecolors=col,
                          facecolors=[pastel], zorder=3)
            if tchest_class not in legend_handles:
                legend_handles[tchest_class] = h

        row = stats_df.loc[stats_df['chrNum'] == ch].iloc[0]
        mean_y = row['mean_val']
        ax.plot([x0 - 0.18, x0 + 0.18], [mean_y, mean_y],
                lw=1.6, color='black', solid_capstyle='butt', zorder=4)

    ax.set_xlim(0.5, n_chr + 0.5)
    ax.set_ylabel('Telomere Length (kbp)', labelpad=4)

    vmax = max(30, float(np.nanmax(values_df['value'].values)) if len(values_df) else 30)
    ax.set_ylim(0, vmax * 1.05)
    ax.set_yticks([5, 10, 15, 20, 25, 30])

    ax.set_xticks(range(1, n_chr + 1))
    ax.set_xticklabels(order, rotation=0)
    for tick, ch in zip(ax.get_xticklabels(), order):
        cls = tchest_map.get(ch, 'Absent')
        tick.set_color(TCHEST_COLORS.get(cls, 'black'))

    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)

    ax.set_xlabel('Chromosome', labelpad=6)

    # Add TCHEST type legend (lines for color only)
    legend_order = [c for c in TCHEST_PRIORITY if c in legend_handles]
    tchest_handles = [Line2D([0], [0], color=TCHEST_COLORS[c], lw=2, label=c)
                      for c in legend_order]
    leg1 = ax.legend(tchest_handles, legend_order, loc='upper left', title='TCHEST type',
                     fontsize=7, title_fontsize=8, bbox_to_anchor=(0.01, 1), framealpha=0.9)
    ax.add_artist(leg1)

    # Add Chr. arm legend (shapes in black) - positioned to the right of TCHEST legend
    arm_handles = [
        Line2D([0], [0], marker='s', color='black', linestyle='None',
               markersize=5, markerfacecolor='lightgray', markeredgewidth=0.8, label='p'),
        Line2D([0], [0], marker='o', color='black', linestyle='None',
               markersize=5, markerfacecolor='lightgray', markeredgewidth=0.8, label='q'),
    ]
    ax.legend(handles=arm_handles, loc='upper left', title='Chr. arm',
              fontsize=7, title_fontsize=8, bbox_to_anchor=(0.13, 1), framealpha=0.9)

    fig.savefig(out_pdf)
    fig.savefig(out_svg)
    fig.savefig(out_png, dpi=600)
    plt.close(fig)
    print(f"Saved: {out_pdf}, {out_svg}, {out_png}")


def plot_chr_canper(values_df: pd.DataFrame,
                    stats_df: pd.DataFrame,
                    order: list[str],
                    out_pdf: str,
                    out_svg: str,
                    out_png: str) -> None:
    """
    Scatter 4 values per chromosome using % canonical repeats.
    Legend in lower right.
    """
    tchest_map = dict(zip(stats_df['chrNum'], stats_df['TCHEST_chr']))
    chr_to_x = {ch: i for i, ch in enumerate(order, start=1)}
    n_chr = len(order)

    width = 0.16 * 41 + 1.2
    height = 3.0
    fig, ax = plt.subplots(figsize=(width, height), constrained_layout=True)

    # Grey horizontal line at overall mean (%)
    overall_mean = values_df['value'].mean()
    ax.axhline(overall_mean, color="0.85", lw=1.0, zorder=1)
    print(f"1B: Overall mean canonical proportion = {overall_mean:.2f}%")

    slots = np.array([-0.15, -0.05, 0.05, 0.15])

    # Track which TCHEST classes are present for legend
    legend_handles = {}

    for ch in order:
        sub = values_df[values_df['chrNum'] == ch].copy()
        x0 = chr_to_x[ch]

        end_order = ['mP', 'mQ', 'pP', 'pQ']
        sub['end_key'] = pd.Categorical(sub['end'], categories=end_order, ordered=True)
        sub = sub.sort_values('end_key').reset_index(drop=True)

        sub['x'] = x0 + slots[:len(sub)]
        tchest_class = tchest_map.get(ch, 'Absent')
        col = TCHEST_COLORS.get(tchest_class, 'black')
        rgb = np.array(mpl.colors.to_rgb(col))
        pastel = 0.6 + 0.4 * rgb

        mp = sub['label'] == 'p'
        mq = sub['label'] == 'q'
        
        if mp.any():
            h = ax.scatter(sub.loc[mp, 'x'], sub.loc[mp, 'value'],
                          s=18, marker='s', linewidths=0.8, edgecolors=col,
                          facecolors=[pastel], zorder=3)
            if tchest_class not in legend_handles:
                legend_handles[tchest_class] = h
        if mq.any():
            h = ax.scatter(sub.loc[mq, 'x'], sub.loc[mq, 'value'],
                          s=18, marker='o', linewidths=0.8, edgecolors=col,
                          facecolors=[pastel], zorder=3)
            if tchest_class not in legend_handles:
                legend_handles[tchest_class] = h

        row = stats_df.loc[stats_df['chrNum'] == ch].iloc[0]
        mean_y = row['mean_val']
        ax.plot([x0 - 0.18, x0 + 0.18], [mean_y, mean_y],
                lw=1.6, color='black', solid_capstyle='butt', zorder=4)

    ax.set_xlim(0.5, n_chr + 0.5)
    ax.set_ylabel('Canonical proportion (%)', labelpad=4)
    ax.set_ylim(75, 100)
    ax.set_yticks([75, 80, 85, 90, 95, 100])

    ax.set_xticks(range(1, n_chr + 1))
    ax.set_xticklabels(order, rotation=0)
    for tick, ch in zip(ax.get_xticklabels(), order):
        cls = tchest_map.get(ch, 'Absent')
        tick.set_color(TCHEST_COLORS.get(cls, 'black'))

    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)

    ax.set_xlabel('Chromosome', labelpad=6)

    # Add TCHEST type legend (lines for color only)
    legend_order = [c for c in TCHEST_PRIORITY if c in legend_handles]
    tchest_handles = [Line2D([0], [0], color=TCHEST_COLORS[c], lw=2, label=c)
                      for c in legend_order]
    leg1 = ax.legend(tchest_handles, legend_order, loc='lower right', title='TCHEST type',
                     fontsize=7, title_fontsize=8, bbox_to_anchor=(0.99, 0), framealpha=0.9)
    ax.add_artist(leg1)

    # Add Chr. arm legend (shapes in black) - positioned to the left of TCHEST legend
    arm_handles = [
        Line2D([0], [0], marker='s', color='black', linestyle='None',
               markersize=5, markerfacecolor='lightgray', markeredgewidth=0.8, label='p'),
        Line2D([0], [0], marker='o', color='black', linestyle='None',
               markersize=5, markerfacecolor='lightgray', markeredgewidth=0.8, label='q'),
    ]
    ax.legend(handles=arm_handles, loc='lower right', title='Chr. arm',
              fontsize=7, title_fontsize=8, bbox_to_anchor=(0.87, 0), framealpha=0.9)

    fig.savefig(out_pdf)
    fig.savefig(out_svg)
    fig.savefig(out_png, dpi=600)
    plt.close(fig)
    print(f"Saved: {out_pdf}, {out_svg}, {out_png}")


# --------------------------------- Main ----------------------------------- #

def main():
    # Set up paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    print(f"Working directory: {script_dir}")

    plot_dir = os.path.join(script_dir, "plots")
    os.makedirs(plot_dir, exist_ok=True)
    print(f"Output directory: {plot_dir}")

    # Load data
    try:
        df_terminal = load_data(script_dir)
    except FileNotFoundError as e:
        print(f"ERROR: {e}")
        return

    # Prepare stats for telomere length (kbp)
    print("\n--- Processing telomere length (kbp) ---")
    values_len, stats_len = prepare_stats(df_terminal, 'lengthKb')
    order_len = order_chromosomes_by_mean(stats_len)
    print(f"Chromosome order (by mean length): {order_len}")

    # Prepare stats for canonical proportion (%)
    print("\n--- Processing canonical proportion (%) ---")
    values_pct, stats_pct = prepare_stats(df_terminal, 'canonicalProp')
    order_pct = order_chromosomes_by_mean(stats_pct)
    print(f"Chromosome order (by mean %): {order_pct}")

    # Generate Figure A: Telomere Length sorted by mean
    print("\n--- Generating Figure A ---")
    plot_chr_telomere_length(
        values_df=values_len,
        stats_df=stats_len,
        order=order_len,
        out_pdf=os.path.join(plot_dir, "supfig16_1A_sorted_mean.pdf"),
        out_svg=os.path.join(plot_dir, "supfig16_1A_sorted_mean.svg"),
        out_png=os.path.join(plot_dir, "supfig16_1A_sorted_mean.png")
    )

    # Generate Figure B: Canonical % sorted by mean
    print("\n--- Generating Figure B ---")
    plot_chr_canper(
        values_df=values_pct,
        stats_df=stats_pct,
        order=order_pct,
        out_pdf=os.path.join(plot_dir, "supfig16_1B_sorted_mean.pdf"),
        out_svg=os.path.join(plot_dir, "supfig16_1B_sorted_mean.svg"),
        out_png=os.path.join(plot_dir, "supfig16_1B_sorted_mean.png")
    )

    print("\n--- Done ---")


if __name__ == "__main__":
    main()
