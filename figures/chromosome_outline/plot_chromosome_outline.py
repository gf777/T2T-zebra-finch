#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import argparse
import re
from collections import defaultdict

plt.rcParams.update({
    'font.size': 28,
    'axes.titlesize': 32,
    'axes.labelsize': 28,
    'xtick.labelsize': 24,
    'ytick.labelsize': 24,
    'legend.fontsize': 24,
})


def natural_sort_key(s):
    return [int(t) if t.isdigit() else t for t in re.split(r'(\d+)', s)]


def plot_chromosome_full(gff_file, fai_file, out_file, flipped_file=None, submetacentric_file=None):
    flipped_chromosomes = set()
    if flipped_file:
        with open(flipped_file) as f:
            flipped_chromosomes = set(line.strip() for line in f if line.strip())

    submetacentric_chroms = set()
    if submetacentric_file:
        with open(submetacentric_file) as f:
            submetacentric_chroms = set(line.strip() for line in f if line.strip())

    df = pd.read_csv(gff_file, sep='\t', comment='#', header=None,
                     names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
    df = df[df['type'] == 'dispersed_repeat']
    df['repeat'] = df['attributes'].str.extract(r'Target "Motif:(Tgut(?:716A|191A))"')
    df = df.dropna(subset=['repeat'])
    df['base_chrom'] = df['seqid'].str.replace("_mat", "").str.replace("_pat", "").str.replace("^chr", "", regex=True)
    df['haplotype'] = df['seqid'].str.extract(r'_(mat|pat)')
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)

    chrom_sizes = pd.read_csv(fai_file, sep='\t', header=None, usecols=[0, 1], names=['seqid', 'length'])
    chrom_sizes = chrom_sizes[chrom_sizes['seqid'].str.contains('_mat|_pat')]
    chrom_sizes['base_chrom'] = chrom_sizes['seqid'].str.replace("_mat", "").str.replace("_pat", "").str.replace("^chr",
                                                                                                                 "",
                                                                                                                 regex=True)

    colors = {
        "Tgut716A": "#d94c8a",
        "Tgut191A": "#5c956f"
    }

    grouped = defaultdict(dict)
    for _, row in chrom_sizes.iterrows():
        grouped[row['base_chrom']][row['seqid']] = row['length']

    base_chromosomes = sorted(grouped.keys(), key=natural_sort_key)
    submetacentric_list = [c for c in base_chromosomes if c in submetacentric_chroms]
    telocentric_list = [c for c in base_chromosomes if c not in submetacentric_chroms]

    cols = 4
    mod = len(telocentric_list) % cols
    telocentric_list.extend([None for _ in range(cols - mod)])
    hap_spacing = 0.7
    chr_height = 0.5

    fig, ax = plt.subplots(figsize=(30, 24))
    y = 0
    label_x_positions = []
    label_y_positions = []
    label_names = []

    # Telocentric layout first (on lower side)
    column_chroms = [[] for _ in range(cols)]
    for idx, chrom in enumerate(reversed(telocentric_list)):
        column_chroms[idx % cols].append(chrom)

    column_chroms.reverse()

    col_offsets = [0]
    col_gap = 0.05 * max(chrom_sizes['length'])
    for i in range(1, cols):
        prev_chroms = column_chroms[i - 1]
        max_len = max([grouped[c][f"chr{c}_{h}"] for c in prev_chroms for h in ['mat', 'pat']
                       if f"chr{c}_{h}" in grouped[c]], default=0)
        col_offsets.append(col_offsets[-1] + max_len + col_gap)

    telocentric_tuples = list(zip(*column_chroms))
    for row_idx, chrom_row in enumerate(telocentric_tuples):
        y_offset = y + row_idx * (2 * hap_spacing + 1)
        for col, chrom in enumerate(chrom_row):
            if chrom is None:
                continue
            x_offset = col_offsets[col]
            chr_count = 0
            for i, hap in enumerate(["mat", "pat"]):
                hap_id = f"chr{chrom}_{hap}"
                if hap_id not in grouped[chrom]:
                    continue
                chrom_len = grouped[chrom][hap_id]
                sub = df[df['seqid'] == hap_id]
                y_pos = y_offset + i * hap_spacing
                ax.add_patch(mpatches.Rectangle((x_offset, y_pos - chr_height / 2), chrom_len, chr_height,
                                                facecolor='gainsboro', edgecolor='black', linewidth=0.5, zorder=0))
                for repeat_type in ["Tgut191A", "Tgut716A"]:
                    for _, row in sub[sub['repeat'] == repeat_type].iterrows():
                        start = row['start']
                        length = max(row['end'] - row['start'], 100_000)
                        if hap_id in flipped_chromosomes:
                            start = chrom_len - row['end']
                        ax.add_patch(mpatches.Rectangle((x_offset + start, y_pos - chr_height / 2), length, chr_height,
                                                        facecolor=colors[row['repeat']], edgecolor='none', zorder=2))
                chr_count += 1
            if chr_count > 0:
                label_x_positions.append(x_offset)
                label_y_positions.append(y_offset + hap_spacing / 2)
                label_names.append(chrom)

    last_col = column_chroms[-1]
    max_len_last_col = max(
        grouped[c][hap]
        for c in last_col if c is not None
        for hap in [f"chr{c}_mat", f"chr{c}_pat"]
        if hap in grouped[c]
    )
    telocentric_width = col_offsets[-1] + max_len_last_col
    ax.text(
        telocentric_width / 2,
        y_offset + (hap_spacing * 2),
        "Telocentric chromosomes",
        fontsize=32,
        ha="center",
        va="bottom"
    )
    y = y_offset + 2 * hap_spacing + 1

    # Plot submetacentric chromosomes higher up
    for chrom in reversed(submetacentric_list):
        plotted_y = []
        for i, hap in enumerate(["mat", "pat"]):
            hap_id = f"chr{chrom}_{hap}"
            if hap_id not in grouped[chrom]:
                continue
            chrom_len = grouped[chrom][hap_id]
            sub = df[df['seqid'] == hap_id]
            y_pos = y + i * hap_spacing
            ax.add_patch(mpatches.Rectangle((0, y_pos - chr_height / 2), chrom_len, chr_height,
                                            facecolor='lightgray', edgecolor='black', linewidth=0.5, zorder=0))
            for repeat_type in ["Tgut191A", "Tgut716A"]:
                for _, row in sub[sub['repeat'] == repeat_type].iterrows():
                    start = row['start']
                    length = max(row['end'] - row['start'], 100_000)
                    if hap_id in flipped_chromosomes:
                        start = chrom_len - row['end']
                    ax.add_patch(mpatches.Rectangle((start, y_pos - chr_height / 2), length, chr_height,
                                                    facecolor=colors[row['repeat']], edgecolor='none', zorder=2))
            plotted_y.append(y_pos)
        if plotted_y:
            label_y = sum(plotted_y) / len(plotted_y)  # center of whatever was actually plotted
            label_x_positions.append(0)
            label_y_positions.append(label_y)
            label_names.append(chrom)
        y += 2 * hap_spacing + 1

    ax.text(
        telocentric_width / 2,  # same horizontal center
        y - 1,  # slightly above first submetacentric
        "Metacentric, submetacentric, and acrocentric chromosomes",
        fontsize=32,
        ha="center",
        va="bottom"
    )

    ax.set_xlim(0, max(chrom_sizes['length']))
    ax.set_ylim(-1, y)

    for x, y, name in zip(label_x_positions, label_y_positions, label_names):
        ax.text(x - 450000, y, name, ha='right', va='center', fontsize=24)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines[:].set_visible(False)
    ax.legend(handles=[mpatches.Patch(color=v, label=k) for k, v in colors.items()], loc='upper right')

    plt.subplots_adjust(left=0.03, right=0.99, top=0.97, bottom=0.03)
    plt.savefig(out_file, dpi=300)
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot full chromosome structure with Tgut repeats.")
    parser.add_argument("--gff", required=True, help="RepeatMasker-style GFF with dispersed_repeat entries.")
    parser.add_argument("--fai", required=True, help="FASTA index file with chromosome sizes.")
    parser.add_argument("--out", default="full_chromosomes.png", help="Output image file (e.g. PNG, PDF)")
    parser.add_argument("--flip", help="Text file with chr_hap (e.g. chr16_pat) to flip vertically, one per line")
    parser.add_argument("--submetacentric",
                        help="Text file with base chromosome names that are submetacentric, one per line")
    args = parser.parse_args()

    plot_chromosome_full(args.gff, args.fai, args.out, flipped_file=args.flip, submetacentric_file=args.submetacentric)
