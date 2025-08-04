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
    'svg.fonttype': 'none',
    'figure.facecolor': 'none',
    'axes.facecolor':   'none',
})

def natural_sort_key(s):
    return [int(t) if t.isdigit() else t for t in re.split(r'(\d+)', s)]

def plot_chromosome_full(gff_file, fai_file, out_file,
                         flipped_file=None, submetacentric_file=None):
    # — read flip/submeta lists —
    flipped_chromosomes = set()
    if flipped_file:
        with open(flipped_file) as f:
            flipped_chromosomes = {l.strip() for l in f if l.strip()}
    submetacentric_chroms = set()
    if submetacentric_file:
        with open(submetacentric_file) as f:
            submetacentric_chroms = {l.strip() for l in f if l.strip()}

    # — load GFF —
    df = pd.read_csv(
        gff_file, sep='\t', comment='#', header=None,
        names=['seqid','source','type','start','end','score','strand','phase','attributes']
    )
    df = df[df['type']=='dispersed_repeat']
    df['repeat'] = df['attributes'].str.extract(
        r'Target "Motif:(Tgut(?:716A|191A))"'
    )
    df = df.dropna(subset=['repeat'])
    df['base_chrom'] = (
        df['seqid']
        .str.replace("_mat","")
        .str.replace("_pat","")
        .str.replace("^chr","",regex=True)
    )
    df['haplotype'] = df['seqid'].str.extract(r'_(mat|pat)')
    df[['start','end']] = df[['start','end']].astype(int)

    # — load .fai —
    chrom_sizes = pd.read_csv(
        fai_file, sep='\t', header=None, usecols=[0,1],
        names=['seqid','length']
    )
    chrom_sizes = chrom_sizes[
        chrom_sizes['seqid'].str.contains('_mat|_pat')
    ]
    chrom_sizes['base_chrom'] = (
        chrom_sizes['seqid']
        .str.replace("_mat","")
        .str.replace("_pat","")
        .str.replace("^chr","",regex=True)
    )

    colors = {"Tgut716A":"#d94c8a","Tgut191A":"#5c956f"}

    # — group lengths by chrom+hap —
    grouped = defaultdict(dict)
    for _, r in chrom_sizes.iterrows():
        grouped[r['base_chrom']][r['seqid']] = r['length']

    # — split lists —
    base_chromosomes    = sorted(grouped.keys(), key=natural_sort_key)
    submetacentric_list = [c for c in base_chromosomes if c in submetacentric_chroms]
    telocentric_list     = [c for c in base_chromosomes if c not in submetacentric_chroms]

    # — pair Z+W as single item —
    paired_sub = []
    paired = False
    for c in submetacentric_list:
        if c in ('Z','W') and not paired:
            paired_sub.append(('W','Z'))
            paired = True
        elif c not in ('Z','W'):
            paired_sub.append(c)
    submetacentric_items = paired_sub

    # — layout params —
    cols        = 4
    mod         = len(telocentric_list) % cols
    telocentric_list.extend([None] * (cols - mod))
    hap_spacing = 0.7
    chr_height  = 0.5

    fig, ax = plt.subplots(figsize=(30,15))
    y = 0.0
    label_x_positions = []
    label_y_positions = []
    label_names      = []

    # — telocentric grid —
    column_chroms = [[] for _ in range(cols)]
    for idx, chrom in enumerate(reversed(telocentric_list)):
        column_chroms[idx % cols].append(chrom)
    column_chroms.reverse()

    col_offsets = [0.0]
    col_gap     = 0.05 * max(chrom_sizes['length'])
    for i in range(1, cols):
        prev = column_chroms[i-1]
        max_len = (
            max(
                grouped[c][f"chr{c}_{h}"]
                for c in prev if c
                for h in ('mat','pat')
                if f"chr{c}_{h}" in grouped[c]
            ) if prev else 0
        )
        col_offsets.append(col_offsets[-1] + max_len + col_gap)

    for row_idx, chrom_row in enumerate(zip(*column_chroms)):
        y_offset = y + row_idx * (2 * hap_spacing + 1)
        for col_idx, chrom in enumerate(chrom_row):
            if chrom is None: continue
            x_off = col_offsets[col_idx]
            drew = 0
            for i, hap in enumerate(('mat','pat')):
                hap_id = f"chr{chrom}_{hap}"
                if hap_id not in grouped[chrom]: continue
                length = grouped[chrom][hap_id]
                sub    = df[df['seqid'] == hap_id]
                y_pos  = y_offset + i * hap_spacing

                ax.add_patch(mpatches.Rectangle(
                    (x_off, y_pos - chr_height/2),
                    length, chr_height,
                    facecolor='gainsboro', edgecolor='black',
                    linewidth=0.5, zorder=0
                ))
                for rt in ('Tgut191A','Tgut716A'):
                    for _, r in sub[sub['repeat']==rt].iterrows():
                        st = r['start']
                        seg = max(r['end'] - r['start'], 100_000)
                        if hap_id in flipped_chromosomes:
                            st = length - r['end']
                        ax.add_patch(mpatches.Rectangle(
                            (x_off + st, y_pos - chr_height/2),
                            seg, chr_height,
                            facecolor=colors[rt], edgecolor='none', zorder=2
                        ))
                drew += 1

            if drew:
                label_x_positions.append(x_off)
                label_y_positions.append(y_offset + hap_spacing/2)
                label_names.append(chrom)

    # telocentric title
    last_col = column_chroms[-1]
    max_last = max(
        grouped[c][hap]
        for c in last_col if c
        for hap in (f"chr{c}_mat",f"chr{c}_pat")
        if hap in grouped[c]
    )
    tel_width = col_offsets[-1] + max_last
    ax.text(
        tel_width/2,
        y_offset + 2 * hap_spacing,
        "Telocentric chromosomes",
        fontsize=32, ha='center', va='bottom'
    )

    # move into submetacentric
    y = y_offset + 2 * hap_spacing

    # — submetacentric, Z+W paired —
    row_gap = 0.5
    for item in reversed(submetacentric_items):
        if isinstance(item, tuple):
            # always plot Z first, then W
            for idx_sub, chrom in enumerate(item):
                if chrom not in item: continue
                avail = [h for h in ('mat','pat')
                         if f"chr{chrom}_{h}" in grouped[chrom]]
                if not avail: continue
                hap = avail[0]
                hap_id = f"chr{chrom}_{hap}"
                length = grouped[chrom][hap_id]
                sub    = df[df['seqid'] == hap_id]
                y_pos  = y + idx_sub * hap_spacing

                ax.add_patch(mpatches.Rectangle(
                    (0, y_pos - chr_height/2),
                    length, chr_height,
                    facecolor='lightgray', edgecolor='black',
                    linewidth=0.5, zorder=0
                ))
                for rt in ('Tgut191A','Tgut716A'):
                    for _, r in sub[sub['repeat']==rt].iterrows():
                        st = r['start']
                        seg = max(r['end'] - r['start'], 100_000)
                        if hap_id in flipped_chromosomes:
                            st = length - r['end']
                        ax.add_patch(mpatches.Rectangle(
                            (st, y_pos - chr_height/2),
                            seg, chr_height,
                            facecolor=colors[rt], edgecolor='none',
                            zorder=2
                        ))
                label_x_positions.append(0)
                label_y_positions.append(y_pos)
                label_names.append(chrom)
            y += 2 * hap_spacing + row_gap
            continue

        # normal submetacentric
        plotted = []
        for i, hap in enumerate(('mat','pat')):
            hap_id = f"chr{item}_{hap}"
            if hap_id not in grouped[item]: continue
            length = grouped[item][hap_id]
            sub    = df[df['seqid'] == hap_id]
            y_pos  = y + i * hap_spacing

            ax.add_patch(mpatches.Rectangle(
                (0, y_pos - chr_height/2),
                length, chr_height,
                facecolor='lightgray', edgecolor='black',
                linewidth=0.5, zorder=0
            ))
            for rt in ('Tgut191A','Tgut716A'):
                for _, r in sub[sub['repeat']==rt].iterrows():
                    st = r['start']
                    seg = max(r['end'] - r['start'], 100_000)
                    if hap_id in flipped_chromosomes:
                        st = length - r['end']
                    ax.add_patch(mpatches.Rectangle(
                        (st, y_pos - chr_height/2),
                        seg, chr_height,
                        facecolor=colors[rt], edgecolor='none',
                        zorder=2
                    ))
            plotted.append(y_pos)

        if plotted:
            label_x_positions.append(0)
            label_y_positions.append(sum(plotted)/len(plotted))
            label_names.append(item)
        y += 2 * hap_spacing + row_gap

    ax.text(
        tel_width/2,
        y - 1,
        "Metacentric, submetacentric, and acrocentric chromosomes",
        fontsize=32, ha='center', va='bottom'
    )

    # — finalize & save —
    ax.set_xlim(0, max(chrom_sizes['length']))
    ax.set_ylim(-1, y)
    for xx, yy, nm in zip(label_x_positions, label_y_positions, label_names):
        ax.text(xx - 450_000, yy, nm, ha='right', va='center', fontsize=24)
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.legend(
        handles=[mpatches.Patch(color=v, label=k) for k,v in colors.items()],
        loc='upper right'
    )
    plt.subplots_adjust(left=0.03, right=0.99, top=0.97, bottom=0.03)

    fig = plt.gcf()
    if out_file.lower().endswith('.svg'):
        fig.savefig(out_file, format='svg', dpi=300,
                    bbox_inches='tight', transparent=True)
    else:
        fig.savefig(out_file, dpi=300,
                    bbox_inches='tight', transparent=True)
    plt.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Plot full chromosome structure with Tgut repeats."
    )
    parser.add_argument('--gff',   required=True,
                        help="RepeatMasker-style GFF with dispersed_repeat entries.")
    parser.add_argument('--fai',   required=True,
                        help="FASTA index file with chromosome sizes.")
    parser.add_argument('--out',   default='full_chromosomes.png',
                        help="Output image file (PNG or SVG).")
    parser.add_argument('--flip',
                        help="File of chr_hap to flip vertically, one per line.")
    parser.add_argument('--submetacentric',
                        help="File of base chromosome names that are submetacentric.")
    args = parser.parse_args()
    plot_chromosome_full(
        args.gff, args.fai, args.out,
        flipped_file=args.flip,
        submetacentric_file=args.submetacentric
    )
