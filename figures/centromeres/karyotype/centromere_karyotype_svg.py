#!/usr/bin/env python3
"""
Generate an SVG "karyotype cartoon" of chromosomes with centromere positions.

- Chromosomes are drawn as vertical rounded bars, scaled by total length.
- The centromere span is shaded gray.
- A short red tick marks the midpoint of the centromere span.
- Layout is in rows, with 8 chromosomes per row by default.
- A blank "image slot" to the right of each chromosome can be reserved to paste real images later.

Inputs:
  1) A GFF3 (or GTF-like) centromere file with one or more "centromere" features per chromosome.
     We use columns 1 (seqid), 4 (start), 5 (end). If multiple per chr, the widest is taken.
  2) A FASTA index (.fai) file providing chromosome lengths (2nd column).

Example:
  python centromere_karyotype_svg.py --gff centromeres.gff3 --fai genome.fa.fai --out karyotype.svg \
    --columns 8 --chrom-height 360 --image-slot-width 120 --show-image-slots
"""
import argparse
import math
import re
import sys
from typing import Dict, Tuple, List

def read_fai(path):
    lens = {}
    with open(path) as f:
        for line in f:
            if not line.strip(): continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 2: continue
            try:
                lens[p[0]] = int(p[1])
            except ValueError:
                pass
    return lens

def read_centromeres(path):
    best = {}
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"): continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 5: continue
            name = p[0]
            try:
                s = int(p[3]); e = int(p[4])
            except ValueError:
                continue
            if e < s: s, e = e, s
            span = e - s
            if name not in best or span > (best[name][1] - best[name][0]):
                best[name] = (s, e)
    return best

def sort_chromosomes(names):
    def key(n):
        if n.startswith("chrZ"): return (10_000_000, 0, 0, n)
        if n.startswith("chrW"): return (10_000_001, 0, 0, n)
        m = re.match(r"^chr(\d+)([A-Za-z]?)(_.*)?$", n)
        if m:
            num = int(m.group(1))
            suff = m.group(2) or ""
            return (num, 0 if suff == "" else 1, ord(suff[0]) if suff else -1, n)
        return (9_999_999, 0, 0, n)
    return sorted(names, key=key)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--columns", type=int, default=8)
    ap.add_argument("--chrom-height", type=float, default=360.0, help="Reference height for the longest chromosome")
    ap.add_argument("--bar-width", type=float, default=14.0)
    ap.add_argument("--tick-width", type=float, default=18.0)
    ap.add_argument("--image-slot-width", type=float, default=140.0)
    ap.add_argument("--col-gap", type=float, default=70.0)
    ap.add_argument("--row-gap", type=float, default=160.0)
    ap.add_argument("--margin-left", type=float, default=60.0)
    ap.add_argument("--margin-right", type=float, default=60.0)
    ap.add_argument("--margin-top", type=float, default=40.0)
    ap.add_argument("--margin-bottom", type=float, default=60.0)
    ap.add_argument("--font-size", type=float, default=14.0)
    ap.add_argument("--label-dy", type=float, default=20.0)
    ap.add_argument("--show-image-slots", action="store_true")
    ap.add_argument("--show-centromere-band", action="store_true", help="Draw gray centromere span (off by default)")
    args = ap.parse_args()

    lengths = read_fai(args.fai)
    centrs = read_centromeres(args.gff)
    names = [n for n in sort_chromosomes(list(centrs.keys())) if n in lengths]
    if not names:
        sys.stderr.write("No overlapping chromosomes between GFF and FAI.\n")
        sys.exit(1)

    n = len(names); cols = max(1, args.columns); rows = math.ceil(n/cols)
    L_max = max(lengths[nm] for nm in names)
    heights = {nm: (args.chrom_height * (lengths[nm] / L_max)) for nm in names}

    # layout
    cell_w = args.bar_width + args.image_slot_width
    row_max_heights = []
    for r in range(rows):
        row_names = names[r*cols : min((r+1)*cols, n)]
        row_max_heights.append(max(heights[nm] for nm in row_names) if row_names else 0.0)

    total_w = args.margin_left + args.margin_right + (cols * cell_w) + ((cols - 1) * args.col_gap)
    total_h = args.margin_top + args.margin_bottom + sum(row_max_heights) + rows*args.label_dy + ((rows - 1) * args.row_gap)

    row_tops = []
    acc = args.margin_top
    for r in range(rows):
        row_tops.append(acc)
        acc += row_max_heights[r] + args.label_dy + (args.row_gap if r < rows-1 else 0.0)

    out = []
    def add(s): out.append(s)
    add(f'<svg xmlns="http://www.w3.org/2000/svg" width="{int(total_w)}" height="{int(total_h)}" viewBox="0 0 {int(total_w)} {int(total_h)}">')
    add('<style> text { font-family: Arial, Helvetica, sans-serif; } </style>')
    add(f'<rect x="0" y="0" width="{int(total_w)}" height="{int(total_h)}" fill="white"/>')

    for idx, name in enumerate(names):
        L = lengths[name]
        cstart, cend = centrs[name]
        h_i = heights[name]

        row = idx // cols; col = idx % cols
        top_y = row_tops[row]
        x0 = args.margin_left + col * (cell_w + args.col_gap) + args.bar_width / 2.0

        # ClipPath with rounded rect
        add(f'<clipPath id="clip_{idx}">')
        add(f'<rect x="{x0 - args.bar_width/2:.2f}" y="{top_y:.2f}" width="{args.bar_width:.2f}" height="{h_i:.2f}" '
            f'rx="{args.bar_width/2:.2f}" />')
        add(f'</clipPath>')

        # Base bar and outline
        add(f'<rect x="{x0 - args.bar_width/2:.2f}" y="{top_y:.2f}" width="{args.bar_width:.2f}" height="{h_i:.2f}" '
            f'rx="{args.bar_width/2:.2f}" fill="#ffffff" stroke="#222" stroke-width="1"/>')

        # Optional centromere band (clipped) with strict mirroring to upper half when centroid > 0.5
        s_raw, e_raw = (cstart, cend) if cstart <= cend else (cend, cstart)
        mid_frac = ((s_raw + e_raw) / 2.0) / L
        s_frac = s_raw / L
        e_frac = e_raw / L
        if mid_frac > 0.5:
            s_frac = 1.0 - s_frac
            e_frac = 1.0 - e_frac
        y_start = top_y + min(s_frac, e_frac) * h_i
        y_end   = top_y + max(s_frac, e_frac) * h_i
        cheight = max(1.0, y_end - y_start)
        if args.show_centromere_band:
            add(f'<rect clip-path="url(#clip_{idx})" x="{x0 - args.bar_width/2:.2f}" y="{y_start:.2f}" width="{args.bar_width:.2f}" height="{cheight:.2f}" '
                f'fill="#b7b7b7" stroke="none"/>')

        # Outline again on top
        add(f'<rect x="{x0 - args.bar_width/2:.2f}" y="{top_y:.2f}" width="{args.bar_width:.2f}" height="{h_i:.2f}" '
            f'rx="{args.bar_width/2:.2f}" fill="none" stroke="#222" stroke-width="1"/>')

        # Pink/red tick mirrored to upper half
        frac = ((cstart + cend) / (2.0 * L))
        if frac > 0.5:
            frac = 1.0 - frac
        y_mid = top_y + frac * h_i
        add(f'<line x1="{x0 - args.tick_width/2:.2f}" y1="{y_mid:.2f}" x2="{x0 + args.tick_width/2:.2f}" y2="{y_mid:.2f}" '
            f'stroke="#d94c8a" stroke-width="2"/>')

        if args.show_image_slots:
            slot_x = x0 + args.bar_width/2 + 6
            add(f'<rect x="{slot_x:.2f}" y="{top_y:.2f}" width="{args.image_slot_width - 6:.2f}" height="{h_i:.2f}" '
                f'fill="none" stroke="#cccccc" stroke-dasharray="4,4" stroke-width="1"/>')

        label_y = top_y + h_i + args.label_dy
        add(f'<text x="{x0:.2f}" y="{label_y:.2f}" font-size="{args.font_size:.1f}" text-anchor="middle">{name}</text>')

    add('</svg>')

    with open(args.out, "w") as fh:
        fh.write("\\n".join(out))

if __name__ == "__main__":
    main()
