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

def read_fai(path: str) -> Dict[str, int]:
    lens = {}
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            name = parts[0]
            try:
                length = int(parts[1])
            except ValueError:
                continue
            lens[name] = length
    return lens

def read_centromeres(path: str) -> Dict[str, Tuple[int, int]]:
    cent = {}
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            chrn = parts[0]
            try:
                start = int(parts[3])
                end = int(parts[4])
            except ValueError:
                continue
            if end < start:
                start, end = end, start
            # Keep the widest span per chromosome
            span = end - start
            if chrn not in cent or span > (cent[chrn][1] - cent[chrn][0]):
                cent[chrn] = (start, end)
    return cent

def sort_chromosomes(names: List[str]) -> List[str]:
    # Place autosomes and lettered variants first (1,1A,2,...), then Z then W at the very end
    def key(n: str):
        if n.startswith("chrZ"):
            return (10_000_000, 0, 0, n)
        if n.startswith("chrW"):
            return (10_000_001, 0, 0, n)
        m = re.match(r"^chr(\d+)([A-Za-z]?)(_.*)?$", n)
        if m:
            num = int(m.group(1))
            suff = m.group(2) or ""
            # '' sorts before 'A' naturally
            return (num, 0 if suff == "" else 1, ord(suff[0]) if suff else -1, n)
        # Fallback
        return (9_999_999, 0, 0, n)
    return sorted(names, key=key)

def clamp(v, lo, hi):
    return max(lo, min(hi, v))

def main():
    ap = argparse.ArgumentParser(description="Draw an SVG karyotype with centromeres.")
    ap.add_argument("--gff", required=True, help="Centromere GFF3 file")
    ap.add_argument("--fai", required=True, help="FASTA index (.fai) file with lengths")
    ap.add_argument("--out", required=True, help="Output SVG path")
    ap.add_argument("--columns", type=int, default=8, help="Chromosomes per row (default: 8)")
    ap.add_argument("--chrom-height", type=float, default=360.0, help="Height of each chromosome bar in pixels")
    ap.add_argument("--bar-width", type=float, default=14.0, help="Width of chromosome bar in pixels")
    ap.add_argument("--tick-width", type=float, default=18.0, help="Width of the red tick mark in pixels")
    ap.add_argument("--image-slot-width", type=float, default=140.0, help="Blank width to the right of each bar for images")
    ap.add_argument("--col-gap", type=float, default=70.0, help="Gap between columns")
    ap.add_argument("--row-gap", type=float, default=160.0, help="Gap between rows")
    ap.add_argument("--margin-left", type=float, default=60.0, help="Left margin")
    ap.add_argument("--margin-right", type=float, default=60.0, help="Right margin")
    ap.add_argument("--margin-top", type=float, default=40.0, help="Top margin")
    ap.add_argument("--margin-bottom", type=float, default=60.0, help="Bottom margin")
    ap.add_argument("--font-size", type=float, default=14.0, help="Label font size")
    ap.add_argument("--label-dy", type=float, default=20.0, help="Label y offset below bars")
    ap.add_argument("--show-image-slots", action="store_true", help="Draw faint rectangles showing the reserved image slot")
    args = ap.parse_args()

    lengths = read_fai(args.fai)
    centrs = read_centromeres(args.gff)

    # Keep only chromosomes that have both a length and a centromere span
    names = [n for n in centrs.keys() if n in lengths]
    if not names:
        sys.stderr.write("No overlapping chromosomes between GFF and FAI.\n")
        sys.exit(1)

    # Order with Z then W at the end
    names = sort_chromosomes(names)

    n = len(names)
    cols = max(1, args.columns)
    rows = math.ceil(n / cols)

    # Each "cell" contains [bar] + [blank image slot]
    cell_w = args.bar_width + args.image_slot_width
    total_w = args.margin_left + args.margin_right + (cols * cell_w) + ((cols - 1) * args.col_gap)
    total_h = args.margin_top + args.margin_bottom + (rows * (args.chrom_height + args.label_dy)) + ((rows - 1) * args.row_gap)

    # SVG header
    out = []
    def add(s): out.append(s)
    add(f'<svg xmlns="http://www.w3.org/2000/svg" width="{int(total_w)}" height="{int(total_h)}" viewBox="0 0 {int(total_w)} {int(total_h)}">')
    add('<style> text { font-family: Arial, Helvetica, sans-serif; } </style>')

    # Background
    add(f'<rect x="0" y="0" width="{int(total_w)}" height="{int(total_h)}" fill="white"/>')

    # Draw each chromosome
    for idx, name in enumerate(names):
        L = lengths[name]
        cstart, cend = centrs[name]
        if L <= 0:  # safety
            continue

        row = idx // cols
        col = idx % cols

        top_y = args.margin_top + row * (args.chrom_height + args.row_gap)
        # x coordinate for the center of the chromosome bar
        x0 = args.margin_left + col * (cell_w + args.col_gap) + args.bar_width / 2.0

        
        # --- Rounded-rect clipPath per chromosome to keep gray inside the capsule ---
        add(f'<clipPath id="clip_{idx}">')
        add(f'<rect x="{x0 - args.bar_width/2:.2f}" y="{top_y:.2f}" width="{args.bar_width:.2f}" height="{args.chrom_height:.2f}" '
            f'rx="{args.bar_width/2:.2f}" />')
        add(f'</clipPath>')

        # Base chromosome (white fill + stroke)
        add(f'<rect x="{x0 - args.bar_width/2:.2f}" y="{top_y:.2f}" width="{args.bar_width:.2f}" height="{args.chrom_height:.2f}" '
            f'rx="{args.bar_width/2:.2f}" fill="#ffffff" stroke="#222" stroke-width="1"/>')

        # Centromere rectangle (gray), clipped to rounded outline
        y_start = top_y + (cstart / L) * args.chrom_height
        y_end   = top_y + (cend   / L) * args.chrom_height
        y_start, y_end = min(y_start, y_end), max(y_start, y_end)
        cheight = max(1.0, y_end - y_start)
        add(f'<rect clip-path="url(#clip_{idx})" x="{x0 - args.bar_width/2:.2f}" y="{y_start:.2f}" width="{args.bar_width:.2f}" height="{cheight:.2f}" '
            f'fill="#b7b7b7" stroke="none"/>')

        # Re-draw a stroke-only outline so the border sits above the gray
        add(f'<rect x="{x0 - args.bar_width/2:.2f}" y="{top_y:.2f}" width="{args.bar_width:.2f}" height="{args.chrom_height:.2f}" '
            f'rx="{args.bar_width/2:.2f}" fill="none" stroke="#222" stroke-width="1"/>')

        # Red tick at the midpoint of the centromere
        y_mid = top_y + ((cstart + cend) / (2.0 * L)) * args.chrom_height
        add(f'<line x1="{x0 - args.tick_width/2:.2f}" y1="{y_mid:.2f}" x2="{x0 + args.tick_width/2:.2f}" y2="{y_mid:.2f}" '
            f'stroke="#d94c8a" stroke-width="2"/>')

        # Optional image slot rectangle (faint, to the right)
        if args.show_image_slots:
            slot_x = x0 + args.bar_width/2 + 6  # small padding from bar
            add(f'<rect x="{slot_x:.2f}" y="{top_y:.2f}" width="{args.image_slot_width - 6:.2f}" height="{args.chrom_height:.2f}" '
                f'fill="none" stroke="#cccccc" stroke-dasharray="4,4" stroke-width="1"/>')

        # Label (below bar)
        label_y = top_y + args.chrom_height + args.label_dy
        add(f'<text x="{x0:.2f}" y="{label_y:.2f}" font-size="{args.font_size:.1f}" text-anchor="middle">{name}</text>')

    add('</svg>')

    with open(args.out, "w") as fh:
        fh.write("\n".join(out))

if __name__ == "__main__":
    main()
