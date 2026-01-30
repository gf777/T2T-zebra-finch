#!/usr/bin/env python3
"""
paf_chr_homology.py

Summarize a PAF file to report chromosome–chromosome homologies.

This version:
  - Legends at the bottom (lane-separated), with more breathing room (less overlap).
  - Saves plot as SVG by default (or by --plot-out extension).
  - Keeps the critical fix: Y-lanes are forced to match the final heatmap box height
    after aspect='equal' is applied.

Notes:
  - DPI is ignored for SVG (vector), but kept for compatibility.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, Iterator, List, Tuple, Optional
import matplotlib as mpl
mpl.rcParams["svg.fonttype"] = "none"   # <-- keep text as <text> in SVG (not paths)
mpl.rcParams["svg.hashsalt"] = "paf_chr_homology"  # optional: stable IDs across runs

# ---------------- I/O ----------------

def open_maybe_gzip(path: str):
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def bp_pretty(n: int) -> str:
    if n < 1000:
        return str(n)
    for unit in ["K", "M", "G", "T", "P"]:
        n_f = n / 1000.0
        if n_f < 1000:
            return f"{n_f:.3g}{unit}"
        n = int(n_f)
    return f"{n}E"


# ---------------- intervals / coverage ----------------

@dataclass(frozen=True)
class Interval:
    start: int
    end: int  # half-open [start, end)

    def length(self) -> int:
        return max(0, self.end - self.start)


def merge_intervals(intervals: List[Interval]) -> List[Interval]:
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: (x.start, x.end))
    merged: List[Interval] = []
    s, e = intervals[0].start, intervals[0].end
    for iv in intervals[1:]:
        if iv.start <= e:
            e = max(e, iv.end)
        else:
            merged.append(Interval(s, e))
            s, e = iv.start, iv.end
    merged.append(Interval(s, e))
    return merged


def parse_paf(handle: Iterable[str]) -> Iterator[Tuple[str, int, int, int, str, str, int, int]]:
    """Yield: (qname, qlen, qstart, qend, strand, tname, tlen, alnlen)"""
    for ln, line in enumerate(handle, 1):
        if not line.strip() or line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            raise ValueError(f"Line {ln}: malformed PAF (expected >=12 fields, got {len(f)})")
        yield (
            f[0],          # qname
            int(f[1]),     # qlen
            int(f[2]),     # qstart
            int(f[3]),     # qend
            f[4],          # strand
            f[5],          # tname
            int(f[6]),     # tlen
            int(f[10]),    # alnlen
        )


# ---------------- chromosome sorting ----------------

_SPECIAL_RANK = {
    "x": 1000,
    "y": 1001,
    "z": 1002,
    "w": 1003,
    "m": 2000,
    "mt": 2000,
    "mitochondrion": 2000,
    "mitochondria": 2000,
    "un": 9000,
    "unknown": 9001,
}


def chrom_sort_key(name: str):
    s = name.strip()
    s_low = s.lower()
    core = s_low[3:] if s_low.startswith("chr") else s_low

    m = re.match(r"^([a-z0-9]+)(.*)$", core)
    head = m.group(1) if m else core
    tail = m.group(2) if m else ""

    if head.isdigit():
        return (0, int(head), tail, s_low)

    if head in _SPECIAL_RANK:
        return (1, _SPECIAL_RANK[head], tail, s_low)

    m2 = re.match(r"^(\d+)([a-z].*)$", head)
    if m2:
        return (0, int(m2.group(1)), m2.group(2) + tail, s_low)

    return (2, head, tail, s_low)


def chrom_head_token(name: str) -> str:
    s = name.strip().lower()
    if s.startswith("chr"):
        s = s[3:]
    m = re.match(r"^([a-z0-9]+)", s)
    return m.group(1) if m else s


def strip_hap_suffix(label: str) -> str:
    return re.sub(r"(_mat|_pat)$", "", label, flags=re.IGNORECASE)


# ---------------- matrix I/O ----------------

def write_matrix(out, title: str, queries: List[str], targets: List[str], values: Dict[Tuple[str, str], int]):
    out.write(f"\n# {title}\n")
    out.write("query\t" + "\t".join(targets) + "\n")
    for q in queries:
        row = [q]
        for t in targets:
            row.append(str(values.get((q, t), 0)))
        out.write("\t".join(row) + "\n")


# ---------------- normalization (plot only) ----------------

def normalize_value(primary_bp: int, qlen: int, tlen: int, mode: str) -> float:
    if qlen <= 0 or tlen <= 0:
        return 0.0
    if mode == "raw":
        return float(primary_bp)
    if mode == "frac_query":
        return primary_bp / qlen
    if mode == "frac_target":
        return primary_bp / tlen
    if mode == "sym_geom":
        return primary_bp / math.sqrt(qlen * tlen)
    if mode == "jaccard":
        denom = (qlen + tlen - primary_bp)
        return (primary_bp / denom) if denom > 0 else 0.0
    raise ValueError(f"Unknown normalization mode: {mode}")


def maybe_log_transform(x: float, mode: str) -> float:
    if mode == "none":
        return x
    if mode == "log1p":
        return math.log1p(max(0.0, x))
    raise ValueError(f"Unknown log mode: {mode}")


# ---------------- annotation loaders ----------------

def load_chicken_chr_classification_tsv(tsv_path: str) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    with open_maybe_gzip(tsv_path) as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip() and not ln.startswith("#")]

    if not lines:
        raise ValueError(f"Empty chicken classification TSV: {tsv_path}")

    def split_row(r: str) -> List[str]:
        if "\t" in r:
            return r.split("\t")
        if "," in r:
            return r.split(",")
        return r.split()

    first = split_row(lines[0])
    headerish = any(x.strip().lower() in ("chr", "chrom", "chromosome", "name", "class", "type") for x in first)
    start_i = 1 if headerish else 0

    for r in lines[start_i:]:
        f = split_row(r)
        if len(f) < 2:
            continue
        tok = chrom_head_token(f[0])
        cat = f[1].strip()
        if tok and cat:
            mapping[tok] = cat

    if not mapping:
        raise ValueError(f"Failed to parse chicken classification TSV: {tsv_path}")
    return mapping


def load_zebrafinch_lanes_csv(csv_path: str) -> Dict[str, Dict[str, str]]:
    required = ["name", "chr_type", "classification_size", "dot", "tchest"]

    for skip in (0, 1):
        with open_maybe_gzip(csv_path) as fh:
            rows = list(csv.reader(fh))
        if not rows:
            break
        rows2 = rows[skip:] if skip else rows
        if not rows2:
            continue

        header = [h.strip() for h in rows2[0]]
        header_lc = [h.strip().lower() for h in header]
        if "name" not in header_lc:
            continue

        idx = {}
        ok = True
        for col in required:
            if col not in header_lc:
                ok = False
                break
            idx[col] = header_lc.index(col)
        if not ok:
            continue

        out: Dict[str, Dict[str, str]] = {}
        for r in rows2[1:]:
            if not r or len(r) <= max(idx.values()):
                continue
            name = r[idx["name"]].strip()
            if not name:
                continue
            out[name] = {
                "Chr_type": (r[idx["chr_type"]] or "").strip(),
                "Classification_size": (r[idx["classification_size"]] or "").strip(),
                "Dot": (r[idx["dot"]] or "").strip(),
                "TCHEST": (r[idx["tchest"]] or "").strip(),
            }

        if out:
            return out

    raise ValueError(
        f"Failed to parse zebra finch CSV lanes: {csv_path}. "
        "Need header with: Name, Chr_type, Classification_size, Dot, TCHEST."
    )


# ---------------- plotting helpers ----------------

def _parse_lane_labels(pairs: List[str]) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for p in pairs:
        if "=" not in p:
            continue
        k, v = p.split("=", 1)
        k = k.strip()
        v = v.strip()
        if k:
            out[k] = v
    return out


def _parse_lane_list(s: str) -> List[str]:
    s = (s or "").strip()
    if not s:
        return []
    return [x.strip() for x in s.split(",") if x.strip()]


def _infer_format_from_path(path: str) -> str:
    ext = os.path.splitext(path)[1].lower().lstrip(".")
    if ext in ("svg", "png", "pdf"):
        return ext
    return "svg"


def plot_matrix_svg(
    out_path: str,
    title: str,
    queries: List[str],
    targets: List[str],
    values: Dict[Tuple[str, str], float],
    *,
    dpi: int = 250,
    cell_in: float = 0.28,
    pad_in: float = 5.6,
    fontsize: int = 7,
    rotate_x: int = 90,
    x_anno: str = "none",
    y_anno: str = "none",
    x_lanes: Optional[List[str]] = None,
    y_lanes: Optional[List[str]] = None,
    lane_label_map: Optional[Dict[str, str]] = None,
    chicken_map: Optional[Dict[str, str]] = None,
    zebra_lanes: Optional[Dict[str, Dict[str, str]]] = None,
    show_anno: bool = True,
    # legend packing controls (more generous defaults)
    legend_fontsize: int = 9,
    legend_title_fontsize: int = 10,
    legend_cols: int = 3,
    legend_row_height: float = 1.55,   # inches per legend row
    legend_pad_x: float = 0.03,
    legend_pad_y: float = 0.12,
    legend_gap_x: float = 0.03,
    legend_gap_y: float = 0.18,
) -> None:
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.colors import to_rgba

    if not queries or not targets:
        return

    lane_label_map = lane_label_map or {}

    x_mode = (x_anno or "none").lower()
    y_mode = (y_anno or "none").lower()
    for m in (x_mode, y_mode):
        if m not in ("none", "chicken", "zebrafinch"):
            raise ValueError(f"Unknown anno mode: {m}")

    # numeric matrix
    mat = np.zeros((len(queries), len(targets)), dtype=float)
    for i, q in enumerate(queries):
        for j, t in enumerate(targets):
            mat[i, j] = values.get((q, t), 0.0)

    zebra_lane_defaults = ["Chr_type", "Classification_size", "Dot", "TCHEST"]

    def zebra_lookup(name: str, lane: str) -> str:
        if zebra_lanes is None:
            return "NA"
        if name in zebra_lanes:
            return (zebra_lanes[name].get(lane, "") or "").strip() or "NA"
        nl = name.lower()
        for k, vdict in zebra_lanes.items():
            if k.lower() == nl:
                return (vdict.get(lane, "") or "").strip() or "NA"
        return "NA"

    def chicken_lookup(name: str) -> str:
        if chicken_map is None:
            return "OTHER"
        tok = chrom_head_token(name)
        v = chicken_map.get(tok, "OTHER")
        return (v or "OTHER").strip().upper()

    # build lanes
    top_lane_names: List[str] = []
    top_lane_vals: Dict[str, List[str]] = {}
    left_lane_names: List[str] = []
    left_lane_vals: Dict[str, List[str]] = {}

    if show_anno and x_mode != "none":
        if x_mode == "chicken":
            if chicken_map is None:
                raise ValueError("x_anno=chicken requires chicken_map")
            top_lane_names = x_lanes if (x_lanes and len(x_lanes) > 0) else ["Chicken_class"]
            if top_lane_names != ["Chicken_class"]:
                raise ValueError("For x_anno=chicken, valid lane list is: Chicken_class")
            top_lane_vals["Chicken_class"] = [chicken_lookup(t) for t in targets]
        else:
            if zebra_lanes is None:
                raise ValueError("x_anno=zebrafinch requires zebra_lanes")
            top_lane_names = x_lanes if (x_lanes and len(x_lanes) > 0) else list(zebra_lane_defaults)
            for lane in top_lane_names:
                top_lane_vals[lane] = [zebra_lookup(t, lane) for t in targets]

    if show_anno and y_mode != "none":
        if y_mode == "chicken":
            if chicken_map is None:
                raise ValueError("y_anno=chicken requires chicken_map")
            left_lane_names = y_lanes if (y_lanes and len(y_lanes) > 0) else ["Chicken_class"]
            if left_lane_names != ["Chicken_class"]:
                raise ValueError("For y_anno=chicken, valid lane list is: Chicken_class")
            left_lane_vals["Chicken_class"] = [chicken_lookup(q) for q in queries]
        else:
            if zebra_lanes is None:
                raise ValueError("y_anno=zebrafinch requires zebra_lanes")
            left_lane_names = y_lanes if (y_lanes and len(y_lanes) > 0) else list(zebra_lane_defaults)
            for lane in left_lane_names:
                left_lane_vals[lane] = [zebra_lookup(q, lane) for q in queries]

    n_top = len(top_lane_names)
    n_left = len(left_lane_names)

    # legend blocks in plotting order
    legend_lane_order: List[Tuple[str, str]] = []
    for ln in top_lane_names:
        legend_lane_order.append(("X", ln))
    for ln in left_lane_names:
        legend_lane_order.append(("Y", ln))

    n_blocks = len(legend_lane_order)
    n_cols = max(1, int(legend_cols))
    n_rows = (n_blocks + n_cols - 1) // n_cols if n_blocks else 0

    # sizing
    w = max(10.0, cell_in * len(targets) + pad_in) + 0.65 * n_left
    h = max(10.0, cell_in * len(queries) + pad_in) + 0.40 * n_top + (legend_row_height * n_rows)

    fig = plt.figure(figsize=(w, h), constrained_layout=False)

    # generous margins to reduce overlap
    bottom_margin = 0.10 + 0.07 * n_rows
    bottom_margin = min(0.40, bottom_margin)
    fig.subplots_adjust(left=0.16, right=0.985, bottom=bottom_margin, top=0.92, wspace=0.02, hspace=0.02)

    top_ratio = max(0.75, float(n_top)) if n_top else 0.75
    left_ratio = max(0.75, float(n_left)) if n_left else 0.75
    legend_ratio = max(1.0, float(n_rows) * 1.4)

    gs = fig.add_gridspec(
        nrows=3, ncols=2,
        height_ratios=[top_ratio, len(queries), legend_ratio],
        width_ratios=[left_ratio, len(targets)],
        wspace=0.02, hspace=0.04,
    )

    ax_top = fig.add_subplot(gs[0, 1])
    ax_left = fig.add_subplot(gs[1, 0])
    ax = fig.add_subplot(gs[1, 1])
    ax_leg = fig.add_subplot(gs[2, :])
    ax_leg.axis("off")

    im = ax.imshow(mat, aspect="equal", interpolation="nearest")
    ax.set_title(title, pad=18)
    ax.set_xlabel("Chicken" if x_mode == "chicken" else "Target", labelpad=18)
    ax.set_ylabel("Zebra finch" if y_mode == "zebrafinch" else "Query", labelpad=18)

    disp_targets = [strip_hap_suffix(t) for t in targets]
    disp_queries = [strip_hap_suffix(q) for q in queries]

    ax.set_xticks(range(len(targets)))
    ax.set_xticklabels(disp_targets, rotation=rotate_x, fontsize=fontsize)
    ax.set_yticks(range(len(queries)))
    ax.set_yticklabels(disp_queries, fontsize=fontsize)
    ax.tick_params(axis="x", pad=8)
    ax.tick_params(axis="y", pad=6)

    ax.set_xlim(-0.5, len(targets) - 0.5)
    ax.set_ylim(len(queries) - 0.5, -0.5)

    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("normalized homology")

    lane_cmaps = ["tab10", "tab20", "Set3", "Dark2"]

    CHICKEN_PALETTE = {
        "MACRO": "lightgray",
        "MICRO": "gray",
        "DOT":   "black",
        "OTHER": "#7f7f7f",
        "NA":    "#cccccc",
    }

    ZEBRA_PALETTES = {
        "Dot": {"TRUE": "black", "FALSE": "gray", "NA": "#cccccc"},
        "Chr_type": {
            "Metacentric":    "#F0F0F0",
            "Submetacentric": "#CFCFCF",
            "Acrocentric":    "#8F8F8F",
            "Telocentric":    "#404040",
            "NA":             "#cccccc",
        },
        "Classification_size": {"Macro": "gray", "Micro": "black", "NA": "#cccccc"},
        "TCHEST": {"Absent": "lightgray", "Relaxed": "gray", "Strict": "black", "NA": "#cccccc"},
    }

    def build_lut(all_vals: List[str], cmap_name: str, fixed: Optional[Dict[str, str]] = None):
        vals = [(v.strip() or "NA") for v in all_vals]
        uniq = sorted(set(vals))
        cmap = plt.get_cmap(cmap_name)
        n = max(1, len(uniq))
        lut: Dict[str, Tuple[float, float, float, float]] = {}
        for i, v in enumerate(uniq):
            if fixed and v in fixed:
                lut[v] = to_rgba(fixed[v])
            else:
                lut[v] = to_rgba(cmap(i / max(1, n - 1)))
        return lut, uniq

    # lane legends
    import matplotlib.patches as mpatches
    lane_legends: Dict[str, Dict[str, object]] = {}

    def add_lane_legend(axis_prefix: str, lane: str, uniq_vals: List[str], lut):
        key = f"{axis_prefix}:{lane}"
        lane_title = lane_label_map.get(lane, lane)
        title_txt = f"{axis_prefix}  {lane_title}"
        handles = [mpatches.Patch(color=lut[v], label=v) for v in uniq_vals]
        lane_legends[key] = {"title": title_txt, "handles": handles}

    # top strip
    import numpy as np
    if n_top:
        top_rgba = np.zeros((n_top, len(targets), 4), dtype=float)
        for li, lane in enumerate(top_lane_names):
            cmap_name = lane_cmaps[li % len(lane_cmaps)]
            fixed = CHICKEN_PALETTE if lane == "Chicken_class" else ZEBRA_PALETTES.get(lane)
            vals = top_lane_vals[lane]
            lut, uniq = build_lut(vals, cmap_name, fixed=fixed)
            for j, v in enumerate(vals):
                top_rgba[li, j, :] = lut[v.strip() or "NA"]
            add_lane_legend("X", lane, uniq, lut)

        ax_top.imshow(top_rgba, interpolation="nearest", aspect="auto", origin="upper")
        ax_top.set_xticks([])
        ax_top.set_yticks(range(n_top))
        ax_top.set_yticklabels([lane_label_map.get(l, l) for l in top_lane_names], fontsize=fontsize)
        ax_top.tick_params(axis="y", pad=12)
        for s in ax_top.spines.values():
            s.set_visible(False)
    else:
        ax_top.axis("off")

    # left strip
    if n_left:
        left_rgba = np.zeros((len(queries), n_left, 4), dtype=float)
        for li, lane in enumerate(left_lane_names):
            cmap_name = lane_cmaps[(li + n_top) % len(lane_cmaps)]
            fixed = CHICKEN_PALETTE if lane == "Chicken_class" else ZEBRA_PALETTES.get(lane)
            vals = left_lane_vals[lane]
            lut, uniq = build_lut(vals, cmap_name, fixed=fixed)
            for i, v in enumerate(vals):
                left_rgba[i, li, :] = lut[v.strip() or "NA"]
            add_lane_legend("Y", lane, uniq, lut)

        ax_left.imshow(left_rgba, interpolation="nearest", aspect="auto", origin="upper")
        ax_left.set_yticks([])
        ax_left.set_xticks([])

        for li, lane in enumerate(left_lane_names):
            lab = lane_label_map.get(lane, lane)
            x = (li + 0.5) / max(1, n_left)
            ax_left.text(
                x, 1.18, lab,
                rotation=90, ha="center", va="bottom",
                fontsize=fontsize,
                transform=ax_left.transAxes,
                clip_on=False,
            )

        for s in ax_left.spines.values():
            s.set_visible(False)
    else:
        ax_left.axis("off")

    # --- CRITICAL ALIGNMENT STEP ---
    fig.canvas.draw()
    pos_main = ax.get_position()
    pos_left = ax_left.get_position()
    pos_top = ax_top.get_position()

    if n_left:
        ax_left.set_position([pos_left.x0, pos_main.y0, pos_left.width, pos_main.height])
        ax_left.set_ylim(ax.get_ylim())
    if n_top:
        ax_top.set_position([pos_main.x0, pos_top.y0, pos_main.width, pos_top.height])
        ax_top.set_xlim(ax.get_xlim())

    fig.canvas.draw()

    # legends bottom: less packed, larger gaps, slightly smaller fonts by default
    legend_keys = [f"{a}:{l}" for (a, l) in legend_lane_order]
    legend_blocks = [lane_legends[k] for k in legend_keys if k in lane_legends]

    if legend_blocks:
        n_blocks = len(legend_blocks)
        n_cols = max(1, int(legend_cols))
        n_rows = (n_blocks + n_cols - 1) // n_cols

        pad_x = float(legend_pad_x)
        pad_y = float(legend_pad_y)
        gap_x = float(legend_gap_x)
        gap_y = float(legend_gap_y)

        cell_w = (1.0 - 2 * pad_x - (n_cols - 1) * gap_x) / n_cols
        cell_h = (1.0 - 2 * pad_y - (n_rows - 1) * gap_y) / n_rows

        for idx, block in enumerate(legend_blocks):
            r = idx // n_cols
            c = idx % n_cols

            x0 = pad_x + c * (cell_w + gap_x)
            y_top = 1.0 - pad_y - r * (cell_h + gap_y)
            y0 = y_top - cell_h

            iax = ax_leg.inset_axes([x0, y0, cell_w, cell_h])
            iax.axis("off")

            iax.legend(
                handles=block["handles"],
                title=block["title"],
                loc="upper left",
                frameon=False,
                fontsize=legend_fontsize,
                title_fontsize=legend_title_fontsize,
                borderaxespad=0.0,
                handlelength=1.2,
                labelspacing=0.40,
                handletextpad=0.6,
                columnspacing=0.9,
            )

    fmt = _infer_format_from_path(out_path)
    fig.savefig(out_path, format=fmt, dpi=dpi)
    plt.close(fig)


def default_plot_path(out_path: str) -> str:
    if out_path and out_path != "-":
        return out_path + ".matrix_primary.svg"
    return "matrix_primary.svg"


# ---------------- main ----------------

def main() -> int:
    ap = argparse.ArgumentParser(
        description="Summarize chromosome–chromosome homology from PAF (long table + strand summary + matrices + SVG plot)."
    )
    ap.add_argument("paf", help="Input PAF (.paf/.paf.gz or '-')")
    ap.add_argument("-o", "--out", default="-", help="Output TSV (default: stdout)")

    ap.add_argument("--min-bp", type=int, default=0, help="Minimum primary bp to report (applies to long table & matrices)")
    ap.add_argument("--top", type=int, default=0, help="Top N targets per query in the long table (0=all)")
    ap.add_argument("--matrix-top", type=int, default=0, help="Top N targets per query to keep in matrices/plot (0=all)")

    ap.add_argument("--no-coverage", action="store_true", help="Disable merged query coverage (covered_bp)")
    ap.add_argument("--no-strand", action="store_true", help="Disable strand-separated summary")
    ap.add_argument("--no-pretty", action="store_true", help="Disable human-readable bp column")
    ap.add_argument("--no-normalize", action="store_true", help="Disable frac_of_query in long table")
    ap.add_argument("--no-matrix", action="store_true", help="Disable matrix text output blocks")
    ap.add_argument("--matrix-only-primary", action="store_true", help="Only output matrix for primary bp (skip sum_alnlen and covered matrices)")
    ap.add_argument("--no-plot", action="store_true", help="Disable plotted matrix output")

    # plot controls
    ap.add_argument("--plot-out", default="", help="Plot path (SVG recommended). Default: <out>.matrix_primary.svg or ./matrix_primary.svg")
    ap.add_argument("--plot-dpi", type=int, default=250, help="DPI (ignored for SVG; kept for PNG/PDF).")
    ap.add_argument("--plot-fontsize", type=int, default=7, help="Tick label fontsize (default: 7)")
    ap.add_argument("--plot-rotate-x", type=int, default=90, help="X label rotation (default: 90)")
    ap.add_argument("--plot-cell-in", type=float, default=0.28, help="Inches per matrix cell (default: 0.28)")
    ap.add_argument("--plot-pad-in", type=float, default=5.6, help="Extra inches padding for labels/colorbar (default: 5.6)")
    ap.add_argument(
        "--plot-norm",
        default="sym_geom",
        choices=["raw", "frac_query", "frac_target", "sym_geom", "jaccard"],
        help="Normalization used ONLY for the plotted matrix (default: sym_geom).",
    )
    ap.add_argument(
        "--plot-log",
        default="log1p",
        choices=["none", "log1p"],
        help="Transform for plotted values (default: log1p).",
    )

    # annotation configuration per axis
    ap.add_argument("--x-anno", default="chicken", choices=["none", "chicken", "zebrafinch"])
    ap.add_argument("--y-anno", default="zebrafinch", choices=["none", "chicken", "zebrafinch"])
    ap.add_argument("--chicken-class", default="chicken_chr_classification.tsv")
    ap.add_argument("--zebra-table", default="")
    ap.add_argument("--no-anno", action="store_true")

    # lane ordering + renaming
    ap.add_argument("--x-lanes", default="", help="Comma-separated lane order for X strip.")
    ap.add_argument("--y-lanes", default="", help="Comma-separated lane order for Y strip.")
    ap.add_argument("--lane-label", action="append", default=[], help="Override lane label: KEY=Pretty label.")

    # legend packing controls
    ap.add_argument("--legend-fontsize", type=int, default=9)
    ap.add_argument("--legend-title-fontsize", type=int, default=10)
    ap.add_argument("--legend-cols", type=int, default=3)
    ap.add_argument("--legend-row-height", type=float, default=1.55)
    ap.add_argument("--legend-pad-x", type=float, default=0.03)
    ap.add_argument("--legend-pad-y", type=float, default=0.12)
    ap.add_argument("--legend-gap-x", type=float, default=0.03)
    ap.add_argument("--legend-gap-y", type=float, default=0.18)

    args = ap.parse_args()

    x_lanes = _parse_lane_list(args.x_lanes)
    y_lanes = _parse_lane_list(args.y_lanes)
    lane_label_map = _parse_lane_labels(args.lane_label)

    need_chicken = (args.x_anno == "chicken") or (args.y_anno == "chicken")
    need_zebra = (args.x_anno == "zebrafinch") or (args.y_anno == "zebrafinch")

    chicken_map: Optional[Dict[str, str]] = None
    zebra_lane_map: Optional[Dict[str, Dict[str, str]]] = None

    if need_chicken:
        chicken_map = load_chicken_chr_classification_tsv(args.chicken_class)

    if need_zebra:
        if not args.zebra_table:
            raise SystemExit("ERROR: zebrafinch annotation requested but --zebra-table <csv> not provided.")
        zebra_lane_map = load_zebrafinch_lanes_csv(args.zebra_table)

    qlen: Dict[str, int] = {}
    tlen: Dict[str, int] = {}
    sum_aln: Dict[Tuple[str, str], int] = defaultdict(int)
    sum_aln_strand: Dict[Tuple[str, str, str], int] = defaultdict(int)
    intervals: Dict[Tuple[str, str], List[Interval]] = defaultdict(list)

    with open_maybe_gzip(args.paf) as fh:
        for q, ql, qs, qe, strand, t, tl, alnlen in parse_paf(fh):
            qlen.setdefault(q, ql)
            tlen.setdefault(t, tl)
            sum_aln[(q, t)] += alnlen
            sum_aln_strand[(q, t, strand)] += alnlen
            if not args.no_coverage and qe > qs:
                intervals[(q, t)].append(Interval(qs, qe))

    covered: Dict[Tuple[str, str], int] = {}
    if not args.no_coverage:
        for key, ivs in intervals.items():
            covered[key] = sum(iv.length() for iv in merge_intervals(ivs))

    primary: Dict[Tuple[str, str], int] = {}
    for key, aln_bp in sum_aln.items():
        primary[key] = covered.get(key, 0) if not args.no_coverage else aln_bp

    # ----- Long table -----
    rows = []
    for (q, t), prim_bp in primary.items():
        if prim_bp < args.min_bp:
            continue
        aln_bp = sum_aln.get((q, t), 0)
        cov_bp = covered.get((q, t), None) if not args.no_coverage else None
        frac = None
        if not args.no_normalize and qlen.get(q, 0) > 0:
            frac = prim_bp / qlen[q]
        rows.append((q, t, prim_bp, aln_bp, cov_bp, frac))

    rows.sort(key=lambda x: (chrom_sort_key(x[0]), -x[2], chrom_sort_key(x[1])))

    if args.top and args.top > 0:
        filtered = []
        seen = defaultdict(int)
        for r in rows:
            q = r[0]
            if seen[q] < args.top:
                filtered.append(r)
                seen[q] += 1
        rows = filtered

    out = sys.stdout if args.out == "-" else open(args.out, "w")

    header = ["query", "target", "bp_primary"]
    if not args.no_pretty:
        header.append("bp_primary_pretty")
    header.append("sum_alnlen_bp")
    if not args.no_coverage:
        header.append("covered_bp")
    if not args.no_normalize:
        header.append("query_len")
        header.append("frac_of_query")
    out.write("\t".join(header) + "\n")

    for q, t, prim_bp, aln_bp, cov_bp, frac in rows:
        fields = [q, t, str(prim_bp)]
        if not args.no_pretty:
            fields.append(bp_pretty(prim_bp))
        fields.append(str(aln_bp))
        if not args.no_coverage:
            fields.append(str(cov_bp if cov_bp is not None else 0))
        if not args.no_normalize:
            fields.append(str(qlen.get(q, 0)))
            fields.append(f"{frac:.6f}" if frac is not None else "")
        out.write("\t".join(fields) + "\n")

    # ----- Strand summary -----
    if not args.no_strand:
        out.write("\n# strand_summary\n")
        out.write("query\ttarget\tstrand\tsum_alnlen_bp\n")
        strand_items = sorted(
            sum_aln_strand.items(),
            key=lambda kv: (chrom_sort_key(kv[0][0]), -kv[1], chrom_sort_key(kv[0][1]), kv[0][2]),
        )
        for (q, t, s), v in strand_items:
            if v >= args.min_bp:
                out.write(f"{q}\t{t}\t{s}\t{v}\n")

    # ----- Matrices -----
    matrix_pairs = [(q, t, bp) for (q, t), bp in primary.items() if bp >= args.min_bp]

    if args.matrix_top and args.matrix_top > 0:
        per_q = defaultdict(list)
        for q, t, bp in matrix_pairs:
            per_q[q].append((t, bp))
        kept = set()
        for q, lst in per_q.items():
            lst.sort(key=lambda x: -x[1])
            for t, _ in lst[: args.matrix_top]:
                kept.add((q, t))
        matrix_pairs = [(q, t, bp) for (q, t, bp) in matrix_pairs if (q, t) in kept]

    queries = sorted({q for q, _, _ in matrix_pairs}, key=chrom_sort_key)
    targets = sorted({t for _, t, _ in matrix_pairs}, key=chrom_sort_key)

    primary_mat_int = {(q, t): bp for q, t, bp in matrix_pairs}

    if not args.no_matrix:
        write_matrix(out, "matrix_primary_bp", queries, targets, primary_mat_int)

        if not args.matrix_only_primary:
            aln_mat = {(q, t): sum_aln.get((q, t), 0) for q, t, _ in matrix_pairs}
            write_matrix(out, "matrix_sum_alnlen_bp", queries, targets, aln_mat)

            if not args.no_coverage:
                cov_mat = {(q, t): covered.get((q, t), 0) for q, t, _ in matrix_pairs}
                write_matrix(out, "matrix_covered_bp", queries, targets, cov_mat)

    # ----- Plot -----
    if not args.no_plot:
        plot_path = args.plot_out.strip() or default_plot_path(args.out)
        out_dir = os.path.dirname(plot_path)
        if out_dir and not os.path.isdir(out_dir):
            os.makedirs(out_dir, exist_ok=True)

        plot_vals: Dict[Tuple[str, str], float] = {}
        for q, t, bp in matrix_pairs:
            nv = normalize_value(bp, qlen.get(q, 0), tlen.get(t, 0), args.plot_norm)
            nv = maybe_log_transform(nv, args.plot_log)
            plot_vals[(q, t)] = nv

        plot_title = f"Homology matrix (norm={args.plot_norm}, log={args.plot_log})"
        plot_matrix_svg(
            out_path=plot_path,
            title=plot_title,
            queries=queries,
            targets=targets,
            values=plot_vals,
            dpi=args.plot_dpi,
            cell_in=args.plot_cell_in,
            pad_in=args.plot_pad_in,
            fontsize=args.plot_fontsize,
            rotate_x=args.plot_rotate_x,
            x_anno=args.x_anno,
            y_anno=args.y_anno,
            x_lanes=(x_lanes if x_lanes else None),
            y_lanes=(y_lanes if y_lanes else None),
            lane_label_map=lane_label_map,
            chicken_map=chicken_map,
            zebra_lanes=zebra_lane_map,
            show_anno=(not args.no_anno),
            legend_fontsize=args.legend_fontsize,
            legend_title_fontsize=args.legend_title_fontsize,
            legend_cols=args.legend_cols,
            legend_row_height=args.legend_row_height,
            legend_pad_x=args.legend_pad_x,
            legend_pad_y=args.legend_pad_y,
            legend_gap_x=args.legend_gap_x,
            legend_gap_y=args.legend_gap_y,
        )
        print(f"Wrote matrix plot: {plot_path}", file=sys.stderr)

    if out is not sys.stdout:
        out.close()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
