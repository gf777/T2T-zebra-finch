#!/usr/bin/env python3
import argparse
import csv
import re
import os
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional, Set
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec

# ---- Global style ----
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.family"] = "Arial"

TISSUES = ["brain", "liver", "spleen", "skin", "testis"]  # columns


# ----------------------------
# Helpers
# ----------------------------
def normalize_gene_id(g: str) -> str:
    if g is None:
        return ""
    g = g.strip()
    if "_" in g:
        g = g.split("_", 1)[0]
    return g


def gene_key(g: str):
    m = re.match(r"^PAK3L-(\d+)$", g)
    return (0, int(m.group(1))) if m else (1, g)


def human_bp(n: int) -> str:
    if n >= 1_000_000:
        return f"{n/1_000_000.0:.2f} Mb"
    if n >= 10_000:
        return f"{n/1_000.0:.1f} Kb"
    if n >= 1_000:
        return f"{n/1_000.0:.2f} Kb"
    return f"{n} bp"


def round_to_kbp(x: int) -> int:
    return int(round(x / 1000.0))


def fmt_int(x: int) -> str:
    return f"{int(x):,}"


def tissues_present(expr: str) -> Dict[str, int]:
    present = {t: 0 for t in TISSUES}
    if not expr or expr in {"-", "N.A.", "NA"}:
        return present
    parts = [p.strip().lower() for p in expr.split(",") if p.strip()]
    for p in parts:
        if p in present:
            present[p] = 1
    return present


# ----------------------------
# I/O
# ----------------------------
def read_counts_tsv(path: str):
    rows = []
    with open(path, "r", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        need = {"cluster_id", "cluster_chrom", "cluster_start", "cluster_end", "gene", "count"}
        miss = need - set(r.fieldnames or [])
        if miss:
            raise SystemExit(f"ERROR: {path} missing columns: {sorted(miss)}")
        for row in r:
            row["cluster_id"] = int(row["cluster_id"])
            row["cluster_start"] = int(row["cluster_start"])
            row["cluster_end"] = int(row["cluster_end"])
            row["count"] = int(row["count"])
            row["gene"] = normalize_gene_id(row["gene"])
            rows.append(row)
    return rows


def read_annos_tsv(path: str):
    ann = {}
    with open(path, "r", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        need = {"gene", "chr", "expr"}
        miss = need - set(r.fieldnames or [])
        if miss:
            raise SystemExit(f"ERROR: {path} missing columns: {sorted(miss)}")
        for row in r:
            gene = normalize_gene_id(row["gene"])
            chr_ = row["chr"].strip()
            expr = row["expr"].strip()
            if expr == ".":
                expr = ""
            ann[gene] = {"chr": chr_, "expr": expr}
    return ann


# ----------------------------
# Minimal Newick parser + pruning + flipping
# ----------------------------
@dataclass
class Node:
    name: str = ""
    length: float = 0.0
    children: Optional[List["Node"]] = None

    def is_leaf(self) -> bool:
        return not self.children


def _parse_newick(s: str) -> Node:
    s = s.strip()
    if s.endswith(";"):
        s = s[:-1]

    i = 0
    n = len(s)

    def parse_name_and_len() -> Tuple[str, float]:
        nonlocal i
        start = i
        while i < n and s[i] not in [":", ",", "(", ")"]:
            i += 1
        name = s[start:i].strip()

        length = 0.0
        if i < n and s[i] == ":":
            i += 1
            start = i
            while i < n and s[i] not in [",", "(", ")"]:
                i += 1
            val = s[start:i].strip()
            try:
                length = float(val)
            except ValueError:
                length = 0.0
        return name, length

    def parse_subtree() -> Node:
        nonlocal i
        if i < n and s[i] == "(":
            i += 1
            children: List[Node] = []
            while True:
                child = parse_subtree()
                children.append(child)
                if i >= n:
                    break
                if s[i] == ",":
                    i += 1
                    continue
                if s[i] == ")":
                    i += 1
                    break
            name, length = parse_name_and_len()
            return Node(name=name, length=length, children=children)
        else:
            name, length = parse_name_and_len()
            return Node(name=name, length=length, children=None)

    return parse_subtree()


def read_newick(path: str) -> Node:
    with open(path, "r") as f:
        s = f.read().strip()
    return _parse_newick(s)


def leaf_names(node: Node) -> List[str]:
    if node.is_leaf():
        return [node.name]
    out: List[str] = []
    for ch in node.children or []:
        out.extend(leaf_names(ch))
    return out


def prune_tree(node: Node, keep: Set[str]) -> Optional[Node]:
    if node.is_leaf():
        nm = normalize_gene_id(node.name)
        return node if nm in keep else None

    new_children: List[Node] = []
    for ch in node.children or []:
        pr = prune_tree(ch, keep)
        if pr is not None:
            new_children.append(pr)

    if not new_children:
        return None
    if len(new_children) == 1:
        only = new_children[0]
        only.length += node.length
        return only

    return Node(name=node.name, length=node.length, children=new_children)


def flip_tree_for_consistency(node: Node):
    if node.is_leaf():
        return
    for ch in node.children or []:
        flip_tree_for_consistency(ch)

    def subtree_min_key(n: Node):
        leaves = [normalize_gene_id(x) for x in leaf_names(n)]
        leaves = [x for x in leaves if x] or [""]
        return min((gene_key(x) for x in leaves))

    node.children = sorted(node.children or [], key=subtree_min_key)


def get_leaf_order(node: Node) -> List[str]:
    return [normalize_gene_id(x) for x in leaf_names(node)]


# ----------------------------
# Tree plotting (aligned to heatmap rows)
# ----------------------------
def layout_tree(node: Node, y_map: Dict[str, float], x0: float = 0.0, scale: float = 1.0):
    coords: Dict[int, Tuple[float, float]] = {}
    max_x = 0.0

    def rec(n: Node, x: float):
        nonlocal max_x
        x2 = x + (n.length or 0.0) * scale
        max_x = max(max_x, x2)

        if n.is_leaf():
            nm = normalize_gene_id(n.name)
            y = y_map[nm]
            coords[id(n)] = (x2, y)
            return y

        ys = []
        for ch in n.children or []:
            ys.append(rec(ch, x2))
        y = float(np.mean(ys)) if ys else 0.0
        coords[id(n)] = (x2, y)
        return y

    rec(node, x0)
    return coords, max_x


def draw_tree(ax, node: Node, coords: Dict[int, Tuple[float, float]]):
    def rec(n: Node):
        x, _y = coords[id(n)]
        if n.is_leaf():
            return

        cys = [coords[id(ch)][1] for ch in (n.children or [])]
        if cys:
            ax.plot([x, x], [min(cys), max(cys)], linewidth=1, color="black")

        for ch in n.children or []:
            cx, cy = coords[id(ch)]
            ax.plot([x, cx], [cy, cy], linewidth=1, color="black")
            rec(ch)

    rec(node)


# ----------------------------
# Output helpers
# ----------------------------
def derive_outpaths(out: str) -> Tuple[str, str]:
    base, ext = os.path.splitext(out)
    ext = ext.lower()
    if ext in {".svg", ".png"}:
        svg = base + ".svg"
        png = base + ".png"
    else:
        svg = out + ".svg"
        png = out + ".png"
    return svg, png


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True)
    ap.add_argument("--annos", required=True)
    ap.add_argument("--tree", required=True)
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--title", default="PAK3L copies per cluster (genes found on chrZ)")
    ap.add_argument("--tree_cluster_prefix", default="chrZ")
    ap.add_argument("--tree-scale", type=float, default=0.6)
    ap.add_argument("--png-dpi", type=int, default=300)

    # Colorbar spacing controls (figure fraction)
    ap.add_argument("--cbar-pad", type=float, default=0.006, help="Gap heatmap->colorbar (figure fraction)")
    ap.add_argument("--cbar-width", type=float, default=0.012, help="Colorbar width (figure fraction)")
    ap.add_argument("--right-pad", type=float, default=0.018, help="Gap colorbar->right panel (figure fraction)")

    # Right panel spacing controls (data units)
    ap.add_argument("--chr-col-width", type=float, default=None,
                    help="Reserved width for Chr column (right-panel data units). If omitted, auto-fit to longest label.")
    ap.add_argument("--chr-tissue-gap", type=float, default=0.35,
                    help="Gap between Chr column and first tissue column (right-panel data units).")
    ap.add_argument("--tissue-dx", type=float, default=0.55,
                    help="Spacing between tissue columns (right-panel data units).")
    ap.add_argument("--tissue-right-pad", type=float, default=0.40,
                    help="Extra right padding after last tissue column (right-panel data units) to avoid clipping.")
    ap.add_argument("--tissue-size", type=float, default=18,
                    help="Marker size (points^2) for tissue presence squares.")
    args = ap.parse_args()

    rows = read_counts_tsv(args.counts)
    ann = read_annos_tsv(args.annos)

    # Genes found on target clusters
    genes_found_on_target = set()
    for row in rows:
        if row["count"] > 0 and row["gene"] != "." and row["cluster_chrom"].startswith(args.tree_cluster_prefix):
            genes_found_on_target.add(row["gene"])
    if not genes_found_on_target:
        raise SystemExit(f"ERROR: no genes with count>0 on clusters starting with {args.tree_cluster_prefix}")

    # Tree prune/order
    tree = read_newick(args.tree)
    pruned = prune_tree(tree, genes_found_on_target)
    if pruned is None:
        raise SystemExit("ERROR: pruning removed all leaves. Check tree tip labels vs gene IDs.")
    flip_tree_for_consistency(pruned)
    genes = [g for g in get_leaf_order(pruned) if g in genes_found_on_target]
    genes_set = set(genes)

    # Cluster ordering
    cluster_meta: Dict[int, Tuple[str, int, int]] = {}
    for row in rows:
        cluster_meta[row["cluster_id"]] = (row["cluster_chrom"], row["cluster_start"], row["cluster_end"])
    cluster_ids = sorted(cluster_meta.keys(), key=lambda cid: (cluster_meta[cid][0], cluster_meta[cid][1], cluster_meta[cid][2], cid))

    # Count matrix
    counts = defaultdict(int)
    for row in rows:
        g = row["gene"]
        if g != "." and g in genes_set:
            counts[(row["cluster_id"], g)] += row["count"]

    M = np.zeros((len(genes), len(cluster_ids)), dtype=int)
    for j, cid in enumerate(cluster_ids):
        for i, g in enumerate(genes):
            M[i, j] = counts.get((cid, g), 0)

    row_tot = M.sum(axis=1)
    col_tot = M.sum(axis=0)

    # Right annotations
    chr_labels = [ann.get(g, {}).get("chr", "unknown") for g in genes]
    T = np.zeros((len(genes), len(TISSUES)), dtype=int)
    for i, g in enumerate(genes):
        expr = ann.get(g, {}).get("expr", "")
        pres = tissues_present(expr)
        for k, t in enumerate(TISSUES):
            T[i, k] = pres[t]

    # X labels: Kbp coords with commas
    xlabels = []
    for cid in cluster_ids:
        chrom, s, e = cluster_meta[cid]
        s_k = round_to_kbp(s)
        e_k = round_to_kbp(e)
        size = human_bp(e - s)
        if chrom.startswith(args.tree_cluster_prefix):
            coord = f"{fmt_int(s_k)}-{fmt_int(e_k)} Kbp"
        else:
            coord = f"{chrom}:{fmt_int(s_k)}-{fmt_int(e_k)} Kbp"
        xlabels.append(f"{coord} ({size})")

    n_genes = len(genes)
    n_clust = len(cluster_ids)

    # Figure
    width = max(12, min(0.55 * n_clust + 9.0, 32))
    height = max(7, min(0.35 * n_genes + 4.5, 24))
    fig = plt.figure(figsize=(width, height))

    gs = GridSpec(nrows=1, ncols=3, width_ratios=[2.8, 10.0, 1.35], wspace=0.01)
    ax_tree = fig.add_subplot(gs[0, 0])
    ax = fig.add_subplot(gs[0, 1])
    ax_right = fig.add_subplot(gs[0, 2])

    fig.suptitle(args.title, y=0.995)
    fig.subplots_adjust(top=0.90)

    # Heatmap (square cells)
    im = ax.imshow(M, interpolation="nearest", aspect="equal")
    ax.set_xlabel("Cluster (span)")
    ax.set_ylabel("")

    ax.set_xticks(range(n_clust))
    ax.set_xticklabels(xlabels, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(n_genes))
    ax.set_yticklabels(genes, fontsize=10)

    ax.set_ylim(n_genes - 0.5, -0.5)
    ax.set_xlim(-0.5, n_clust - 0.5)

    # Column totals (fixed band)
    xform = ax.get_xaxis_transform()
    for j in range(n_clust):
        ax.text(j, 1.01, fmt_int(col_tot[j]), transform=xform,
                ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax.text(1.005, 1.01, "Total", transform=ax.transAxes,
            ha="left", va="bottom", fontsize=9, fontweight="bold")

    # Row totals
    yform = ax.get_yaxis_transform()
    for i in range(n_genes):
        ax.text(1.01, i, fmt_int(row_tot[i]), transform=yform,
                ha="left", va="center", fontsize=9, fontweight="bold")

    # Per-cell counts
    cell_fs = 10 if n_clust <= 18 else 9 if n_clust <= 28 else 8
    for i in range(n_genes):
        for j in range(n_clust):
            v = int(M[i, j])
            if v > 0:
                ax.text(j, i, fmt_int(v), ha="center", va="center",
                        fontsize=cell_fs, fontweight="bold", color="white")

    # Tree aligned to heatmap
    y_map = {g: float(i) for i, g in enumerate(genes)}
    coords, max_x = layout_tree(pruned, y_map=y_map, x0=0.0, scale=args.tree_scale)
    draw_tree(ax_tree, pruned, coords)
    ax_tree.set_ylim(ax.get_ylim())
    ax_tree.set_yticks([])
    ax_tree.set_xlim(-0.02 * max_x, 1.02 * max_x)
    ax_tree.set_title("")
    ax_tree.set_xlabel("")
    for sp in ax_tree.spines.values():
        sp.set_visible(False)
    ax_tree.tick_params(axis="x", bottom=False, labelbottom=False)

    # ---- Manual colorbar + padding controls (figure fraction) ----
    pos_hm = ax.get_position()
    pos_rt = ax_right.get_position()

    cbar_x0 = pos_hm.x1 + args.cbar_pad
    cbar_w = args.cbar_width

    new_rt_x0 = cbar_x0 + cbar_w + args.right_pad
    rt_w = max(0.05, pos_rt.x1 - new_rt_x0)
    ax_right.set_position([new_rt_x0, pos_rt.y0, rt_w, pos_rt.height])

    cax = fig.add_axes([cbar_x0, pos_hm.y0, cbar_w, pos_hm.height])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("Copies")
    cax.tick_params(labelsize=9)

    # ---- Right panel: tissues extend fully + direct chr<->tissue spacing control ----
    ax_right.set_ylim(ax.get_ylim())
    ax_right.set_yticks([])
    ax_right.set_xlabel("")
    for sp in ax_right.spines.values():
        sp.set_visible(False)
    ax_right.tick_params(axis="x", length=0)

    tissue_dx = args.tissue_dx

    # Reserve Chr column width (either fixed by user, or auto-fit to longest label)
    chr_col_w = args.chr_col_width
    if chr_col_w is None:
        # auto-fit based on text extents (convert pixel width to data units)
        tmp = [ax_right.text(0.0, i, chr_labels[i], ha="left", va="center", fontsize=9, alpha=0.0)
               for i in range(n_genes)]
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        max_px = 0.0
        for t in tmp:
            bb = t.get_window_extent(renderer=renderer)
            max_px = max(max_px, bb.width)
            t.remove()

        # establish pixels-per-data-unit on x
        ax_right.set_xlim(0, 1)  # temporary
        fig.canvas.draw()
        p0 = ax_right.transData.transform((0, 0))[0]
        p1 = ax_right.transData.transform((1, 0))[0]
        px_per_data = max(1e-6, (p1 - p0))
        chr_col_w = max_px / px_per_data

    x_chr = 0.0
    x_t0 = chr_col_w + args.chr_tissue_gap
    x_tissues = [x_t0 + i * tissue_dx for i in range(len(TISSUES))]

    # IMPORTANT: pad the xlim so last column markers are not clipped (fixes "rectangles")
    x_max = x_tissues[-1] + args.tissue_right_pad
    ax_right.set_xlim(-0.05, x_max)

    # ticks/labels
    ax_right.set_xticks([0.15] + x_tissues)
    ax_right.set_xticklabels(["Chr"] + TISSUES, rotation=45, ha="right", fontsize=9)

    # draw chr labels
    for i in range(n_genes):
        ax_right.text(x_chr, i, chr_labels[i], ha="left", va="center", fontsize=9)

    # tissues
    ys, xs = np.where(T == 1)
    if len(xs) > 0:
        ax_right.scatter([x_tissues[x] for x in xs], ys, marker="s",
                         s=args.tissue_size, c="black", linewidths=0)

    # Save
    out_svg, out_png = derive_outpaths(args.out)
    fig.savefig(out_svg)
    fig.savefig(out_png, dpi=args.png_dpi)
    plt.close(fig)

    print(f"Wrote {out_svg}")
    print(f"Wrote {out_png}")


if __name__ == "__main__":
    main()
