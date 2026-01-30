#!/usr/bin/env python3

import argparse
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pybedtools
from adjustText import adjust_text
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter

plt.rcParams.update({
    "svg.fonttype": "none",           # keep text as <text>, not paths
    "text.usetex": False,             # avoid TeX -> paths
    "axes.formatter.use_mathtext": False,
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans", "Nimbus Sans"]
})

def simplify_chrom_name(chrom):
    if chrom.startswith("chr"):
        chrom = chrom[3:]  # Remove "chr"
    if chrom.endswith("_pat") or chrom.endswith("_mat"):
        return chrom[:-4]
    return chrom

def merge_primers(df, max_dist=1000, verbose=False):
    df = df.copy()

    # Extract base name (strip trailing -F/-R or _F/_R or similar)
    df["pair_id"] = df["name"].str.extract(r"^(.+?)[_-][FR]$", expand=False)
    df["pair_id"] = df["pair_id"].fillna(df["name"])

    merged_rows = []
    grouped = df.groupby(["chrom", "pair_id"])

    for (chrom, pid), group in grouped:
        if len(group) == 2:
            s1, e1 = group["start"].values[0], group["end"].values[0]
            s2, e2 = group["start"].values[1], group["end"].values[1]

            # Check if the primers are too far apart or if they overlap
            if abs(s1 - s2) <= max_dist or (s1 <= e2 and s2 <= e1):  # Overlap condition
                # Merge the primers
                merged = group.iloc[0].copy()
                merged["start"] = min(s1, s2)  # Merge by taking the min start position
                merged["end"] = max(e1, e2)    # Merge by taking the max end position
                merged["name"] = pid
                merged_rows.append(merged.to_dict())  # Convert Series to dict
            else:
                # Keep both primers if they can't be merged
                merged_rows.extend(group.to_dict("records"))
        else:
            # If there's only one primer, just add it as is
            merged_rows.extend(group.to_dict("records"))

    # Ensure that merged_rows is a list of dictionaries
    merged_df = pd.DataFrame(merged_rows)

    if verbose:
        print(f"Merged primers: {len(df)} → {len(merged_df)} after collapsing close pairs")

    return merged_df

def load_primers(bed_file):
    df = pd.read_csv(bed_file, sep="\t", header=None, names=["chrom", "start", "end", "name"])
    df["type"] = df["name"].str.extract(r"^(Centromere|Distal|Unknown)", expand=False).str.lower()
    return df

def load_repeats(gff_file, repeat_name, shared_chroms=None, min_length=0, largest_only=False, verbose=False):
    gff = pybedtools.BedTool(gff_file)
    filtered = gff.filter(lambda f: f[2] == "dispersed_repeat" and repeat_name in f[8])

    # Convert to BED-like format with length info
    bed_like = filtered.each(lambda f: pybedtools.create_interval_from_list([f[0], str(int(f[3]) - 1), f[4], repeat_name])).saveas()

    bed_like_df = bed_like.to_dataframe(names=["chrom", "start", "end", "name"])
    bed_like_df["length"] = bed_like_df["end"] - bed_like_df["start"]
    bed_like_df = bed_like_df[bed_like_df["length"] >= min_length]

    if largest_only:
        bed_like_df = bed_like_df.sort_values("length", ascending=False)
        bed_like_df = bed_like_df.drop_duplicates("chrom", keep="first")
        if verbose:
            print(f"Repeat '{repeat_name}': retained only longest per chromosome")

    if shared_chroms:
        before = len(bed_like_df)
        bed_like_df = bed_like_df[bed_like_df["chrom"].isin(shared_chroms)]
        after = len(bed_like_df)
        if verbose:
            print(f"Repeat '{repeat_name}': {before} entries before chrom filter, {after} after")
    else:
        if verbose:
            print(f"Repeat '{repeat_name}': {len(bed_like_df)} entries total")

    return pybedtools.BedTool.from_dataframe(bed_like_df[["chrom", "start", "end", "name"]]).saveas()

def compute_distance_dict(primers_df, repeat_bed, verbose=False):
    primer_chroms = set(primers_df["chrom"].unique())
    repeat_chroms = set([f.chrom for f in repeat_bed])
    shared_chroms = primer_chroms & repeat_chroms

    if verbose:
        print("Primer chromosomes:", sorted(primer_chroms))
        print("Repeat chromosomes:", sorted(repeat_chroms))
        print("Shared chromosomes:", sorted(shared_chroms))

    if not shared_chroms:
        return {}

    filtered_primers = primers_df[primers_df["chrom"].isin(shared_chroms)].copy()
    filtered_primers = filtered_primers.sort_values(["chrom", "start", "end"])
    primers_bt = pybedtools.BedTool.from_dataframe(filtered_primers).saveas()
    repeat_bed = repeat_bed.filter(lambda f: f[0].strip() in shared_chroms).sort().saveas()

    closest = primers_bt.closest(repeat_bed, d=True, io=True)
    distances = {}
    for f in closest:
        key = (f[0], int(f[1]), int(f[2]), f[3])
        distances[key] = int(f[-1])

    if verbose:
        print(f"Returning {len(distances)} distances")
    return distances

def plot_subset(df, out_file, label_chrom=False, title=""):
    fig, ax = plt.subplots(figsize=(6, 6))

    # Color encodes "type"; shape encodes haplotype (maternal = circle, paternal = square)
    type_palette = {"centromere": "#1f77b4", "distal": "#ff7f0e", "unknown": "#d3d3d3"}
    hue_order = ["centromere", "distal", "unknown"]
    style_order = ["maternal", "paternal"]
    marker_map = {"maternal": "o", "paternal": "s"}

    # Match x/y limits to ensure equal visual scale in log space
    xlims = df["dist_716A"].agg(["min", "max"])
    ylims = df["dist_191A"].agg(["min", "max"])
    lower = min(xlims["min"], ylims["min"])
    upper = max(xlims["max"], ylims["max"])

    # Grid lines at 10^6 and 10^7 (behind data)
    tick_vals = [1e6, 1e7]
    for v in tick_vals:
        ax.axvline(v, color="lightgrey", lw=0.5, zorder=0)
        ax.axhline(v, color="lightgrey", lw=0.5, zorder=0)

    # y = x reference (behind data, same color as grid)
    ax.plot([lower, upper], [lower, upper], ls="--", color="lightgrey", lw=0.5, zorder=0)

    # Scatter plot (on top of grid)
    sns.scatterplot(
        data=df,
        x="dist_716A",
        y="dist_191A",
        hue="type",
        hue_order=hue_order,
        palette=type_palette,
        style="hap",
        style_order=style_order,
        markers=marker_map,
        s=50,
        alpha=0.55,
        legend=False,
        ax=ax,
        zorder=2
    )

    # Log scales
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Distance to Tgut716A (bp)", fontsize=9, fontfamily="Arial")
    ax.set_ylabel("Distance to Tgut191A (bp)", fontsize=9, fontfamily="Arial")

    # Add small padding in log space (shrink lower, expand upper by ~5%)
    pad_factor = 0.05
    log_lower = lower * (10 ** (-pad_factor))
    log_upper = upper * (10 ** pad_factor)
    ax.set_xlim(log_lower, log_upper)
    ax.set_ylim(log_lower, log_upper)

    # Chromosome labels with adjust_text to prevent overlaps
    texts = [
        ax.text(row["dist_716A"], row["dist_191A"], row["simple_chrom"],
                fontsize=7, fontfamily="Arial", ha="center", va="center", zorder=3,
                color='dimgrey')
        for _, row in df.iterrows()
    ]
    adjust_text(texts,
                arrowprops=dict(arrowstyle="-", color='lightgrey', lw=0.3),
                expand_points=(3, 3),
                force_points=(1.0, 1.0),
                force_text=(0.5, 0.5))

    # Tick label font and spine/tick thickness
    ax.tick_params(axis='both', labelsize=9, width=0.5)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontfamily("Arial")

    # Set spine (box) thickness
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)

    # Build legend handles for "type" and "hap"
    type_handles = [Line2D([0], [0], marker='o', linestyle='', color=type_palette[t],
                           label=t, markersize=6) for t in hue_order]
    hap_handles = [Line2D([0], [0], marker=marker_map[h], linestyle='', color="black",
                          label=h, markersize=6) for h in style_order]

    # --- Single box, vertical two-column (Type left, Haplotype right) ---
    from matplotlib.legend_handler import HandlerBase
    # Add invisible spacer for column titles
    spacer = Line2D([], [], linestyle='', marker='', label='')
    # Left column: Type title + items; Right column: Haplotype title + items
    all_handles = [spacer] + type_handles + [spacer] + hap_handles
    all_labels = ["Type"] + [h.get_label() for h in type_handles] + ["Haplotype"] + [h.get_label() for h in hap_handles]
    leg = ax.legend(handles=all_handles, labels=all_labels, loc="upper center",
                    ncol=2, frameon=True, fontsize=9, columnspacing=1.2,
                    handletextpad=0.5)
    leg.get_frame().set_linewidth(0.5)
    for t in leg.get_texts():
        t.set_fontfamily("Arial")

    plt.tight_layout()

    # Save in multiple formats
    base, _ = os.path.splitext(out_file)
    plt.savefig(base + ".svg")
    plt.savefig(base + ".pdf")
    plt.savefig(base + ".png", dpi=600)
    plt.close()

def main(primers_file, gff_file, out_plot, verbose=False, min_repeat_len=0,
         label_chrom=False, macrochr_file=None, largest_repeat_only=False):
    
    # Ensure plots go under ./plots
    out_dir = "./plots"
    os.makedirs(out_dir, exist_ok=True)
    out_plot = os.path.join(out_dir, os.path.basename(out_plot))

    primers = load_primers(primers_file)
    primers = merge_primers(primers, max_dist=1000, verbose=verbose)

    if verbose:
        print(f"Loaded {len(primers)} primers")
    primers["key"] = list(zip(
        primers["chrom"],
        primers["start"].astype(int),
        primers["end"].astype(int),
        primers["name"]
    ))

    shared_chroms = set(primers["chrom"].unique())
    repeats_716A = load_repeats(gff_file, "Tgut716A", shared_chroms, min_length=min_repeat_len, largest_only=largest_repeat_only, verbose=verbose)
    repeats_191A = load_repeats(gff_file, "Tgut191A", shared_chroms, min_length=min_repeat_len, largest_only=largest_repeat_only, verbose=verbose)

    dist_716A = compute_distance_dict(primers, repeats_716A, verbose)
    dist_191A = compute_distance_dict(primers, repeats_191A, verbose)

    primers["dist_716A"] = primers["key"].map(dist_716A)
    primers["dist_191A"] = primers["key"].map(dist_191A)

    primers_clean = primers.dropna(subset=["dist_716A", "dist_191A"]).copy()
    primers_clean["simple_chrom"] = primers_clean["chrom"].apply(simplify_chrom_name)
    primers_clean[["start", "end", "dist_716A", "dist_191A"]] = primers_clean[["start", "end", "dist_716A", "dist_191A"]].astype(int)

    # Infer haplotype from chromosome suffix (rely on original chrom only)
    def infer_hap(chrom):
        if chrom.endswith("_mat"):
            return "maternal"
        if chrom.endswith("_pat"):
            return "paternal"
        return "maternal"  # default fallback
    primers_clean["hap"] = primers_clean["chrom"].apply(infer_hap)

    if verbose:
        print(f"{len(primers_clean)} primers have distances to both repeats")

    # Save distance table alongside plots in ./plots
    out_base = os.path.splitext(os.path.basename(out_plot))[0]
    tsv_path = os.path.join(out_dir, out_base + ".tsv")
    primers_clean[["chrom", "start", "end", "name", "type", "dist_716A", "dist_191A"]].to_csv(
        tsv_path, sep="\t", index=False)
    if verbose:
        print(f"Saved distance table to {tsv_path}")

    if macrochr_file:
        with open(macrochr_file) as f:
            macrochrs = set(line.strip() for line in f if line.strip())
        primers_clean["chr_class"] = primers_clean["chrom"].apply(lambda c: "macro" if c in macrochrs else "micro")
        primers_macro = primers_clean[primers_clean["chr_class"] == "macro"].copy()
        primers_micro = primers_clean[primers_clean["chr_class"] == "micro"].copy()
        plot_subset(primers_macro, os.path.join(out_dir, out_base + "_macro.png"), label_chrom, title="Macrochromosomes")
        plot_subset(primers_micro, os.path.join(out_dir, out_base + "_micro.png"), label_chrom, title="Microchromosomes")
    else:
        # Unified plot for all chromosomes — no title per request
        plot_subset(primers_clean, out_plot, label_chrom, title="")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot distances of primers to centromeric repeats.")
    parser.add_argument("--primers", required=True, help="Tab-delimited BED file with primer coordinates and names.")
    parser.add_argument("--gff", required=True, help="GFF file with centromeric repeat annotations.")
    parser.add_argument("--out", default="primer_distances.svg", help="Output base name for plots.")
    parser.add_argument("--min-repeat-len", type=int, default=0, help="Minimum repeat length (in bp) to include.")
    parser.add_argument("--macrochrs", help="File listing macrochromosome names (one per line)")
    parser.add_argument("--verbose", action="store_true", help="Enable debug output.")
    parser.add_argument("--label-chrom", action="store_true", help="Label each point with its chromosome name")
    parser.add_argument("--largest-repeat-only", action="store_true", help="Only retain the largest repeat unit per chromosome"
)
    args = parser.parse_args()

    main(
        primers_file=args.primers,
        gff_file=args.gff,
        out_plot=args.out,
        verbose=args.verbose,
        min_repeat_len=args.min_repeat_len,
        label_chrom=args.label_chrom,
        macrochr_file=args.macrochrs,
        largest_repeat_only=args.largest_repeat_only
    )
