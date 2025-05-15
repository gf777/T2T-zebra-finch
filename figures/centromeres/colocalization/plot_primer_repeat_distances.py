#!/usr/bin/env python3

import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pybedtools
from adjustText import adjust_text


def simplify_chrom_name(chrom):
    if chrom.startswith("chr"):
        chrom = chrom[3:]  # Remove "chr"
    if chrom.endswith("_mat"):
        return chrom[:-4] + "m"
    elif chrom.endswith("_pat"):
        return chrom[:-4] + "p"
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
            s1, s2 = group["start"].values
            if abs(s1 - s2) <= max_dist:
                merged = group.iloc[0].copy()
                merged["start"] = int(group["start"].mean())
                merged["end"] = int(group["end"].mean())
                merged["name"] = pid
                merged_rows.append(merged)
            else:
                merged_rows.extend(group.to_dict("records"))
        else:
            merged_rows.extend(group.to_dict("records"))

    merged_df = pd.DataFrame(merged_rows)
    if verbose:
        print(f"Merged primers: {len(df)} → {len(merged_df)} after collapsing close pairs")

    return merged_df


def load_primers(bed_file):
    df = pd.read_csv(bed_file, sep="\t", header=None, names=["chrom", "start", "end", "name"])
    df["type"] = df["name"].str.extract(r"^(Centromere|Distal)", expand=False).str.lower()
    return df


def load_repeats(gff_file, repeat_name, shared_chroms=None, min_length=0, largest_only=False, verbose=False):
    gff = pybedtools.BedTool(gff_file)
    filtered = gff.filter(lambda f: f[2] == "dispersed_repeat" and repeat_name in f[8])

    # Convert to BED-like format with length info
    bed_like = filtered.each(
        lambda f: pybedtools.create_interval_from_list([f[0], str(int(f[3]) - 1), f[4], repeat_name])).saveas()

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
    plt.figure(figsize=(6, 6))
    style_order = ["centromere", "distal"]
    style_dict = {"centromere": "D", "distal": "o"}
    ax = sns.scatterplot(
        data=df,
        x="dist_716A",
        y="dist_191A",
        hue="type",
        style="type",
        s=80,
        style_order=style_order, markers=style_dict
    )

    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Distance to Tgut716A (bp)")
    plt.ylabel("Distance to Tgut191A (bp)")
    plt.title(title)
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    # Match x/y limits to ensure equal visual scale in log space
    xlims = df["dist_716A"].agg(["min", "max"])
    ylims = df["dist_191A"].agg(["min", "max"])
    lower = min(xlims["min"], ylims["min"])
    upper = max(xlims["max"], ylims["max"])
    plt.xlim(lower, upper)
    plt.ylim(lower, upper)
    plt.plot([lower, upper], [lower, upper], ls="--", color="gray", lw=0.75, zorder=0, label="y = x")

    if label_chrom:
        texts = [
            plt.text(row["dist_716A"], row["dist_191A"], row["simple_chrom"],
                     fontsize=7, alpha=0.6, ha="center", va="center")
            for _, row in df.iterrows()
        ]
        adjust_text(texts, arrowprops=dict(arrowstyle="-", color='gray', lw=0.5))

    plt.tight_layout()
    plt.savefig(out_file, dpi=300)
    plt.close()


def main(primers_file, gff_file, out_plot, verbose=False, min_repeat_len=0, label_chrom=False, macrochr_file=None,
         largest_repeat_only=False):
    with open(macrochr_file) as f:
        macrochrs = set(line.strip() for line in f if line.strip())

    primers = load_primers(primers_file)
    primers = merge_primers(primers, max_dist=1000, verbose=verbose)
    primers["chr_class"] = primers["chrom"].apply(lambda c: "macro" if c in macrochrs else "micro")

    if verbose:
        print(f"Loaded {len(primers)} primers")
    primers["key"] = list(zip(
        primers["chrom"],
        primers["start"].astype(int),
        primers["end"].astype(int),
        primers["name"]
    ))

    shared_chroms = set(primers["chrom"].unique())
    repeats_716A = load_repeats(gff_file, "Tgut716A", shared_chroms, min_length=min_repeat_len,
                                largest_only=largest_repeat_only, verbose=verbose)
    repeats_191A = load_repeats(gff_file, "Tgut191A", shared_chroms, min_length=min_repeat_len,
                                largest_only=largest_repeat_only, verbose=verbose)

    dist_716A = compute_distance_dict(primers, repeats_716A, verbose)
    dist_191A = compute_distance_dict(primers, repeats_191A, verbose)

    primers["dist_716A"] = primers["key"].map(dist_716A)
    primers["dist_191A"] = primers["key"].map(dist_191A)

    primers_clean = primers.dropna(subset=["dist_716A", "dist_191A"]).copy()
    primers_clean["simple_chrom"] = primers_clean["chrom"].apply(simplify_chrom_name)
    primers_clean[["start", "end", "dist_716A", "dist_191A"]] = primers_clean[
        ["start", "end", "dist_716A", "dist_191A"]].astype(int)

    if verbose:
        print(f"{len(primers_clean)} primers have distances to both repeats")

    primers_macro = primers_clean[primers_clean["chr_class"] == "macro"].copy()
    primers_micro = primers_clean[primers_clean["chr_class"] == "micro"].copy()

    # Save distance table
    out_base = out_plot.rsplit(".", 1)[0]
    primers_clean[["chrom", "start", "end", "name", "type", "chr_class", "dist_716A", "dist_191A"]].to_csv(
        out_base + ".tsv", sep="\t", index=False)
    if verbose:
        print(f"Saved distance table to {out_base + '.tsv'}")

    # Plot macro and micro separately
    plot_subset(primers_macro, out_base + "_macro.png", label_chrom, title="Macrochromosomes")
    plot_subset(primers_micro, out_base + "_micro.png", label_chrom, title="Microchromosomes")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot distances of primers to centromeric repeats.")
    parser.add_argument("--primers", required=True, help="Tab-delimited BED file with primer coordinates and names.")
    parser.add_argument("--gff", required=True, help="GFF file with centromeric repeat annotations.")
    parser.add_argument("--out", default="primer_distances.png", help="Output base name for plots.")
    parser.add_argument("--min-repeat-len", type=int, default=0, help="Minimum repeat length (in bp) to include.")
    parser.add_argument("--macrochrs", required=True, help="File listing macrochromosome names (one per line)")
    parser.add_argument("--verbose", action="store_true", help="Enable debug output.")
    parser.add_argument("--label-chrom", action="store_true", help="Label each point with its chromosome name")
    parser.add_argument("--largest-repeat-only", action="store_true",
                        help="Only retain the largest repeat unit per chromosome"
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
