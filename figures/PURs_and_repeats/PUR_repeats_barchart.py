#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import pysam
from collections import defaultdict


def parse_gtf_attributes(attr_str):
    attrs = {}
    for part in attr_str.strip().split(";"):
        if part.strip() == "":
            continue
        key_value = part.strip().split(" ", 1)
        if len(key_value) == 2:
            key, val = key_value
            attrs[key] = val.strip('"')
    return attrs


def classify_repeat(feature_type):
    cls = feature_type.lower()
    if "transposon" in cls:
        return "Transposon"
    elif "retrotransposon" in cls:
        return "Retrotransposon"
    elif "satellite" in cls:
        return "Satellite"
    elif "simple" in cls:
        return "Simple"
    elif "ltr" in cls:
        return "LTR"
    elif "line" in cls:
        return "LINE"
    elif "sine" in cls:
        return "SINE"
    elif "repeat" in cls:
        return "Repeat"
    return "Other"


def simplify_chrom(chrom):
    chrom = chrom.replace("chr", "")
    if chrom.endswith("_mat"):
        return chrom.replace("_mat", ""), "mat"
    elif chrom.endswith("_pat"):
        return chrom.replace("_pat", ""), "pat"
    else:
        return chrom, "unphased"


def calculate_overlap_tabix(pur_bed, gtf_file):
    tbx = pysam.TabixFile(gtf_file)
    data = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))  # chrom -> hap -> class -> bp

    pur_df = pd.read_csv(pur_bed, sep="\t", header=None, names=["chrom", "start", "end"])

    for row in pur_df.itertuples(index=False):
        chrom, start, end = row.chrom, row.start, row.end
        base, hap = simplify_chrom(chrom)
        pur_len = end - start

        try:
            records = tbx.fetch(chrom, start, end)
        except ValueError:
            data[base][hap]["Unannotated"] += pur_len
            continue

        hit = False
        for rec in records:
            fields = rec.strip().split("\t")
            if len(fields) < 5:
                continue
            r_start, r_end = int(fields[3]), int(fields[4])
            ov_start = max(start, r_start)
            ov_end = min(end, r_end)
            ov_len = max(0, ov_end - ov_start)
            if ov_len > 0:
                cls = classify_repeat(fields[2])
                data[base][hap][cls] += ov_len
                hit = True

        if not hit:
            data[base][hap]["Unannotated"] += pur_len

    return data


def plot_stacked_two_bars_per_chromosome(data, out_file):
    records = []
    all_classes = set()

    for chrom in data:
        for hap in data[chrom]:
            for cls in data[chrom][hap]:
                mb = data[chrom][hap][cls] / 1e6
                records.append({
                    "Chromosome": chrom,
                    "Haplotype": hap,
                    "RepeatClass": cls,
                    "Mb": mb
                })
                all_classes.add(cls)

    df = pd.DataFrame(records)
    all_classes = sorted(all_classes)
    hap_order = ["mat", "pat"]

    # Sort chromosomes: numeric first, then Z/W
    def chrom_sort_key(x):
        try:
            return (0, int(x))
        except ValueError:
            return (1, x)

    chrom_order = sorted(set(df["Chromosome"]), key=chrom_sort_key)

    fig, ax = plt.subplots(figsize=(10, 7))

    # Assign one color per repeat class
    from matplotlib.cm import get_cmap
    cmap = get_cmap("tab20")
    class_colors = {cls: cmap(i % 20) for i, cls in enumerate(all_classes)}

    bar_width = 0.1
    bar_spacing = 0.1
    tick_positions = []
    tick_labels = []

    current_x = 0
    for chrom in chrom_order:
        chrom_df = df[df["Chromosome"] == chrom]
        plotted_bars = 0

        for i, hap in enumerate(hap_order):
            hap_df = chrom_df[chrom_df["Haplotype"] == hap]
            if hap_df.empty:
                continue

            bottom = 0
            for cls in all_classes:
                h = hap_df[hap_df["RepeatClass"] == cls]["Mb"].sum()
                if h > 0:
                    ax.bar(
                        current_x,
                        h,
                        bar_width,
                        bottom=bottom,
                        color=class_colors[cls],
                        label=cls if current_x == 0 else ""
                    )

                    bottom += h
            current_x += bar_width
            plotted_bars += 1

        # Add label in the center of the pair
        if plotted_bars > 0:
            center = current_x - (bar_width * plotted_bars) / 2
            tick_positions.append(center)
            tick_labels.append(chrom)

        current_x += bar_spacing  # space between chromosomes

    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, rotation=45, ha='right')
    ax.set_ylabel("Mb of PUR bases")
    ax.set_title("Stacked Repeat Content in PURs by Chromosome (mat/pat)")
    ax.legend(title="Repeat Class", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(out_file)
    print(f"Saved stacked bar chart to: {out_file}")


def main():
    parser = argparse.ArgumentParser(description="Stacked bar chart of PUR repeat content by haplotype.")
    parser.add_argument("--pur-bed", required=True, help="BED file with PUR coordinates")
    parser.add_argument("--repeat-gtf", required=True, help="bgzipped + tabix-indexed GTF with repeat annotations")
    parser.add_argument("--out", default="pur_repeat_stacked.png", help="Output PNG file")
    args = parser.parse_args()

    data = calculate_overlap_tabix(args.pur_bed, args.repeat_gtf)
    plot_stacked_two_bars_per_chromosome(data, args.out)


if __name__ == "__main__":
    main()
