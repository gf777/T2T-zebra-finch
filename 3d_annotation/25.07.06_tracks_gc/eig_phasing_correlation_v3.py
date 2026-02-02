#!/usr/bin/env python3
"""
Generate per-chromosome GC cov vs. E1 scatter plots and genome-wide correlation summaries
for dip, mat+Z, and pat+WZ datasets.

Run from the directory containing your *.cis.vecs.tsv and *.gc.bedgraph files.  
Outputs will be written under ./plots/{dip,mat+Z,pat+WZ}/.
"""

import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

def make_plots():
    wd = os.getcwd()
    plot_root = os.path.join(wd, "plots")
    os.makedirs(plot_root, exist_ok=True)

    # find all cis.vecs.tsv files
    vecs_files = glob.glob("*.cis.vecs.tsv")

    for vecs_path in vecs_files:
        # infer haplotype key from filename: bTaeGut7.{hap}.cur....cis.vecs.tsv
        hap = os.path.basename(vecs_path).split(".")[1]

        # corresponding gc file
        gc_path = f"bTaeGut7.{hap}.cur.20250313.bwa.200000.gc.bedgraph"
        if not os.path.exists(gc_path):
            print(f"[!] Missing GC track: {gc_path}, skipping {hap}")
            continue

        # prepare output directory
        out_dir = os.path.join(plot_root, hap)
        os.makedirs(out_dir, exist_ok=True)

        # load data
        df_eigs = pd.read_csv(vecs_path, sep="\t")
        df_gc   = pd.read_csv(gc_path,   sep="\t")

        # merge on chrom/start/end
        df = pd.merge(
            df_eigs[["chrom", "start", "end", "E1"]],
            df_gc  [["chrom", "start", "end", "GC"]],
            on=["chrom", "start", "end"],
            how="inner",
        ).dropna(subset=["E1", "GC"])

        # compute per-chrom correlations and scatter plots
        corrs = {}
        for chrom, grp in df.groupby("chrom"):
            r = grp["E1"].corr(grp["GC"])
            corrs[chrom] = r

            fig, ax = plt.subplots(figsize=(4, 3))
            ax.scatter(grp["GC"] * 100, grp["E1"], s=5, alpha=0.5)
            ax.axhline(0, color="gray", linewidth=1)
            ax.set_title(f"{hap} · {chrom} (r={r:.2f})")
            ax.set_xlabel("GC content (%)")
            ax.set_ylabel("E1")
            fig.tight_layout()
            fig.savefig(os.path.join(out_dir, f"{chrom}.png"), dpi=600)
            plt.close(fig)

        # summary barplot of all correlations
        chs = sorted(corrs)
        vals = [corrs[c] for c in chs]
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.bar(chs, vals)
        ax.axhline(0, color="gray", linewidth=1)
        ax.set_xticks(range(len(chs)))
        ax.set_xticklabels(chs, rotation=90, fontsize=6)
        ax.set_ylabel("Pearson r (E1 vs GC)")
        ax.set_title(f"{hap} genome-wide E1/GC correlation")
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, f"{hap}_correlation_summary.png"), dpi=600)
        plt.close(fig)

if __name__ == "__main__":
    make_plots()
