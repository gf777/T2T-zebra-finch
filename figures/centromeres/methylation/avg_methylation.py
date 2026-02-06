import argparse
import subprocess
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.stats import gaussian_kde  # for half-violin KDE
from matplotlib.lines import Line2D   # legend proxies

# Keep fonts as text in SVGs
plt.rcParams['svg.fonttype'] = 'none'


def bigwig_to_bedgraph(bigwig_file, chrom, start, end, bedgraph_file):
    """Extracts specific region from BigWig to BedGraph using bigwigToBedGraph."""
    try:
        subprocess.run([
            "bigWigToBedGraph",
            bigwig_file,
            bedgraph_file,
            f"-chrom={chrom}",
            f"-start={start}",
            f"-end={end}"
        ], check=True)
        print(f"Extracted region {chrom}:{start}-{end} to {bedgraph_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error during BigWig extraction for region {chrom}:{start}-{end}: {e}")
        return None
    if not os.path.exists(bedgraph_file):
        print(f"Error: BedGraph file {bedgraph_file} was not created.")
        return None
    return bedgraph_file


def get_methylation_from_bigwig(bigwig_file, gff_file, min_length=0, select_best=False, verbose=False):
    """
    Get methylation levels from a BigWig file based on regions in a GFF file.
    Extracts only the relevant regions from BigWig to a temporary BedGraph file.
    """
    regions_per_chromosome = {"Tgut716A": {}, "Tgut191A": {}}  # Store regions per chromosome for filtering
    generated_bedgraph_files = []  # List to store names of generated BedGraph files

    if verbose:
        print(f"BigWig file: {bigwig_file}")

    methylation_data = []  # Temporary list to hold methylation data for the DataFrame

    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue  # Skip comment lines
            fields = line.strip().split("\t")

            # Extract the repeat motif name from the Target attribute in the GFF
            target = None
            if "Target" in fields[8]:
                target = fields[8].split("Target ")[1].split()[0].split(":")[1].strip('"')  # Remove any double quotes

                if verbose:
                    print(f"Extracted Target: '{target}'")

            # Only collect Tgut716A and Tgut191A if the target matches
            if target in regions_per_chromosome:
                chrom = fields[0]
                start = int(fields[3])
                end = int(fields[4])

                # Apply the minimum length filter
                if (end - start) >= min_length:
                    # Convert the BigWig region to BedGraph for the specified region
                    bedgraph_file = f"{chrom}_{start}_{end}.bedgraph"
                    bedgraph_file = bigwig_to_bedgraph(bigwig_file, chrom, start, end, bedgraph_file)

                    if bedgraph_file and os.path.exists(bedgraph_file):
                        # Add BedGraph file to list for later removal
                        generated_bedgraph_files.append(bedgraph_file)

                        # Load the extracted BedGraph file into a DataFrame
                        df = pd.read_csv(bedgraph_file, sep="\t", header=None, names=["chrom", "start", "end", "value"])

                        # Filter the DataFrame to get the valid methylation values
                        valid_values = df["value"].dropna()

                        if verbose:
                            print(f"Valid methylation values for {chrom}:{start}-{end}: {valid_values}")

                        if not valid_values.empty:
                            # Calculate the average methylation for valid values
                            avg_methylation = valid_values.mean()

                            if chrom not in regions_per_chromosome[target]:
                                regions_per_chromosome[target][chrom] = []
                            regions_per_chromosome[target][chrom].append((start, end, avg_methylation))
                        else:
                            print(f"No valid methylation values for {chrom}, {target}, {start}, {end}")

                    else:
                        print(
                            f"Error with {chrom}, {target}, {start}, {end}: BedGraph file not found or extraction failed.")

                    if verbose and not valid_values.empty:
                        print(f"Region for {target}: {chrom}, {start}, {end} --> Avg methylation: {avg_methylation}")
                else:
                    if verbose:
                        print(f"Skipping region for {target}: {chrom}, {start}, {end} (too short)")

    # After gathering all regions, apply the "select best" logic if needed
    if select_best:
        print("Selecting the best (lowest methylation) region for each chromosome.")
        # Add a 'Match_Type' field in methylation_data to label best and other regions
        for repeat in regions_per_chromosome:
            for chrom, regions in regions_per_chromosome[repeat].items():
                # Sort the regions for each chromosome based on methylation (lowest first)
                sorted_regions = sorted(regions, key=lambda x: x[2])
                # Keep only the region with the lowest methylation
                best_region = sorted_regions[0]
                methylation_data.append((repeat, best_region[2], chrom, "best"))  # Tag best region
                # Add the rest of the regions (non-best) with "other" label
                for region in sorted_regions[1:]:
                    methylation_data.append((repeat, region[2], chrom, "other"))

    # Build the DataFrame from methylation_data
    if methylation_data:
        df_result = pd.DataFrame(methylation_data, columns=["Repeat", "Methylation", "Chromosome", "Match_Type"])

        # Cleanup: Remove all generated BedGraph files after use
        for bedgraph_file in generated_bedgraph_files:
            os.remove(bedgraph_file)
            print(f"Removed temporary file: {bedgraph_file}")

        return df_result
    else:
        print("No methylation data found.")
        return pd.DataFrame()  # Return an empty DataFrame


# Raincloud helpers (half-violin + thin box + jitter)
def _draw_half_violin(ax, center_x, vals, color, width_violin=0.32):
    if len(vals) > 1 and np.max(vals) > np.min(vals):
        kde = gaussian_kde(vals, bw_method=0.3)
        y = np.linspace(vals.min(), vals.max(), 300)
        dens = kde(y)
        if dens.max() > 0:
            scale = width_violin / dens.max()
            x_right = center_x + dens * scale
            ax.fill_betweenx(y, center_x, x_right, alpha=0.6, linewidth=0, color=color, zorder=1)
            ax.plot(x_right, y, linewidth=1.0, color=color, zorder=2)


def _draw_box(ax, center_x, vals, color, box_width=0.1):
    sns.boxplot(
        x=np.full_like(vals, center_x, dtype=float),
        y=vals,
        width=box_width,
        showcaps=True,
        boxprops={'facecolor': color, 'edgecolor': 'black', 'linewidth': 1.0},
        whiskerprops={'color': 'black'},
        capprops={'color': 'black'},
        medianprops={'color': 'black', 'linewidth': 1.0},
        showfliers=False,
        orient='v',
        ax=ax,
        zorder=3
    )


def main():
    # Set up the command line argument parser
    parser = argparse.ArgumentParser(description="Compare methylation levels at CpG sites for two repeat regions.")
    parser.add_argument("bigwig_file", help="Path to the BigWig file containing methylation data")
    parser.add_argument("gff_file", help="Path to the GFF file containing repeat annotations")
    parser.add_argument("--output", help="Output base name or filename (PNG & SVG will be written)",
                        default="methylation_comparison.png")
    parser.add_argument("--verbose", help="Enable verbose output for debugging", action="store_true")
    parser.add_argument("--min-length", type=int, default=0, help="Minimum length for repeat regions to be included")
    parser.add_argument("--select-best", help="Select the single best (lowest) candidate region per chromosome",
                        action="store_true")

    # Parse the command line arguments
    args = parser.parse_args()

    # Get methylation data
    methylation_df = get_methylation_from_bigwig(args.bigwig_file, args.gff_file,
                                                 min_length=args.min_length,
                                                 select_best=args.select_best,
                                                 verbose=args.verbose)

    if methylation_df.empty:
        print("No methylation data found. Please check your input files and regions.")
        return

    fig, ax = plt.subplots(figsize=(4, 6), layout='constrained')

    # Fixed order & palette like v10
    fam_order = ["Tgut716A", "Tgut191A"]
    palette = {"Tgut191A": "#009E73", "Tgut716A": "#CC79A7"}

    # Geometry tuned for 2-group panels
    rng = np.random.default_rng(6)
    jitter_sd = 0.04
    x_offset = 0.25
    violin_w = 0.32
    point_size = 15

    # Match-Type colors/styles
    dark_gray = "#4D4D4D"
    light_gray = "#C0C0C0"

    pos = np.arange(len(fam_order), dtype=float)
    for i, fam in enumerate(fam_order):
        sub = methylation_df[methylation_df["Repeat"].eq(fam)]
        vals = sub["Methylation"].dropna().values
        if len(vals) == 0:
            continue
        deep = palette[fam]

        _draw_half_violin(ax, pos[i], vals, deep, width_violin=violin_w)
        _draw_box(ax, pos[i], vals, deep)

        # Jittered points with "Match type" styling
        if "Match_Type" in sub.columns:
            sub_best = sub[sub["Match_Type"].eq("best")]
            sub_other = sub[sub["Match_Type"].eq("other")]

            if len(sub_other):
                x_jit_o = rng.normal(loc=pos[i] - x_offset, scale=jitter_sd, size=len(sub_other))
                ax.scatter(x_jit_o, sub_other["Methylation"].values, s=point_size,
                           linewidths=0.7, facecolors='none', edgecolors=light_gray, zorder=4)
            if len(sub_best):
                x_jit_b = rng.normal(loc=pos[i] - x_offset, scale=jitter_sd, size=len(sub_best))
                ax.scatter(x_jit_b, sub_best["Methylation"].values, s=point_size,
                           linewidths=0.7, facecolors=dark_gray, edgecolors=dark_gray, zorder=5)
        else:
            # Fallback: hollow points with family-colored edges
            x_jit = rng.normal(loc=pos[i] - x_offset, scale=jitter_sd, size=len(vals))
            ax.scatter(x_jit, vals, s=point_size, linewidths=0.6, facecolors='none', edgecolors=deep, zorder=4)

    ax.set_xticks(pos, fam_order, fontsize=10)

    # Compact Y label with % symbol
    ax.set_ylabel("Average CpG methylation %")
    ax.set_xlabel("")  # drop X-axis title

    # Clean spines
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    # Add a little vertical padding so nothing gets clipped
    yvals = methylation_df["Methylation"].to_numpy()
    yvals = yvals[np.isfinite(yvals)]
    if yvals.size:
        y_min, y_max = np.min(yvals), np.max(yvals)
        if y_max == y_min:
            pad = max(1.0, abs(y_max) if y_max != 0 else 1.0) * 0.05
        else:
            pad = (y_max - y_min) * 0.05
        ax.set_ylim(y_min - pad, y_max + pad)
    ax.margins(y=0.05)

    # Match type legend (only if present)
    if "Match_Type" in methylation_df.columns:
        handles = [
            Line2D([0], [0], marker='o', linestyle='None', markerfacecolor='none',
                   markeredgecolor=light_gray, markeredgewidth=1.0, label='other'),
            Line2D([0], [0], marker='o', linestyle='None', markerfacecolor=dark_gray,
                   markeredgecolor=dark_gray, label='best'),
        ]
        ax.legend(handles=handles, title="Match type",
                  loc="lower right", frameon=True, borderaxespad=0.5)

    # ---- Save PNG and SVG ----
    base, _ext = os.path.splitext(args.output)
    png_path = base + ".png"
    svg_path = base + ".svg"
    pdf_path = base + ".pdf"
    fig.savefig(png_path, dpi=600, bbox_inches='tight', pad_inches=0.02)
    fig.savefig(svg_path,           bbox_inches='tight', pad_inches=0.02)
    fig.savefig(pdf_path,           bbox_inches='tight', pad_inches=0.02)
    print(f"Plots saved to {png_path}, {svg_path} and {pdf_path}")

    # Optional interactive view
    plt.show()


if __name__ == "__main__":
    main()
