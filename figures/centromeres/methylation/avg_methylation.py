import argparse
import subprocess
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os

# Keep fonts as text in SVGs
plt.rcParams['svg.fonttype'] = 'none'


def bigwig_to_bedgraph(bigwig_file, chrom, start, end, bedgraph_file):
    """Extracts specific region from BigWig to BedGraph using bigwigToBedGraph."""
    try:
        subprocess.run([
            "./bigWigToBedGraph",
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

    # ---- Plot ----
    fig, ax = plt.subplots(figsize=(4, 6), layout='constrained')

    # Violin plot: cut=0 prevents KDE extending beyond data, which can push the y-axis oddly
    sns.violinplot(x="Repeat", y="Methylation", data=methylation_df, ax=ax, cut=0)

    # Jittered points; highlight best in red, others black
    if "Match_Type" in methylation_df.columns:
        sns.stripplot(x="Repeat", y="Methylation", data=methylation_df, hue="Match_Type",
                      palette={"best": "red", "other": "black"}, dodge=False, alpha=0.6, jitter=True, ax=ax)
        # Place legend outside to avoid crowding/cropping
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(handles, labels, title="Match type",
                      loc="lower right", frameon=True, borderaxespad=0.5)
    else:
        sns.stripplot(x="Repeat", y="Methylation", data=methylation_df,
                      color="black", alpha=0.6, jitter=True, ax=ax)

    ax.set_ylabel("Average methylation at CpG sites")
    ax.set_xlabel("Repeat region")

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

    # ---- Save PNG and SVG ----
    base, _ext = os.path.splitext(args.output)
    png_path = base + ".png"
    svg_path = base + ".svg"
    fig.savefig(png_path, dpi=300, bbox_inches='tight', pad_inches=0.02)
    fig.savefig(svg_path,           bbox_inches='tight', pad_inches=0.02)
    print(f"Plots saved to {png_path} and {svg_path}")

    # Optional interactive view
    plt.show()


if __name__ == "__main__":
    main()

