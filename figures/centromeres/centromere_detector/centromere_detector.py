import argparse
import subprocess
import os
import numpy as np
from pyBigWig import open as open_bw


def parse_gff(gff_file):
    """Parse GFF file and extract all repeat blocks."""
    blocks = []
    with open(gff_file) as f:
        for line in f:
            if not line.startswith("#"):  # Ignore comment lines
                parts = line.strip().split("\t")

                # Ensure we have at least 9 columns before proceeding
                if len(parts) >= 9:
                    chrom, source, feature, start, end, score, strand, phase, attributes = parts[:9]

                    # Add all repeat blocks to the list (no filtering at this point)
                    blocks.append((chrom, int(start), int(end), attributes))

                    # Check if Tgut716A (or any relevant motif) is in the attributes for centromere candidate
                    if "Tgut716A" in attributes:
                        print(f"Found Tgut716A block in {chrom}:{start}-{end}")
                        blocks[-1] = (chrom, int(start), int(end), attributes, "centromere_candidate")
                else:
                    print(f"Skipping line with insufficient columns: {line.strip()}")
    return blocks


def bigwig_to_bedgraph(bigwig_file, chrom, start, end, bedgraph_file):
    """Extract a specific region from BigWig to BedGraph using bigwigToBedGraph."""
    try:
        # Run bigwigToBedGraph for the specified region and ensure the BedGraph file is created
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


def calculate_average_methylation(bedgraph_file):
    """Calculate the average methylation level from a BedGraph file."""
    methylation_values = []
    with open(bedgraph_file) as f:
        for line in f:
            values = line.strip().split()
            methylation_values.append(float(values[3]))  # Methylation value is in the 4th column
    return np.mean(methylation_values)


def find_centromere(blocks, bw_file, min_length):
    """Identify the Tgut716A block with the lowest average methylation level per chromosome."""
    chromosome_centromeres = {}

    for block in blocks:
        # Ensure we correctly unpack 4 or 5 element tuples
        if len(block) == 4:  # Normal repeat block
            chrom, start, end, attributes = block
            centromere_candidate = None
        elif len(block) == 5:  # Centromere candidate block
            chrom, start, end, attributes, centromere_candidate = block
        else:
            print(f"Skipping block with unexpected format: {block}")
            continue

        # Check if the block length is within the specified range
        region_length = end - start
        if region_length < min_length:
            continue  # Skip blocks that do not meet the length criteria

        # Only process centromere candidate blocks
        if centromere_candidate == "centromere_candidate":
            bedgraph_file = f"temp_{chrom}_{start}_{end}.bedgraph"
            bigwig_to_bedgraph(bw_file, chrom, start, end, bedgraph_file)
            avg_methylation = calculate_average_methylation(bedgraph_file)

            # Ensure avg_methylation is a float
            if isinstance(avg_methylation, str):
                try:
                    avg_methylation = float(avg_methylation)
                except ValueError:
                    avg_methylation = 0.0  # Default value in case of an error

            # Check if this is the block with the lowest methylation for this chromosome
            if chrom not in chromosome_centromeres:
                # First centromere candidate for the chromosome
                chromosome_centromeres[chrom] = (chrom, start, end, attributes, avg_methylation)
            else:
                # Compare with existing methylation value and select the one with the lowest methylation
                current_best = chromosome_centromeres[chrom]
                if avg_methylation < current_best[4]:  # Compare methylation values
                    chromosome_centromeres[chrom] = (chrom, start, end, attributes, avg_methylation)

            # Remove the temporary BedGraph file after use
            if os.path.exists(bedgraph_file):
                os.remove(bedgraph_file)
                print(f"Removed temporary file: {bedgraph_file}")

    # Return the best centromere candidates per chromosome with their methylation value
    return [chromosome_centromeres[chrom] for chrom in chromosome_centromeres]

def annotate_kinetochore_binding_site(chrom, start, end, bw_file, avg_methylation):
    """Identify regions with significantly lower methylation within the centromere."""

    # Ensure avg_methylation is a float
    if isinstance(avg_methylation, str):
        try:
            avg_methylation = float(avg_methylation)
        except ValueError:
            avg_methylation = 0.0  # Default value in case of an error

    bedgraph_file = f"temp_{chrom}_{start}_{end}.bedgraph"
    bigwig_to_bedgraph(bw_file, chrom, start, end, bedgraph_file)
    methylation_values = []
    positions = []

    with open(bedgraph_file) as f:
        for line in f:
            values = line.strip().split()
            methylation_values.append(float(values[3]))
            positions.append(int(values[1]))  # Start position of the region

    # Determine significant low-methylation regions (e.g., 30% below average)
    threshold = avg_methylation * 0.7
    low_methylation_regions = []

    for i in range(len(methylation_values)):
        if methylation_values[i] < threshold:
            low_methylation_regions.append((positions[i], methylation_values[i]))

    # Remove the temporary BedGraph file after use
    if os.path.exists(bedgraph_file):
        os.remove(bedgraph_file)
        print(f"Removed temporary file: {bedgraph_file}")

    return low_methylation_regions


def output_gff(candidate_centromeres, kinetochore_regions, output_file):
    """Output the GFF file with annotated regions."""
    with open(output_file, "w") as f:
        for candidate_centromere in candidate_centromeres:
            # Unpack the 5-element tuple
            chrom, start, end, attributes, avg_methylation = candidate_centromere
            # Write the candidate centromere with methylation information
            f.write(
                f"{chrom}\tcentromere_detector\tcentromere\t{start}\t{end}\t.\t.\t.\t{attributes};average_methylation={avg_methylation}\n")

        # Write the kinetochore binding site annotations nested inside the centromere
        for pos, methylation in kinetochore_regions:
            f.write(
                f"{chrom}\tcentromere_detector\tkinetochore_binding_site\t{pos}\t{pos + 1}\t.\t.\t.\tmethylation={methylation}\n")


def main():
    parser = argparse.ArgumentParser(description="Identify putative centromere locations.")
    parser.add_argument("gff_file", help="Input GFF file with merged repeat blocks")
    parser.add_argument("bw_file", help="Input BigWig file for methylation data")
    parser.add_argument("output_file", help="Output GFF file with centromere and kinetochore annotations")
    parser.add_argument("--min_length", type=int, default=1000,
                        help="Minimum length of the centromere candidate region")
    args = parser.parse_args()

    # Parse the GFF file to get the repeat blocks
    blocks = parse_gff(args.gff_file)

    # Find the candidate centromeres with the lowest methylation and length filter
    candidate_centromeres = find_centromere(blocks, args.bw_file, args.min_length)

    if candidate_centromeres:
        # Annotate the kinetochore binding sites for each centromere
        kinetochore_regions = []
        for candidate_centromere in candidate_centromeres:
            # Unpack the 5-element tuple: (chrom, start, end, attributes, avg_methylation)
            chrom, start, end, attributes, avg_methylation = candidate_centromere
            kinetochore_regions.extend(
                annotate_kinetochore_binding_site(chrom, start, end, args.bw_file, avg_methylation))

        # Output the result to a GFF file
        output_gff(candidate_centromeres, kinetochore_regions, args.output_file)
        print(f"Centromere candidates and kinetochore binding sites have been annotated in {args.output_file}")
    else:
        print("No Tgut716A blocks found in the GFF file.")


if __name__ == "__main__":
    main()
