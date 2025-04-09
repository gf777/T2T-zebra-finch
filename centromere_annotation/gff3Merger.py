#!/usr/bin/env python3
import sys

def parse_target(attributes):
    """Extract the first word of the Target attribute and the third numerical value."""
    for attr in attributes.split(';'):
        if attr.startswith("Target="):
            target_parts = attr.split('=')[1].split()
            if len(target_parts) >= 3:
                return target_parts[0], int(target_parts[2])  # Extract first word and length
    return None, 0

def update_target_length(attributes, new_length):
    """Update the third column of the Target attribute to reflect the merged feature length."""
    parts = attributes.split(';')
    for i, attr in enumerate(parts):
        if attr.startswith("Target="):
            target_parts = attr.split('=')
            target_values = target_parts[1].split()
            if len(target_values) >= 3:
                target_values[2] = str(new_length)  # Update length as a string
                parts[i] = f"Target={' '.join(target_values)}"
                break
    return ";".join(parts)

if len(sys.argv) != 4:
    print(f"Usage: {sys.argv[0]} input.gff3 output.gff3 distance")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]
distance = int(sys.argv[3])

merged_entries = []
with open(input_file) as fin:
    header = []
    prev_entry = None

    for line in fin:
        if line.startswith("#"):
            header.append(line)
            continue

        fields = line.strip().split("\t")
        if len(fields) < 9:
            continue

        chrom, source, feature, start, end, score, strand, phase, attributes = fields
        start, end = int(start), int(end)
        target_id, target_length = parse_target(attributes)

        if prev_entry:
            (prev_chrom, prev_source, prev_feature, prev_start, prev_end, prev_score, 
             prev_strand, prev_phase, prev_attributes, prev_target_id, prev_target_length) = prev_entry

            if (chrom == prev_chrom and target_id == prev_target_id
                    and start < int(prev_end) + distance):  # Ensure they're consecutive

                # Extend the feature and accumulate length
                prev_entry[4] = str(end)  # Ensure end remains string
                prev_entry[10] = int(prev_entry[10]) + target_length  # Accumulate target length as int
                prev_entry[8] = update_target_length(prev_attributes, prev_entry[10])  # Update target length in attributes
                continue  # Continue merging

            # If not mergeable, save the previous entry and move to the next
            merged_entries.append("\t".join(map(str, prev_entry[:-1])))

        prev_entry = [chrom, source, feature, start, end, score, strand, phase, attributes, target_id, target_length]

    if prev_entry:
        merged_entries.append("\t".join(map(str, prev_entry[:-1])))

with open(output_file, "w") as fout:
    fout.writelines(header)
    fout.writelines("\n".join(merged_entries) + "\n")

print(f"Processed file saved to {output_file}")

