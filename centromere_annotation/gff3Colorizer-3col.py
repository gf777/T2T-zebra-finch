#!/usr/bin/env python3
import sys
import random
import re
import hashlib

def hash_color(feature_type):
    """Generate a deterministic color based on the feature type string using a hash."""
    hash_object = hashlib.md5(feature_type.encode())  # Using MD5 to hash the feature type
    hash_value = int(hash_object.hexdigest(), 16)  # Get a large integer from the hash
    return "#{:06x}".format(hash_value & 0xFFFFFF)  # Limit to 6 hex digits for color

if len(sys.argv) != 3:
    print(f"Usage: {sys.argv[0]} input.gff3 output.gff3")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file) as fin, open(output_file, 'w') as fout:
    for line in fin:
        if line.startswith('#') or not line.strip():
            fout.write(line)
            continue

        fields = line.strip().split('\t')
        if len(fields) < 9:
            fout.write(line)
            continue

        feature_type = fields[2]  # The third column (feature type) of the GFF

        # Use hash-based color generation for consistency based on the feature type
        color = hash_color(feature_type)

        color_tag = f"color={color}"
        attributes = fields[8]  # Attribute column
        fields[8] = f"{attributes};{color_tag}"

        fout.write('\t'.join(fields) + '\n')

