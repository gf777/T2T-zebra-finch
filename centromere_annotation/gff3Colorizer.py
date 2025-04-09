#!/usr/bin/env python3
import sys
import random
import re
import hashlib

def hash_color(target):
    """Generate a deterministic color based on the target string using a hash."""
    hash_object = hashlib.md5(target.encode())  # Using MD5 to hash the target
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

        attributes = fields[8]
        # Search for Target attribute (format: Target=target_value)
        match = re.search(r'Target=([^;\s]+)', attributes)
        target = match.group(1) if match else 'NO_TARGET'

        # Use hash-based color generation for consistency
        color = hash_color(target)

        color_tag = f"color={color}"
        fields[8] = f"{attributes};{color_tag}"

        fout.write('\t'.join(fields) + '\n')

