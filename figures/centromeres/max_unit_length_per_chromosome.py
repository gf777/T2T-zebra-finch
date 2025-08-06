#!/usr/bin/env python3

import sys
from collections import defaultdict

def parse_gff_compute_lengths(file_path):
    max_len = defaultdict(int)

    with open(file_path) as f:
        for line in f:
            if line.strip() == "" or line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue  # skip malformed lines
            chrom = fields[0]
            try:
                start = int(fields[3])
                end = int(fields[4])
                unit_len = end - start + 1
                if unit_len > max_len[chrom]:
                    max_len[chrom] = unit_len
            except ValueError:
                continue  # skip lines with bad coordinates

    for chrom in sorted(max_len):
        print(f"{chrom}\t{max_len[chrom]}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python max_unit_len_per_chrom.py input.gff")
        sys.exit(1)

    parse_gff_compute_lengths(sys.argv[1])