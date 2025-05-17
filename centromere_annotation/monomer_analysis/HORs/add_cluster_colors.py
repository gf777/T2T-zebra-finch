#!/usr/bin/env python3
import sys
import random

import re

def parse_uc(uc_file):
    """Map each sequence ID to its cluster number, parsing ID=...; from full headers."""
    id_to_cluster = {}
    with open(uc_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if parts[0] in ('S', 'H'):  # seed or hit
                cluster = int(parts[1])
                full_id = parts[8]  # e.g. ID=246437;Target=...
                m = re.search(r'ID=([^\s;]+)', full_id)
                if m:
                    seq_id = m.group(1)
                    id_to_cluster[seq_id] = cluster
    return id_to_cluster

def random_color(seed=None):
    """Generate a truly random hex color (avoids pastel bias)."""
    return "#{:06x}".format(random.randint(0, 0xFFFFFF))

def extract_id(header):
    m = re.search(r'ID=([^\s;]+)', header)
    if m:
        return m.group(1)
    return None

def main(gff_file, uc_file, fasta_file, output_file):
    id_to_cluster = parse_uc(uc_file)

    # Map FASTA headers to cluster IDs
    fasta_id_to_cluster = {}
    with open(fasta_file) as f:
        for line in f:
            if line.startswith('>'):
                seq_id = extract_id(line[1:])
                if seq_id in id_to_cluster:
                    fasta_id_to_cluster[seq_id] = id_to_cluster[seq_id]

    cluster_to_color = {cl: random_color(cl) for cl in set(id_to_cluster.values())}

    with open(gff_file) as fin, open(output_file, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                fout.write(line)
                continue
            parts = line.strip().split('\t')
            attr_field = parts[8]
            id_val = None
            for field in attr_field.split(';'):
                if field.startswith('ID='):
                    id_val = field[3:]
                    break
            if id_val and id_val in fasta_id_to_cluster:
                cluster = fasta_id_to_cluster[id_val]
                color = cluster_to_color[cluster]
                attr_field += f";color={color}"
                parts[8] = attr_field
            fout.write('\t'.join(parts) + '\n')

if __name__ == '__main__':
    if len(sys.argv) != 5:
        print("Usage: python add_cluster_colors.py input.gff clusters.uc clustered.fa output.gff")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

