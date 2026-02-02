#!/usr/bin/env bash
set -euo pipefail

# Wrapper script to extract gene biotypes from a GTF
# Usage: ./extract_gene_biotype.sh

# Input GTF
GTF="bTaeGut7v0.4_MT_rDNA.EGAPx.v0.1.gtf"

# Derive filenames
GENE_GTF="${GTF%.gtf}.genes.gtf"
BED="${GTF%.gtf}.genes.BIOTYPE.bed"
CODING="${GTF%.gtf}.genes.CODING.bed"

# # Step 1: Subset only 'gene' features
# awk -F$'\t' '$3 == "gene"' "$GTF" > "$GENE_GTF"

# echo "Created subset GTF: $GENE_GTF"

# # Step 2: Extract gene_biotype into BED-like file
# awk -F$'\t' '{
#     match($9, /gene_biotype "([^\"]+)"/, arr)
#     print $1, $4, $5, arr[1]
# }' OFS=$'\t' "$GENE_GTF" > "$BED"

# echo "Generated biotype BED: $BED"

# Step 3: Extract coding genes into a separate BED file
awk -F$'\t' '$4 == "protein_coding"' "$BED" > "$CODING"

echo "Generated coding genes BED: $CODING"