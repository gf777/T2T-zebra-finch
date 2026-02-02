#!/usr/bin/env bash
set -euo pipefail

# generate_coding_cov.sh
# ---------------------
# For each ${PREFIX}.${RES}.bins.tsv in the working directory,
# computes coding coverage per bin (fraction of bin covered by coding genes)
# and writes ${PREFIX}.${RES}.coding_cov.bedgraph with columns:
#   chrom  start  end  weight  coding
#
# Prerequisites:
#  - bedtools must be installed and in PATH
#  - coding genes BED: bTaeGut7v0.4_MT_rDNA.EGAPx.v0.1.genes.CODING.bed

# Path to the coding-genes BED file (contains all chromosomes)
CODING_BED="bTaeGut7v0.4_MT_rDNA.EGAPx.v0.1.genes.CODING.bed"

# Check that the coding BED exists
if [[ ! -f "$CODING_BED" ]]; then
  echo "Error: coding genes BED file not found: $CODING_BED" >&2
  exit 1
fi

# Process each bins TSV
for TSV in *.bins.tsv; do
  # Base name (PREFIX.RES)
  BASE="${TSV%.bins.tsv}"
  OUT="${BASE}.coding_cov.bedgraph"
  TMP_BED="${BASE}.bins.bed"
  TMP_COV="${BASE}.cov.tmp"
  TMP_FRAC="${BASE}.frac.tmp"
  TMP_DATA="${BASE}.data.tmp"

  echo "Processing $TSV → $OUT"

  # Convert bins TSV (chrom, start, end, weight[,GC]) to BED (first 3 cols)
  cut -f1-3 "$TSV" > "$TMP_BED"

  # Compute coverage; bedtools outputs 7 columns with fraction in $7
  bedtools coverage -a "$TMP_BED" -b "$CODING_BED" > "$TMP_COV"

  # Extract coding fraction (column 7)
  awk '{ print $7 }' "$TMP_COV" > "$TMP_FRAC"

  # Extract original columns 1-4 (skip header)
  tail -n +2 "$TSV" | cut -f1-4 > "$TMP_DATA"

  # Write header, updating GC → coding
  echo -e "chrom	start	end	weight	coding" > "$OUT"

  # Combine original data with coding fraction
  paste "$TMP_DATA" "$TMP_FRAC" >> "$OUT"

  # Clean up temporary files
  rm "$TMP_BED" "$TMP_COV" "$TMP_FRAC" "$TMP_DATA"

done

echo "All coding coverage bedgraphs generated with preserved structure."
