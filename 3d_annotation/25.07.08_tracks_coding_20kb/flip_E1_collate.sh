#!/usr/bin/env bash
set -euo pipefail

# flip_and_collate_tracks.sh
# -------------------------
# Flips compartment E1 bigWigs (scale -1), converts to bedGraph + bigWig,
# and collates maternal/paternal into a diploid track.
# Additionally creates an A/B compartment BED track (A: E1>0, B: E1<0).
#
# Loops over HAPLOTYPES and RESOLUTIONS.
# For HAP = dip, mat+Z, pat+WZ:
#   - scale -1 the cis.bw → .flipped.bg
#   - convert .flipped.bg → .flipped.bw
# For mat+Z & pat+WZ:
#   - remove sex-chromosome duplicates
#   - concatenate filtered .flipped.bg → dip.collated.bg
#   - convert dip.collated.bg → dip.collated.bw

# 1) Define arrays and constants
HAPS=("dip" "mat+Z" "pat+WZ")
RESOLUTIONS=(20000)        # list of resolutions to process
DATE="20250313"          # date string in filenames
METHOD="bwa"              # aligner suffix in prefixes
BASE_PREFIX="bTaeGut7"    # common prefix
COLL_PREFIX="bTaeGut7v0.4_MT_rDNA.Cooltools.v0.2.E1.20Kb.flipped.dip.collated"
COLL_BG="${COLL_PREFIX}.bg"
COLL_BW="${COLL_PREFIX}.bw"

# 2) Flip each haplotype's compartments track
for HAP in "${HAPS[@]}"; do
  ASM="${BASE_PREFIX}.${HAP}.cur.${DATE}"
  SIZES_FILE="${ASM}.chrom.sizes"

  for RES in "${RESOLUTIONS[@]}"; do
    # Paths for input and outputs
    IN_BW="${ASM}.${METHOD}.${RES}.pairs.sorted.dedup.cool.compartments.cis.bw"
    OUT_BG="${ASM}.${METHOD}.${RES}.pairs.sorted.dedup.cool.compartments.cis.flipped.bg"
    OUT_BW="${ASM}.${METHOD}.${RES}.pairs.sorted.dedup.cool.compartments.cis.flipped.bw"

    echo "[${HAP} @ ${RES}] Inverting E1 sign → ${OUT_BG} + ${OUT_BW}"
    wiggletools write_bg "$OUT_BG" scale -1 "$IN_BW"
    bedGraphToBigWig "$OUT_BG" "$SIZES_FILE" "$OUT_BW"
  done
done

# 3) Collate maternal & paternal into diploid track, filtering sex chromosomes
for RES in "${RESOLUTIONS[@]}"; do
  PAT_BG_RAW="${BASE_PREFIX}.pat+WZ.cur.${DATE}.${METHOD}.${RES}.pairs.sorted.dedup.cool.compartments.cis.flipped.bg"
  MAT_BG_RAW="${BASE_PREFIX}.mat+Z.cur.${DATE}.${METHOD}.${RES}.pairs.sorted.dedup.cool.compartments.cis.flipped.bg"
  PAT_BG_FILTERED="${BASE_PREFIX}.pat.cur.${DATE}.${METHOD}.${RES}.filtered.bg"
  MAT_BG_FILTERED="${BASE_PREFIX}.mat.cur.${DATE}.${METHOD}.${RES}.filtered.bg"
  DIP_SIZES="${BASE_PREFIX}.dip.cur.${DATE}.chrom.sizes"

  echo "[collate @ ${RES}] Removing sex chromosome duplicates"
  grep -v "chrW_mat" "$PAT_BG_RAW" > "$PAT_BG_FILTERED"
  grep -v "chrZ_pat" "$MAT_BG_RAW" > "$MAT_BG_FILTERED"

  echo "[collate @ ${RES}] Merging pat+WZ + mat+Z → ${COLL_BG}"
  cat "$PAT_BG_FILTERED" "$MAT_BG_FILTERED" > "$COLL_BG"

  echo "[collate @ ${RES}] Converting collated bedGraph → bigWig"
  bedGraphToBigWig "$COLL_BG" "$DIP_SIZES" "$COLL_BW"

done

# 4) Consolidate A/B compartments and merge adjacent bins
AB_BED="${COLL_PREFIX}.AB.bed"
echo "Generating consolidated A/B BED → $AB_BED"

awk 'BEGIN{
  OFS="\t";
  last_chr=""; last_start=0; last_end=0; last_comp=""
}
# Skip NaN lines and treat >0 as A, <0 as B
$4!="nan" {
  comp=($4>0)?"A":"B";
  # If same chromosome, same compartment, contiguous, extend
  if ($1==last_chr && comp==last_comp && $2==last_end) {
    last_end=$3;
  } else {
    # flush previous block
    if (last_comp!="") {
      color=(last_comp=="A"?"#E38AAA":"#6F9DD0");
      print last_chr, last_start, last_end, last_comp, 0, ".", last_start, last_end, color;
    }
    # start new block
    last_chr=$1; last_start=$2; last_end=$3; last_comp=comp;
  }
}
END {
  # flush final block
  if (last_comp!="") {
    color=(last_comp=="A"?"#E38AAA":"#6F9DD0");
    print last_chr, last_start, last_end, last_comp, 0, ".", last_start, last_end, color;
  }
}' "$COLL_BG" > "$AB_BED"

echo "All tracks flipped, collated, and A/B BED generated."
