#!/bin/bash

# Default flanking size
flank_size=0

# Parse command-line arguments
while getopts "f:" opt; do
  case $opt in
    f) flank_size=$OPTARG ;;
    *) echo "Usage: $0 [-f flank_size] input.gff"; exit 1 ;;
  esac
done
shift $((OPTIND - 1))

# Check for input file
if [ $# -ne 1 ]; then
  echo "Usage: $0 [-f flank_size] input.gff"
  exit 1
fi
input_gff=$1

# Process GFF file with added flanking sequence
awk -v flank="$flank_size" 'BEGIN {OFS="\t"} 
  /^[^#]/ {
    start = $4 - flank;
    if (start < 1) start = 1;  # Ensure start is at least 1
    end = $5 + flank;
    $4 = start;
    $5 = end;
    print $0;
  }' "$input_gff"

