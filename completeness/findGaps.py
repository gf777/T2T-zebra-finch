import sys
from collections import defaultdict

def load_reference_lengths(fai_file):
    lengths = {}
    with open(fai_file) as f:
        for line in f:
            chrom, length, *_ = line.strip().split("\t")
            lengths[chrom] = int(length)
    return lengths

def load_alignments(paf_file):
    aligned_intervals = defaultdict(list)
    with open(paf_file) as f:
        for line in f:
            fields = line.strip().split("\t")
            target = fields[5]
            target_start = int(fields[7])
            target_end = int(fields[8])
            aligned_intervals[target].append((target_start, target_end))
    return aligned_intervals

def merge_intervals(intervals):
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged

def find_gaps(reference_lengths, aligned_intervals):
    gaps = []
    for chrom in reference_lengths:
        intervals = merge_intervals(aligned_intervals.get(chrom, []))
        last_end = 0
        for start, end in intervals:
            if start > last_end:
                gaps.append((chrom, last_end, start))
            last_end = end
        if last_end < reference_lengths[chrom]:
            gaps.append((chrom, last_end, reference_lengths[chrom]))
    return gaps

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python find_gaps.py alignment.paf ref.fai", file=sys.stderr)
        sys.exit(1)

    paf_file = sys.argv[1]
    fai_file = sys.argv[2]

    ref_lengths = load_reference_lengths(fai_file)
    alignments = load_alignments(paf_file)
    gaps = find_gaps(ref_lengths, alignments)

    for chrom, start, end in gaps:
        print(f"{chrom}\t{start}\t{end}")

