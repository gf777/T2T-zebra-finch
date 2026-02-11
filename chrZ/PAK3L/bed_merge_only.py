#!/usr/bin/env python3
import argparse
from typing import List, Tuple, Iterator
import sys

Interval = Tuple[str, int, int]

def read_bed3(path: str) -> Iterator[Interval]:
    with open(path, "r") as f:
        for ln, line in enumerate(f, start=1):
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                parts = line.split()  # allow space-separated
            if len(parts) < 3:
                raise ValueError(f"{path}:{ln}: expected >=3 columns (chrom, start, end)")
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            if start < 0 or end < 0 or end < start:
                raise ValueError(f"{path}:{ln}: invalid interval start={start} end={end}")
            yield (chrom, start, end)

def merge_intervals(intervals: List[Interval], merge_bp: int) -> List[Interval]:
    if not intervals:
        return []
    intervals.sort(key=lambda x: (x[0], x[1], x[2]))
    out: List[Interval] = []
    cur_chrom, cur_s, cur_e = intervals[0]

    for chrom, s, e in intervals[1:]:
        if chrom == cur_chrom and s <= cur_e + merge_bp:
            if e > cur_e:
                cur_e = e
        else:
            out.append((cur_chrom, cur_s, cur_e))
            cur_chrom, cur_s, cur_e = chrom, s, e

    out.append((cur_chrom, cur_s, cur_e))
    return out

def apply_flank(intervals: List[Interval], flank: int) -> List[Interval]:
    if flank <= 0:
        return intervals
    out: List[Interval] = []
    for chrom, s, e in intervals:
        ns = max(0, s - flank)
        ne = e + flank
        out.append((chrom, ns, ne))
    return out

def main():
    ap = argparse.ArgumentParser(description="Merge nearby BED3 intervals (merge-only mode).")
    ap.add_argument("in_bed", help="Input BED (>=3 columns: chrom start end).")
    ap.add_argument("-o", "--out-bed", required=True, help="Output BED3.")
    ap.add_argument("--merge-bp", type=int, default=50000, help="Merge if next.start <= cur.end + merge_bp (default: 50000).")
    ap.add_argument("--flank", type=int, default=0, help="Expand merged intervals by this many bp on both sides (default: 0).")
    args = ap.parse_args()

    intervals = list(read_bed3(args.in_bed))
    merged = merge_intervals(intervals, merge_bp=args.merge_bp)
    merged = apply_flank(merged, flank=args.flank)

    with open(args.out_bed, "w") as out:
        for chrom, s, e in merged:
            out.write(f"{chrom}\t{s}\t{e}\n")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(2)

