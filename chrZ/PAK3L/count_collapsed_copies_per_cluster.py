#!/usr/bin/env python3
import argparse
from dataclasses import dataclass
from typing import Dict, List, Tuple
import sys

@dataclass
class Cluster:
    chrom: str
    start: int
    end: int
    idx: int

@dataclass
class Locus:
    chrom: str
    start: int
    end: int
    gene: str

def read_bed3(path: str) -> List[Cluster]:
    clusters: List[Cluster] = []
    with open(path) as f:
        for ln, line in enumerate(f, start=1):
            line = line.strip()
            if not line or line.startswith(("#","track","browser")):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                parts = line.split()
            if len(parts) < 3:
                raise ValueError(f"{path}:{ln}: expected BED3")
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            if end < start:
                raise ValueError(f"{path}:{ln}: invalid interval {start}-{end}")
            clusters.append(Cluster(chrom, start, end, len(clusters)+1))
    return clusters

def read_collapsed_bed(path: str) -> List[Locus]:
    loci: List[Locus] = []
    with open(path) as f:
        for ln, line in enumerate(f, start=1):
            line = line.strip()
            if not line or line.startswith(("#","track","browser")):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                parts = line.split()
            if len(parts) < 4:
                raise ValueError(f"{path}:{ln}: expected BED4+ (chrom start end name)")
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            gene = parts[3]
            loci.append(Locus(chrom, start, end, gene))
    return loci

def overlap_bp(a0:int,a1:int,b0:int,b1:int) -> int:
    s = max(a0,b0)
    e = min(a1,b1)
    return max(0, e-s)

def main():
    ap = argparse.ArgumentParser(description="Count collapsed miniprot loci (copies) per gene within each cluster BED interval.")
    ap.add_argument("--clusters-bed", required=True, help="BED3 of clusters")
    ap.add_argument("--collapsed-bed", required=True, help="BED from collapse_miniprot_isoforms.py (gene in 4th column)")
    ap.add_argument("-o", "--out", required=True, help="Output TSV (for plotting)")
    ap.add_argument("--min-overlap-bp", type=int, default=1, help="Minimum bp overlap to count (default 1)")
    args = ap.parse_args()

    clusters = read_bed3(args.clusters_bed)
    loci = read_collapsed_bed(args.collapsed_bed)

    # index loci by chrom, sorted
    by_chr: Dict[str, List[Locus]] = {}
    for L in loci:
        by_chr.setdefault(L.chrom, []).append(L)
    for chrom in by_chr:
        by_chr[chrom].sort(key=lambda x: x.start)

    with open(args.out, "w") as out:
        out.write("\t".join(["cluster_id","cluster_chrom","cluster_start","cluster_end","gene","count"]) + "\n")

        for c in clusters:
            lst = by_chr.get(c.chrom, [])
            counts: Dict[str, int] = {}

            # scan overlaps (clusters are few; loci likely not huge)
            for L in lst:
                if L.end <= c.start:
                    continue
                if L.start >= c.end:
                    break
                if overlap_bp(L.start, L.end, c.start, c.end) >= args.min_overlap_bp:
                    counts[L.gene] = counts.get(L.gene, 0) + 1

            if not counts:
                out.write(f"{c.idx}\t{c.chrom}\t{c.start}\t{c.end}\t.\t0\n")
            else:
                for gene in sorted(counts.keys()):
                    out.write(f"{c.idx}\t{c.chrom}\t{c.start}\t{c.end}\t{gene}\t{counts[gene]}\n")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(2)

