#!/usr/bin/env python3
import argparse
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
import sys

@dataclass
class Mrna:
    seqid: str
    start: int   # 0-based inclusive
    end: int     # 0-based exclusive
    strand: str
    score: float
    mrna_id: str
    gene: str
    rank: Optional[int]
    identity: Optional[float]
    target_span: Optional[Tuple[int,int]]  # query aa coords if present

def parse_attrs(attr_str: str) -> Dict[str, str]:
    d: Dict[str, str] = {}
    for part in attr_str.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = v
        else:
            d[part] = ""
    return d

def parse_target(target: str) -> Tuple[str, Optional[Tuple[int,int]]]:
    """
    miniprot Target: "<query> <qstart> <qend>"
    """
    if not target:
        return ("NA", None)
    toks = target.strip().split()
    gene = toks[0]
    span = None
    if len(toks) >= 3:
        try:
            span = (int(toks[1]), int(toks[2]))
        except ValueError:
            span = None
    return (gene, span)

def normalize_gene_id(g: str) -> str:
    # PAK3L-3_finch_ENST... -> PAK3L-3
    g = g.strip()
    return g.split("_", 1)[0] if "_" in g else g

def read_mrnas(gff: str) -> List[Mrna]:
    out: List[Mrna] = []
    with open(gff) as f:
        for ln, line in enumerate(f, start=1):
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            seqid, source, ftype, s, e, score_s, strand, phase, attrs_s = parts
            if ftype != "mRNA":
                continue
            try:
                s1 = int(s)   # 1-based inclusive
                e1 = int(e)
            except ValueError:
                continue
            start0 = s1 - 1
            end0 = e1

            try:
                score = float(score_s) if score_s != "." else 0.0
            except ValueError:
                score = 0.0

            attrs = parse_attrs(attrs_s)
            mrna_id = attrs.get("ID", f"mRNA_line_{ln}")

            gene, qspan = parse_target(attrs.get("Target", ""))
            gene = normalize_gene_id(gene)

            rank = None
            if "Rank" in attrs:
                try:
                    rank = int(attrs["Rank"])
                except ValueError:
                    pass

            identity = None
            if "Identity" in attrs:
                try:
                    identity = float(attrs["Identity"])
                except ValueError:
                    pass

            out.append(Mrna(
                seqid=seqid, start=start0, end=end0, strand=strand,
                score=score, mrna_id=mrna_id, gene=gene,
                rank=rank, identity=identity, target_span=qspan
            ))
    return out

def merge_loci(mrnas: List[Mrna], merge_bp: int) -> List[List[Mrna]]:
    """
    Cluster mRNAs into loci by coordinate overlap/proximity on the same seqid.
    """
    if not mrnas:
        return []
    mrnas.sort(key=lambda m: (m.seqid, m.start, m.end))

    loci: List[List[Mrna]] = []
    cur: List[Mrna] = [mrnas[0]]
    cur_seq = mrnas[0].seqid
    cur_s = mrnas[0].start
    cur_e = mrnas[0].end

    for m in mrnas[1:]:
        if m.seqid != cur_seq:
            loci.append(cur)
            cur = [m]
            cur_seq = m.seqid
            cur_s, cur_e = m.start, m.end
            continue

        if m.start <= cur_e + merge_bp:
            cur.append(m)
            if m.end > cur_e:
                cur_e = m.end
        else:
            loci.append(cur)
            cur = [m]
            cur_s, cur_e = m.start, m.end

    loci.append(cur)
    return loci

def pick_best(locus: List[Mrna]) -> Mrna:
    """
    Best representative = best gene assignment for the locus.
    Ranking:
      1) score (desc)
      2) identity (desc)
      3) rank (asc; missing rank treated as very large)
      4) longer genomic span (desc)
    """
    def key(m: Mrna):
        ident = m.identity if m.identity is not None else -1.0
        r = m.rank if m.rank is not None else 10**9
        span = m.end - m.start
        return (-m.score, -ident, r, -span, m.start, m.end)
    return sorted(locus, key=key)[0]

def main():
    ap = argparse.ArgumentParser(description="Collapse miniprot 'isoforms' to one best gene assignment per locus.")
    ap.add_argument("--gff", required=True, help="miniprot GFF3 output")
    ap.add_argument("--merge-bp", type=int, default=2000,
                    help="Merge mRNAs into one locus if within this distance (default 2000)")
    ap.add_argument("--out-prefix", required=True, help="Output prefix")
    args = ap.parse_args()

    mrnas = read_mrnas(args.gff)
    if not mrnas:
        raise SystemExit("ERROR: no mRNA records found in GFF")

    loci = merge_loci(mrnas, merge_bp=args.merge_bp)

    # output TSV + BED
    tsv_path = f"{args.out_prefix}.collapsed.tsv"
    bed_path = f"{args.out_prefix}.collapsed.bed"

    with open(tsv_path, "w") as tsv, open(bed_path, "w") as bed:
        tsv.write("\t".join([
            "locus_id","seqid","start","end","strand",
            "best_gene","best_mrna_id","best_score","best_identity","best_rank",
            "n_isoforms","all_genes","all_mrna_ids"
        ]) + "\n")

        for i, locus in enumerate(loci, start=1):
            best = pick_best(locus)
            locus_seq = best.seqid
            ls = min(m.start for m in locus)
            le = max(m.end for m in locus)
            strand = best.strand

            all_genes = sorted({m.gene for m in locus})
            all_ids = [m.mrna_id for m in sorted(locus, key=lambda x: (-x.score, x.start, x.end))]

            locus_id = f"LOCUS{i:06d}"
            # BED for IGV
            bed.write(f"{locus_seq}\t{ls}\t{le}\t{best.gene}\t0\t{strand}\n")

            tsv.write("\t".join(map(str, [
                locus_id, locus_seq, ls, le, strand,
                best.gene, best.mrna_id, best.score,
                best.identity if best.identity is not None else ".",
                best.rank if best.rank is not None else ".",
                len(locus),
                ",".join(all_genes),
                ",".join(all_ids)
            ])) + "\n")

    print(f"Read mRNAs: {len(mrnas)}")
    print(f"Loci after collapsing: {len(loci)} (merge_bp={args.merge_bp})")
    print(f"Wrote: {tsv_path}")
    print(f"Wrote: {bed_path}")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(2)

