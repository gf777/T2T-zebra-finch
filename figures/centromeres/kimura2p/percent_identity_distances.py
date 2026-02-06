#!/usr/bin/env python3
"""
percent_identity_distances.py

Compute per-sequence identities and Kimura 2-parameter distances for:
  (1) global consensus vs globally aligned sequences,
  (2) Takki (aligned first row) vs globally aligned sequences,
  (3) per-chromosome consensus vs per-chromosome aligned sequences.

Two identity flavors are reported:
  - *_maskN  : excludes columns where the reference is ambiguous (N),
               penalizes base-vs-gap as mismatch, ignores gap-gap.
  - *_countN : same as above but ALSO penalizes base-vs-N as mismatch
               (still ignores N-vs-gap and gap-gap).

For Kimura 2P, we use:
  P = (# transition differences) / (total positions scored)
  Q = (# transversion differences) / (total positions scored)
  d = -1/2 * ln(1 - 2P - Q) - 1/4 * ln(1 - 2Q)
where the denominator "total positions scored" follows the same inclusion
rules as the corresponding identity flavor.

Input tree (per family = Tgut716A, Tgut191A):
  {wd}/{fam}_asm_alignments/{fam}_aln.fasta        (global MSA; Takki is first row)
  {wd}/{fam}_asm_alignments/{fam}_consensus.fa     (aligned-length consensus; A/C/G/T/N)
  {wd}/{fam}_chr_alignments/{chr}/{fam}_{chr}_aln.fa
  {wd}/{fam}_chr_alignments/{chr}/{fam}_{chr}_consensus.fa

Annotation:
  {wd}/inputs/PAR.bed
  {wd}/inputs/bTaeGut7v0.4_MT_rDNA.centromere_detector.v0.1.gff

Output:
  {wd}/{fam}_pid_kimura.tsv

Notes:
- Query sequences are those whose headers contain coordinates after '::',
  e.g., 'dispersed_repeat::chrX:123-456' (old) or 'Tgut191A::chrX:123-456|strand=+' (new).
- Consensus strings contain A/C/G/T/N (no '-'); Takki row may contain '-'.
- All strings are uppercased; comparisons are column-wise in the MSA frame.

"""

import os
import math
import csv
from collections import defaultdict
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from Bio import AlignIO, SeqIO

# ───────────────────────────────────────────── Config
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
wd = os.path.normpath(os.path.join(SCRIPT_DIR, ".."))
families = ["Tgut716A", "Tgut191A"]

PAR_BED = os.path.join(wd, "inputs", "PAR.bed")
CEN_GFF = os.path.join(wd, "inputs", "bTaeGut7v0.4_MT_rDNA.centromere_detector.v0.1.gff")

# ───────────────────────────────────────────── Helpers: IO

def read_alignment_dict(path: str) -> Dict[str, str]:
    """Read a FASTA alignment into {id: aligned_seq (upper)}."""
    aln = AlignIO.read(path, "fasta")
    return {rec.id: str(rec.seq).upper() for rec in aln}

def read_alignment_ordered(path: str) -> List[Tuple[str, str]]:
    """Read FASTA alignment into a list of (id, seq) preserving order, uppercased."""
    aln = AlignIO.read(path, "fasta")
    return [(rec.id, str(rec.seq).upper()) for rec in aln]

def read_single_fasta_seq(path: str) -> Tuple[str, str]:
    """Read a single-record FASTA, return (id, seq upper)."""
    rec = next(SeqIO.parse(path, "fasta"))
    return rec.id, str(rec.seq).upper()

# ───────────────────────────────────────────── Helpers: parsing coordinates

def parse_seqname(header: str) -> Tuple[str, int, int]:
    """
    Parse headers into (chr, start, end).

    Supports:
      - 'dispersed_repeat::chr1_pat:19669813-19670024'
      - 'Tgut191A::chr1_pat:23262-23438|strand=+'

    Strategy:
      - Split once on '::', take the RHS.
      - If '|' is present (e.g., '|strand=+'), strip it off.
      - Then split 'chr:start-end'.
    """
    try:
        body = header.split("::", 1)[1]
        if "|" in body:
            body = body.split("|", 1)[0]
        chrom, coord = body.split(":")
        start_s, end_s = coord.split("-")
        return chrom, int(start_s), int(end_s)
    except Exception as e:
        raise ValueError(f"Unrecognized seqname format: {header}") from e

def is_query_header(header: str) -> bool:
    """Return True iff the header parses as a coordinate-bearing instance."""
    try:
        parse_seqname(header)
        return True
    except Exception:
        return False

# ───────────────────────────────────────────── Helpers: interval annotation

def load_bed_intervals(path: str) -> Dict[str, List[Tuple[int, int]]]:
    """
    Load BED (chr, start, end, ...). Returns dict: chr -> [(start, end), ...].
    Treats intervals as half-open [start, end); overlap check uses (a_start < b_end and a_end > b_start).
    """
    iv = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom, s, e = parts[0], int(parts[1]), int(parts[2])
            iv[chrom].append((s, e))
    return iv

def load_gff_intervals(path: str) -> Dict[str, List[Tuple[int, int]]]:
    """
    Load GFF3 and return dict: chr -> [(start0, end0), ...] in 0-based half-open.
    (GFF is 1-based inclusive; convert to 0-based half-open by start-1, end).
    Include all feature rows (non-comment).
    """
    iv = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            chrom = parts[0]
            start1 = int(parts[3])  # 1-based inclusive
            end1 = int(parts[4])    # 1-based inclusive
            start0 = start1 - 1
            end0 = end1
            if start0 < end0:
                iv[chrom].append((start0, end0))
    return iv

def overlaps(anylist: List[Tuple[int, int]], start: int, end: int) -> bool:
    """Return True if [start,end) overlaps any interval in anylist using half-open overlap test."""
    if not anylist:
        return False
    for s, e in anylist:
        if start < e and end > s:
            return True
    return False

# ───────────────────────────────────────────── Helpers: scoring

ACGT = set("ACGT")

def is_base(x: str) -> bool:
    return x in ACGT

def is_gap(x: str) -> bool:
    return x == "-"

def is_amb(x: str) -> bool:
    # Ambiguous (e.g., N) but not gap
    return (not is_base(x)) and (not is_gap(x))

# transitions: A<->G, C<->T
TRANSITION_PAIRS = {("A","G"), ("G","A"), ("C","T"), ("T","C")}

def is_transition(a: str, b: str) -> bool:
    return (a, b) in TRANSITION_PAIRS

def column_included_maskN(ref: str, qry: str) -> bool:
    """
    Include in denominator for PID_maskN / P / Q if:
      - not gap-gap,
      - reference is not ambiguous (N etc).
    Gap vs base IS included (penalized as mismatch).
    """
    if is_gap(ref) and is_gap(qry):
        return False
    if is_amb(ref):
        return False
    return True

def column_included_countN(ref: str, qry: str) -> bool:
    """
    Include in denominator for PID_countN if:
      - not gap-gap,
      - not (ref is ambiguous AND qry is gap)  [N vs '-' excluded],
    Otherwise include (so base vs N is counted as mismatch).
    """
    if is_gap(ref) and is_gap(qry):
        return False
    if is_amb(ref) and is_gap(qry):
        return False
    return True

def score_pairwise(ref_aln: str, qry_aln: str, flavor: str) -> Tuple[float, float, float]:
    """
    Compute (PID, P, Q) for aligned reference vs aligned query.
    flavor: "maskN" or "countN" (see rules above).
    - PID = matches / denom * 100
    - P   = transitions / denom
    - Q   = transversions / denom
    where denom is the count of included columns for the chosen flavor.
    """
    assert len(ref_aln) == len(qry_aln), "Aligned strings must be same length"
    include = column_included_maskN if flavor == "maskN" else column_included_countN

    denom = 0
    matches = 0
    trans = 0
    tvans = 0

    for r, q in zip(ref_aln, qry_aln):
        if not include(r, q):
            continue
        denom += 1

        if is_base(r) and is_base(q):
            if r == q:
                matches += 1
            else:
                if is_transition(r, q):
                    trans += 1
                else:
                    tvans += 1
        else:
            # Any case with gap vs base → mismatch counted in PID (no P/Q increment).
            # Base vs N (countN flavor) also falls here (mismatch, no P/Q).
            # N never appears in query by our assumptions.
            pass

    if denom == 0:
        return (float("nan"), float("nan"), float("nan"))

    pid = (matches / denom) * 100.0
    P = trans / denom
    Q = tvans / denom
    return (pid, P, Q)

def kimura2p(P: float, Q: float) -> float:
    """
    Kimura 2-parameter distance:
      d = -1/2 ln(1 - 2P - Q) - 1/4 ln(1 - 2Q)
    Returns NaN if arguments to ln are non-positive or P/Q are NaN.
    """
    if any(map(lambda x: x is None or (isinstance(x, float) and math.isnan(x)), (P, Q))):
        return float("nan")
    a = 1.0 - 2.0*P - Q
    b = 1.0 - 2.0*Q
    if a <= 0.0 or b <= 0.0:
        return float("nan")
    return -0.5 * math.log(a) - 0.25 * math.log(b)

# ───────────────────────────────────────────── Core computation

def compute_family_table(fam: str,
                         par_iv: Dict[str, List[Tuple[int,int]]],
                         cen_iv: Dict[str, List[Tuple[int,int]]]) -> pd.DataFrame:
    """
    Build the per-sequence table for one family.
    """
    # Paths
    global_dir = os.path.join(wd, f"{fam}_asm_alignments")
    global_aln_path = os.path.join(global_dir, f"{fam}_aln.fasta")
    global_cons_path = os.path.join(global_dir, f"{fam}_consensus.fa")

    # Read global alignment (ordered to find first row = Takki)
    global_rows = read_alignment_ordered(global_aln_path)
    if not global_rows:
        raise RuntimeError(f"No rows in alignment: {global_aln_path}")

    # Reference strings
    takki_id, takki_aln = global_rows[0]  # first row
    _, global_cons_aln = read_single_fasta_seq(global_cons_path)
    # sanity
    L = len(global_rows[0][1])
    if len(global_cons_aln) != L:
        raise RuntimeError(f"Global consensus length != alignment length for {fam}")

    # Collect query rows (instances): any header we can parse as coordinates
    queries_global = [(sid, s) for sid, s in global_rows if is_query_header(sid)]

    # Initialize records per seqname
    recs = {}
    def ensure_rec(seqname: str):
        if seqname not in recs:
            recs[seqname] = {
                "seqname": seqname,
                "chr": "",
                "repeat_length": np.nan,
                "PAR": False,
                "CEN": False,
                # Pre-fill all metrics with NaN
                "global_PID_maskN": np.nan, "global_P_maskN": np.nan, "global_Q_maskN": np.nan, "global_kimura2p_maskN": np.nan,
                "global_PID_countN": np.nan, "global_P_countN": np.nan, "global_Q_countN": np.nan, "global_kimura2p_countN": np.nan,
                "takki_PID_maskN": np.nan, "takki_P_maskN": np.nan, "takki_Q_maskN": np.nan, "takki_kimura2p_maskN": np.nan,
                "takki_PID_countN": np.nan, "takki_P_countN": np.nan, "takki_Q_countN": np.nan, "takki_kimura2p_countN": np.nan,
                "chr_PID_maskN": np.nan, "chr_P_maskN": np.nan, "chr_Q_maskN": np.nan, "chr_kimura2p_maskN": np.nan,
                "chr_PID_countN": np.nan, "chr_P_countN": np.nan, "chr_Q_countN": np.nan, "chr_kimura2p_countN": np.nan,
            }

    # --- 2.1 global consensus vs sequences
    for sid, seq in queries_global:
        ensure_rec(sid)
        pid, P, Q = score_pairwise(global_cons_aln, seq, "maskN")
        recs[sid]["global_PID_maskN"] = pid
        recs[sid]["global_P_maskN"] = P
        recs[sid]["global_Q_maskN"] = Q
        recs[sid]["global_kimura2p_maskN"] = kimura2p(P, Q)

        pid, P, Q = score_pairwise(global_cons_aln, seq, "countN")
        recs[sid]["global_PID_countN"] = pid
        recs[sid]["global_P_countN"] = P
        recs[sid]["global_Q_countN"] = Q
        recs[sid]["global_kimura2p_countN"] = kimura2p(P, Q)

    # --- 2.2 Takki vs sequences (Takki is first aligned row)
    for sid, seq in queries_global:
        ensure_rec(sid)
        pid, P, Q = score_pairwise(takki_aln, seq, "maskN")
        recs[sid]["takki_PID_maskN"] = pid
        recs[sid]["takki_P_maskN"] = P
        recs[sid]["takki_Q_maskN"] = Q
        recs[sid]["takki_kimura2p_maskN"] = kimura2p(P, Q)

        pid, P, Q = score_pairwise(takki_aln, seq, "countN")
        recs[sid]["takki_PID_countN"] = pid
        recs[sid]["takki_P_countN"] = P
        recs[sid]["takki_Q_countN"] = Q
        recs[sid]["takki_kimura2p_countN"] = kimura2p(P, Q)

    # --- 2.3 per-chromosome consensus vs per-chrom aligned sequences
    chr_root = os.path.join(wd, f"{fam}_chr_alignments")
    if os.path.isdir(chr_root):
        for chr_name in sorted(os.listdir(chr_root)):
            chr_dir = os.path.join(chr_root, chr_name)
            if not os.path.isdir(chr_dir):
                continue
            aln_path = os.path.join(chr_dir, f"{fam}_{chr_name}_aln.fa")
            cons_path = os.path.join(chr_dir, f"{fam}_{chr_name}_consensus.fa")
            if not (os.path.isfile(aln_path) and os.path.isfile(cons_path)):
                continue

            rows = read_alignment_ordered(aln_path)
            _, chr_cons_aln = read_single_fasta_seq(cons_path)
            # sanity: aligned-length consensus must match
            if len(rows) == 0:
                continue
            if len(chr_cons_aln) != len(rows[0][1]):
                raise RuntimeError(f"Chr consensus length != alignment length for {fam} {chr_name}")

            for sid, seq in rows:
                if not is_query_header(sid):
                    continue
                ensure_rec(sid)

                pid, P, Q = score_pairwise(chr_cons_aln, seq, "maskN")
                recs[sid]["chr_PID_maskN"] = pid
                recs[sid]["chr_P_maskN"] = P
                recs[sid]["chr_Q_maskN"] = Q
                recs[sid]["chr_kimura2p_maskN"] = kimura2p(P, Q)

                pid, P, Q = score_pairwise(chr_cons_aln, seq, "countN")
                recs[sid]["chr_PID_countN"] = pid
                recs[sid]["chr_P_countN"] = P
                recs[sid]["chr_Q_countN"] = Q
                recs[sid]["chr_kimura2p_countN"] = kimura2p(P, Q)

    # --- annotate PAR / CEN and chr / repeat_length
    for sid in recs.keys():
        try:
            chrom, start, end = parse_seqname(sid)
        except ValueError:
            # If the row is unexpected, leave defaults
            continue
        # record chr and 1-based inclusive length
        recs[sid]["chr"] = chrom
        recs[sid]["repeat_length"] = (end - start + 1)

        # Convert to 0-based half-open for overlap tests
        start0, end0 = (start - 1), end
        par_hit = overlaps(par_iv.get(chrom, []), start0, end0)
        cen_hit = overlaps(cen_iv.get(chrom, []), start0, end0)
        recs[sid]["PAR"] = bool(par_hit)
        recs[sid]["CEN"] = bool(cen_hit)

    # DataFrame
    df = pd.DataFrame.from_records(list(recs.values()))
    # sort columns in the requested order explicitly
    cols = [
        "seqname", "chr", "repeat_length", "PAR", "CEN",
        "global_PID_maskN", "global_P_maskN", "global_Q_maskN", "global_kimura2p_maskN",
        "global_PID_countN", "global_P_countN", "global_Q_countN", "global_kimura2p_countN",
        "takki_PID_maskN", "takki_P_maskN", "takki_Q_maskN", "takki_kimura2p_maskN",
        "takki_PID_countN", "takki_P_countN", "takki_Q_countN", "takki_kimura2p_countN",
        "chr_PID_maskN", "chr_P_maskN", "chr_Q_maskN", "chr_kimura2p_maskN",
        "chr_PID_countN", "chr_P_countN", "chr_Q_countN", "chr_kimura2p_countN",
    ]
    df = df.reindex(columns=cols)
    return df

# ───────────────────────────────────────────── Main

def main():
    # Load annotation intervals once
    par_iv = load_bed_intervals(PAR_BED)
    cen_iv = load_gff_intervals(CEN_GFF)

    for fam in families:
        df = compute_family_table(fam, par_iv, cen_iv)
        out_tsv = os.path.join(wd, f"{fam}_pid_kimura.tsv")
        df.to_csv(out_tsv, sep="\t", index=False, quoting=csv.QUOTE_NONE)
        print(f"Wrote {out_tsv} with {len(df)} rows.")

if __name__ == "__main__":
    main()
