#!/usr/bin/env python3
"""
alignment_consensus.py — Minimal MAFFT alignment + consensus builder

What it does
------------
Provides *two* independent, comment-to-disable blocks:

A) GLOBAL (per family):
   - Input:  wd/{fam}_asm_alignments/{fam}_for_aln.fa  (already includes Takki)
   - Output: wd/{fam}_asm_alignments/{fam}_aln.fasta
             wd/{fam}_asm_alignments/{fam}_consensus.fa
             wd/{fam}_asm_alignments/{fam}_nonN_consensus.fa

B) PER-CHROMOSOME (per family × chromosome):
   - Input:  wd/{fam}_chr_alignments/{chr}/{fam}_{chr}_for_aln.fa
   - Output: wd/{fam}_chr_alignments/{chr}/{fam}_{chr}_aln.fa
             wd/{fam}_chr_alignments/{chr}/{fam}_{chr}_consensus.fa
             wd/{fam}_chr_alignments/{chr}/{fam}_{chr}_nonN_consensus.fa

This script performs *only* alignment + consensus. No plotting, no identities.

Consensus rule
-----------------------------------
• For each alignment column:
  - Ignore gaps; count only A/C/G/T.
  - Select the modal base (plurality). If none present, write 'N'.
• Reason: gap-aware plurality avoids bias from indels and keeps consensus usable
  for downstream identity calculations without over-complication.

Non-N consensus
---------------
• After computing the consensus, create a second sequence with all 'N' removed
  (i.e., keeping only A/C/G/T). This supports downstream alignment of consensus
  sequences that must be free of ambiguous bases.


"""

import os
import sys
import subprocess
from Bio import AlignIO, SeqIO

# ───────────────────────────────────── Config
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
WD = os.path.normpath(os.path.join(SCRIPT_DIR, ".."))

# Set this to the absolute path of your indexed genome FASTA (outside repo)
GENOME_FASTA = "/mnt/d/research/zebrafinch_t2t_article/split_genome/v0.4/bTaeGut7v0.4.fa"
FAMILIES = ["Tgut191A", "Tgut716A"]
MAFFT = "mafft"                  # change if needed
MAFFT_ARGS = ["--auto", "--thread", "16"]          # or adjust threads via env (e.g., MAFFT_NTHREADS)

OVERWRITE = True                # if False, skip work when outputs already exist

# Optional consensus knobs (kept simple by default)
ALLOW_TIES_AS_N = True           # if two bases tie, write 'N'
MIN_NON_GAP_FRACTION = 0.5       # if non-gap fraction <= this, write 'N' (0.0 disables)

# ───────────────────────────────────── Utilities (kept minimal)
ACGT = set("ACGT")

def compute_consensus_from_alignment(aln,
                                     allow_ties_as_n=ALLOW_TIES_AS_N,
                                     min_non_gap_frac=MIN_NON_GAP_FRACTION):
    """
    Gap-aware plurality consensus:
      - ignore '-' columns in counts;
      - choose modal A/C/G/T; ties -> 'N' if allow_ties_as_n;
      - if non-gap coverage fraction <= min_non_gap_frac, output 'N'.
    """
    # alignment length in columns
    L = aln.get_alignment_length()
    nseq = len(aln)
    cols = []
    for i in range(L):
        col = [rec.seq[i].upper() for rec in aln]
        nongap = [b for b in col if b != "-"]
        frac = len(nongap) / nseq if nseq else 0.0
        if len(nongap) == 0 or frac <= min_non_gap_frac:
            cols.append("N")
            continue
        # count only A/C/G/T
        counts = {}
        for b in nongap:
            if b in ACGT:
                counts[b] = counts.get(b, 0) + 1
        if not counts:
            cols.append("N")
            continue
        # choose modal base
        maxc = max(counts.values())
        winners = [b for b, c in counts.items() if c == maxc]
        if len(winners) > 1 and allow_ties_as_n:
            cols.append("N")
        else:
            cols.append(winners[0])
    return "".join(cols)

def run_mafft(input_fa, output_aln, mafft_bin=MAFFT, extra_args=MAFFT_ARGS):
    with open(output_aln, "w") as out:
        subprocess.run([mafft_bin, *extra_args, input_fa],
                       stdout=out, stderr=subprocess.DEVNULL, check=True)

def write_fasta(path, header, seq):
    with open(path, "w") as h:
        h.write(f">{header}\n")
        # wrap to 80 cols for readability
        for i in range(0, len(seq), 80):
            h.write(seq[i:i+80] + "\n")

def read_first_seq(path):
    rec = next(SeqIO.parse(path, "fasta"), None)
    return str(rec.seq).upper() if rec is not None else ""

def strip_N(seq):
    # Keep only A/C/G/T; consensus only produces A/C/G/T/N, so removing N is sufficient.
    return seq.replace("N", "")

# ───────────────────────────────────── A) GLOBAL (per family)
for fam in FAMILIES:
    fam_dir = os.path.join(WD, f"{fam}_asm_alignments")
    in_fa   = os.path.join(fam_dir, f"{fam}_for_aln.fa")
    aln_fa  = os.path.join(fam_dir, f"{fam}_aln.fasta")
    cons_fa = os.path.join(fam_dir, f"{fam}_consensus.fa")
    nonn_fa = os.path.join(fam_dir, f"{fam}_nonN_consensus.fa")

    if not os.path.isfile(in_fa):
        print(f"[GLOBAL] Input not found, skipping: {in_fa}")
        continue

    # align
    if OVERWRITE or not os.path.isfile(aln_fa):
        print(f"[GLOBAL] MAFFT {fam} …")
        run_mafft(in_fa, aln_fa)
    else:
        print(f"[GLOBAL] Alignment exists, skipping: {aln_fa}")

    # consensus (+ nonN consensus)
    cons_seq = None
    if OVERWRITE or not os.path.isfile(cons_fa):
        print(f"[GLOBAL] Consensus {fam} …")
        aln = AlignIO.read(aln_fa, "fasta")
        cons_seq = compute_consensus_from_alignment(aln)
        write_fasta(cons_fa, f"{fam}_consensus", cons_seq)
    else:
        print(f"[GLOBAL] Consensus exists, skipping: {cons_fa}")

    if OVERWRITE or not os.path.isfile(nonn_fa):
        if cons_seq is None:
            cons_seq = read_first_seq(cons_fa)
        nonn_seq = strip_N(cons_seq)
        write_fasta(nonn_fa, f"{fam}_nonN_consensus", nonn_seq)
        print(f"[GLOBAL] nonN consensus {fam} written.")
    else:
        print(f"[GLOBAL] nonN consensus exists, skipping: {nonn_fa}")

# ───────────────────────────────────── B) PER-CHROMOSOME (per family × chr)
for fam in FAMILIES:
    fam_root = os.path.join(WD, f"{fam}_chr_alignments")
    if not os.path.isdir(fam_root):
        print(f"[PER-CHR] Directory not found, skipping family: {fam_root}")
        continue

    # iterate chromosome directories present on disk
    for chr_name in sorted(os.listdir(fam_root)):
        chr_dir = os.path.join(fam_root, chr_name)
        if not os.path.isdir(chr_dir):
            continue

        in_fa    = os.path.join(chr_dir, f"{fam}_{chr_name}_for_aln.fa")
        aln_fa   = os.path.join(chr_dir, f"{fam}_{chr_name}_aln.fa")
        cons_fa  = os.path.join(chr_dir, f"{fam}_{chr_name}_consensus.fa")
        nonn_fa  = os.path.join(chr_dir, f"{fam}_{chr_name}_nonN_consensus.fa")  # NEW

        if not os.path.isfile(in_fa):
            print(f"[PER-CHR] Input not found, skipping: {in_fa}")
            continue

        # align
        if OVERWRITE or not os.path.isfile(aln_fa):
            print(f"[PER-CHR] MAFFT {fam} {chr_name} …")
            run_mafft(in_fa, aln_fa)
        else:
            print(f"[PER-CHR] Alignment exists, skipping: {aln_fa}")

        # consensus (+ nonN consensus)
        cons_seq = None
        if OVERWRITE or not os.path.isfile(cons_fa):
            print(f"[PER-CHR] Consensus {fam} {chr_name} …")
            aln = AlignIO.read(aln_fa, "fasta")
            cons_seq = compute_consensus_from_alignment(aln)
            write_fasta(cons_fa, f"{fam}_{chr_name}_consensus", cons_seq)
        else:
            print(f"[PER-CHR] Consensus exists, skipping: {cons_fa}")

        if OVERWRITE or not os.path.isfile(nonn_fa):
            if cons_seq is None:
                cons_seq = read_first_seq(cons_fa)
            nonn_seq = strip_N(cons_seq)
            write_fasta(nonn_fa, f"{fam}_{chr_name}_nonN_consensus", nonn_seq)
            print(f"[PER-CHR] nonN consensus {fam} {chr_name} written.")
        else:
            print(f"[PER-CHR] nonN consensus exists, skipping: {nonn_fa}")

print("Done.")
