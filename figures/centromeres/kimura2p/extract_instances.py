#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
extract_instances.py

Goal: Extract repeat instances from GFFs (one per family) using the indexed genome FASTA,
      writing (a) a single global instances FASTA and (b) per-chromosome instance FASTAs,
      in one pass to minimize runtime. If a feature is on the '-' strand, we reverse-complement
      after extracting the genomic interval.

Additionally:
  - Create both {fam}_instances.fa (instances only) and {fam}_for_aln.fa.
  - Prepend the corresponding Takki consensus (WD/inputs/{fam}.fasta) as the FIRST record of {fam}_for_aln.fa.

Inputs/Assumptions:
  - Genome FASTA is indexed (fai present) and readable.
  - WD/inputs contains {fam}.colored.gff and {fam}.fasta for each fam in FAMILIES.
  - Output directories exist (we will create per-chromosome subdirs if needed).

Outputs:
  - WD/{fam}_asm_alignments/{fam}_instances.fa
  - WD/{fam}_asm_alignments/{fam}_for_aln.fa   (Takki first, then instances)
  - WD/{fam}_chr_alignments/{chr}/{fam}_{chr}_for_aln.fa

Rules honored:
  - Use simple, Biopython-based extraction (SeqIO) and reverse-complement logic.
  - No extra validation/file checkers beyond minimal overwrite handling.
"""

import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ---- Config (required) -------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
WD = os.path.normpath(os.path.join(SCRIPT_DIR, ".."))

# Set this to the absolute path of your indexed genome FASTA (outside repo)
GENOME_FASTA = "/mnt/d/research/zebrafinch_t2t_article/split_genome/v0.4/bTaeGut7v0.4.fa"
FAMILIES = ["Tgut191A", "Tgut716A"]
OVERWRITE = True  # if False, skip work when outputs already exist
# -----------------------------------------------------------------------------


# Load genome sequences once (keyed by chromosome/contig id)
# Simplicity > memory: this keeps full sequences in memory for fast slicing.
genome = SeqIO.to_dict(SeqIO.parse(GENOME_FASTA, "fasta"))

for fam in FAMILIES:
    gff_path = os.path.join(WD, "inputs", f"{fam}.colored.gff")
    takki_path = os.path.join(WD, "inputs", f"{fam}.fasta")

    # Global output paths for this family
    fam_asm_dir = os.path.join(WD, f"{fam}_asm_alignments")
    os.makedirs(fam_asm_dir, exist_ok=True)
    instances_fa = os.path.join(fam_asm_dir, f"{fam}_instances.fa")
    global_fa = os.path.join(fam_asm_dir, f"{fam}_for_aln.fa")

    # Skip if outputs already exist and OVERWRITE is False
    if (not OVERWRITE) and os.path.isfile(global_fa) and os.path.isfile(instances_fa):
        print(f"[skip] {fam}: outputs exist and OVERWRITE=False")
        continue

    # Prepare per-chromosome root dir
    fam_chr_root = os.path.join(WD, f"{fam}_chr_alignments")
    os.makedirs(fam_chr_root, exist_ok=True)

    # We'll keep one open handle per chromosome and close them at the end
    chr_handles = {}

    # Helper to get an open handle for a chr-specific file
    def get_chr_handle(chr_name: str):
        if chr_name not in chr_handles:
            chr_dir = os.path.join(fam_chr_root, chr_name)
            os.makedirs(chr_dir, exist_ok=True)
            chr_path = os.path.join(chr_dir, f"{fam}_{chr_name}_for_aln.fa")
            # Always (over)write in this run
            chr_handles[chr_name] = open(chr_path, "w")
        return chr_handles[chr_name]

    # Open the global outputs (always overwrite in this run)
    with open(global_fa, "w") as global_out, open(instances_fa, "w") as instances_out:
        # 1) Write Takki consensus FIRST into the global_for_aln file
        takki_rec = next(SeqIO.parse(takki_path, "fasta"))
        global_out.write(f">{takki_rec.description}\n{str(takki_rec.seq)}\n")

        n_written = 0

        # 2) Stream GFF and append instances to (a) instances.fa, (b) global_for_aln.fa, and (c) per-chr files
        with open(gff_path, "r") as gff:
            for line in gff:
                if not line or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 8:
                    continue  # minimal robustness; assume valid GFF otherwise

                seqid = fields[0]
                # fields[1] = source
                # fields[2] = type
                try:
                    start = int(fields[3])  # 1-based inclusive
                    end = int(fields[4])    # 1-based inclusive
                except ValueError:
                    continue

                strand = fields[6]  # '+' or '-'
                # fields[7] = phase
                # fields[8] = attributes (optional)

                if seqid not in genome:
                    # Unknown contig; skip
                    continue

                # Extract subsequence (convert 1-based inclusive to 0-based slicing)
                # Python slice end is exclusive, so use end as-is to include GFF end.
                subseq = genome[seqid].seq[start - 1:end]

                # Reverse-complement if on '-' strand
                if strand == "-":
                    subseq = subseq.reverse_complement()

                # Build a stable, informative header compatible with chr-splitting logic
                # Pattern: {fam}::{chr}:{start}-{end}|strand={+/-}
                header = f"{fam}::{seqid}:{start}-{end}|strand={strand}"

                # Write to instances-only global file
                instances_out.write(f">{header}\n{str(subseq)}\n")

                # Write to global_for_aln (Takki already written first)
                global_out.write(f">{header}\n{str(subseq)}\n")

                # Write to per-chromosome
                chr_out = get_chr_handle(seqid)
                chr_out.write(f">{header}\n{str(subseq)}\n")

                n_written += 1

        print(f"[done] {fam}: wrote {n_written} instances "
              f"(global: {global_fa} + {instances_fa}; per-chr in {fam_chr_root}/<chr>/)")

    # Close all per-chromosome handles
    for h in chr_handles.values():
        h.close()
