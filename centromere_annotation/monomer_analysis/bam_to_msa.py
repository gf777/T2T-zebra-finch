import sys
import pysam
from Bio import SeqIO

bam_file = sys.argv[1]
ref_fasta = sys.argv[2]
output_fasta = sys.argv[3]

# Load reference sequences
ref_seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(ref_fasta, "fasta")}

# Open BAM file
bam = pysam.AlignmentFile(bam_file, "rb")

aligned_seqs = {}

# Parse alignments and store gapped sequences
for read in bam.fetch(until_eof=True):
    if read.is_unmapped:
        continue

    ref_name = bam.get_reference_name(read.reference_id)
    ref_seq = ref_seqs[ref_name]
    query_seq = read.query_sequence

    ref_start = read.reference_start
    q_pos = 0
    r_pos = ref_start
    q_aln = []

    # Add leading gaps if alignment doesn't start at position 0
    if ref_start > 0:
        q_aln.append('-' * ref_start)

    for op, length in read.cigartuples:
        if op == 0:  # match/mismatch
            q_aln.append(query_seq[q_pos:q_pos + length])
            q_pos += length
            r_pos += length
        elif op == 1:  # insertion to ref → skip
            q_pos += length
        elif op == 2:  # deletion from ref → pad query with gaps
            q_aln.append('-' * length)
            r_pos += length
        elif op == 4:  # soft clip
            q_pos += length
        elif op == 5:  # hard clip → ignore
            pass

    aligned_query = ''.join(q_aln)
    aligned_seqs[read.query_name] = aligned_query

# Find the max alignment length
max_len = max(len(seq) for seq in aligned_seqs.values())

# Write output with trailing padding
with open(output_fasta, "w") as out:
    for name, seq in aligned_seqs.items():
        padded_seq = seq.ljust(max_len, '-')  # pad with trailing '-'
        out.write(f">{name}\n{padded_seq}\n")

