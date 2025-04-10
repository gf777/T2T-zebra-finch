import sys
from Bio import SeqIO

query_fasta = sys.argv[1]
paf_file = sys.argv[2]
aligned_output = sys.argv[3]          # e.g. output_oriented.fa
unaligned_output = sys.argv[4]        # e.g. output_unaligned.fa

# Step 1: Parse strand info from PAF
strand_map = {}
with open(paf_file) as f:
    for line in f:
        fields = line.strip().split('\t')
        if len(fields) > 4:
            qname = fields[0]
            strand = fields[4]
            if qname not in strand_map:
                strand_map[qname] = strand  # take first valid alignment

# Step 2: Process sequences
aligned_seqs = []
unaligned_seqs = []

for rec in SeqIO.parse(query_fasta, "fasta"):
    if rec.id in strand_map:
        strand = strand_map[rec.id]
        if strand == "-":
            rec.seq = rec.seq.reverse_complement()
            rec.id += "_revcomp"
            rec.description += " (reverse complemented)"
        aligned_seqs.append(rec)
    else:
        unaligned_seqs.append(rec)

# Step 3: Write output
SeqIO.write(aligned_seqs, aligned_output, "fasta")
SeqIO.write(unaligned_seqs, unaligned_output, "fasta")

print(f"✓ {len(aligned_seqs)} sequences written to {aligned_output}")
print(f"✓ {len(unaligned_seqs)} sequences written to {unaligned_output}")

