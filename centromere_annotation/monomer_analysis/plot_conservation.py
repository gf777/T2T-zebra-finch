from Bio import AlignIO
import matplotlib.pyplot as plt
import csv
import sys

alignment_file = sys.argv[1]
output_plot = sys.argv[2]
output_csv = sys.argv[3]

# Load alignment
alignment = AlignIO.read(alignment_file, "fasta")
alignment_length = alignment.get_alignment_length()
num_sequences = len(alignment)

# Compute conservation and most common base
conservation = []
most_common_bases = []

for i in range(alignment_length):
    column = alignment[:, i]
    counts = {}
    for base in column:
        counts[base] = counts.get(base, 0) + 1
    most_common_base = max(counts, key=counts.get)
    max_count = counts[most_common_base]

    # Set conservation to 0 if the most common base is a gap
    score = 0.0 if most_common_base == '-' else max_count / num_sequences

    conservation.append(score)
    most_common_bases.append(most_common_base)

# Write CSV: position, most common base, conservation score
with open(output_csv, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Position", "Most_Common_Base", "Conservation"])
    for pos, (base, score) in enumerate(zip(most_common_bases, conservation), 1):
        writer.writerow([pos, base, round(score, 6)])

# Plot barplot
plt.figure(figsize=(15, 4))
plt.bar(range(1, alignment_length + 1), conservation, width=1.0, color="steelblue")
plt.xlabel("Alignment Position")
plt.ylabel("Conservation")
plt.ylim(0, 1.05)
plt.title("Per-Base Sequence Conservation")
plt.tight_layout()
plt.savefig(output_plot, dpi=300)
plt.show()

