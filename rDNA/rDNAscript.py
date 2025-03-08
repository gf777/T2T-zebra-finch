import sys
import random
import re
from Bio import SeqIO


def parse_fasta(fasta_file, unit_coverage):
    morphs = []
    total_expected_copies = 0

    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            match = re.search(r'coverage(\d+)', record.description)  # Updated regex
            if match:
                observed_copies = int(match.group(1))
                expected_copies = round(observed_copies / unit_coverage)
                morphs.append((record.id, str(record.seq), expected_copies, observed_copies))
                total_expected_copies += expected_copies

    return morphs, total_expected_copies


def generate_random_morph_sequence(morphs, total_expected_copies):
    random.seed(42)  # Added seed for reproducibility
    weighted_morphs = []
    actual_counts = {morph_id: 0 for morph_id, _, _, _ in morphs}

    for morph_id, sequence, count, _ in morphs:
        for _ in range(count):
            weighted_morphs.append((morph_id, sequence))

    random.shuffle(weighted_morphs)

    mixed_sequence = ""
    for morph_id, sequence in weighted_morphs:
        mixed_sequence += sequence
        actual_counts[morph_id] += 1

    return mixed_sequence, actual_counts, total_expected_copies


def main():
    if len(sys.argv) < 6:
        print("Usage: python script.py <input.fasta> -o <output.fasta> -c <unit_coverage>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    output_file = sys.argv[3] if len(sys.argv) > 3 and sys.argv[2] == "-o" else None
    unit_coverage = float(sys.argv[5]) if len(sys.argv) > 5 and sys.argv[4] == "-c" else None

    if unit_coverage is None:
        print("Error: unit_coverage is a required parameter.")
        sys.exit(1)

    morphs, total_expected_copies = parse_fasta(fasta_file, unit_coverage)
    mixed_sequence, actual_counts, total_expected_copies = generate_random_morph_sequence(morphs, total_expected_copies)

    if output_file:
        with open(output_file, "w") as f:
            f.write(f">randomized_morphs\n{mixed_sequence}\n")

    print("\n# Morph Copy Statistics")
    for morph_id, _, expected_copies, observed_copies in morphs:
        actual_copies = actual_counts[morph_id]
        frequency = actual_copies / total_expected_copies if total_expected_copies else 0
        print(
            f"{morph_id}: Observed={observed_copies}, Expected={expected_copies}, Actual={actual_copies}, Frequency={frequency:.4f}")


if __name__ == "__main__":
    main()
