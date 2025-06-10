import sys

def read_lookup_file(filename):
    with open(filename) as f:
        return [line.strip() for line in f if line.strip()]

def normalize(chrom):
    return chrom.replace("chr", "").replace("_mat", "").replace("_pat", "").upper()

# Manual override
manual_map = {
    "1B": "chr1_mat",
}

def best_match(chrom, lookup_values):
    if chrom in manual_map:
        return manual_map[chrom]

    norm = normalize(chrom)
    exact_matches = [v for v in lookup_values if normalize(v) == norm]
    if len(exact_matches) == 1:
        return exact_matches[0]

    prefix_matches = [v for v in lookup_values if normalize(v).startswith(norm)]
    if len(prefix_matches) == 1:
        return prefix_matches[0]

    partial_matches = [v for v in lookup_values if norm in normalize(v)]
    if len(partial_matches) == 1:
        return partial_matches[0]

    return None  # ambiguous or not found

def prompt_user(key, lookup_values):
    print(f"Cannot determine match for '{key}'. Options include:", file=sys.stderr)
    for i, opt in enumerate(lookup_values):
        print(f"  {i}: {opt}", file=sys.stderr)
    while True:
        user_input = input(f"Enter index, name, or 'undefined' for '{key}': ").strip()
        if user_input.lower() in {"undefined", ""}:
            return ""
        elif user_input.isdigit() and int(user_input) < len(lookup_values):
            return lookup_values[int(user_input)]
        elif user_input in lookup_values:
            return user_input
        else:
            print("Invalid input. Try again.", file=sys.stderr)

def main(input_file, lookup_file, output_file):
    lookup_values = read_lookup_file(lookup_file)
    resolved = {}

    with open(input_file) as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            line = line.rstrip()
            if not line.strip():
                continue
            fields = line.split('\t')
            key = fields[0]

            if key in resolved:
                match = resolved[key]
            else:
                match = best_match(key, lookup_values)
                if match is None:
                    match = prompt_user(key, lookup_values)
                resolved[key] = match

            f_out.write('\t'.join(fields + [match]) + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python match_chromosomes.py input.tsv lookup.txt output.tsv", file=sys.stderr)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])