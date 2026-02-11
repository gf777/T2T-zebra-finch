#!/usr/bin/env python3
import argparse
import re
import sys

def norm_chr(token: str) -> str:
    # token like chrZ:586... or chrUn:...
    token = token.strip()
    if ":" in token:
        token = token.split(":", 1)[0]
    return token

def norm_expr(expr: str) -> str:
    expr = expr.strip()
    if expr in {"-", "N.A.", "NA", ""}:
        return ""
    # normalize separators/spaces
    expr = expr.replace(" ", "")
    return expr

def expr_category(expr: str) -> str:
    """
    Categories for side-bar:
      - none
      - testis-only
      - multi (>=2 tissues)
      - single:<tissue> (brain/liver/spleen/skin)
    """
    if not expr:
        return "none"
    parts = [p for p in expr.split(",") if p]
    if len(parts) >= 2:
        # special case: if it's just "testis,..." still multi
        return "multi"
    # single tissue
    t = parts[0].lower()
    if t == "testis":
        return "testis-only"
    return f"single:{t}"

def chr_category(ch: str) -> str:
    ch = ch.strip()
    if ch.startswith("chrZ_random"):
        return "chrZ_random"
    if ch == "chrZ":
        return "chrZ"
    if ch == "chrUn":
        return "chrUn"
    return "other"

def main():
    ap = argparse.ArgumentParser(description="Parse PAK3L supplementary table text into annotation TSV.")
    ap.add_argument("supp_table_txt", help="Text file with lines like: 'PAK3L-1 ENST... chrZ:... - testis'")
    ap.add_argument("-o", "--out", default="pak3l_annotations.tsv", help="Output TSV (default: pak3l_annotations.tsv)")
    args = ap.parse_args()

    out_lines = []
    with open(args.supp_table_txt, "r") as f:
        for ln, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue

            # Split into: gene, enst, locus, strand, expr(optional...)
            # Example: PAK3L-4 ENST... chr4A:... - brain,testis
            # Some have "-" as expr placeholder.
            parts = line.split()
            if len(parts) < 4:
                print(f"WARNING: skipping line {ln} (too few fields): {line}", file=sys.stderr)
                continue

            gene = parts[0]                    # PAK3L-#
            locus = parts[2]                   # chrX:...
            strand = parts[3]                  # +/-
            expr = " ".join(parts[4:]) if len(parts) > 4 else ""
            expr = norm_expr(expr)

            ch = norm_chr(locus)

            out_lines.append((gene, ch, chr_category(ch), strand, expr, expr_category(expr)))

    if not out_lines:
        raise SystemExit("No rows parsed. Check input file formatting.")

    with open(args.out, "w") as out:
        out.write("gene\tchr\tchr_cat\tstrand\texpr\texpr_cat\n")
        for gene, ch, chcat, strand, expr, ecat in out_lines:
            out.write(f"{gene}\t{ch}\t{chcat}\t{strand}\t{expr if expr else '.'}\t{ecat}\n")

    print(f"Wrote {args.out} ({len(out_lines)} genes)")

if __name__ == "__main__":
    main()

