#!/usr/bin/env python3
import re
import sys
from collections import OrderedDict

PDF_PATH = sys.argv[1]
OUT_FASTA = sys.argv[2]

# Regex for alignment lines: "NAME  SEQCHUNK"
# We keep '-' gaps. We ignore "Consensus/80%" and other annotation lines.
line_re = re.compile(r"^(\S+)\s+([A-Za-z\-\(\)0-9]+)\s*$")

# Tokens to ignore as "names"
IGNORE_NAMES = {
    "Consensus/80%",
    "Consensus/80",  # just in case
}

def clean_chunk(s: str) -> str:
    # remove parenthetical residue counts like "(23)" seen in the PDF alignment
    s = re.sub(r"\(\d+\)", "", s)
    # keep only letters and gaps
    s = re.sub(r"[^A-Za-z\-]", "", s)
    return s.upper()

def main():
    try:
        import pdfplumber
    except ImportError:
        print("ERROR: pdfplumber not installed. Try: pip install pdfplumber", file=sys.stderr)
        sys.exit(2)

    seqs = OrderedDict()  # name -> concatenated alignment string

    with pdfplumber.open(PDF_PATH) as pdf:
        for page in pdf.pages:
            text = page.extract_text() or ""
            for raw in text.splitlines():
                line = raw.rstrip()
                m = line_re.match(line)
                if not m:
                    continue
                name, chunk = m.group(1), m.group(2)

                # filter non-sequence lines
                if name in IGNORE_NAMES:
                    continue
                # skip exon/structure labels that sometimes get parsed as "names"
                if name.startswith("Exon") or name in {"0", "1", "2", "3"}:
                    continue
                if name.startswith(("α", "β", "", "")):
                    continue

                chunk = clean_chunk(chunk)
                if not chunk:
                    continue

                if name not in seqs:
                    seqs[name] = ""
                seqs[name] += chunk

    if not seqs:
        print("ERROR: no sequences extracted. Are you pointing to the right PDF?", file=sys.stderr)
        sys.exit(2)

    # sanity: all same length?
    lengths = {len(v) for v in seqs.values()}
    if len(lengths) != 1:
        # still write output, but warn loudly (usually means you missed some chunks)
        lens_preview = sorted((k, len(v)) for k, v in seqs.items())
        print("WARNING: extracted sequences have differing lengths.", file=sys.stderr)
        for k, L in lens_preview[:10]:
            print(f"  {k}\t{L}", file=sys.stderr)
        print("  ...", file=sys.stderr)
        print("This often means the PDF wrapped/omitted some lines; check output carefully.", file=sys.stderr)

    with open(OUT_FASTA, "w") as out:
        for name, aln in seqs.items():
            out.write(f">{name}\n")
            for i in range(0, len(aln), 80):
                out.write(aln[i:i+80] + "\n")

    print(f"Wrote {OUT_FASTA} with {len(seqs)} sequences.")
    print(f"Alignment length(s): {sorted(lengths)}")

if __name__ == "__main__":
    main()

