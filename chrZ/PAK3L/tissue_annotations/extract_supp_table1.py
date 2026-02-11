#!/usr/bin/env python3
import pdfplumber
import re
import sys

pdf_path = sys.argv[1]
out_tsv = sys.argv[2]

rows = []
in_table = False

with pdfplumber.open(pdf_path) as pdf:
    for page in pdf.pages:
        text = page.extract_text() or ""
        lines = text.splitlines()

        for line in lines:
            if "Supplementary Table 1" in line:
                in_table = True
                continue

            if in_table:
                # stop if next supplementary section begins
                if line.startswith("Supplementary"):
                    in_table = False
                    continue

                # heuristic: table rows start with gene name like PAK3L-1_finch or similar
                if re.match(r"PAK3L", line):
                    parts = re.split(r"\s{2,}", line.strip())
                    rows.append(parts)

# write raw table
with open(out_tsv, "w") as out:
    for r in rows:
        out.write("\t".join(r) + "\n")

print(f"Wrote {out_tsv}")

