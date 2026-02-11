#!/usr/bin/env python3
import sys
import re
from collections import OrderedDict

inp = sys.argv[1]
prefix = sys.argv[2] if len(sys.argv) > 2 else "pak_alignment_fixed"

def read_fasta(path):
    seqs = OrderedDict()
    name = None
    buf = []
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(buf)
            name = line[1:].split()[0]
            buf = []
        else:
            buf.append(line)
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs

def clean_seq(s):
    s = s.strip()
    s = re.sub(r"\s+", "", s)
    s = s.upper()
    # keep AA letters + gap + unknowns
    s = re.sub(r"[^ACDEFGHIKLMNPQRSTVWYXBZJUO\-]", "", s)
    return s

seqs = read_fasta(inp)
seqs = OrderedDict((k, clean_seq(v)) for k, v in seqs.items())

lengths = {len(v) for v in seqs.values()}
Lmax = max(lengths)
Lmin = min(lengths)

# report
with open(f"{prefix}.lengths.tsv", "w") as rep:
    rep.write("id\tlen\n")
    for k, v in seqs.items():
        rep.write(f"{k}\t{len(v)}\n")

if Lmax != Lmin:
    # pad shorter sequences with '-' at end
    for k in list(seqs.keys()):
        v = seqs[k]
        if len(v) < Lmax:
            seqs[k] = v + ("-" * (Lmax - len(v)))

# write fixed FASTA
with open(f"{prefix}.faa", "w") as out:
    for k, v in seqs.items():
        out.write(f">{k}\n")
        for i in range(0, len(v), 80):
            out.write(v[i:i+80] + "\n")

# write strict PHYLIP with stable short IDs + mapping
mapping = []
short_ids = []
for i, full in enumerate(seqs.keys(), start=1):
    short = f"S{i:07d}"  # 8 chars
    mapping.append((short, full))
    short_ids.append(short)

with open(f"{prefix}.name_map.tsv", "w") as m:
    for short, full in mapping:
        m.write(f"{short}\t{full}\n")

with open(f"{prefix}.phy", "w") as f:
    f.write(f"{len(seqs)} {Lmax}\n")
    for short, (full, seq) in zip(short_ids, seqs.items()):
        f.write(f"{short} {seq}\n")

# write NEXUS with full IDs (MrBayes handles longer names fine)
with open(f"{prefix}.nex", "w") as f:
    f.write("#NEXUS\n\nbegin data;\n")
    f.write(f"  dimensions ntax={len(seqs)} nchar={Lmax};\n")
    f.write("  format datatype=protein gap=- missing=?;\n")
    f.write("  matrix\n")
    for k, v in seqs.items():
        f.write(f"{k}\t{v}\n")
    f.write("  ;\nend;\n")

print(f"Read {len(seqs)} seqs. Lengths were {sorted(lengths)}. Now all length = {Lmax}.")
print(f"Wrote: {prefix}.faa  {prefix}.phy  {prefix}.nex  {prefix}.name_map.tsv  {prefix}.lengths.tsv")

