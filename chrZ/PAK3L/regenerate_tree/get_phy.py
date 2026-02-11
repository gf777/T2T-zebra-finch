# Convert FASTA alignment to PHYLIP deterministically (no external tools)
from collections import OrderedDict
import sys

inp="pak_alignment_from_pdf.faa"
out="pak_alignment_from_pdf.phy"

seqs=OrderedDict()
name=None
buf=[]
for line in open(inp):
    line=line.rstrip("\n")
    if not line: continue
    if line.startswith(">"):
        if name is not None:
            seqs[name]="".join(buf)
        name=line[1:].split()[0]
        buf=[]
    else:
        buf.append(line.strip())
if name is not None:
    seqs[name]="".join(buf)

lengths={len(s) for s in seqs.values()}
assert len(lengths)==1, f"Non-uniform lengths: {sorted((k,len(v)) for k,v in list(seqs.items())[:10])}"
L=next(iter(lengths))

# PHYLIP strict: names <=10 chars; make stable mapping
# We'll write relaxed PHYLIP-like with full names if your PhyML supports it;
# to be safe, we do strict + mapping file.
mapping=[]
names=[]
for i,full in enumerate(seqs.keys(), start=1):
    short=f"S{i:07d}"  # 8 chars
    names.append(short)
    mapping.append((short, full))

with open(out,"w") as f:
    f.write(f"{len(seqs)} {L}\n")
    for short,(full,seq) in zip(names, seqs.items()):
        f.write(f"{short} {seq}\n")

with open("pak_alignment_from_pdf.name_map.tsv","w") as m:
    for short,full in mapping:
        m.write(f"{short}\t{full}\n")

print("Wrote", out, "and pak_alignment_from_pdf.name_map.tsv")

