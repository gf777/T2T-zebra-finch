import re, sys

map_tsv = "pak_alignment_fixed.name_map.tsv"
tree_fn = "pak_alignment_fixed.phy_phyml_tree.txt"
out_fn  = "pak_alignment_fixed.phy_phyml_tree.fixed_names.nwk"

# read map
mp = {}
with open(map_tsv) as f:
    for line in f:
        line=line.rstrip("\n")
        if not line or line.startswith("#"): 
            continue
        old,new = line.split("\t", 1)
        mp[old]=new

# read tree
tree = open(tree_fn).read().strip()

# Newick labels are sequences not containing: whitespace, : , ( ) ; [ ]
# We'll replace only full label tokens.
pat = re.compile(r'(?<![\w.-])(' + "|".join(map(re.escape, sorted(mp, key=len, reverse=True))) + r')(?![\w.-])')

def repl(m):
    return mp[m.group(1)]

new_tree, n = pat.subn(repl, tree)

# sanity: report unmapped S000... labels remaining
remaining = sorted(set(re.findall(r'\bS\d{7}\b', new_tree)))
sys.stderr.write(f"Replaced {n} labels\n")
if remaining:
    sys.stderr.write(f"WARNING: still present {len(remaining)} unmapped IDs, e.g. {remaining[:10]}\n")

with open(out_fn, "w") as f:
    f.write(new_tree + ("\n" if not new_tree.endswith("\n") else ""))
print(out_fn)

