# PAK3L (finch) array analysis: miniprot → collapsed loci → cluster copy heatmap + chr/expr + subtree

This directory contains a reproducible pipeline to:

1. Map a set of query proteins (PAK3L family) to a target genome with **miniprot**.
2. Collapse miniprot “isoforms” (multiple mRNA models per genomic placement) into **one locus (copy)** with a **single best gene assignment**.
3. Define genomic **clusters** (arrays) via BED intervals and count **copies per cluster per gene**.
4. Plot a **cluster×gene** copy heatmap with:
   - a **left phylogenetic subtree** (subset to genes actually found on chrZ clusters),
   - a **right annotation panel** with chromosome label and qualitative tissue-expression ticks.
5. Perform a final **sanity check** by extracting miniprot-predicted proteins, aligning, and visually verifying that sequences assigned to the same “best gene” are more similar to each other.

---

## Contents

- `collapse_miniprot_isoforms.py`  
  Collapse miniprot mRNA models into loci; assign each locus to a single best gene.
- `count_collapsed_copies_per_cluster.py`  
  Count collapsed loci per gene within each BED-defined cluster interval.
- `plot_all.py` (or your final plot script name)  
  Plot: tree (left) + heatmap (center) + colorbar + (chr + expression ticks on right).  
  **Tree is subset to genes actually found on chrZ clusters** (from counts table), not from annotations.
- `extract_miniprot_proteins.py`  
  Translate proteins for each miniprot mRNA from the genome using CDS features.
- `label_proteins_by_best_gene.py`  
  Rename protein headers to encode best-gene assignments from the collapsed loci table.

Input files typically used:
- `queries.faa` (original protein sequences you map; e.g., PAK3L-1_finch_ENST..., etc.)
- `genome.fa` (target genome used for miniprot; must match the GFF coordinates)
- `clusters.bed` (BED3 intervals defining genomic clusters/arrays)
- `pak3l_annotations.tsv` (gene annotations: gene, chr, expr)
- `pak_alignment_fixed.phy_phyml_tree.fixed_names.nwk` (Newick tree for the family)

---

## Software requirements

### Core tools
- `miniprot` (protein→genome spliced alignment)
- `python3` (>=3.9 recommended) with:
  - `numpy`
  - `matplotlib`
- `mafft` (sanity-check multiple sequence alignment)
- `seaview` (optional, for interactive alignment visualization)

Optional:
- `FastTree` (quick tree sanity check from aligned proteins)

### Installation (examples)

Conda/mamba:
```bash
micromamba create -n pak3l -c conda-forge -c bioconda python=3.11 numpy matplotlib mafft fasttree -y
micromamba activate pak3l
# miniprot may be available via bioconda on your platform; otherwise install from release
````

---

## 0) Input preparation

### 0.1 Query proteins

Your query FASTA should contain the protein sequences you want to map and/or classify. Example:

* `queries.faa`

In our case we are using the protein sequences available in the Supplementary Material of Kong et al. 2010. Headers can be anything, but downstream we normalize IDs by taking everything before the first `_`:

* `PAK3L-3_finch_ENSTGUG...` → `PAK3L-3`

### 0.2 Target genome

The genome FASTA used for miniprot, e.g. `genome.fa`, in this case zebra finch T2T chrZ.

> Important: `genome.fa` must match the assembly used to interpret the miniprot GFF coordinates.

### 0.3 Cluster intervals (arrays)

A BED3 file defining the array clusters you want to summarize (0-based, half-open):

* `clusters.bed`

Example:

```text
chrZ_pat    5036113   5863399
chrZ_pat    9792418   10190885
...
```

This can be generated using [bed_merge_only.py](clusters/bed_merge_only.py) (adjust parameters to get the desired number of clusters.

### 0.4 Gene annotation table for plotting

A TSV file:

* `pak3l_annotations.tsv`

Columns:

* `gene`  (e.g., `PAK3L-1` or `PAK3L-1_finch_ENST...`; will be normalized)
* `chr`   (e.g., `chrZ`, `chrUn`, `chr4A`, etc.)
* `expr`  (comma-separated tissues or `-`)

Example:

```text
gene    chr     expr
PAK3L-1 chrZ    testis
PAK3L-4 chr4A   brain,testis
PAK3L-9 chr3_random spleen,skin,testis
```

To generate this follow the instructions in the [tissue_annotations](tissue_annotations) folder.

---

## 1) Map query proteins to the genome with miniprot

Run miniprot with sensitive settings (example; adjust to your preference):

```bash
miniprot -t 8 -N 500 -p 0.1 --outs 0.6 --outc 0.03 \
  genome.fa queries.faa > pak3l.sensitive.gff3
```

Notes:

* This often finds **all copies**, but may output many **isoforms** per genomic locus.
* The next step collapses these isoforms.

---

## 2) Collapse isoforms into loci + single best-gene assignment

This step groups overlapping/nearby miniprot mRNA models into a single **locus** and assigns the locus to **one best gene** (based on best scoring mRNA model within the locus).

```bash
python3 collapse_miniprot_isoforms.py \
  --gff pak3l.sensitive.gff3 \
  --merge-bp 2000 \
  --out-prefix pak3l
```

Outputs:

* `pak3l.collapsed.bed`
  BED6: `chrom start end best_gene 0 strand`
  (good for IGV; one interval per inferred copy/locus)
* `pak3l.collapsed.tsv`
  Per-locus trace table including:

  * `best_gene`, `best_mrna_id`
  * `all_mrna_ids` that were collapsed into the locus

Tuning:

* `--merge-bp` controls isoform collapsing:

  * Increase if obvious single copies still have multiple loci.
  * Decrease if adjacent tandem copies get merged.

---

## 3) Count copies per cluster per gene (using collapsed loci)

Now treat each collapsed locus as **one copy**. Intersect loci with cluster BED intervals and count copies per gene.

```bash
python3 count_collapsed_copies_per_cluster.py \
  --clusters-bed clusters.bed \
  --collapsed-bed pak3l.collapsed.bed \
  --min-overlap-bp 1 \
  -o cluster_gene_counts.tsv
```

Output schema (used by plotting):

* `cluster_id`
* `cluster_chrom`
* `cluster_start`
* `cluster_end`
* `gene`
* `count`  (copies)

---

## 4) Plot: phylogeny + copy heatmap + chr/expr annotations

You have a Newick tree for the family, e.g.:

* `pak_alignment_fixed.phy_phyml_tree.fixed_names.nwk`

For Kong, 2010, this can be generated as described in [regenerate_tree](regenerate_tree).

Run the final plotting script:

```bash
python3 plot_all.py \
  --counts cluster_gene_counts.tsv \
  --annos pak3l_annotations.tsv \
  --tree pak_alignment_fixed.phy_phyml_tree.fixed_names.nwk \
  -o pak3l_chrZ_tree_heatmap.svg
```

### What the plot does (important details)

* **Heatmap rows** (genes) are sorted by the **tree leaf order**.
* The **tree is pruned** to include only genes that were *actually found on chrZ clusters*, defined as:

  * rows where `count > 0` and `cluster_chrom` starts with `chrZ`
* The right-side annotation panel shows:

  * `Chr` label (from `pak3l_annotations.tsv`, for context only)
  * tissue tick marks (presence/absence, qualitative)

If you see multicolored branches in the tree:

* it was Matplotlib’s default color cycle
* fix by forcing `color="black"` in tree plotting calls (already recommended)

---

## 5) Final sanity check: extract miniprot proteins → align → validate gene assignments visually

Goal:

* Show that sequences assigned to the same `best_gene` are more similar to each other.

### 5.1 Extract proteins for each miniprot mRNA

This translates each miniprot mRNA’s CDS from the genome.

```bash
python3 extract_miniprot_proteins.py \
  --gff pak3l.sensitive.gff3 \
  --genome genome.fa \
  -o miniprot_models.faa
```

Output:

* `miniprot_models.faa` with headers as miniprot mRNA IDs (e.g., `MP000123`)

### 5.2 Relabel proteins by best-gene assignment

Use `pak3l.collapsed.tsv` to map each isoform mRNA ID to its locus’s `best_gene`.

```bash
python3 label_proteins_by_best_gene.py \
  --proteins miniprot_models.faa \
  --collapsed-tsv pak3l.collapsed.tsv \
  -o miniprot_models.labeled.faa
```

FASTA headers become:

* `PAK3L-14|MP000123`

### 5.3 Align proteins with MAFFT

```bash
mafft --auto miniprot_models.labeled.faa > miniprot_models.labeled.aln.faa
```

### 5.4 Visual inspection in SeaView

```bash
seaview miniprot_models.labeled.aln.faa &
```

What to check:

* Sequences with the same `PAK3L-#|...` label should show higher similarity and shared motifs.
* Between-label divergence should be larger.

Optional: quick clustering tree from proteins

```bash
FastTree -wag miniprot_models.labeled.aln.faa > miniprot_models.labeled.aln.tree.nwk
```

---

## Outputs summary

Primary analysis products:

* `pak3l.collapsed.bed`  (IGV-ready loci/copies, labeled by best gene)
* `cluster_gene_counts.tsv`  (copy counts per cluster per gene)
* `pak3l_chrZ_tree_heatmap.svg`  (final figure)

Sanity-check products:

* `miniprot_models.faa` (proteins per miniprot mRNA)
* `miniprot_models.labeled.aln.faa` (alignment with best-gene labels)

---

## Troubleshooting

### Tree leaf names don’t match gene IDs

* The plotting and pruning normalize leaf names by splitting at `_`, so:

  * `PAK3L-14_finch` → `PAK3L-14`
* If your tree uses different leaf naming, rename leaves or add a preprocessing map.

### Too many loci after collapsing isoforms

* Increase:

  * `collapse_miniprot_isoforms.py --merge-bp 5000` (or 10000)

### Adjacent tandem copies merged incorrectly

* Decrease:

  * `--merge-bp 1000` or `500`

### Extracted proteins look wrong (many X or many stops)

* Confirm:

  * `genome.fa` matches the miniprot GFF coordinates
* If still wrong, we may need to adjust phase handling or use a different reconstruction strategy.

---

## Reproducibility notes

* Record tool versions:

  ```bash
  miniprot --version || true
  mafft --version
  python3 --version
  ```
* Keep the exact command lines (above) in your run logs.
* Keep `queries.faa`, `genome.fa`, and `pak3l.sensitive.gff3` as immutable inputs for reruns.

---
