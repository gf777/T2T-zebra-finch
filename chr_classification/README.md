## Generating the chromosome-level figures

This repository contains scripts to generate three chromosome-level figures comparing zebra finch and chicken genomes, based on whole-genome alignments and chromosome metadata tables.

### Overview of the workflow

1. **Generate whole-genome alignments**
2. **Convert alignments to PAF**
3. **Run plotting scripts to generate the figures**

Each step is described below.

---

## 1. Whole-genome alignment (FastGA)

Genome-wide alignments between the zebra finch maternal assembly and the chicken reference are generated using **FastGA**:

```bash
FastGA bTaeGut7v0.4_MT_rDNA.mat.wZ.fa chicken.v23.woM.fa -1:mat_vs_chicken.1aln
```

This produces a FastGA alignment file (`.1aln`) containing all primary alignments.

---

## 2. Convert alignments to PAF

The FastGA output is converted to PAF format using `ALNtoPAF`:

```bash
ALNtoPAF mat_vs_chicken.1aln > mat_vs_chicken.paf
```

The resulting PAF file is the main input for the heatmap-style chromosome comparison.

---

## 3. Chromosome heatmap / dot-chromosome annotation

The chromosome–chromosome correspondence heatmap is generated with `chr_heatmap.py`, combining the PAF alignment with zebra finch and chicken chromosome annotations:

```bash
python chr_heatmap.py mat_vs_chicken.paf \
  -o mat_vs_chicken.tsv \
  --zebra-table PCA/Supplementary_Table_9.csv \
  --chicken-class chicken_chr_classification.tsv \
  --x-anno chicken \
  --y-anno zebrafinch \
  --y-lanes "Classification_size,Chr_type,TCHEST,Dot" \
  --lane-label Chr_type="Centromere position" \
  --lane-label Dot="Dot chromosome" \
  --lane-label TCHEST="TCHEST" \
  --lane-label Classification_size="Size" \
  --lane-label Chicken_class="Chicken classification" \
  --legend-fontsize 11 \
  --legend-title-fontsize 12
```

This command produces:

```text
mat_vs_chicken.tsv.matrix_primary.svg
```

which is the main chromosome-level heatmap figure.

---

## Other figures

Additional scripts in this repository generate complementary figures using the same metadata tables:

* **GC% vs chromosome size** (`size_vs_gc.py`)
* **Chromosome PCA based on size-corrected genomic features** (`PCA.py`)

All scripts use `Supplementary_Table_9.csv` as the primary zebra finch chromosome annotation table and are designed to produce publication-ready SVG figures.

---

### Notes

* All plots are generated as **SVG** to preserve vector quality.
* No external PCA or repel-label libraries are required; all analyses rely on NumPy, pandas, and matplotlib.
* Dot chromosomes, macro/microchromosomes, and haplotypes are encoded consistently across figures.