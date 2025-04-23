# T2T zebra finch (bTaeGut7) assembly project
Description of folders:
- [graph_annotation](graph_annotation) contains one liners and scripts to annotate the Verkko graph prior to assembly graph manual curation
- [rDNA](rDNA) contains one liners and scripts to generate haplotype-specific rDNA models
- [centromere analysis](centromere_analysis) contains one liners and scripts to generate to:
  - merge, colorize, manipulate repeat annotations
  - extract, align, and analyze the putative centromere monomer sequence ([monomer analysis](centromere_analysis/monomer_analysis))
  - extract specific sequences and their protein counterparts, e.g. CENPA ([protein finder](centromere_analysis/protein_finder))
- [completeness](completeness) contains one liners and scripts to annotate gaps in the previous reference relative to the T2T assembly
- - [TEs](TEs) contains one liners and scripts to annotate TEs