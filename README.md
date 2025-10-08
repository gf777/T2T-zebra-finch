# T2T zebra finch (bTaeGut7) assembly project
Description of folders:
- [centromere_annotation](centromere_annotation) contains one liners and scripts to generate to:
  - merge, colorize, manipulate repeat annotations
  - extract, align, and analyze the putative centromere monomer sequence ([monomer analysis](centromere_analysis/monomer_analysis))
  - extract specific sequences and their protein counterparts, e.g. CENPA ([protein finder](centromere_analysis/protein_finder))
- [completeness](completeness) contains one liners and scripts to annotate gaps in the previous reference relative to the T2T assembly
- [coverage](coverage) contains scripts to compute read coverage statistics
- [figures](figures) contains scripts used to generate panels in the associated manuscript
- [graph_annotation](graph_annotation) contains one liners and scripts to annotate the Verkko graph prior to assembly graph manual curation
- [hpc_to_ref](hpc_to_ref) contains commands to generate alignments from homopolymer-compressed space to the lineare reference
- [rDNA](rDNA) contains one liners and scripts to generate haplotype-specific rDNA models
- [TEs](TEs) contains one liners and scripts to annotate TEs
- [verkko_consensus](verkko_consensus) contains commands and scripts helpful when iteratively improving the consensus sequence towards T2T