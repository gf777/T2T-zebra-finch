# Map homopolymer-compressed utigs to final sequence
First extract, the utigs as hpc fasta:
```
gfastats assembly.homopolymer-compressed.gfa --discover-p -o assembly.homopolymer-compressed.fasta # all
gfastats  assembly.homopolymer-compressed.fasta utig4-23_path -outig4-23.fasta # extract utigs of interest
```
Next extract/compress/index the reference (one chr at a time) and align:
```
gfastats bTaeGut7v0.4_MT_rDNA.fa.gz chrZ_pat -o chrZ_pat.fasta --homo 1
minimap2 -t 32 -H -d chrZ_pat.uncompressed.mmi chrZ_pat.uncompressed.fasta
minimap2 -t 32 -a chrZ_pat.uncompressed.mmi utigs.fasta | samtools sort -@ 8 -o asm_chrZ-utigs_to_chrZ_pat.uncompressed.bam
```
Postprocess to merge the alignments:
```
bedtools bamtobed -i asm_chrZ-utigs_to_chrZ_pat.uncompressed.bam > asm_chrZ-utigs_to_chrZ_pat.uncompressed.bed
at asm_chrZ-utigs_to_chrZ_pat.uncompressed.bed | sort -k4,4 | groupBy -g 1,4,6 -c 2,3 -o min,max, | awk -v OFS='\t' '{print $1, $4, $5, $2, "0", $3, $4, $5, 255,0,0}' | sort -k1,1 -k2,2n > asm_chrZ-utigs_to_chrZ_pat.uncompressed.merge.bed
```