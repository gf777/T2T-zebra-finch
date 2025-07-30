# Compute average read length per base
To achieve this, a strategy modified from https://github.com/arangrhie/T2T-Polish can be used.
Reads should be first aligned (see [winnowmap_array.sh](winnowmap_array.sh) or [align.sh](align.sh) in this folder).
Next, the alignment needs to be filtered to remove secondary alignments and unmapped reads:
```
samtools view -F0x104 -@32 -hb v0.4_dip_hifi.bam > v0.4_dip_hifi.filtered.bam
```
The alignment can then be converted to paf and the relevant wig files generated (it needs the binaries in this folder):
```
sam2paf.sh v0.4_dip_hifi.filtered.bam v0.4_dip_hifi.filtered.paf bTaeGut7.v0.4_dip_hifi
bash collect_summary.sh v0.4_dip_hifi.filtered.paf v0.4_dip_hifi.filtered full
```
Wig files can be further converted to bigWig for visualization:
```
wigToBigWig v0.4_dip_hifi.filtered.filtered.idy.med.wig bTaeGut7v0.4_MT_rDNA.fa.chrom.sizes bTaeGut1.4_CLR.idy.med.bw -clip
```
To view them:
```
bigWigToWig bTaeGut1.4_CLR.cov.wig -chrom=chr1_mat bTaeGut1.4_CLR.cov.bw
bigWigToBedGraph bTaeGut7.mat+Z.cur.20250313.bwa.10000.pairs.sorted.dedup.cool.compartments.cis.bw
```