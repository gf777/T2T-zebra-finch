# HOR identification
To identify candidate High Order Repeats (HORs) in centromere, we can extract the monomer sequences, cluster them, and paint them back on the sequence to see if they show periodicity. Here is an example using only chr37:
```
bedtools intersect -a bTaeGut7.mat+Z.cur.20250313.RM.Takki2022.v0.1.lib-only.gff -b chr37.mat.Tgut716A.centromere-coordinates.bed > bTaeGut7.mat+Z.cur.20250313.RM.Takki2022.v0.1.lib-only.chr37.gff
awk '$0 !~ /^#/ { print $1"\t"($4-1)"\t"$5"\t"$9"\t.\t"$7 }' bTaeGut7.mat+Z.cur.20250313.RM.Takki2022.v0.1.lib-only.chr37.gff > bTaeGut7.mat+Z.cur.20250313.RM.Takki2022.v0.1.lib-only.chr37.bed # convert to bed
bedtools getfasta -fi bTaeGut7v0.4_MT_rDNA.fa -bed bTaeGut7.mat+Z.cur.20250313.RM.Takki2022.v0.1.lib-only.chr37.bed -s -name > bTaeGut7.mat+Z.cur.20250313.RM.Takki2022.v0.1.lib-only.chr37.fasta
```
Once we have a FASTA file with the individual units we can use `vsearch` to do the clustering:
```
vsearch --cluster_fast bTaeGut7.mat+Z.cur.20250313.RM.Takki2022.v0.1.lib-only.chr37.fasta --id 0.95 --centroids centroids.fa --uc clusters.uc
```
Clusters may be well-defined when there are not too many singletons or too few clusters.
Then we can manually label/color the units according to their cluster:
```
python add_cluster_colors.py bTaeGut7.mat+Z.cur.20250313.RM.Takki2022.v0.1.lib-only.chr37.gff clusters.uc bTaeGut7.mat+Z.cur.20250313.RM.Takki2022.v0.1.lib-only.chr37.fasta bTaeGut7.mat+Z.cur.20250313.RM.Takki2022.v0.1.lib-only.chr37.clustered.gff
```
We can also run a PCA to get a sense of the distribution:
```
python PCA.py  bTaeGut7.mat+Z.cur.20250313.RM.Takki2022.v0.1.lib-only.chr1.Tgut191A.fasta --kmer 8
```