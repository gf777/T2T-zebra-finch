# extract and align centromeric tandem repeats
First, we extract the monomer sequences from an existing annotation:
```
bedtools getfasta -fi bTaeGut7.mat+Z.cur.20250313.fasta -bed <(grep Tgut716A units/bTaeGut7.mat+Z.cur.20250313.RM.Takki2022.v0.1.lib-only.gff) -fo bTaeGut7.mat+Z.cur.20250313.RM.Takki2022.v0.1.lib-only.Tgut716A.fasta
```
Next, we should select the longest and most frequent copy as reference, so that we can reorient all copies relative to the reference to make sure they have the same orientation:
```
# get length distribution and manually extract one of the sequences with the most common length
gfastats bTaeGut7.mat+Z.cur.20250313.RM.Takki2022.v0.1.lib-only.Tgut716A.fasta --path-r | cut -f4 | sort | uniq -c | sort -nk1 | tail
# use relaxed mm2 parameters for short sequences
minimap2 -x sr -t 4 -k7 -w1 -A2 -B4 -O6,24 -E2,1 --secondary=no --score-N=0 --no-end-flt ref.fa bTaeGut7.mat+Z.cur.20250313.RM.Takki2022.v0.1.lib-only.Tgut716A.fasta > alignments.paf
# use the alignment to reorient the sequences
python rvcp.py bTaeGut7.mat+Z.cur.20250313.RM.Takki2022.v0.1.lib-only.Tgut716A.fasta alignments.paf bTaeGut7.mat+Z.cur.20250313.RM.Takki2022.v0.1.lib-only.Tgut716A.oriented.fasta bTaeGut7.mat+Z.cur.20250313.RM.Takki2022.v0.1.lib-only.Tgut716A.unaligned.fasta
```
Finally, we can generate new alignment to convert to msa:
```
minimap2 -ax sr -t 4 -k7 -w1 -A2 -B4 -O6,24 -E2,1 --secondary=no --score-N=0 --no-end-flt ref.fa bTaeGut7.mat+Z.cur.20250313.RM.Takki2022.v0.1.lib-only.Tgut716A.oriented.fasta | samtools view -bS - > alignments.bam
bam_to_msa.py alignments.bam bTaeGut7.mat+Z.cur.20250313.RM.Takki2022.v0.1.lib-only.Tgut716A.oriented.fasta bTaeGut7.mat+Z.cur.20250313.RM.Takki2022.v0.1.lib-only.Tgut716A.aligned.fasta
```
We can then plot alignment conservation:
```
python plot_conservation.py bTaeGut7.mat+Z.cur.20250313.RM.Takki2022.v0.1.lib-only.Tgut716A.aligned.fasta bTaeGut7.mat+Z.cur.20250313.RM.Takki2022.v0.1.lib-only.Tgut716A.aligned.png bTaeGut7.mat+Z.cur.20250313.RM.Takki2022.v0.1.lib-only.Tgut716A.aligned.csv
```