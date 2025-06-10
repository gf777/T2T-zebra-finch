# Compute gaps in previous reference bTaeGut1.4
We can use an aligner to assess the gaps in the previous reference relative to the T2T reference. Input to the next steps should be a PAF file. One way to generate it is to use minimap2:
```
minimap2 -x asm5 bTaeGut7.mat+Z.cur.20250313.fasta.gz GCF_003957565.2_bTaeGut1.4.pri_genomic.fna.gz > GCF_003957565_vs_bTaeGut7.mat+Z.cur.20250313.paf
```
An alternative tool is FastGA:
```
FastGA GCF_003957565.2_bTaeGut1.4.pri_genomic.fna.gz bTaeGut7.mat+Z.cur.20250313.fasta.gz -1:GCF_003957565_vs_bTaeGut7.mat+Z.cur.20250313.1aln
ALNtoPAF GCF_003957565_vs_bTaeGut7.mat+Z.cur.20250313.1aln > GCF_003957565_vs_bTaeGut7.mat+Z.cur.20250313.paf
```
First, we need a lookup table with the chromosome names from the old and new reference. We can fetch chromosome names and assign them to NCBI accessions using the [assembly report](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/957/565/GCF_003957565.2_bTaeGut1.4.pri/GCF_003957565.2_bTaeGut1.4.pri_assembly_report.txt):
```
grep -v "^#" GCF_003957565.2_bTaeGut1.4.pri_assembly_report.txt | cut -f3,5 | grep -v na > GCF_003957565.chrTable.tsv
```
In case this genome is not available in NCBI we can use the curation file:
```
awk -F',' '{print $2"\t"$1}' bTaeGut2.hap1.W.cur.20220905.csv > bTaeGut2.hap1.W.cur.20220905.chrTable.tsv
```
Then, using the chr name in the genome as query, we need to add the corresponding chr name in the genome we are using as reference to assess completeness:
```
python match_chromosomes.py GCF_003957565.chrTable.tsv <(cut -f1 bTaeGut7.mat+Z.cur.20250313.fasta.fai) GCF_003957565.chrTable.with_mat.tsv
```
The script adds the lookup value as the 3rd column of `GCF_003957565.chrTable.with_mat.tsv`.
The following command extracts from a paf file (either the one from minimap or from FastGA) only the alignments from the same chromosome, otherwise another chromosome sequence can mask the gap in the chromosome being considered:
```
bash query_uniq.sh GCF_003957565_vs_bTaeGut7.mat+Z.cur.20250313.paf GCF_003957565.chrTable.with_mat.tsv
```
With this alignment, we can extract the gaps in the sequence as a bed file:
```
python findGaps.py GCF_003957565_vs_bTaeGut7.mat+Z.cur.20250313.one2one.paf bTaeGut7.mat+Z.cur.20250313.fasta.fai > GCF_003957565.gaps.bed
```
This gives us the total base pairs of gap sequence:
```
awk '{sum+=$3-$2}END{print sum}' GCF_003957565.gaps.bed
```
This gives us chromosome-specific statistics:
```
bash gap_stats.sh GCF_003957565.gaps.bed bTaeGut7.mat+Z.cur.20250313.fasta.fai
```