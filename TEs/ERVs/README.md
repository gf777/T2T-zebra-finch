# Check TE annotation
This can be done either on individual units, or extracting entire regions. In the latter case, first we can merge multiple annotations in an array. This basically concatenate annotations of a specific element if they are next/close to each other, and there is nothing in between. The example was used to identify an issue with EDTA2, annotating a satellite repeat (Tgut716A) as an ERV:
```
python ../gff3Merger.py units/bTaeGut7.mat+Z.cur.20250313.Takki2022.v0.1.gff merged/bTaeGut7.mat+Z.cur.20250313.Takki2022.v0.1.merged.gff
grep Tgut716A merged/bTaeGut7.mat+Z.cur.20250313.Takki2022.v0.1.merged.gff > bTaeGut7.mat+Z.cur.20250313.candidate_centromeres.v0.1.gff
```
Then, when dealing with multiple annotations (whether merged or not), after having extracted their sequence and doing some filtering, we can run gmsuite to generate ab initio protein annotations and then interproscan to annotate them functionally:
```
while read line; do grep $line -A1 Tgut716A.repeats-only.fasta.masked.nonewlines.singleNs.gt100.notelo.fasta > out/$line.fa; perl ../../gmsuite/gmsn.pl out/$line.fa --faa --fnn; done<scaffold.ls
./../interProScan/interproscan-5.73-104.0/interproscan.sh -i mat.fasta -f tsv
```
In this case, the presence of GAG, POL or ENV proteins would be indicative of an ERV (but look for consistency accross units and that the protein-coding sequences are not flanking the elements but part of them):
```
grep -i "gag\|pol\|env" mat.fasta.tsv | cut -f1,6
```
