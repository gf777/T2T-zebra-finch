repeatMasker
python ../gff3Merger.py units/bTaeGut7.mat+Z.cur.20250313.Takki2022.v0.1.gff merged/bTaeGut7.mat+Z.cur.20250313.Takki2022.v0.1.merged.gff
grep Tgut716A merged/bTaeGut7.mat+Z.cur.20250313.Takki2022.v0.1.merged.gff > bTaeGut7.mat+Z.cur.20250313.candidate_centromeres.v0.1.gff
while read line; do grep $line -A1 Tgut716A.repeats-only.fasta.masked.nonewlines.singleNs.gt100.notelo.fasta > out/$line.fa; perl ../../gmsuite/gmsn.pl out/$line.fa --faa --fnn; done<scaffold.ls
./../interProScan/interproscan-5.73-104.0/interproscan.sh -i mat.fasta -f tsv
grep -i "gag\|pol\|env" mat.fasta.tsv | cut -f1,6
replace with this pattern in full header \|G.*>
#bash group_interprot.sh <(grep -i "gag\|pol\|env" mat.fasta.interprot.tsv) mat.fasta.interprot.grouped.ERVonly.tsv
bash annotate_gm_w_interpro.sh mat.fasta.interprot.ERVonly.tsv out/all.faa out/all.nocomments.fnn out/all.fnn.annotated.fa