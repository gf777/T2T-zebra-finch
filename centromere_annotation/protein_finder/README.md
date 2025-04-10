# Extract DNA and protein sequence for specific proteins
We can use the genome annotation to extract useful sequences, e.g. to design centromere-specific probes.
First extract the protein annotation (in this case for the centromere-associated protein CENP-C):
```
CENPC bTaeGut7.mat+Z.cur.20250313.EGAPx.v0.1.gtf > bTaeGut7.mat+Z.cur.20250313.EGAPx.v0.1.CENPC.gtf
```
`gffread` has a known issue with missing transcript ID present in genes from NCBI. Therefore we remove the gene annotation, leaving only transcript and CDS features:
```
awk '{if($3!="gene") print}' bTaeGut7.mat+Z.cur.20250313.EGAPx.v0.1.CENPC.gtf > bTaeGut7.mat+Z.cur.20250313.EGAPx.v0.1.CENPC.nogene.gtf
```
Then we get the mRNA and coding sequence:
```
gffread bTaeGut7.mat+Z.cur.20250313.EGAPx.v0.1.CENPC.nogene.gtf -g bTaeGut7.mat+Z.cur.20250313.fasta -x CENPC.cds.fa -w CENPC.transcript.fa -y CENPC.protein.fa
```