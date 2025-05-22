# Identify kinetochore binding region
To precisely pin down the putative centromere location once a candidate satellite (in our case Tgut716A) has been identified, we can use methylation information:
```
python centromere_detector.py bTaeGut7v0.4_MT_rDNA.RM.Takki2022.v0.1.colored.merged.gff bTaeGut7v0.4_MT_rDNA.PBmethylation.v0.1.bw out.gff --min_length 5000
grep -v kinetochore_binding_site out2.gff > bTaeGut7v0.4_MT_rDNA.centromere_detector.v0.1.gff
```
The function to identify the putative kinetochore binding site based on methylation is not working yet, so we remove that from the annotation.
Note: in the specific case of chr5_mat, because of the inversion, the centromere sequence will be too small. So need a cutoff of 1000 bp only for this one chromosome