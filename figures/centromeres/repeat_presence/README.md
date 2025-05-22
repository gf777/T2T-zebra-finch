# compute repeat abundance in genome
To show the relative abudance of two repeat classes per chromosome we can first count how many bases are annotated for each repeat class:
```
grep Tgut716A bTaeGut7v0.4_MT_rDNA.RM.Takki2022.v0.1.colored.merged.gff | awk '!/^#/ {sum[$1] += $5 - $4} END {for (chr in sum) print chr","sum[chr]}' > Tgut716A.csv
grep Tgut191A bTaeGut7v0.4_MT_rDNA.RM.Takki2022.v0.1.colored.merged.gff | awk '!/^#/ {sum[$1] += $5 - $4} END {for (chr in sum) print chr","sum[chr]}' > Tgut191A.csv 
```
Then combine them into a plot:
```
python presence.py Tgut716A.csv Tgut191A.csv 
```