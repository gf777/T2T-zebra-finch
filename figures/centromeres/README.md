# Annotate satellite repeats
```
RepeatMasker -lib Takki2022.library.fasta bTaeGut7v0.4_MT_rDNA.fa.gz -nolow -no_is -pa 8
```
Or in slurm:
```
sbatch ... --array=1-84 repeatMasker_parallel.sh
```