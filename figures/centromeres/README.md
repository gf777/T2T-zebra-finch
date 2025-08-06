# Annotate satellite repeats
```
RepeatMasker -lib Takki2022.library.fasta bTaeGut7v0.4_MT_rDNA.fa.gz -nolow -no_is -pa 8
```
Or in slurm:
```
sbatch ... --array=1-84 repeatMasker_parallel.sh
```
# Get max unit length per chromosome
```
python max_unit_length_per_chromosome.py <(grep Tgut716A bTaeGut7v0.4_MT_rDNA.RM.Takki2022.v0.1.colored.merged.gff) | sort -nk2
```