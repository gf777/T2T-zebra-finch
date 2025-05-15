#!/bin/bash
name=$(sed -n "${SLURM_ARRAY_TASK_ID}p" scaffolds.ls)
mkdir -p out
echo $name
gfastats ../../bTaeGut7v0.4_MT_rDNA.fa $name -o out/$name.fasta
RepeatMasker -lib Takki2022.library.fasta out/$name.fasta -nolow -no_is -pa 8 -gff