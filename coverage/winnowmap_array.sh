#!/bin/bash
#SBATCH --job-name=winnowmap_array
#SBATCH --cpus-per-task=32
#SBATCH --mem=350G
#SBATCH --time=24:00:00
#SBATCH --output=log/log_%x_%A_%a.out
#SBATCH --error=log/log_%x_%A_%a.err
#SBATCH --partition=vgl_a
#SBATCH --account=vgl_condo_bank 

cpus=$SLURM_CPUS_PER_TASK
mkdir -p bams

# Reference and options
ref=bTaeGut7v0.4_MT_rDNA.fa
map=map-pb
opt=""
repeats=repetitive_k15.txt

# Get FASTQ file (or pair base) from list
read_line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" fastq_list.txt)
base=$(basename "$read_line" | sed 's/_1\.fastq//;s/\.fastq.gz//')

echo "[INFO] Mapping and sorting: $reads"

# Map
winnowmap --MD -W "$repeats" -ax $map $opt -t$cpus "$ref" $read_line | samtools sort -@4 -m2G -T "$base.tmp" -O bam -o "bams/$base.sort.bam"
