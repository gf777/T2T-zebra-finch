#!/bin/bash
minimap2 -t 32 -ax map-pb --secondary=no bTaeGut7.mat+Z.cur.20250313.fasta.gz CLR/*gz \
  | samtools view -@ 8 -b -o aln.bam
