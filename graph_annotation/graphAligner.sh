#!/bin/bash

prefix=$(echo graph.index | sed s/.index//)
fileN=$(printf "%03d" $SLURM_ARRAY_TASK_ID)
mkdir -p alignments

GraphAligner -t 32 -g assembly.homopolymer-compressed.gfa -f split/ont${fileN}.fasta.gz -a alignments/aligned${fileN}.WORKING.gaf \
  --diploid-heuristic 21 31 --diploid-heuristic-cache diploid.index \
  --seeds-mxm-cache-prefix $prefix \
  --seeds-mxm-windowsize 5000 \
  --seeds-mxm-length 30 \
  --seeds-mem-count 10000 \
  --bandwidth 15 \
  --multimap-score-fraction 0.99 \
  --precise-clipping 0.85 \
  --min-alignment-score 5000 \
  --hpc-collapse-reads \
  --discard-cigar \
  --clip-ambiguous-ends 100 \
  --overlap-incompatible-cutoff 0.15 \
  --max-trace-count 5 \
  --mem-index-no-wavelet-tree

