#!/bin/bash
BED=$1
FAI=$2
awk '
FILENAME == ARGV[1] {
    seq_order[NR] = $1         # preserve .fai order
    seq_len[$1] = $2
    next
}
{
    gap_len[$1] += $3 - $2
}
END {
    for (i = 1; i <= length(seq_order); i++) {
        seq = seq_order[i]
        g = gap_len[seq] + 0
        t = seq_len[seq]
        frac = g / t
        printf "%s\t%d\t%d\t%.6f\n", seq, g, t, frac
    }
}' $FAI $BED
