
        set terminal png size 600,400 truecolor
        set output "pat+WZ_bt2/pat+WZ_bt2-indel-dist.png"
        set grid xtics ytics y2tics back lc rgb "#cccccc"
        set style line 1 linetype 1  linecolor rgb "red"
        set style line 2 linetype 2  linecolor rgb "black"
        set style line 3 linetype 3  linecolor rgb "green"
        set style increment user
        set ylabel "Indel count [log]"
        set xlabel "Indel length"
        set y2label "Insertions/Deletions ratio"
        set log y
        set y2tics nomirror
        set ytics nomirror
        set title "bTaeGut7.pat+WZ.cur.20250313.bt2.sorted.bam.stats" noenhanced
        plot '-' w l ti 'Insertions', '-' w l ti 'Deletions', '-' axes x1y2 w l ti "Ins/Dels ratio"
    1	37723851
2	6811602
3	4046885
4	896428
5	202397
6	19752
7	1541
8	857
end
1	19410482
2	3398589
3	1761080
4	320311
5	78971
6	7678
7	110
8	94
end
1	1.943478
2	2.004244
3	2.297956
4	2.798618
5	2.562928
6	2.572545
7	14.009091
8	9.117021
end
