
        set terminal png size 600,400 truecolor
        set output "mat+Z_bt2/mat+Z_bt2-indel-dist.png"
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
        set title "bTaeGut7.mat+Z.cur.20250313.bt2.sorted.bam.stats" noenhanced
        plot '-' w l ti 'Insertions', '-' w l ti 'Deletions', '-' axes x1y2 w l ti "Ins/Dels ratio"
    1	37723915
2	6815314
3	4049505
4	898391
5	202232
6	19713
7	1509
8	920
end
1	19409681
2	3399517
3	1760570
4	320893
5	78709
6	7569
7	135
8	77
end
1	1.943562
2	2.004789
3	2.300110
4	2.799659
5	2.569363
6	2.604439
7	11.177778
8	11.948052
end
