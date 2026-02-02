
        set terminal png size 600,400 truecolor
        set output "dip_bt2-indel-dist.png"
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
        set title "bTaeGut7.dip.cur.20250313.bt2.sorted.bam.stats" noenhanced
        plot '-' w l ti 'Insertions', '-' w l ti 'Deletions', '-' axes x1y2 w l ti "Ins/Dels ratio"
    1	42464129
2	7734419
3	4594823
4	1055850
5	212675
6	19070
7	1267
8	825
end
1	20985723
2	3748816
3	1977436
4	360641
5	83653
6	8088
7	39
8	31
end
1	2.023477
2	2.063163
3	2.323627
4	2.927704
5	2.542348
6	2.357814
7	32.487179
8	26.612903
end
