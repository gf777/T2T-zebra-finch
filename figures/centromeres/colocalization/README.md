# Plot colocalization of centromere/distal probes
Primers have been previously associated in Knief et al. 2016 with either centromeres or distal from centromeres in acrocentric chromosomes. The distance from the largest satellite unit for Tgut716A and Tgut191A can be computed and plotted for macro and microchromosomes.
```
python plot_primer_repeat_distances.py --primers bTaeGut7v0.4_MT_rDNA.BLAST.Knief2016.v0.2.bed --gff bTaeGut7v0.4_MT_rDNA.RM.Takki2022.v0.1.colored.merged.gff --out plot.png --macrochrs macros.ls --label-chrom --largest-repeat-only
```
