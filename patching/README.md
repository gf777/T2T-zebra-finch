
# Patching of chrZ
Remaining gaps (GCF_048771995.1_bTaeGut7.mat_genomic.fna.gz):
```
NC_133057.1	269085	269118	33
NC_133062.1	5270652	5370737	100085
NC_133063.1	10173891	10273966	100075
NC_133063.1	10273967	10274062	95
NC_133063.1	22867643	22967723	100080
NC_133063.1	22967724	22967809	85
NC_133063.1	63253513	63261013	7500
NC_133063.1	76526006	76533594	7588
NC_133063.1	76533648	76533673	25
```

Workflow:
1. fix gap1 in verkko
2. fix gap2,3,4 using lastal + gfastats

## 1. gap1 fix
First gap resolution is done using verkko and comparing the hifi-only assembly vs the hifi+correctedONT assembly:
```
# chrZ path in asm5-redo:
>utig4-1327[N5000N:ambig_path]>utig4-1173>utig4-1241<utig4-1861<utig4-1857<utig4-1847[N5000N:ambig_path]<utig4-3505[N100000N:scaffold]>utig4-1030>utig4-1036>utig4-1037[N100000N:scaffold]<utig4-2418
# gap1 (utig4-2418_uti4-1037) patch:
>utig4-1037>utig4-3697<utig4-56>utig4-53<utig4-825>utig4-824>utig4-827<utig4-3284<utig4-1833<utig4-651>utig4-652>utig4-2224<utig4-265<utig4-262<utig4-227>utig4-229>utig4-230<utig4-3020<utig4-3018<utig4-263<utig4-259>utig4-260<utig4-3773<utig4-1660>utig4-1659<utig4-2418
# patched path:
>utig4-1327[N5000N:ambig_path]>utig4-1173>utig4-1241<utig4-1861<utig4-1857<utig4-1847[N5000N:ambig_path]<utig4-3505[N100000N:scaffold]>utig4-1030>utig4-1036>utig4-1037>utig4-3697<utig4-56>utig4-53<utig4-825>utig4-824>utig4-827<utig4-3284<utig4-1833<utig4-651>utig4-652>utig4-2224<utig4-265<utig4-262<utig4-227>utig4-229>utig4-230<utig4-3020<utig4-3018<utig4-263<utig4-259>utig4-260<utig4-3773<utig4-1660>utig4-1659<utig4-2418
# rerun verkko with patch:
sbatch --mem 350G --cpus-per-task=32 --time=5-0 --partition=vgl_a --account=vgl_condo_bank regenerate.sh traversal.test1.gaf asm4/ genomic_data/hifi/\*fastq.gz genomic_data/ont/raw/\*fastq.gz asm4_edited
```

## 2. gap2-4 fix
We want to cut as much into the model as possible, removing all overhangs from the model not from the consensus sequence. But, as going in the tangle, when the consensus and the model disagree, we prefer the model. We noticed the presence of ONT adaptor sequences (identified by FCS) flanking the tangles (not the short masking, which probably doesn't cover it all). This may at least partially explain why they diverge. They likely also diverge because the tips of utigs are imperfect and the unresolved path may be incorrect.

First get base-level alignments:
```
# RefDB
srun --mem 350G --cpus-per-task=32 --time=5-0 --partition=vgl_a --account=vgl_condo_bank lastdb -P32 assembly.rvcp.gap1_resolved assembly.rvcp.gap1_resolved.fasta

srun --mem 350G --cpus-per-task=32 --time=5-0 --partition=vgl_a --account=vgl_condo_bank lastal -P32 assembly.rvcp.gap1_resolved contig-0000001.rvcp.fasta | last-split > gap2.aln.maf
srun --mem 350G --cpus-per-task=32 --time=5-0 --partition=vgl_b --account=vgl_condo_bank lastal -P32 assembly.rvcp.gap1_resolved contig-0000002.fasta | last-split > gap3.aln.maf

gfastats haplotype2-0000040.rvcp.fasta haplotype2-0000040:75049783-87539355 -ohaplotype2-0000040.rvcp.75049783-87539355.fasta # subset at a convenient region to speed up
srun --mem 350G --cpus-per-task=32 --time=5-0 --partition=vgl_a --account=vgl_condo_bank lastal -P32 assembly.rvcp.gap1_resolved haplotype2-0000040.rvcp.75049783-87539355.fasta | last-split > gap4.aln.maf
```

Relevant alignments:
```
gap2 (lastal maf output, usually the last and first alignments are the useful ones):
asm				 22498965 2030847 + 83927455
tangle1					0 2031099 +  7917070

asm			 	 24630063 3279781 + 83927455
tangle1			  4636877 3280192 +  7917070

gap3:
asm			 	 62837139 2105041 + 83927455
tangle2                 1 2105237 +  4220777

asm			 	  64970029 338171 + 83927455
tangle2            3882572 338205 +  4220777

gap4:
asm			 	   75079251 3085522 + 83927455
seq1			          0 3085951 + 12489572

asm			 	    78250110 676851 + 83927455
seq1			    11812660 676912 + 12489572

asm=asm5_redo chrZ
seq1=fragment of asm chrZ (haplotype2-0000040)
```

Example for gap2
The first 2031099 bases of the tangle are aligning from position 22498965 for 2030847 bp. The assembly gap is at:
asm		24530034	24630034
The alignment to the assembly cover 22498965+2030847=24529812 bases, meaning that 24530034-24529812=222 bases don't align an need to be trimmed.

The final sequence should be (see also proper AGP):
```
gaps_resolved			   1	24529812	asm	   			1	24529812
gaps_resolved		24529813	24529814	gap	1
gaps_resolved		24529815	27135591	tangle1	  2031100	 4636876
gaps_resolved		27135592	27135593	gap	1
gaps_resolved		27135594    64942180	asm		 24630035	62436621
gaps_resolved		64942181	64942182	gap	1
gaps_resolved		64942183	66719517	tangle2	  2105238	 3882572
gaps_resolved		66719518	66719519	gap	1
gaps_resolved		66719520    79576093	asm		 65308200	78164773
gaps_resolved		79576094	79576095	gap	1
gaps_resolved		79576096	88302804	seq1	  3085952	11812660
gaps_resolved		88302805	88302806	gap	1
gaps_resolved		88302807    93980152	asm		 78250110	83927455
```

Actual AGP (gaps_resolved.agp):
```
##agp-version	2.0
# ORGANISM: zebra finch (bTaeGut7)
# ASSEMBLY NAME: chrZ T2T gaps_resolved
# DESCRIPTION: In this version the 4 gaps were closed by 1 tangle traversal using the hifi-only assembly in the hifi+correctedONT graph (gap1), tangle modelling using TTT (gaps2,3), alignment to the hifi-only assembly and stitching (gap4; the tangle traversal approach did fully extend through the tangle in the hifi+correctedONT assembly).
# --------------------------------------------------------------
gaps_resolved	1	24529812	1	W	contig-0000001	1	24529812	+
gaps_resolved	24529813	24529814	2	N	1	scaffold	yes	paired-ends
gaps_resolved	24529815	27135591	3	W	tangle1	2031100	4636876	+
gaps_resolved	27135592	27135593	4	N	1	scaffold	yes	paired-ends
gaps_resolved	27135594	64942180	5	W	contig-0000001	24630035	62436621	+
gaps_resolved	64942181	64942182	6	N	1	scaffold	yes	paired-ends
gaps_resolved	64942183	66719517	7	W	tangle2	2105238	3882572	+
gaps_resolved	66719518	66719519	8	N	1	scaffold	yes	paired-ends
gaps_resolved	66719520	79576093	9	W	contig-0000001	65308200	78164773	+
gaps_resolved	79576094	79576095	10	N	1	scaffold	yes	paired-ends
gaps_resolved	79576096	88302804	11	W	haplotype2-0000040	3085952	11812660	+
gaps_resolved	88302805	88302806	12	N	1	scaffold	yes	paired-ends
gaps_resolved	88302807	93980152	13	W	contig-0000001	78250110	83927455	+
```

Generate consensus:
```
cat assembly.rvcp.gap1_resolved.fasta contig-0000001.rvcp.fasta contig-0000002.fasta haplotype2-0000040.rvcp.75049783-87539355.fasta > combined.fasta
gfastats -a gaps_resolved.agp combined.fasta -bg -oassembly.rvcp.gaps_resolved.fasta
gfastats assembly.rvcp.gaps_resolved.fasta -bg | awk '{print $0, $3-$2}' # sanity check
```

# Adaptors removal on chr32:
chr32_mat:268,981 <- first base to be removed
chr32_mat:269,120 <- last base to be removed
resulting sequence is
CTCTTTT[]CCCTGGG