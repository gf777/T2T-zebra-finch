#activate verkko env
micromamba activate verkko

######## telo and rDNA #########

# first link these files and folder in a new folder (just to avoid messing up with the assembly folder)
ln -s [path/to/asm] asm

ln -s asm/assembly.fasta
ln -s asm/assembly.colors.csv
ln -s asm/assembly.paths.tsv
ln -s asm/assembly.scfmap
ln -s asm/assembly.homopolymer-compressed.gfa
ln -s asm/assembly.homopolymer-compressed.noseq.gfa
ln -s asm/8-hicPipeline

# generate telomere bed track
seqtk telo -d 50000 assembly.fasta > assembly.telomere.bed
# simplify/remove rDNA
bash removeRDNA.sh human-rdna-KY962518.1.hpc.fasta

# extract telo colors for metadata table (probably same as assembly.colors.telo_rdna.csv)
grep telo assembly.homopolymer-compressed.noseq.telo_rdna.gfa | grep ^S | awk 'BEGIN{print "utg,assignment,node_color"}{print $2",telo,#008000"}' > telo_utgs.csv

######## hapmers ########

# first count kmers in father mother and child in compressed mode, something like:
meryl count compress k=21 *fastq.gz output child.meryl
meryl count compress k=21 *fastq.gz output mat.meryl
meryl count compress k=21 *fastq.gz output pat.meryl

# then generate hapmers with merqury:
$MERQURY/_submit_hapmers.sh mat.meryl pat.meryl child.meryl

# use hapmers to assign unitigs (make sure to use compressed kmers)
gfastats assembly.homopolymer-compressed.gfa -o assembly.homopolymer-compressed.fasta --discover-paths # first get the unitig sequence that meryl can work with
meryl-lookup -existence -sequence assembly.homopolymer-compressed.fasta -mers mat.hapmer.meryl pat.hapmer.meryl -output hapmers.out

# extract hapmer assignments
awk 'BEGIN{print "utg,assignment,node_color,hap1,hap2"}{IFS="\t"; gsub(/\_path/, "", $1); printf $1","; if ($4 == $6) {printf "NA,#88FF88"} else if ($4 > $6) {printf "mat,#FF8888"} else {printf "pat,#8888FF"}; printf","; {if ($2!=0) {print $4/$2*100","$6/$2*100} else {print 0","0}}}' hapmers.out > hapmers.csv

######## chr assignment ########

bash getChrNames.sh GCF_003957565.2_bTaeGut1.4.pri_genomic.fna.gz

# then we can map refSeq IDs to numbers (report from RefSeq)
grep -v "^\#" GCF_003957565.2_bTaeGut1.4.pri_assembly_report.txt | cut -f3,7 > refseq_chr_assignments.ls

# or chr .csv to chr numbers:
cut -d',' -f1,2 bTaeGut2.hap1.W.cur.20220905.csv | awk -F',' '{print $2"\t"$1}' > bTaeGut2.hap1.W.cur.20220905_chr_assignments.ls

# verkko
awk 'NR == FNR {chrs[$2]=$1;next}{print $0"\t"chrs[$2]}' bTaeGut2.hap1.W.cur.20220905_chr_assignments.ls translation_hap1 > translation_hap1_with_numbers.tsv
awk 'NR == FNR {chrs[$2]=$1;next}{print $0"\t"chrs[$2]}' bTaeGut2.hap1.W.cur.20220905_chr_assignments.ls translation_hap2 > translation_hap2_with_numbers.tsv

# alternative approach (e.g. for hifiasm)
mashmap -r ../GCF_003957565.2_bTaeGut1.4.pri_genomic.fna.gz -q bTaeGut7.asm.hic.hap1.p_ctg.fasta.gz --pi 95 -s 10000 -t 8 -o hap1.mashmap.out
mashmap -r ../GCF_003957565.2_bTaeGut1.4.pri_genomic.fna.gz -q bTaeGut7.asm.hic.hap2.p_ctg.fasta.gz --pi 95 -s 10000 -t 8 -o hap2.mashmap.out
awk 'NR == FNR {chrs[$2]=$1;next}{print $0"\t"chrs[$6]}' bTaeGut2.hap1.W.cur.20220905_chr_assignments.ls hap1.mashmap.out > translation_hap1_with_numbers.tsv
awk 'NR == FNR {chrs[$2]=$1;next}{print $0"\t"chrs[$6]}' bTaeGut2.hap1.W.cur.20220905_chr_assignments.ls hap2.mashmap.out > translation_hap2_with_numbers.tsv

# project on the unitigs
grep path assembly.scfmap | cut -d' ' -f2,3 | tr ' ' ',' > hap_to_path.csv

######## GraphAligner ########

# first simlink the following:
ln -s asm/3-align/split/
ln -s asm/3-align/diploid.index
ln -s asm/assembly.homopolymer-compressed.gfa

# then run something like and cat the alignments together
sbatch -pvgl -c32 --array=1-154%5 graphAligner.sh
# you can add edge support to the graph with
gfalign eval -g aln.gaf -f assembly.homopolymer-compressed.rDNA_telo.gfa -o assembly.homopolymer-compressed.rDNA_telo.edge_support.gfa

######## combine ########
awk 'NR==1 {print "utg"} NR>1 {split($2, a, ","); for (i in a) print substr(a[i], 1, length(a[i])-1)}' assembly.paths.tsv > utg.ls # list of unitigs
awk -F',' 'NR==1 {print "utg,hap"} NR==FNR{p[$2]=$1; next}; {FS="\t"; split($2, a, ","); for (i in a) print a[i]","p[$1]}' hap_to_path.csv assembly.paths.tsv > utg_to_hap.csv # paths to haplotypes
awk -F'\t' 'NR==1 {print "utg,chr"} NR==FNR{p[$1]=$2; next}; {FS=","; print substr($1, 1, length($1)-1)","p[$2]}' <(cut -f1,5 translation_hap1_with_numbers.tsv translation_hap2_with_numbers.tsv) utg_to_hap.csv > utg_to_chr.nosign.csv # haplotypes to chromosomes, beware 2nd row may need to be removed manually
python combine_annotations.py utg.ls utg_to_chr.nosign.csv hapmers.csv telo_utgs.csv

######## generate verkko consensus ########
python path_to_gaf.py ${gaf_output} assembly.paths.tsv
regenerate.sh assembly.paths.curated.v0.0.1.gaf asm/ genomic_data/hifi/\*fastq.gz genomic_data/ont/raw/\*fastq.gz asm_edited
# get original scaffold to path mapping (do for both paternal and maternal)
awk 'NR == FNR {scaff[$1]=$2;next} {if(scaff[$1]) print $2"\t"scaff[$1]}' scaffold_to_chr.paternal.tsv <(grep path old_assembly.scfmap | sed -e "s/path //") > path_to_chr.paternal.tsv
awk 'NR == FNR {scaff[$1]=$2;next} {if(scaff[$2]) print $1"\t"scaff[$2]}' path_to_chr.paternal.tsv <(grep path new_assembly.scfmap | sed -e "s/path //") > new_scaffold_to_chr.paternal.tsv
gfastats --include <(cut -f1 new_scaffold_to_chr.paternal.tsv) new_assembly.fasta -k <(awk '{print "COMMENT\t"$0}' new_scaffold_to_chr.paternal.tsv) -o paternal.fasta
gfastats -k <(awk '{print "RENAME\t"$1"\t"$1"_asm5-redo"}' new_scaffold_to_chr.paternal.tsv) paternal.fasta -o paternal.renamed.fasta

######## extra ########
# basic dotplot
mashmap -r ref.fa -q qry.fa --pi 95 -s 10000 -t 32 -o mashmap.out
perl generateDotPlot png large mashmap.out

# to better visualize sequence direction in bandage remember to:
- load graph into bandageNG
- Select BandageNG-> Preferences
- Under graph appearence -> "Arrowheads in sigle node style" switch to on.
- Draw graph as usual