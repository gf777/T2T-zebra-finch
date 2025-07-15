######## generate verkko consensus ########
python path_to_gaf.py ${gaf_output} assembly.paths.tsv
regenerate.sh assembly.paths.curated.v0.0.1.gaf asm/ genomic_data/hifi/\*fastq.gz genomic_data/ont/raw/\*fastq.gz asm_edited
# get original scaffold to path mapping (do for both paternal and maternal)
awk 'NR == FNR {scaff[$1]=$2;next} {if(scaff[$1]) print $2"\t"scaff[$1]}' scaffold_to_chr.paternal.tsv <(grep path old_assembly.scfmap | sed -e "s/path //") > path_to_chr.paternal.tsv
awk 'NR == FNR {scaff[$1]=$2;next} {if(scaff[$2]) print $1"\t"scaff[$2]}' path_to_chr.paternal.tsv <(grep path new_assembly.scfmap | sed -e "s/path //") > new_scaffold_to_chr.paternal.tsv
gfastats --include <(cut -f1 new_scaffold_to_chr.paternal.tsv) new_assembly.fasta -k <(awk '{print "COMMENT\t"$0}' new_scaffold_to_chr.paternal.tsv) -o paternal.fasta
gfastats -k <(awk '{print "RENAME\t"$1"\t"$1"_asm5-redo"}' new_scaffold_to_chr.paternal.tsv) paternal.fasta -o paternal.renamed.fasta
