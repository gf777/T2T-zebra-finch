#!/bin/bash
PAF=$1
PAIRS=$2
awk 'NR==FNR {keep[$2" "$3]; next} {if (($1" "$6) in keep) print}' $PAIRS $PAF > ${PAF%.*}.one2one.paf
