#!/bin/bash

#Maintainer: Ian Rambo :: ian.rambo@utexas.edu
#Thirteen... that's a mighty unlucky number... for somebody!

#Script to extract VIBRANT phage FASTA by viral cycle (lytic, lysogenic) and genome quality (high,medium,low)

has_command () {
    command -v "$1" >/dev/null 2>&1 || { echo "Requires $1. Ensure that $1 is in your \$PATH."; exit 1; }
}

has_command parallel && has_command pullseq

sfile=/home/rambo/projects/CRISPRCas_Sediment/GuaymasC/guaymas_sampleids_withorf.txt
basein=/home/rambo/projects/CRISPRCas_Sediment/GuaymasC
baseout=/home/rambo/projects/CRISPRCas_Sediment/GuaymasC/GuaymasC_VIBRANT_extracted_phages

ext=faa

find ${basein}/VIBRANT* -type f -not -empty -name "*phages_combined.${ext}" | \
    rev | \
    cut -f1 -d '/' | \
    rev | \
    cut -f1 -d'.' | \
    sort -u > $sfile


declare -a qual=("low" "medium" "high")
declare -a vcycle=("lytic" "lysogenic")


for sid in $(cat $sfile); do
    for vc in "${vcycle[@]}"; do
        outdir=${baseout}/${vc}; test -d $outdir || mkdir -p $outdir
        for q in "${qual[@]}"; do
            echo "$sid --- $vc --- $q quality: extracting"
            grep $vc ${baseout}/VIBRANT_${sid}/VIBRANT_results_${sid}/VIBRANT_genome_quality_${sid}.tsv | \
            grep "$q quality" | \
            cut -f1 | \
            sort -u | \
            parallel pullseq --input $baseout/VIBRANT_${sid}/VIBRANT_phages_${sid}/${sid}.phages_${vc}.${ext} --regex {}_.* '>' ${outdir}/${sid}_{}_${q}.${ext}
            echo "$sid --- $vc --- $q quality: finished"
        done
    done
done
