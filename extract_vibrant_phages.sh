#!/bin/bash

#Maintainer: Ian Rambo :: ian.rambo@utexas.edu
#Thirteen... that's a mighty unlucky number... for somebody!

#Script to extract VIBRANT phage FASTA by viral cycle (lytic, lysogenic) and genome quality (high,medium,low)

has_command () {
    #Make sure you can execute certain commands
    command -v "$1" >/dev/null 2>&1 || { echo "Requires $1. Ensure that $1 is in your \$PATH."; exit 1; }
}

check_failed () {
    #Check a GNU parallel joblog for any jobs that had a non-zero exit code
    tail -n +2 $1 | awk '$7 != 0'
}

has_command parallel && has_command pullseq


usage="$(basename "$0"): extract FASTA for VIBRANT phages, split up by virus cycle and genome quality

where:

   -h --- show this help message
   -i --- directory containing VIBRANT output directories
   -o --- base output directory
   -e --- file extenstion to pull from: fna (nucleotide), faa (amino)
   -s --- text file containing sample IDs of VIBRANT output folders you want to pull from
   -l --- path to GNU parallel joblog
   -j --- number of parallel jobs

    "

novar=1

while getopts ':hi:o:e:s:l:j:' option; do
    case "${option}" in
    h) echo "$usage"
       exit ;;
    i) basein=${OPTARG};;
    o) baseout=${OPTARG};;
    e) ext=${OPTARG};;
    s) sfile=${OPTARG};;
    l) joblog=${OPTARG};;
    j) njobs=${OPTARG};;


    :) printf "missing argument for -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1 ;;
    \?) printf "illegal option: -%s\n" "$OPTARG" >&2
        echo "$usage" >&2
        exit 1 ;;
    esac
done

shift $((OPTIND - 1))

#basein=/home/rambo/projects/CRISPRCas_Sediment/GuaymasC
#baseout=/home/rambo/projects/CRISPRCas_Sediment/GuaymasC/GuaymasC_VIBRANT_extracted_phages
#ext=faa
#sfile=${basein}/sampleids_withphage.txt


# find ${basein}/VIBRANT* -type f -not -empty -name "*phages_combined.${ext}" | \
#     rev | \
#     cut -f1 -d '/' | \
#     rev | \
#     cut -f1 -d'.' | \
#     sort -u > $sfile


declare -a qual=("low" "medium" "high")
declare -a vcycle=("lytic" "lysogenic")


for sid in $(cat $sfile); do
    for vc in "${vcycle[@]}"; do
        outdir=${baseout}/${vc}; test -d $outdir || mkdir -p $outdir
        for q in "${qual[@]}"; do
            echo "$sid --- $vc --- $q quality: extracting"
            grep $vc ${basein}/VIBRANT_${sid}/VIBRANT_results_${sid}/VIBRANT_genome_quality_${sid}.tsv | \
            grep "$q quality" | \
            cut -f1 | \
            sort -u | \
            parallel --joblog $joblog --jobs $njobs pullseq --input ${basein}/VIBRANT_${sid}/VIBRANT_phages_${sid}/${sid}.phages_${vc}.${ext} --regex {}_.* '>' ${outdir}/${sid}_{}_${q}.${ext}
            echo "$sid --- $vc --- $q quality: finished"
        done
    done
done

check_failed $joblog
