#!/bin/bash
### Maintainer: Ian Rambo :: ian.rambo@utexas.edu
### Thirteen... that's a mighty unlucky number... for somebody!


### Purpose: create spacepharer target databases for VIBRANT phage FASTA extracted with extract_vibrant_phages.sh
### Databases are split up by phage genome cycle (lytic, lysogenic) and quality (low, medium, high)

has_command () {
    command -v "$1" >/dev/null 2>&1 || { echo "Requires $1. Ensure that $1 is in your \$PATH."; exit 1; }
}


has_command parallel && has_command spacepharer

seqdir=/home/rambo/projects/CRISPRCas_Sediment/GuaymasC/GuaymasC_VIBRANT_extracted_phages
outbase=/home/rambo/projects/CRISPRCas_Sediment/GuaymasC/GuaymasC_VIBRANT_extracted_phages_spacepharerDB
tmpdir=/home/rambo/tmp

prefix=GuaymasC_VIBRANT_extracted_phages

tstamp=$(date +'%Y-%m-%d_%H-%M-%S')
joblog=~/joblogs/spacepharer_${prefix}_targetDB_${tstamp}.joblog
njobs=4
threads=2
ext=faa

#Spacepharer parameters
#createsetdb database type - 0 == auto, 1 == amino, 2 == nuc
dbtype=0
trans_table=11

# Spacepharer database - base directory
dbdir=${outbase}/${prefix}

test -d $dbdir || mkdir -p $dbdir

spacepharer_cmd="spacepharer createsetdb ${seqdir}/{1}/*_{2}.${ext} ${dbdir}/${prefix}_{1}_{2}{4} $tmpdir --translation-table $trans_table --threads $threads --reverse-fragments {3} --dbtype $dbtype"

parallel --joblog $joblog --jobs $njobs $spacepharer_cmd ::: lytic lysogenic ::: high medium low ::: 0 1 :::+ _targetDB _targetDB_rev
