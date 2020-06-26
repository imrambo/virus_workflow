#/bin/bash
# Get the size of a FASTA file in bp
# Maintainer: Ian Rambo

fa_size=$(grep -v '>' $1 | grep -oe "[ACTG]" | wc -l)
echo $fa_size
