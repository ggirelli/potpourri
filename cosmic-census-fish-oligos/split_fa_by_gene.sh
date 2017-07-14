#!/usr/bin/env bash

# Gabriele Girelli - 20170619
# Aim: split fasta file based on gene name
# 
# Usage: ./split_fa_by_gene.sh all_gene.fa outdir
# 
# Note:
# 	- works only if no space is present in the fasta headers
# 	- as sequences are appended, remove previous run outputs when re-running
# 

outdir=$2
mkdir -p $outdir

while IFS='' read -r line || [[ -n "$line" ]]; do
	if [ ${line:0:1} == '>' ]; then
    	name=$(echo "$line" | tr -d '>' | tr '_' ' ' | cut -d ' ' -f 1)
    	echo $line | tr ' ' '\n' >> $outdir"/"$name".fa"
    fi
done < <(cat "$1" | paste - -)
