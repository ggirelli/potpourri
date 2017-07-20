#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.0
# Date: 20170719
# Project: 680 genes
# Description: Filter fasta files based on filterd and rearranged BLASTN output.
# 
# ------------------------------------------------------------------------------



# ENV VAR ======================================================================

export LC_ALL=C

# PARAMS =======================================================================


# Help string
helps="
usage: ./blast_filter.sh [-h][-k k][-g gt][-s st]
                         -i stg -y gs -f fasta_in -o fasta_out

 Description:
  Filter BLASTN output based on: number of OT, and common (saturated) OTs.
  BLASTN output is expected to be already filtered based on homology (#PM/k).
  If BLASTN was run with --outfmt 6, then such a filter can be easily run
  with a single line AWK command: awk '(\$4-\$5)/k >= threshold'.

 Notes:
  PM: Perfect Match.
  OT: Off Target.
  k: oligomer length
  threshold: homology threshold in percentage

 Mandatory arguments:
  -i stg	Table with oligo_sequence|transcript_ID|Gene_Symbol columns.
  -y gs	Gene Symbol.
  -f fasta_in	Gene fasta file.
  -o fasta_out	Output fasta of sequences surviving the filter.

 Optional arguments:
  -h	Show this help page.
  -k k	Oligo length. Default: 30
  -g gt	Threshold on the number of off-targets gene, for a single oligo.
    	(included) Default: 20
  -s st	Threshold on the number of oligos off-targeting a gene, for the
     	selections of 'saturated' off-target genes. Oligos targeting a
     	saturated off-target are filtered out. (included) Default: 5
"

# Default values
k=30
gene_ot_thr=20
saturation_level=5

# Parse options
while getopts hk:g:s:i:y:f:o: opt; do
	case $opt in
		h)
			# Help page
			echo -e "$helps"
			exit 0
		;;
		k)
			# Oligo length
			if [ $OPTARG -gt 0 ]; then
				k=$OPTARG
			else
				msg="Invalid -k option, k must be greater than 0.\n"
				echo -e "$helps\n!!! ERROR! $msg"
				exit 1
			fi
		;;
		g)
			# OT number threshold
			if [ $OPTARG -gt 0 ]; then
				gene_ot_thr=$OPTARG
			else
				msg="Invalid -g option, gt must be greater than 0.\n"
				echo -e "$helps\n!!! ERROR! $msg"
				exit 1
			fi
		;;
		s)
			# Saturated OT definition
			if [ $OPTARG -gt 0 ]; then
				saturation_level=$OPTARG
			else
				msg="Invalid -s option, st must be greater than 0.\n"
				echo -e "$helps\n!!! ERROR! $msg"
				exit 1
			fi
		;;
		i)
			# Input STG table
			if [ -e $OPTARG ]; then
				stg_path=$OPTARG
			else
				msg="Invalid -i option, file not found.\n File: $OPTARG\n"
				echo -e "$helps\n!!! ERROR! $msg"
				exit 1
			fi
		;;
		y)
			# Gene Symbol
			gene_symbol=$OPTARG
		;;
		f)
			# Input Gene fasta file
			if [ -e $OPTARG ]; then
				fain_path=$OPTARG
			else
				msg="Invalid -f option, file not found.\n File: $OPTARG\n"
				echo -e "$helps\n!!! ERROR! $msg"
				exit 1
			fi
		;;
		o)
			# Output fasta file
			if [ -e $OPTARG ]; then
				msg="Invalid -o option, file already exists.\n File: $OPTARG\n"
				echo -e "$helps\n!!! ERROR! $msg"
				exit 1
			else
				faout_path=$OPTARG
			fi
		;;
	esac
done

# Check mandatory options
if [ -z "$stg_path" ]; then
	echo -e "$helps\n!!! ERROR! Missing mandatory -i option.\n"
	exit 1
fi
if [ -z "$gene_symbol" ]; then
	echo -e "$helps\n!!! ERROR! Missing mandatory -y option.\n"
	exit 1
fi
if [ -z "$fain_path" ]; then
	echo -e "$helps\n!!! ERROR! Missing mandatory -f option.\n"
	exit 1
fi
if [ -z "$faout_path" ]; then
	echo -e "$helps\n!!! ERROR! Missing mandatory -o option.\n"
	exit 1
fi

# Additional checks
# ...

# Print settings
echo -e "
 SETTINGS:

        stg_file : $stg_path
        fasta_in : $fain_path
          output : $faout_path
     gene symbol : $gene_symbol
       
               k : $k nt
         #OT thr : $gene_ot_thr
  Saturation lvl : $saturation_level
"

# RUN ==========================================================================

# Read fasta input
fa_out=$(cat $fain_path | paste - - | sed 's/..//' | sort -k2)

# Count oligos
n_oligo=$(cat "$fain_path" | paste - - | wc -l | cut -d ' ' -f 1)
echo -e " · Found $n_oligo $k-mers..."

# Unique sequences from fasta input
echo -e " · Extracting fasta uniqued sequences ..."
useq_in=$(cat $fain_path | paste - - | cut -f 2 | sort | uniq | sort -k1)
useq=$useq_in
n_uniq=$(echo -e "$useq_in" | wc -l)
echo -e " >>> Found $n_uniq unique $k-mers."

# Extract from STG
echo -e " · Extracting STG lines ..."
stg=$(cat $stg_path | join -j 1 -t $'\t' - <(echo -e "$useq_in"))

# Remove correct targets
echo -e " · Focusing on off-targets..."
stg_ot=$(echo -e "$stg" | grep -v "$gene_symbol")


# OFF-TARGET FILTER #1 ---------------------------------------------------------
# #OT filter

# Count off-targets per sequence
echo -e " · Counting off-targets per sequence..."
awkprg='
{
	if ( $1 in a ) {
		a[$1] = a[$1] + 1
	} else {
		a[$1] = 1
	}
}

END{
	for ( k in a ) {
		if ( a[k] >= ot ) {
			print k
		}
	}
}
'
torm_seq=$(echo -e "$stg_ot" | awk -v ot=$gene_ot_thr "$awkprg" | sort)

if [ -n "$torm_seq" ]; then
	echo -e " >>> Removing sequences..."

	# Remove
	fa_out=$(echo -e "$fa_out" | join -12 -21 -t $'\t' -v1 \
		- <(echo -e "$torm_seq") -o 1.1,1.2 | sort -k2)
	stg_ot=$(echo -e "$stg_ot" | join -v1 -j1 -t$'\t' - <(echo -e "$torm_seq"))

	# Count
	if [ -z "$fa_out" ]; then
		n_kept_oligo=0
	else
		n_kept_oligo=$(echo -e "$fa_out" | wc -l)
	fi
	n_rm_oligo=$(bc <<< "$n_oligo - $n_kept_oligo")
	n_oligo=$n_kept_oligo
	n_rm_seq=$(echo -e "$torm_seq" | wc -l)

	# Log
	echo -e " >>> $n_rm_seq removed sequences had too many off-targets."
	echo -e " >>> $n_rm_oligo removed oligos had too many off-targets."
else
	echo -e " >>> 0 removed sequences had too many off-targets."
fi

# OFF-TARGET FILTER #2 ---------------------------------------------------------
# Saturation filter

if [ -z "$fa_out" ]; then
	echo -e " · Skipping saturated off-target transcripts filter..."
else

	# Count hits per off-target transcript
	echo -e " · Identifying saturated off-target transcripts..."
	awkprg='
	{
		if ( $2 in a ) {
			a[$2] = a[$2] + 1
		} else {
			a[$2] = 1
		}
	}

	END {
		for ( k in a ) {
			if ( a[k] >= st ) {
				print k
			}
		}
	}
	'
	torm_trans=$(echo -e "$stg_ot" | \
		awk -v st=$saturation_level "$awkprg" | sort)
	n_rm_trans=$(echo -e "$torm_trans" | wc -l)
	echo -e " >>> Found $n_rm_trans saturated transcripts."

	# Identify sequences hitting on saturated transcripts
	torm_seq=$(echo -e "$stg_ot" | sort -k2 | \
		join -12 -21 -t$'\t' - <(echo -e "$torm_trans") -o 1.1 | \
		cut -f 1 | sort)

	if [ -n "$torm_seq" ]; then
		echo -e " >>> Removing sequences..."

		# Remove
		fa_out=$(echo "$fa_out" | join -12 -21 -t $'\t' -v1 \
			- <(echo -e "$torm_seq") -o 1.1,1.2)

		# Count
		if [ -z "$fa_out"]; then
			n_kept_oligo=0
		else
			n_kept_oligo=$(echo -e "$fa_out" | wc -l)
		fi
		n_rm_oligo=$(bc <<< "$n_oligo - $n_kept_oligo")
		n_rm_seq=$(echo -e "$torm_seq" | wc -l)

		# Log
		echo -e " >>> Found $n_rm_seq sequences hitting saturated transcripts."
		echo -e " >>> Removed $n_rm_oligo oligos hitting saturated transcripts."
	else
		echo -e " >>> Found 0 sequences hitting saturated transcripts."
	fi
fi

# OUTPUT =======================================================================

if [ -n "$fa_out" ]; then
	# Sort output
	awkprg='
	{
		split($1, ff, "_");
		split(ff[3], oo, ":");
		oi=substr(oo[1], 2);
		print ff[2]"\t"oi"\t"ff[1]"_"ff[2]":("oi + 1"-"oi + k - 1"):"oo[2]":"oo[3]":"oo[4]"\t"$2;
	}
	'
	fa_out=$(echo -e "$fa_out" | awk -v k=$k "$awkprg" | \
		sort -k1,1 -k2,2n | cut -f 3,4)

	# Write output
	n_out_oligo=$(echo -e "$fa_out" | wc -l)
	echo -e " · Writing output ($n_out_oligo oligos)..."
	echo -e "$fa_out" | sed 's/^/>/' | tr '\t' '\n' > $faout_path
else
	echo -e " · No output."
	touch $faout_path
fi

# END ==========================================================================

################################################################################
