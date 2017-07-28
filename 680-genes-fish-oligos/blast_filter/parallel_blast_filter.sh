#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: ...
# Email: ...
# Version: X.X.X
# Date: YYYYMMDD
# Project: ...
# Description: ...
# 
# ------------------------------------------------------------------------------



# ENV VAR ======================================================================

export LC_ALL=C

# PARAMS =======================================================================


# Help string
helps="
usage: ./parallel_blast_filter.sh [-h][-k k][-g gt][-s st]
                                  -i stg -y gs -f indir -o outdir

 Description:
  Run blast_filter.sh in parallel on every fasta of indir.

 Mandatory arguments:
  -i stg	Table with oligo_sequence|transcript_ID|Gene_Symbol columns.
  -y gs	Table with Gene_Symbol|gene_ID columns.
  -f indir	Folder with gene fasta file.
  -o outdir	Output folder.

 Optional arguments:
  -h	Show this help page.
  -k k	Oligo length. Default: 30
  -g gt	Threshold on the number of off-targets gene, for a single oligo.
    	Default: 20
  -s st	Threshold on the number of oligos off-targeting a gene, for the
     	selections of 'saturated' off-target genes. Oligos targeting a
     	saturated off-target are filtered out. Default: 5
  -t threads	Number of threads for parallelization. Default: 1
"

# Default values
k=30
gene_ot_thr=20
saturation_level=5
threads=1

# Parse options
while getopts hk:g:s:i:y:f:o:t: opt; do
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
			# Gene_ID|Gene_Symbol table
			if [ -e $OPTARG ]; then
				gst_path=$OPTARG
			else
				msg="Invalid -y option, file not found.\n File: $OPTARG\n"
				echo -e "$helps\n!!! ERROR! $msg"
				exit 1
			fi
		;;
		f)
			# Input folder
			if [ -e $OPTARG ]; then
				fin_path=$OPTARG
			else
				msg="Invalid -f option, file not found.\n File: $OPTARG\n"
				echo -e "$helps\n!!! ERROR! $msg"
				exit 1
			fi
		;;
		o)
			# Output folder
			if [ ! -e $OPTARG ]; then
				msg="Output folder created."
				mkdir -p $OPTARG
			else
				if [ ! -d $OPTARG ]; then
					msg="Invalid -o option, a file exists with the folder name."
					echo -e "$helps\n!!! ERROR! $msg\n"
					exit 1
				fi
			fi
			fout_path=$OPTARG
		;;
		t)
			# Number of threads for parallelization
			if [[ 0 -lt $OPTARG ]]; then
				threads=$OPTARG
			fi
		;;
	esac
done

# Check mandatory options
if [ -z "$stg_path" ]; then
	echo -e "$helps\n!!! ERROR! Missing mandatory -i option.\n"
	exit 1
fi
if [ -z "$gst_path" ]; then
	echo -e "$helps\n!!! ERROR! Missing mandatory -y option.\n"
	exit 1
fi
if [ -z "$fin_path" ]; then
	echo -e "$helps\n!!! ERROR! Missing mandatory -f option.\n"
	exit 1
fi
if [ -z "$fout_path" ]; then
	echo -e "$helps\n!!! ERROR! Missing mandatory -o option.\n"
	exit 1
fi


# Additional checks
# ...

# Print settings
echo -e "
 SETTINGS:

         stg_file : $stg_path
            indir : $fin_path
           outdir : $fout_path
gene symbol table : $gst_path
       
                k : $k nt
          #OT thr : $gene_ot_thr
   Saturation lvl : $saturation_level

          threads : $threads
"

# FUNCTIONS ====================================================================

function blast_filter() {
	# Usage:
	# 	blast_filter fain_path faout_path gene_symbol stg_path
	fain_path=$1
	faout_path=$2
	gene_symbol=$3
	stg_path=$4
	logpath=$2".log"
	gene_ot_thr=$5
	saturation_level=$6
	k=$7

	# Keep track
	awkprg='{ split($NF, ff, "."); print ff[1]; }'
	echo -e "$fain_path" | tr '/' '\t' | awk "$awkprg"

	# Log Gene name
	echo -e " · $gene_symbol" > $logpath

	# Read fasta input
	fa_out=$(cat $fain_path | paste - - | sed 's/..//' | sort -k2)

	# Count oligos
	n_oligo=$(cat "$fain_path" | paste - - | wc -l | cut -d ' ' -f 1)
	echo -e " · Found $n_oligo $k-mers..." >> $logpath

	# Unique sequences from fasta input
	echo -e " · Extracting fasta uniqued sequences ..." >> $logpath
	useq_in=$(cat $fain_path | paste - - | cut -f 2 | sort | uniq | sort -k1)
	useq=$useq_in
	n_uniq=$(echo -e "$useq_in" | wc -l)
	echo -e " >>> Found $n_uniq unique $k-mers." >> $logpath

	# Extract from STG
	echo -e " · Extracting STG lines ..." >> $logpath
	stg=$(cat $stg_path | join -j 1 -t $'\t' - <(echo -e "$useq_in"))

	# Remove correct targets
	echo -e " · Focusing on off-targets..." >> $logpath
	stg_ot=$(echo -e "$stg" | grep -v "$gene_symbol")

	# OFF-TARGET FILTER #1 -----------------------------------------------------
	# #OT filter

	# Count off-targets per sequence
	echo -e " · Counting off-targets per sequence..." >> $logpath
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
		echo -e " >>> Removing sequences..." >> $logpath

		# Remove
		fa_out=$(echo -e "$fa_out" | join -12 -21 -t $'\t' -v1 \
			- <(echo -e "$torm_seq") -o 1.1,1.2 | sort -k2)
		stg_ot=$(echo -e "$stg_ot" | \
			join -v1 -j1 -t$'\t' - <(echo -e "$torm_seq"))

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
		msg=" >>> $n_rm_seq removed sequences had too many off-targets."
		echo -e "$msg" >> $logpath
		msg=" >>> $n_rm_oligo removed oligos had too many off-targets."
		echo -e "$msg" >> $logpath
	else
		echo -e " >>> 0 removed sequences had too many off-targets." >> $logpath
	fi

	# OFF-TARGET FILTER #2 -----------------------------------------------------
	# Saturation filter

	if [ -z "$fa_out" ]; then
		echo -e " · Skipping saturated off-target transcripts filter..." >> $logpath
	else
		# Count hits per off-target transcript
		echo -e " · Identifying saturated off-target transcripts..." >> $logpath
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
		echo -e " >>> Found $n_rm_trans saturated transcripts." >> $logpath

		# Identify sequences hitting on saturated transcripts
		torm_seq=$(echo -e "$stg_ot" | sort -k2 | \
			join -12 -21 -t$'\t' - <(echo -e "$torm_trans") -o 1.1 | \
				cut -f 1 | sort)

		if [ -n "$torm_seq" ]; then
			echo -e " >>> Removing sequences..." >> $logpath

			# Remove
			fa_out=$(echo "$fa_out" | join -12 -21 -t $'\t' -v1 \
				- <(echo -e "$torm_seq") -o 1.1,1.2)

			# Count
			n_kept_oligo=$(echo -e "$fa_out" | wc -l)
			n_rm_oligo=$(bc <<< "$n_oligo - $n_kept_oligo")
			n_rm_seq=$(echo -e "$torm_seq" | wc -l)

			# Log
			msg=" >>> Found $n_rm_seq sequences hitting saturated transcripts."
			echo -e "$msg" >> $logpath
			msg=" >>> Removed $n_rm_oligo oligos hitting saturated transcripts."
			echo -e "$msg" >> $logpath
		else
			msg=" >>> Found 0 sequences hitting saturated transcripts."
			echo -e "$msg" >> $logpath
		fi
	fi
	
	# OUTPUT ===================================================================

	if [ -n "$fa_out" ]; then

		# Sort output
		awkprg='
		{
			split($1, ff, "_");
			split(ff[3], oo, ":");
			oi=substr(oo[1], 2);
			head=ff[2]"\t"oi"\t"ff[1]"_"ff[2]":("oi + 1"-"oi + k"):"oo[2];
			head=head":"oo[3]":"oo[4]"\t"$2;
			print head;
		}
		'
		fa_out=$(echo -e "$fa_out" | awk -v k=$k "$awkprg" | \
			sort -k1,1 -k2,2n | cut -f 3,4)
		#fa_out=$(echo -e "$fa_out" | sort -k1,1 -k2,2n)

		# Write output
		n_out_oligo=$(echo -e "$fa_out" | wc -l)
		echo -e " · Writing output ($n_out_oligo oligos)..." >> $logpath
		echo -e "$fa_out" | sed 's/^/>/' | tr '\t' '\n' > $faout_path
	else
		echo -e " · No output." >> $logpath
		touch $faout_path
	fi
	
}
export -f blast_filter

# RUN ==========================================================================

echo -e " · Preparing jobs..."
args=()
for f in $(ls "$fin_path"/*.fa); do
	f=$(echo "$f" | tr '/' '\t' | awk '{ print $NF; }')
	# Add input fasta
	args+=($fin_path"/"$f)
	# Add output fasta
	args+=($fout_path"/"$f)
	# Gene symbol
	ensg=$(echo "$f" | cut -d '.' -f 1)
	gs=$(cat $gst_path | grep "$ensg" | cut -f 1)
	args+=($gs)
	# STG path
	args+=($stg_path)
	args+=($gene_ot_thr)
	args+=($saturation_level)
	args+=($k)
done

echo -e " · Submitting jobs..."
echo ${args[@]} | tr ' ' '\t' | \
	parallel --jobs $threads -d $'\t' -n 7 blast_filter

echo -e " ~ DONE ~"
# END ==========================================================================

################################################################################
