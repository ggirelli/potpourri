#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.1
# Description: 	convert IDs in the specified field of a GTF attribute column.
# 
# ------------------------------------------------------------------------------



# ENV VAR ======================================================================

export LC_ALL=C

# DEPENDENCIES =================================================================

function join_by { local IFS="$1"; shift; echo "$*"; }

# INPUT ========================================================================

# Help string
helps="
 usage: ./gtf_convert_id.sh [-hvr][-f field] -g gtf -c conv

 Description:
  Convert IDs in the specified field of a GTF attribute column. Requires a GTF
  file and a conversion file obtained from bioDBnet.
  bioDBnet: https://biodbnet-abcc.ncifcrf.gov/db/dbOrg.php

 Mandatory arguments:
  -g gtf	GTF file.
  -c conv	Conversion file, obtained from bioDBnet.

 Optional arguments:
  -h		Show this help page.
  -v		Skip rows without the specified attribute field.
  -r		Skip rows with no conversion ID.
  -f field	Label of the attribute column field containing the ID to convert.
  		Default: ref_gene_id
"

# Default options
field="ref_gene_id"
skipKey=false
skipVal=false

# Parse options
while getopts hvrf:g:c: opt "${bedfiles[@]}"; do
	case $opt in
		h)
			echo -e "$helps\n"
			exit 0
		;;
		v)
			skipKey=true
		;;
		r)
			skipVal=true
		;;
		f)
			field="$OPTARG"
		;;
		g)
			if [ -e $OPTARG ]; then
				gtf=$OPTARG
			else
				msg="!!! Invalid GTF file, file not found.\n    File: $bf"
				echo -e " $helps\n$msg"
				exit 1
			fi
		;;
		c)
			if [ -e $OPTARG ]; then
				conv=$OPTARG
			else
				msg="!!! Invalid conversion file, file not found."
				msg=$msg"\n    File: $bf"
				echo -e " $helps\n$msg"
				exit 1
			fi
		;;
	esac
done

# Check mandatory options
if [ -z "$gtf" ]; then
	msg="!!! No GTF file provided."
	echo -e "$helps\n$msg"
	exit 1
fi
if [ -z "$conv" ]; then
	msg="!!! No conversion file provided."
	echo -e "$helps\n$msg"
	exit 1
fi

# TEST =========================================================================

# RUN ==========================================================================

conversion_program='@include "join"
BEGIN {
	# Set OFS and FS
	OFS=FS="\t";
}

# Make a conversion array (ca) based on conversion file
(FNR == NR) {
	# Remove leading space from conversion target
	gsub(/^ /, "", $2);

	# Store key,val couple
	ca[$1]=$2;

	# Skip rest of the program
	next;
}

# Navigate  GTF file
{
	# Split attributes
	nattrs=split($9, attributes, "; ");

	# Define empty key and value counters
	k=0;
	v=0;

	# Check the attributes
	for ( i = 0; i < nattrs; i++ ) {

		# Identify attribute key and value
		split(attributes[i], cattr, " ");

		# If this is the specified attribute (-f)
		if ( cattr[1] == lab ) {

			# Increase key counter
			k++;

			# Remove double quotes from value
			gsub(/\"/, "", cattr[2]);

			# Check if value can be converted
			if ( cattr[2] in ca ) {
				# If conversion target is not empty, increase value counter
				if ( "-" != ca[cattr[2]] )
					v++;

				# Convert
				cattr[2]="\""ca[cattr[2]]"\"";
			} else {
				# No target conversion
				cattr[2]="\"-\"";
			}
		}

		# Merge updated attributes
		attributes[i]=join(cattr, 1, 2, " ");
	}

	# Update attributes in output line
	$9=join(attributes, 1, nattrs, "; ");

	# Output based on skip options
	if ( "false" == skipKey || 0 != k )
		if ( "false" == skipVal || 0 != v )
			print $0;
}
'

# Run
awk -v lab="$field" -v skipKey="$skipKey" -v skipVal="$skipVal" \
	"$conversion_program" <(cat $conv) <(cat $gtf)

# End --------------------------------------------------------------------------

################################################################################
