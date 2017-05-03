#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.1
# Description: 	copy an attribute value over another
# 
# ------------------------------------------------------------------------------



# ENV VAR ======================================================================

export LC_ALL=C

# DEPENDENCIES =================================================================

function join_by { local IFS="$1"; shift; echo "$*"; }

# INPUT ========================================================================

# Help string
helps="
 usage: ./gtf_copy_attr.sh [-h] -g gtf -1 from -2 to

 Description:
  Overwrites the value of attribute (-2) with the value of attribute (-1)

 Mandatory arguments:
  -g gtf	GTF file.
  -1 from	Attribute value to copy from.
  -2 to	Attribute value to copy to.

 Optional arguments:
  -h		Show this help page.
"

# Default options
field="ref_gene_id"

# Parse options
while getopts hg:1:2: opt "${bedfiles[@]}"; do
	case $opt in
		h)
			echo -e "$helps\n"
			exit 0
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
		1)
			from=$OPTARG
		;;
		2)
			to=$OPTARG
		;;
	esac
done

# Check mandatory options
if [ -z "$gtf" ]; then
	msg="!!! No GTF file provided: missing mandatory (-g) option."
	echo -e "$helps\n$msg"
	exit 1
fi
if [ -z "$from" ]; then
	msg="!!! Missing mandatory (-1) option."
	echo -e "$helps\n$msg"
	exit 1
fi
if [ -z "$to" ]; then
	msg="!!! Missing mandatory (-2) option."
	echo -e "$helps\n$msg"
	exit 1
fi

# TEST =========================================================================

# RUN ==========================================================================

overwrite_program='@include "join"
BEGIN {
	# Set OFS and FS
	OFS=FS="\t";
}

# Navigate  GTF file
{
	# Split attributes
	nattrs=split($9, attributes, "; ");

	# Set default ID fields and counter
	from_id = -1;
	to_id = -1;
	c = 0;

	# Identify from and to attributes
	for ( i = 0; i < nattrs; i++ ) {

		# Identify attribute key and value
		split(attributes[i], cattr, " ");

		# If this is the specified attribute (-1)
		if ( cattr[1] == from ) {
			from_id = i;
			c++;
		}

		# If this is the specified attribute (-2)
		if ( cattr[1] == to ) {
			to_id = i;
			c++;
		}
	}

	# Overwrite value
	if ( 2 == c ) {
		# Identify attribute key and value (-1)
		split(attributes[from_id], from_attr, " ");

		# Identify attribute key and value (-2)
		split(attributes[to_id], to_attr, " ");

		to_attr[2] = from_attr[2];
		attributes[to_id] = join(to_attr, 1, 2, " ");
	}

	# Update attributes in output line
	$9=join(attributes, 1, nattrs, "; ");

	# Output
	print $0;
}
'

# Run
awk -v from="$from" -v to="$to" "$overwrite_program" <(cat $gtf)

# End --------------------------------------------------------------------------

################################################################################
