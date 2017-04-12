#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.0
# Description:
# 	Repeat a bedfile row as many times as specified in the score column.
# 
# ------------------------------------------------------------------------------



# ENV VAR ======================================================================

export LC_ALL=C

# DEPENDENCIES =================================================================

function join_by { local IFS="$1"; shift; echo "$*"; }

# INPUT ========================================================================

# Help string
helps="
 usage: ./bed_rep.sh [-h][-c colID][-d del] -b bedfile

 Description:
  Repeat a bedfile row as many times as specified in the score column.
  Following the bed format, the score column should be the 5th column.
  If not, specify the index of the column using the -c option.

 Mandatory arguments:
  -b bedfile	Bed file.

 Optional arguments:
  -h		Show this help page.
  -c colID	Column index (1-indexed). [Default: 5]
  -d del	Column delimiter. [Default: TAB]
"

# Default options
colID=5
delimiter="\t"

# Parse options
while getopts hc:d:b: opt "${bedfiles[@]}"; do
	case $opt in
		h)
			echo -e "$helps\n"
			exit 0
		;;
		c)
			colID=$OPTARG
		;;
		b)
			if [ -e $OPTARG ]; then
				bedfile=$OPTARG
			else
				msg="!!! Invalid -b option, file not found.\n    File: $OPTARG"
				echo -e "$helps\n$msg"
				exit 1
			fi
		;;
	esac
done

# Check mandatory options
if [ -z "$bedfile" ]; then
	msg="!!! Missing mandatory -b option."
	echo -e "$helps\n$msg"
	exit
fi

# RUN ==========================================================================

awkprogram='{
	OFS=FS=del
	for (i=1; i <= $coli; i++) {
		print $1 OFS $2 OFS $3
	}
}'
cat $bedfile | awk -v del="$delimiter" -v coli="$colID" "$awkprogram"

# End --------------------------------------------------------------------------

################################################################################
