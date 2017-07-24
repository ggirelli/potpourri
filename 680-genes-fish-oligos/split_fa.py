#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.0
# Date: 20170722
# Project: 680 genes
# Description: split fasta based on header pattern
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

import argparse
import os
import progressbar

# PARAMETERS ===================================================================


# Add script description
parser = argparse.ArgumentParser(
	description = 'Split fasta file based on header pattern'
)

# Add mandatory arguments
parser.add_argument('inFasta', type = str, nargs = 1,
	help = 'Input fasta file.')
parser.add_argument('outdir', type = str, nargs = 1,
	help = 'Output folder.')

# Add arguments with default value
parser.add_argument('-d', type = str, nargs = 1,
	metavar = 'delim', help = """
	Delimiter. Default: '_'""", default = ["_"])
parser.add_argument('-f', type = int, nargs = 1,
	metavar = 'field', help = """
	0-indexed field ID. Default: 0""", default = [0])


# Add flags
parser.add_argument('-o',
	action = 'store_const', dest = 'only_once',
	const = True, default = False,
	help = 'Write output once, instead of appending.')

# Parse arguments
args = parser.parse_args()

# Assign to in-script variables
fain_path = args.inFasta[0]
outdir = args.outdir[0]
delim = args.d[0]
field = args.f[0]
only_once = args.only_once

# Create outdir if it does not exist
if not os.path.isdir(outdir):
	os.mkdir(outdir)

# FUNCTIONS ====================================================================

def file_nrow(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return(i + 1)

def split_fa(head, seq, d, only_once):
	k = head.split(delim)[field][1:].strip()
	if k in d.keys():
		d[k] += "%s\n%s\n" % (head, seq,)

		if not only_once:
			fout = open("%s/%s.fa" % (outdir, k), 'a+')
			fout.write("%s\n%s\n" % (head, seq,))
			fout.close()
	else:
		d[k] = "%s\n%s\n" % (head, seq,)

		if not only_once:
			fout = open("%s/%s.fa" % (outdir, k), 'w+')
			fout.write("%s\n%s\n" % (head, seq,))
			fout.close()

	return(d)

# RUN ==========================================================================

d = {}
curr_head = ""
curr_seq = ""

bar = progressbar.ProgressBar(max_value = file_nrow(fain_path))
with open(fain_path, 'r') as fain:
	i = 0
	for row in fain:
		bar.update(i)
		i += 1

		if '>' == row[0]:
			# Run split
			d = split_fa(curr_head, curr_seq, d, only_once)

			# Reset seq
			curr_seq = ""

			# Read head
			curr_head = row.strip()
		else:
			curr_seq += row.strip()
fain.close()

# Run last item
d = split_fa(curr_head, curr_seq, d, only_once)

if only_once:
	# Write final output
	for k in d.keys():
		fout = open("%s/%s.fa" % (outdir, k), 'w+')
		fout.write(d[k])
		fout.close()

# END ==========================================================================

################################################################################
