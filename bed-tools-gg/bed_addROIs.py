#!/usr/bin/python
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.1
# Description: 	Assigns rows in a bed file to a given list of regions of
# 				interest (ROIs). The ROIs can be overlapping. A new column is
# 				added to the end of the bed file, with all the regions
# 				containing it, in the chr:start-end format, space-separated.
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

import argparse
import numpy as np
import pandas as pd

pd.options.mode.chained_assignment = None  # default='warn'

# INPUT ========================================================================

# Add script description
parser = argparse.ArgumentParser(
	description = """Assigns rows in a bed file to a given list of regions of
		interest (ROIs). The ROIs can be overlapping. A new column is
		added to the end of the bed file, with all the regions
		containing it, in the chr:start-end format, space-separated."""
)

# Add params
parser.add_argument('regfile', type = str, nargs = 1,
	help = 'Path to bedfile, containing regions to be assigned to.')
parser.add_argument('bedfile', type = str, nargs = 1,
	help = 'Path to bedfile, containing rows to be assigned.')

# Add flags
parser.add_argument('-u',
	action = 'store_const', const = True, default = False,
	help = 'Keep bedfile rows that do not match any region.')
parser.add_argument('-m',
	action = 'store_const', const = True, default = False,
	help = 'Assign to bedfile rows that partially match a region.')
parser.add_argument('-l',
	action = 'store_const', const = True, default = False,
	help = 'Assign to bedfile rows that include a region.')
parser.add_argument('-o', metavar = 'outfile', type = str, nargs = 1,
	default = [False],
	help = 'Output file (not a bed). Output to stdout if not specified.')

# Parse arguments
args = parser.parse_args()

# Retrieve arguments
regfile = args.regfile[0]
bedfile = args.bedfile[0]
keep_unassigned_rows = args.u
keep_marginal_overlaps = args.m
keep_including = args.l
outfile = args.o[0]

# Default variables
bedcolnames = ['chr', 'start', 'end', 'name', 'score']

# FUNCTIONS ====================================================================

def join_trim(l):
	return ' '.join(' '.join(l).split())

# RUN ==========================================================================

# Read regions file
rois = pd.read_csv(regfile, '\t', names = bedcolnames)

# Read bed file
bed = pd.read_csv(bedfile, '\t', names = bedcolnames, skiprows = [0])

# Add regions column
bed['rois'] = pd.Series(np.array(['' for i in range(bed.shape[0])]),
	index = bed.index)

# Assign reads to rows
chr_set = set(bed['chr'])
for chri in chr_set:
	# Select rois and rows
	chr_bed = bed.iloc[np.where(bed['chr'] == chri)[0], :]
	chr_roi = rois.iloc[np.where(rois['chr'] == chri)[0], :]

	if 0 == chr_bed.shape[0] or 0 == chr_roi.shape[0]:
		continue

	# Prepare roi label
	chr_roi_labels = np.array(chr_roi['chr']).astype('str').tolist()
	chr_roi_labels = np.core.defchararray.add(chr_roi_labels, ':')
	chr_roi_labels = np.core.defchararray.add(chr_roi_labels,
		np.array(chr_roi['start']).astype('str').tolist())
	chr_roi_labels = np.core.defchararray.add(chr_roi_labels, '-')
	chr_roi_labels = np.core.defchararray.add(chr_roi_labels,
		np.array(chr_roi['end']).astype('str').tolist())

	# Build matrix -------------------------------------------------------------
	bed_hash_start = [hash(i) for i in chr_bed['start']]
	bed_hash_end = [hash(i) for i in chr_bed['end']]
	roi_hash_start = [hash(i) for i in chr_roi['start']]
	roi_hash_end = [hash(i) for i in chr_roi['end']]

	# Start should be higher than the region start
	condition_start = np.greater_equal.outer(bed_hash_start, roi_hash_start)

	# End should be lower than the region end
	condition_end = np.logical_not(np.greater.outer(bed_hash_end, roi_hash_end))

	# Perfectly contained (in)
	condition_in = np.logical_and(condition_start, condition_end)

	if keep_marginal_overlaps or keep_including:
		# Start should be lower than the region start
		condition_left_start = np.logical_not(condition_start)

		# End should be higher than the region end
		condition_right_end = np.logical_not(condition_end)

	if keep_marginal_overlaps:
		# End should be lower than the region end and higher than its start
		condition_left_end = np.logical_and(condition_end,
			np.greater_equal.outer(bed_hash_end, roi_hash_start))

		# Start should be higher than the region start and lower than its end
		condition_right_start = np.logical_and(condition_start,
			np.logical_not(np.greater.outer(bed_hash_start, roi_hash_end)))

		# Partial overlap on left margin
		condition_left = np.logical_and(
			condition_left_start, condition_left_end)

		# Partial overlap on right margin
		condition_right = np.logical_and(
			condition_right_start, condition_right_end)

		# Partial overla on a margin
		condition_margins = np.logical_or(condition_left, condition_right)

		# Partial overlap or inside
		condition_in = np.logical_or(condition_in, condition_margins)

	if keep_including:
		# Rows that include the region
		condition_larger = np.logical_and(
			condition_left_start, condition_right_end)

		# Included (inside) or including
		condition_in = np.logical_or(condition_in, condition_larger)

	if 0 == condition_in.sum():
		continue

	# Add rois per row
	labels = np.tile(chr_roi_labels, (condition_in.shape[0], 1))
	labels[np.logical_not(condition_in)] = ''
	bed['rois'][chr_bed.index.values] = [join_trim(labels[rowi, :])
		for rowi in range(labels.shape[0])]

# Remove rows without regions
if not keep_unassigned_rows:
	bed = bed[bed['rois'] != '']

# Output
if False == outfile:
	for i in range(bed.shape[0]):
		print('\t'.join(bed.iloc[i, :].astype('str').tolist()))
else:
	bed.to_csv(outfile, sep = '\t', header = False, index = False)

# END --------------------------------------------------------------------------

################################################################################
