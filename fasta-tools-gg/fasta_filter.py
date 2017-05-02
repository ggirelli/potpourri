#!/usr/bin/python
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.1
# Description: 	Select sequences that match a header regexp.
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

import argparse
import numpy as np
import pandas as pd
import re

pd.options.mode.chained_assignment = None  # default='warn'

# INPUT ========================================================================

# Add script description
parser = argparse.ArgumentParser(
	description = """
	Select sequences from the input FASTA file that match a given regular
	expression in their header.
	"""
)

# Add params
parser.add_argument('FASTA', type = str, nargs = 1,
	help = 'Path to a fasta file.')
parser.add_argument('regexp', type = str, nargs = 1,
	help = 'Regular expression to match to the header.')

# Parse arguments
args = parser.parse_args()

# Retrieve arguments
fasta_path = args.FASTA[0]
regexp = re.compile(args.regexp[0])

# FUNCTIONS ====================================================================

# RUN ==========================================================================

skip = False
current_key = np.nan

with open(fasta_path, 'r') as fasta_file:
	for row in fasta_file:
		if row.startswith('>'):
			current_key = row
			if regexp.search(current_key) is None:
				skip = True
				continue
			else:
				skip = False
				print row.strip()
				continue
		else:
			if not skip:
				print row.strip()
				continue

# END --------------------------------------------------------------------------

################################################################################
