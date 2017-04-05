#!/usr/bin/python
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.1
# Description: Shuffle a certain percentage of reads in a bed file.
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

import argparse
import numpy as np
import os
import pandas as pd
import pickle
import sys

pd.options.mode.chained_assignment = None  # default='warn'

# INPUT ========================================================================

# Add script description
parser = argparse.ArgumentParser(
	description = 'Shuffle bed file read counts.'
)

# Add params
parser.add_argument('seed', type = int, nargs = 1,
	help = 'Seed for random number generation.')
parser.add_argument('bedfile', type = str, nargs = 1,
	help = 'Path to bedfile.')

# Add flags
parser.add_argument('-n', metavar = 'nIter', type = int, nargs = 1,
	default = [100], help = 'Number of iterations.')
parser.add_argument('-p', metavar = 'perc', type = int, nargs = 1,
	default = [10], help = 'Percentage of reads to shuffle.')
parser.add_argument('-o', metavar = 'outDir', type = str, nargs = 1,
	default = ['./shuffled/'], help = 'Output directory.')
parser.add_argument('-k', metavar = 'keepSeed', type = bool, nargs = 1,
	default = [True],
	help = 'Reload previous seed state. Use on subsequent runs.')

# Parse arguments
args = parser.parse_args()

# Retrieve arguments
seed = args.seed[0]
bedfile = args.bedfile[0]
nIter = args.n[0]
perc = args.p[0]
outDir = args.o[0]
keepSeed = args.k[0]

# Output file name prefix
outName = '.'.join(bedfile.split('/')[-1].split('.')[:-1])

# Make output directory if missing
if not os.path.isdir(outDir):
	os.makedirs(outDir)

# RUN ==========================================================================

# Log info
print(' · Shuffling x'+str(nIter)+' '+str(perc)+'% of '+bedfile)

# Set seed
seed = np.random.RandomState(seed)

# Load seed if available
fname = outDir + '/.seed_state.pickle'
if os.path.isfile(fname):
	f = open(fname, 'r')
	seed.set_state(pickle.load(f))
	f.close()

# Read bedfile -----------------------------------------------------------------

bf = pd.read_csv(bedfile, '\t', skiprows = [0],
	names = ['chr', 'start', 'end', 'name', 'score'])

# Shuffle ----------------------------------------------------------------------

# Count reads
nreads = sum(bf['score'])
toShuffle = nreads * perc / 100

print(' >>> Found ' + str(nreads) + ' reads.')


print(' >>> Pre-shuffling...')
preshuffle = np.repeat(np.arange(len(bf['score'])), bf['score'])

for i in range(nIter):
	print(' >>># Iteration #' + str(i+1))
	print(' >>># Shuffle...')
	pos_from = seed.randint(0, len(preshuffle), toShuffle)
	pos_to = seed.randint(0, len(bf['score']), toShuffle)

	print(' >>># Moving reads...')
	shuffled = np.array(preshuffle)
	shuffled[pos_from] = pos_to

	print(' >>># Counting shuffled reads')
	counts = np.unique(shuffled, return_counts = True)
	shuffled = bf.copy()
	shuffled['score'] = np.zeros(shuffled['score'].shape)
	shuffled['score'][counts[0]] = counts[1]

	# Output
	shuffled.to_csv(outDir+outName+'.iter'+str(i+1)+'.'+str(perc)+'perc.bed',
		sep = '\t', header = False, index = False)

# Saving seed state
seed_state = seed.get_state()
f = open(outDir + '/.seed_state.pickle', 'w+')
pickle.dump(seed_state, f)
f.close()

# END --------------------------------------------------------------------------

################################################################################
