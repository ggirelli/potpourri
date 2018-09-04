#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 0.0.1
# Date: 20180904
# Project: GPSeq
# Description: correct translocations on bed files.
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

import argparse
import numpy as np
import os
import pandas as pd
import re
import sys
from tqdm import tqdm

# PARAMETERS ===================================================================

# Add script description
parser = argparse.ArgumentParser(description = '''
Corrects chromosomal rearrangements in bed files. Corrected chromosome sites
format checked with assertions. Only first three columns of the bed files are
expected and manipulated (i.e., BED3). If preasent, header track line is
preserved and translocation site information are added as "transSite1" and
"transSite2". Corrected bed files are exported with suffix ".transCorrected".
''', formatter_class = argparse.RawDescriptionHelpFormatter)

# Add mandatory arguments
parser.add_argument('bedfile',
	type = str, nargs = '+',
	help = 'Bed file(s).')

parser.add_argument('-1', '--first-site',
	type = str, metavar = 'suffix', required = True,
	help = """First translocation site, format chrA:NNNNNN.""")
parser.add_argument('-2', '--second-site',
	type = str, metavar = 'suffix', required = True,
	help = """Second translocation site, format chrB:NNNNNN.""")

# Add arguments with default value
parser.add_argument('-S', '--suffix',
	type = str, metavar = 'suffix', default = "transCorrected",
	help = """Suffix for output bed files. Default: '.transCorrected'""")

# Version flag
version = "0.0.1"
parser.add_argument('--version', action = 'version',
	version = '%s v%s' % (sys.argv[0], version,))

# Parse arguments
args = parser.parse_args()

# FUNCTIONS ====================================================================

# RUN ==========================================================================

if not 0 == len(args.suffix[0]):
	if not "." == args.suffix[0]:
		args.suffix = ".%s" % (args.suffix)

siteRE = "chr[0-9]+:[0-9]+"
assert_msg = "provided first site does not match format."
assert type(None) != type(re.fullmatch(siteRE, args.first_site)), assert_msg
assert_msg = "provided second site does not match format."
assert type(None) != type(re.fullmatch(siteRE, args.second_site)), assert_msg

(chrom1, pos1) = args.first_site.split(":")
(chrom2, pos2) = args.second_site.split(":")
pos1, pos2 = [int(x) for x in [pos1, pos2]]

for bedfile in tqdm(args.bedfile):
	bed_check = open(bedfile, "r")
	headline = bed_check.readlines()[0].strip()
	if headline.startswith("track "):
		bed = pd.read_csv(bedfile, "\t", header = None, skiprows = 1)
	else:
		headline = None
		bed = pd.read_csv(bedfile, "\t", header = None)
	bed_check.close()

	cond  = np.logical_and(bed.iloc[:,0] == chrom1, bed.iloc[:,1] >= pos1)
	bed.loc[cond, 0] = "%s:%s" % (chrom2, chrom1[3:])
	bed.loc[cond, (1, 2)] = bed.loc[cond, (1, 2)] - pos1 + pos2

	cond  = np.logical_and(bed.iloc[:,0] == chrom2, bed.iloc[:,1] >= pos2)
	bed.loc[cond, 0] = "%s:%s" % (chrom1, chrom2[3:])
	bed.loc[cond, (1, 2)] = bed.loc[cond, (1, 2)] - pos2 + pos1

	bed.loc[bed.iloc[:,0] == chrom1, 0] = "%s:%s" % (chrom1, chrom2[3:])
	bed.loc[bed.iloc[:,0] == chrom2, 0] = "%s:%s" % (chrom2, chrom1[3:])

	(basename, ext) = os.path.splitext(bedfile)

	outpath = "%s%s%s" % (basename, args.suffix, ext)
	if type(None) != type(headline):
		headline += ' transSite1="%s" transSite2="%s"\n' % (
			args.first_site, args.second_site)

		with open(outpath, "w") as OH:
			OH.write(headline)

		bed.to_csv(outpath, sep = "\t",
			header = False, index = False, mode = "a")
	else:
		bed.to_csv(outpath, sep = "\t",
			header = False, index = False)

# END ==========================================================================

################################################################################
