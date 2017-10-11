#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.0
# Date: 20171011
# Project: Oligo analysis
# Description: calculate possible secondary structure of provided nucleic acid.
# 
# Rationale:
# 	- Convert provided sequence to numeric (A=0, U|T=1, C=2, G=3).
# 	- Use sums to identify correct matches.
# 	- Identify locations (indices) of each differnt base.
# 	- Save matching base-pairs.
# 	- Identify matching base-pairs that open-close a stem based on its
# 	  neighbours.
# 	- Aggregate matching base-pairs into stems, do not allow for sharp turns.
# 	- Provide all possible combinations of stems.
# 	- Calculate Tm and dG for all secondary structure combinations.
# 	- Produce 'graphical' representation.
# 	- Output.
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

import argparse

# PARAMETERS ===================================================================


# Add script description
parser = argparse.ArgumentParser(description = '''
	Calculate melting temperature and delta free energy of possible secondary
	structure conformation of the provided nucleic acid.
''')

# Add mandatory arguments
parser.add_argument('seq', type = str, nargs = 1,
	help = '''ssDNA or ssRNA sequence.''')

# Add arguments with default value
parser.add_argument('-o', '--oconc', metavar = "oligo_conc",
	type = float, nargs = 1,
	help = '''Oligonucleotide concentration [M].
	Default: 0.25e-6 M''',
	default = [0.25e-6])
parser.add_argument('-n', '--naconc', metavar = "na_conc",
	type = float, nargs = 1,
	help = '''Na+ concentration [M].
	Default: 50e-3 M''',
	default = [50e-3])
parser.add_argument('-m', '--mgconc', metavar = "mg_conc",
	type = float, nargs = 1,
	help = '''Mg2+ concentration [M].
	Default: 0 M''',
	default = [0])
parser.add_argument('-t', '--temperature', metavar = "temp",
	type = float, nargs = 1,
	help = '''Temperature in Celsius.
	Default: 37 degC''',
	default = [37.])
parser.add_argument('-T', '--thermo', type = str, nargs = 1,
	help = '''Thermodynamic table to use in the calculations.
	Possible values: allawi (based on ref.2, default) or 
	freier (based on ref.1).''',
	choices = ['allawi', 'freier'], default = ['allawi'])

# Add flags
parser.add_argument('-C', '--celsius',
	dest = 'celsius', action = 'store_const',
	const = True, default = False,
	help = 'Output temperature in Celsius degrees. Default: Kelvin')
parser.add_argument('-v', '--verbose',
	dest = 'verbose', action = 'store_const',
	const = True, default = False,
	help = 'Verbose output.')

# Parse arguments
args = parser.parse_args()

# Assign to in-script variables ------------------------------------------------

# Input oligo
seq = args.seq[0]

# Concentrations
oligo_conc = args.oconc[0]
na_conc = args.naconc[0]
mg_conc = args.mgconc[0]
temperature = args.temperature[0]

# Thermodynamic table
tt_mode = args.thermo[0]

# Temperature units
celsius = args.celsius

# Verbose mode
is_verbose = args.verbose


# Constants --------------------------------------------------------------------

# Gas constant
R = 1.987 / 1000	# kcal / (K mol)

# Thermodynamic tables ---------------------------------------------------------

# Table from Freier et al, PNAS(83), 1986 - in 1 M NaCl [RNA]
# dH0: kcal / mol
# dS0: eu = cal / (K mol)
# dG0: kcal / mol
freier = {
	#				 dH0	 dS0 	 dG0
	'AA'		:	(-6.6,	-18.4,	-0.9),
	'UU'		:	(-6.6,	-18.4,	-0.9),
	'TT'		:	(-6.6,	-18.4,	-0.9),
	'AU'		:	(-5.7,	-15.5,	-0.9),
	'AT'		:	(-5.7,	-15.5,	-0.9),
	'UA'		:	(-8.1,	-22.6,	-1.1),
	'TA'		:	(-8.1,	-22.6,	-1.1),
	'CA'		:	(-10.5,	-27.8,	-1.8),
	'UG'		:	(-10.5,	-27.8,	-1.8),
	'TG'		:	(-10.5,	-27.8,	-1.8),
	'CU'		:	(-7.6,	-19.2,	-1.7),
	'CT'		:	(-7.6,	-19.2,	-1.7),
	'AG'		:	(-7.6,	-19.2,	-1.7),
	'GA'		:	(-13.3,	-35.5,	-2.3),
	'UC'		:	(-13.3,	-35.5,	-2.3),
	'TC'		:	(-13.3,	-35.5,	-2.3),
	'GU'		:	(-10.2,	-26.2,	-2.1),
	'GT'		:	(-10.2,	-26.2,	-2.1),
	'AC'		:	(-10.2,	-26.2,	-2.1),
	'CG'		:	(-8.0,	-19.4,	-2.0),
	'GC'		:	(-14.2,	-34.9,	-3.4),
	'GG'		:	(-12.2,	-29.7,	-2.9),
	'CC'		:	(-12.2,	-29.7,	-2.9),
	'endA'		:	(0,		-10.8,	3.4),
	'endC'		:	(0,		-10.8,	3.4),
	'endT'		:	(0,		-10.8,	3.4),
	'endG'		:	(0,		-10.8,	3.4),
	'endU'		:	(0,		-10.8,	3.4),
	'has_end'	:	True,
	'sym'		:	(0,		-1.4,	0.4),
	'has_sym'	:	True
}

# Table from Allawi&Santalucia, Biochemistry(36), 1997 - in 1 M NaCl [DNA]
# dH0: kcal / mol
# dS0: eu = cal / (K mol)
# dG0: kcal / mol
allawi = {
	#				 dH0	 dS0 	 dG0
	'AA'		:	(-7.9,	-22.2,	-1.0),
	'TT'		:	(-7.9,	-22.2,	-1.0),
	'AT'		:	(-7.2,	-20.4,	-0.88),
	'TA'		:	(-7.2,	-21.3,	-0.58),
	'CA'		:	(-8.5,	-22.7,	-1.45),
	'TG'		:	(-8.5,	-22.7,	-1.45),
	'GT'		:	(-8.4,	-22.4,	-1.44),
	'AC'		:	(-8.4,	-22.4,	-1.44),
	'CT'		:	(-7.8,	-21.0,	-1.28),
	'AG'		:	(-7.8,	-21.0,	-1.28),
	'GA'		:	(-8.2,	-22.2,	-1.3),
	'TC'		:	(-8.2,	-22.2,	-1.3),
	'CG'		:	(-10.6,	-27.2,	-2.17),
	'GC'		:	(-9.8,	-24.4,	-2.24),
	'GG'		:	(-8.0,	-19.9,	-1.84),
	'CC'		:	(-8.0,	-19.9,	-1.84),
	'endG'		:	(.1,	-2.8,	.98),
	'endC'		:	(.1,	-2.8,	.98),
	'endA'		:	(2.3,	4.1,	1.03),
	'endT'		:	(2.3,	4.1,	1.03),
	'has_end'	:	True,
	'sym'		:	(2.3,	4.1,	1.03),
	'has_sym'	:	True
}

# Save tables
tt = {
	'freier' 	: freier,
	'allawi' 	: allawi
}

# FUNCTIONS ====================================================================

def rc(na, t):
	'''
	Generates the reverse complement of a nucleic acid string.

	Args:
		na (string): nucleic acid sequence.
		t (string): nucleic acid type, either 'dna' or 'rna'.

	Return:
		string: reverse complement of na.
	'''

	# Identify type
	t = t.lower()

	# Select alphabet
	if t == 'dna':
		ab = ["ATCG", "TAGC"]
	elif t == 'rna':
		ab = ["AUCG", "UAGC"]
	else:
		print('ERROR: unknown na type.')
		return()

	rab = ab[1].strip().lower()
	ab = ab[0].strip().lower()

	# Check provided string
	na = na.lower()

	for c in na:
		if not c in ab:
			print('ERROR: provided string conflicts with selected alphabet.')
			return()

	# Calculate reverse
	r = na[::-1]

	# Calculate reverse complement
	rc = []
	for c in r:
		rc.append(rab[ab.index(c)])

	rc=''.join([str(c) for c in rc]).upper()

	return(rc)

def base2num(b):
	'''
	Converts nucleic acid base character to a number.

	Args:
		b (string): nucleic acid base character.

	Returns:
		int: nucleic acid base number.
	'''

def nbase_is_match(b1, b2):
	'''
	Checks if two nucleic acid bases in numeric form are complementary.

	Args:
		b1 (string): nucleic acid base in numeric form.
		b2 (string): nucleic acid base in numeric form.

	Returns:
		bool.
	'''

def find_loops(seq):
	'''
	Identifies all possible loops in a single-strand nucleic acid.

	Args:
		seq (string): nucleic acid string
	'''

# RUN ==========================================================================

# END ==========================================================================

################################################################################
