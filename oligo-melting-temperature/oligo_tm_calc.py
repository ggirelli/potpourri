#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.0
# Date: 20170711
# Project: oligo characterization
# Description:	calculate melting temperature of a provide DNA duplex
# 		
# @TODO:
# 	Correction for Mg2+ concentration does not match IDT OligoAnalyzer output.
# 
# Changelog:
# 		1.0.0: first implementation.
# 
# References:
# 	[1] Freier et al, PNAS(83), 1986
# 	[2] Allawi & Santalucia, Biochemistry(36), 1997
# 	[3] SantaLucia, PNAS(95), 1998
# 	[4] Owczarzy et al, Biochemistry(43), 2004
# 	[5] Owczarzy et al, Biochemistry(47), 2008
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

import argparse
import math

# PARAMETERS ===================================================================


# Add script description
parser = argparse.ArgumentParser(
	description = '''
Calculate melting temeprature of a DNA duplex at provided [oligo],
[Na+], [Mg2+]. References:
 [1] Freier et al, PNAS(83), 1986;
 [2] Allawi & Santalucia, Biochemistry(36), 1997;
 [3] SantaLucia, PNAS(95), 1998;
 [4] Owczarzy et al, Biochemistry(43), 2004;
 [5] Owczarzy et al, Biochemistry(47), 2008.
''')

# Add mandatory arguments
parser.add_argument('seq', type = str, nargs = 1,
	help = 'DNA duplex sequence (one strand only).')

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
parser.add_argument('-t', '--thermo', type = str, nargs = 1,
	help = '''Thermodynamic table to use in the calculations.
	Possible values: allawi (based on ref.2, default) or 
	freier (based on ref.1).''',
	choices = ['allawi', 'freier'], default = ['allawi'])
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

# RUN ==========================================================================

# Make string uppercase
seq = seq.upper()

# Make NN couples
couples = [seq[i:(i+2)] for i in range(len(seq) - 1)]

# Calculate GC content
fgc = (seq.count('G') + seq.count('C')) / len(seq)

# 1 M NaCl case
# Based on SantaLucia, PNAS(95), 1998
# -----------------------------------
na_conc_0 = 1.

# Calculate dH0(N-N)
h = sum([tt[tt_mode][c][0] for c in couples])

# Add initiation enthalpy
if tt[tt_mode]['has_end']:
	h += tt[tt_mode]['end%s' % (seq[0],)][0]
	h += tt[tt_mode]['end%s' % (seq[-1],)][0]

# Correct enthalpy for symmetry
if tt[tt_mode]['has_sym'] and seq == rc(seq, 'dna'):
	h += tt[tt_mode]['sym'][0]

# Calculate dS0(N-N) in kcal / (K mol)
s = sum([tt[tt_mode][c][1] for c in couples])

# Add initiation enthalpy
if tt[tt_mode]['has_end']:
	s += tt[tt_mode]['end%s' % (seq[0],)][1]
	s += tt[tt_mode]['end%s' % (seq[-1],)][1]

# Correct enthalpy for symmetry
if tt[tt_mode]['has_sym'] and seq == rc(seq, 'dna'):
	s += tt[tt_mode]['sym'][1]

s /= 1e3

# Calculate melting temperature in Celsius
Tm1 = h / (s + R * math.log(oligo_conc))

# Adjusted for [Na]
# Based on Owczarzy et al, Biochemistry(43), 2004
# -----------------------------------------------
Tm2 = Tm1

# Parameters from paper
na_a = 4.29e-5
na_b = 3.95e-5
na_c = 9.4e-6

# Adjust
if not 0  == na_conc:
	Tm2r = (1. / Tm1)
	Tm2r += (na_a * fgc - na_b) * math.log(na_conc / na_conc_0)
	Tm2r += na_c * (math.log(na_conc / na_conc_0) ** 2)
	Tm2 = 1. / Tm2r

# Adjusted for Mg
# Based on Owczarzy et al, Biochemistry(47), 2008
# -----------------------------------------------
Tm3 = Tm2

# Parameters from paper
mg_a = 3.92e-5
mg_b = -9.11e-6
mg_c = 6.26e-5
mg_d = 1.42e-5
mg_e = -4.82e-4
mg_f = 5.25e-4
mg_g = 8.31e-5

# Adjust
if not 0  == mg_conc:
	mg_conc_log = math.log(mg_conc)
	Tm3r = 1./Tm2 + mg_a + mg_b*mg_conc_log + fgc*(mg_c + mg_d*mg_conc_log)
	Tm3r_factor = mg_e + mg_f*mg_conc_log + mg_g*(mg_conc_log)**2
	Tm3r_factor *= (1./(2*(len(seq) - 1)))
	Tm3r += Tm3r_factor
	Tm3 = 1./Tm3r

# Log output
# ----------

if not is_verbose:
	if celsius:
		print(Tm3 - 273.15)
	else:
		print(Tm3)
else:
	print("""
	  Oligo sequence : %s
	         [oligo] : %f M
	           [NA+] : %f M
	          [Mg2+] : %f M
	  Thermod. table : %s

	             dH0 : %f kcal/mol 
	             dS0 : %f kcal/(k·mol)
	             dG0 : %f kcal/mol

	  [Na+] = 1 M    : Tm = %f K (= %f degC)

	  [Na+] = %f M   : Tm = %f K (= %f degC)

	  [Na+]  = %f M
	  [Mg2+] = %f M  : Tm = %f K (= %f degC)
	""" % (
		seq, oligo_conc, na_conc, mg_conc, tt_mode,
		h, s, h - (37 + 273.15) * s,
		Tm1, Tm1 - 273.15,
		na_conc, Tm2, Tm2 - 273.15,
		na_conc, mg_conc, Tm3, Tm3 - 273.15,
	))

# END ==========================================================================

################################################################################