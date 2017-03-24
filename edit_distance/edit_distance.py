#!/usr/bin/python3
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 0.0.1
# Description: 	Set of functions to calculate the edit distance between two
# 				lists of elements.
# 
# ------------------------------------------------------------------------------


# DEPENDENCIES =================================================================

import numpy as np

# PARAMETERS ===================================================================

# FUNCTIONS ====================================================================

def dHamming(l1, l2, relative = None):
	'''Calculate Hamming distance between two lists of elements.
	The lists MUST have the same number of elements.
	if relative : (identical) 0 <= d <= 1 (totally different)
	else : (identical) 0 <= d <= len(l1) (totally different)

	Args:
		l1 (list): first list.
		l2 (list): second list.
		relative (boolean): whether to calculate relative Hamming distance.

	Returns:
		int: the Hamming distance between l1 and l2.
		float: the Hamming distance between l1 and l2 relative to their length.
	'''

	# Convert strings to lists
	if type('') == type(l1):
		l1 = [c for c in l1]
	if type('') == type(l2):
		l2 = [c for c in l2]

	# Lists MUST be lists
	if not type(l1) == type([]) or not type(l2) == type([]):
		return None

	# The lists MUST have the same length
	if len(l1) != len(l2):
		return None

	# Default relative value: False
	if relative is None:
		relative = False

	# Calculate Hamming distance
	d = sum([1 if l1[i] != l2[i] else 0 for i in range(len(l1))])

	# Output
	if relative:
		return d / float(len(l1))
	return d

def dJaro(l1, l2):
	'''Calculate Jaro distance between two lists of elements.
	(identical) 0 <= d <= 1 (totally different)

	Args:
		l1 (list): first list.
		l2 (list): second list.

	Returns:
		float: the Jaro distance between l1 and l2.
	'''

	# Convert strings to lists
	if type('') == type(l1):
		l1 = [c for c in l1]
	if type('') == type(l2):
		l2 = [c for c in l2]

	# Lists MUST be lists
	if not type(l1) == type([]) or not type(l2) == type([]):
		return None

	# String lengths
	l1_len = float(len(l1))
	l2_len = float(len(l2))

	# Calculate threshold
	dthr = abs(max(l1_len, l2_len) / 2. - 1)

	# Empty matching lists
	l1_matches = np.array([False] * int(l1_len))
	l2_matches = np.array([False] * int(l2_len))

	# Count matching characters
	matches = 0
	for i in range(int(l1_len)):

		# Check if every char in l2 matched already
		condition = np.logical_not(l2_matches)
		if 0 == sum(condition):
			break

		# Add matching condition to unmatched characters
		condition = np.logical_and(condition, np.array(l2) == l1[i])

		# Distance of matching characters from current position
		dmatch = abs(np.where(condition)[0])

		# Select only those not farther than dthr as matching
		dpass = dmatch[np.where(abs(dmatch - i) <= dthr)[0]]

		# Update match 
		if 0 != len(dpass):
			l1_matches[i] = True
			l2_matches[min(dpass)] = True
			matches += 1

	if 0 == matches:
		return 0

	# Count transpositions
	transpositions = 0
	k = 0
	for i in np.where(l1_matches)[0].tolist():
		# Stop if exceeding l2 length
		if k >= l2_len:
			break

		# Check transpositions
		js = np.where(l2_matches[k:])[0] + k
		if l1[i] != l2[min(js)]:
			transpositions += 1

		# Update k
		ks = np.where(l2_matches[k:])[0] + k
		if 0 != len(ks):
			k = ks.tolist()[0]
		k += 1

	return ((matches / l1_len) +
			(matches / l2_len) +
			((matches - transpositions/2.) / matches)) / 3.

def dWinkler(l1, l2, p = None, boost_thr = None):
	'''Calculate Jaro-Winkler distance between two lists of elements.
	(identical) 0 <= d <= 1 (totally different)

	Args:
		l1 (list): first list.
		l2 (list): second list.
		p (float): significance.
		boost_thr (float): apply prefix bonus only if dJaro >= boost_thr.

	Returns:
		float: the Jaro-Winkler distance between l1 and l2.
	'''

	# Default p
	if type(None) == type(p):
		p = 0.1

	# Default boost_thr
	if type(None) == type(boost_thr):
		boost_thr = 0.7

	# Convert strings to lists
	if type('') == type(l1):
		l1 = [c for c in l1]
	if type('') == type(l2):
		l2 = [c for c in l2]

	# Lists MUST be lists
	if not type(l1) == type([]) or not type(l2) == type([]):
		return None

	# Calculate dJaro
	dj = dJaro(l1, l2)

	# Check prefix
	prefix = 0
	for i in range(min(len(l1), len(l2))):
		if l1[i] != l2[i]:
			break
		prefix = min(4, prefix + 1)

	if dj >= boost_thr:
		return dj + (prefix * p * (1 - dj))
	else:
		return dj


# TEST =========================================================================

# END ==========================================================================

################################################################################
