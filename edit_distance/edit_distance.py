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

def str2list(s):
	'''Convert string to list of characters.

	Args:
		s (strings)

	Returns:
		list: if the input was a string, its conversion to list.
		any: the original input if it was not a string.
	'''
	if type('') == type(s):
		return [c for c in s]
	return s

def unique(l):
	'''Remove duplicated elements from a list.

	Args:
		l (list): list to be de-duplicated.

	Returns:
		list: l without duplicated elements.
	'''
	return [e for e in set(l)]

def dHamming(l1, l2, relative = None):
	'''Calculate Hamming distance between two lists of elements.
	The lists MUST have the same number of elements.
	if relative : (identical) 0 <= d <= 1 (totally different)
	else : (identical) 0 <= d <= len(l1) (totally different)

	Args:
		l1 (list, string): first list.
		l2 (list, string): second list.
		relative (boolean): whether to calculate relative Hamming distance.

	Returns:
		int: the Hamming distance between l1 and l2.
		float: the Hamming distance between l1 and l2 relative to their length.
	'''

	# Default relative value: False
	if type(None) == type(relative):
		relative = False

	# Convert strings to lists
	l1 = str2list(l1)
	l2 = str2list(l2)

	# Lists MUST be lists
	if any([not type(l) == type([]) for l in [l1, l2]]):
		return None

	# The lists MUST have the same length
	if len(l1) != len(l2):
		msg = 'The Hamming distance is not defined '
		msg += 'for sets of different lengths.'
		print(msg)
		return None

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
		l1 (list, string): first list.
		l2 (list, string): second list.

	Returns:
		float: the Jaro distance between l1 and l2.
	'''

	# Convert strings to lists
	l1 = str2list(l1)
	l2 = str2list(l2)

	# Lists MUST be lists
	if any([not type(l) == type([]) for l in [l1, l2]]):
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
		l1 (list, string): first list.
		l2 (list, string): second list.
		p (float): significance.
		boost_thr (float): apply prefix bonus only if dJaro >= boost_thr.

	Returns:
		float: the Jaro-Winkler distance between l1 and l2.
	'''

	# Default p : 0.1
	if type(None) == type(p):
		p = 0.1

	# Default boost_thr : 0.7
	if type(None) == type(boost_thr):
		boost_thr = 0.7

	# Convert strings to lists
	l1 = str2list(l1)
	l2 = str2list(l2)

	# Lists MUST be lists
	if any([not type(l) == type([]) for l in [l1, l2]]):
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

def dOSA(l1, l2, relative = None):
	'''Calculate the optimal string alignment distance between two lists of
	elements, as the number of edit operations (insertions, deletions,
	substitutions or transpositions) needed to make them identical,
	under the condition that no substring is edited more than once.
	if relative: (identical) 0 <= d <= 1 (totally different)
	else: (identical) 0 <= d <= max(len(l1), len(l2)) (totally different)

	Args:
		l1 (list, string): first list.
		l2 (list, string): second list.
		relative (boolean): whether to calculate relative Hamming distance.

	Returns:
		int: the optimal string alignment distance between l1 and l2.
		float: the optimal string alignment distance between l1 and l2, relative
				to the maximum possible distance.
	'''

	# Default relative value: False
	if type(None) == type(relative):
		relative = False

	# Convert strings to lists
	l1 = str2list(l1)
	l2 = str2list(l2)

	# Lists MUST be lists
	if any([not type(l) == type([]) for l in [l1, l2]]):
		return None

	# Lists length
	l1_len = len(l1)
	l2_len = len(l2)

	# Create matrix
	d = np.zeros((l1_len, l2_len))

	# Starting distances
	d[:, 0] = np.arange(l1_len)
	d[0, :] = np.arange(l2_len)

	# Fill the matrix
	for i in range(1, l1_len):
		for j in range(1, l2_len):
			cost = 0
			if not l1[i] == l2[j]:
				cost = 1
			d[i, j] = min([
				d[i-1, j] +1,		# deletion
				d[i, j-1] + 1,		# insertion
				d[i-1, j-1] + cost	# substitution
			])
			if i > 1 and j > 1 and l1[i] == l2[j-1] and l1[i-1] == l2[j]:
				d[i, j] = min([d[i, j], d[i-2, j-2] + cost]) # transposition

	# Output bottom right corner
	if relative:
		return int(d[-1, -1]) / float(max([l1_len, l2_len]))
	return int(d[-1, -1])

def dLevenshtein(l1, l2, relative = None):
	'''Calculate the Levenshtein distance between two lists of
	elements, as the number of edit operations (insertions, deletions or
	substitutions) needed to make them identical.
	if relative: (identical) 0 <= d <= 1 (totally different)
	else: (identical) 0 <= d <= max(len(l1), len(l2)) (totally different)

	Args:
		l1 (list, string): first list.
		l2 (list, string): second list.
		relative (boolean): whether to calculate relative Hamming distance.

	Returns:
		int: the Levenshtein distance between l1 and l2.
		float: the Levenshtein distance between l1 and l2, relative to the
				maximum possible value.
	'''

	# Default relative value: False
	if type(None) == type(relative):
		relative = False

	# Convert strings to lists
	l1 = str2list(l1)
	l2 = str2list(l2)

	# Lists MUST be lists
	if any([not type(l) == type([]) for l in [l1, l2]]):
		return None

	# Lists length
	l1_len = len(l1)
	l2_len = len(l2)

	# Create matrix
	d = np.zeros((l1_len, l2_len))

	# Starting distances
	d[:, 0] = np.arange(l1_len)
	d[0, :] = np.arange(l2_len)

	# Fill the matrix
	for i in range(1, l1_len):
		for j in range(1, l2_len):
			cost = 0
			if not l1[i] == l2[j]:
				cost = 1
			d[i, j] = min([
				d[i-1, j] +1,		# deletion
				d[i, j-1] + 1,		# insertion
				d[i-1, j-1] + cost	# substitution
			])

	# Output bottom right corner
	if relative:
		return int(d[-1, -1]) / float(max([l1_len, l2_len]))
	return int(d[-1, -1])

def dDamerau(l1, l2, absize = None, relative = None):
	'''Calculate the Damerau-Levenshtein distance between two lists of
	elements, as the number of edit operations (insertions, deletions or
	substitutions, or transpositions) needed to make them identical.
	if relative: (identical) 0 <= d <= 1 (totally different)
	else: (identical) 0 <= d <= max(len(l1), len(l2)) (totally different)

	Args:
		l1 (list, string): first list.
		l2 (list, string): second list.
		absize (int): alphabet size.
		relative (boolean): whether to calculate relative Hamming distance.

	Returns:
		int: the Damerau-Levenshtein distance between l1 and l2.
		float: the Damerau-Levenshtein distance between l1 and l2, relative to
				the maximum possible value.
	'''

	# Default relative value: False
	if type(None) == type(relative):
		relative = False

	# Convert strings to lists
	l1 = str2list(l1)
	l2 = str2list(l2)

	# Default alphabet size value
	lu = l1
	lu.extend(l2)
	lu = unique(lu)
	lu_len = len(lu)
	if type(None) == type(absize):
		absize = lu_len
	elif absize < lu_len:
		print('The provided alphabet size does not match the input sets.')
		absize = lu_len
	print('Using an alphabet of up to ' + str(absize) + ' elements.')
	lu.sort()
	print(lu)

	# Lists MUST be lists
	if any([not type(l) == type([]) for l in [l1, l2]]):
		return None

	# Alphabet-sized vector
	da = {}
	for e in lu:
		da[e] = 1

	# Lists length
	l1_len = len(l1)
	l2_len = len(l2)

	# Create matrix
	d = np.zeros((l1_len + 1, l2_len + 1))

	# Starting distances
	maxdist = l1_len + l2_len
	d[:, 0] = maxdist
	d[:, 1] = np.arange(l1_len + 1)
	d[0, :] = maxdist
	d[1, :] = np.arange(l2_len + 1)

	# Fill the matrix
	for i in range(2, l1_len + 1):
		db = 1
		for j in range(2, l2_len + 1):
			k = da[l2[j-1]]
			l = db
			cost = 1
			if l1[i-1] == l2[j-1]:
				cost = 0
				db = j
			d[i, j] = min([
				d[i-1, j-1] + cost,						# substitution
				d[i, j-1] + 1,							# insertion
				d[i-1, j] + 1,							# deletion
				d[k-1, l-1] + (i-k-1) + 1 + (j-l-1)		# transposition
			])
		da[l1[i-1]] = i

	# Output bottom right corner
	if relative:
		return int(d[-1, -1]) / float(max([l1_len, l2_len]))
	return int(d[-1, -1])

def longestCommonSubseqs(s1, s2):
	'''Find the longest common subsequence of two strings.

	Args:
		s1 (string): first string.
		s2 (string): second string.

	Returns:
		string: the longest common subsequence between s1 and s2.
		list: list of multiple longest common subsequences.
	'''

	# Strings MUST be strings
	if any([not type('') == type(e) for e in [s1, s2]]):
		return none

	# Strings length
	s1_len = len(s1)
	s2_len = len(s2)

	# Set a set with the longer string first
	if s1_len > s2_len:
		ss = [s1, s2]
		ss_len = [s1_len, s2_len]
	else:
		ss = [s2, s1]
		ss_len = [s2_len, s1_len]

	# Cycle through the shortest string substrings, from the longest.
	lcs = ['']
	lcs_len = 0
	for s in range(1, ss_len[1] + 1)[::-1]:
		for i in range(0, ss_len[1] - s + 1):
			substr = ss[1][i:(i + s)]
			if substr in ss[0] and s > lcs_len:
				lcs = [substr]
				lcs_len = s
			elif substr in ss[0] and s == lcs_len:
				lcs.append(substr)
	return lcs

def dLCS(s1, s2):
	'''Uses longest common subsequence as a distance between two strings.
	(identical) 0 <= d <= 1 (totally different)

	Args:
		s1 (string): first string.
		s2 (string): second string.

	Returns:
		float: ration between LCS length and min(len(s1), len(s2)).
	'''


	# Strings MUST be strings
	if any([not type('') == type(e) for e in [s1, s2]]):
		return none

	LCS = longestCommonSubseqs(s1, s2)[0]

	if 0 == len(LCS):
		return 1

	return 1 - len(LCS) / float(min(len(s1), len(s2)))

# TEST =========================================================================

# END ==========================================================================

################################################################################
