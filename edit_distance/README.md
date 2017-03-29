Edit distance
===

A set of functions to calculate the edit distance between two lists of elements.

## Edit distances

### Hamming distance

> The **Hamming distance** measures the minimum number of *substitutions* required to change one string into the other, or the minimum number of *errors* that could have transformed one string into the other. ~Wikipedia

Implemented in:

* `edit_distance.py` : `dHamming`

### Jaro & Jaro-Winkler distances

> The **Jaro distance** between two words is the minimum number of single-character transpositions required to change one word into the other. ~Wikipedia

> The **Jaro-Winkler distance** is a variation of the Jaro metric, that gives more favourable ratings to strings that match from the beginning. ~Wikipedia

Implemented in:

* `edit_distance.py` : `dJaro`
* `edit_distance.py` : `dWinkler`

### Optimal string alignment distance

> The optimal string alignment algorithm computes the number of edit operations (consisting of *insertions*, *deletions* or *substitutions* of a single character, or *transposition* of two adjacent characters) needed to make the strings equal under the condition that no substring is edited more than once. ~Wikipedia

Implemented in:

* `edit_distance.py` : `dOSA`

### Levenshtein & Damerau-Levenshtein distances

> The **Levenshtein distance** between two words is the minimum number of single-character edits (*insertions*, *deletions* or *substitutions*) required to change one word into the other. ~Wikipedia

> The **Damerau–Levenshtein distance** between two words is the minimum number of operations (consisting of *insertions*, *deletions* or *substitutions* of a single character, or *transposition* of two adjacent characters) required to change one word into the other. ~Wikipedia

Implemented in:

* `edit_distance.py` : `dLevenshtein`
* `edit_distance.py` : `dDamerau`

### Longest common subsequence size

The size (or length) of the **longest common subsequence** can be used as a measure of similarity.

Implemented in:

* `edit_distance.py` : `dLCS`

### Kendall tau distance

> The Kendall tau distance is equivalent to the number of swaps that the bubble sort algorithm would make to place one list in the same order as the other list. ~Wikipedia

Implemented in:

* `edit_distance.py` : `dKendall`