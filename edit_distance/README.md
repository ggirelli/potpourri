Edit distance
===

A set of functions to calculate the edit distance between two lists of elements.

## Edit distances

### Hamming distance

> The **Hamming distance** measures the minimum number of *substitutions* required to change one string into the other, or the minimum number of *errors* that could have transformed one string into the other. ~Wikipedia

Implemented in:

* `edit_distance.py` : `dHamming`

### Jaro & Jaro-Winkler distance

> The **Jaro distance** between two words is the minimum number of single-character transpositions required to change one word into the other. ~Wikipedia

> The **Jaro-Winkler distance** is a variation of the Jaro metric, that gives more favourable ratings to strings that match from the beginning. ~Wikipedia

Implemented in:

* `edit_distance.py` : `dJaro`
* `edit_distance.py` : `dWinkler`

### Damerau-Levenshtein distance

### Longest common subsequence
