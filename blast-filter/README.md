blast-filter
===

The script is designed to analyze the output of BLASTing oligos for RNA FISH probe design. It filters BLAST output based on homology percentage (as number of perfect matches over query length). Then check for off-targets and saturated off-targets (i.e., transcripts off-targeted by a sufficient number of oligos to generate a false positive).