cosmic-census-fish-oligos
===

Scripts to design RNA FISH probes for COSMIC CENSUS genes.

Of the original 616 genes in the list:

- 13 were discarded from the analysis as missing proper coordinates and Ensembl
    stable gene ID.
- 28 were discarded as missing Ensembl stable gene ID.
- 9 had a non-existeng Ensembl gene ID.
- 1 had no CDS information in its transcripts.

## mk_trans_db.py

Designs databases comprising, for every census gene: all exon, CDS or CDS+UTR. Also, retrieves data on genes and transcripts and genes sequences.

## mk_oligos.py

Generates fasta file with all k-mer from input fasta (database).

## characterize_oligos.py

Calculates melting temperature, GC-content and homopolymer presence of all sequences in input. Input: a file with one oligo per sequence (e.g., fasta without headers).