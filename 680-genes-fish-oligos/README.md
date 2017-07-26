680-genes-fish-oligos
===

Scripts to design RNA FISH probes for 680 genes.

Of the original 616+64 genes in the list some genes were annotated by manually picking the corresponding Ensembl stable gene ID.

## mk_trans_db.py

Designs databases comprising, for every census gene: all exon, CDS or CDS+UTR. Also, retrieves data on genes and transcripts and genes sequences.

## mk_oligos.py

Generates fasta file with all k-mer from input fasta (database).

## characterize_oligos.py

Calculates melting temperature, GC-content and homopolymer presence of all sequences in input. Input: a file with one oligo per sequence (e.g., fasta without headers).

## split_fa_by_gene.sh and split_fa.py

Splits a fasta by gene (based on header pattern).

## blast_filter

Contains scripts for BLASTN output preparation and filtering, in a parallel fashion. The code strictly resembles the `blast_filter.py` script.

## 01_prep.sh

contains the step to be performed before blast_filter can be run.