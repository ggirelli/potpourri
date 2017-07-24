# Filter for GC
cat cdna_seq.30mer.filter.30_70_gc.noHpol.tsv | sed 's/^..//' | awk '($3 >= 0.4){ print "> "$0; }' > cdna_seq.30mer.filter.40_70_gc.noHpol.tsv
cat cds_seq.30mer.filter.30_70_gc.noHpol.tsv | sed 's/^..//' | awk '($3 >= 0.4){ print "> "$0; }' > cds_seq.30mer.filter.40_70_gc.noHpol.tsv
cat exon_seq.30mer.filter.30_70_gc.noHpol.tsv | sed 's/^..//' | awk '($3 >= 0.4){ print "> "$0; }' > exon_seq.30mer.filter.40_70_gc.noHpol.tsv
cat cdna_seq.30mer.filter.30_70_gc.noHpol.tsv | sed 's/^..//' | awk '($3 >= 0.4 && $3 <= 0.6){ print "> "$0; }' > cdna_seq.30mer.filter.40_60_gc.noHpol.tsv
cat cds_seq.30mer.filter.30_70_gc.noHpol.tsv | sed 's/^..//' | awk '($3 >= 0.4 && $3 <= 0.6){ print "> "$0; }' > cds_seq.30mer.filter.40_60_gc.noHpol.tsv
cat exon_seq.30mer.filter.30_70_gc.noHpol.tsv | sed 's/^..//' | awk '($3 >= 0.4 && $3 <= 0.6){ print "> "$0; }' > exon_seq.30mer.filter.40_60_gc.noHpol.tsv

# Generate FASTA file
cat exon_seq.30mer.filter.30_70_gc.noHpol.tsv | sed 's/^..//' | awk '{ print "> "$1":"$3":"$4":"$5"\n"$2 }' > exon_seq.30mer.filter.30_70_gc.noHpol.fa
cat cdna_seq.30mer.filter.40_60_gc.noHpol.tsv | sed 's/^..//' | awk '{ print "> "$1":"$3":"$4":"$5"\n"$2 }' > cdna_seq.30mer.filter.40_60_gc.noHpol.fa
cat cdna_seq.30mer.filter.40_70_gc.noHpol.tsv | sed 's/^..//' | awk '{ print "> "$1":"$3":"$4":"$5"\n"$2 }' > cdna_seq.30mer.filter.40_70_gc.noHpol.fa
cat cds_seq.30mer.filter.40_60_gc.noHpol.tsv | sed 's/^..//' | awk '{ print "> "$1":"$3":"$4":"$5"\n"$2 }' > cds_seq.30mer.filter.40_60_gc.noHpol.fa
cat cds_seq.30mer.filter.40_70_gc.noHpol.tsv | sed 's/^..//' | awk '{ print "> "$1":"$3":"$4":"$5"\n"$2 }' > cds_seq.30mer.filter.40_70_gc.noHpol.fa
cat exon_seq.30mer.filter.40_60_gc.noHpol.tsv | sed 's/^..//' | awk '{ print "> "$1":"$3":"$4":"$5"\n"$2 }' > exon_seq.30mer.filter.40_60_gc.noHpol.fa
cat exon_seq.30mer.filter.40_70_gc.noHpol.tsv | sed 's/^..//' | awk '{ print "> "$1":"$3":"$4":"$5"\n"$2 }' > exon_seq.30mer.filter.40_70_gc.noHpol.fa

# Split FASTA file into single-gene fasta
split_fa.py -o cds_seq.30mer.filter.30_70_gc.noHpol.fa cds.30mer.30_70_gc.noHpol
split_fa.py -o cds_seq.30mer.filter.40_60_gc.noHpol.fa cds.30mer.40_60_gc.noHpol
split_fa.py -o cds_seq.30mer.filter.40_70_gc.noHpol.fa cds.30mer.40_70_gc.noHpol
split_fa.py -o cdna_seq.30mer.filter.30_70_gc.noHpol.fa cdna.30mer.30_70_gc.noHpol
split_fa.py -o cdna_seq.30mer.filter.40_60_gc.noHpol.fa cdna.30mer.40_60_gc.noHpol
split_fa.py -o cdna_seq.30mer.filter.40_70_gc.noHpol.fa cdna.30mer.40_70_gc.noHpol
split_fa.py -o exon_seq.30mer.filter.30_70_gc.noHpol.fa exon.30mer.30_70_gc.noHpol
split_fa.py -o exon_seq.30mer.filter.40_60_gc.noHpol.fa exon.30mer.40_60_gc.noHpol
split_fa.py -o exon_seq.30mer.filter.40_70_gc.noHpol.fa exon.30mer.40_70_gc.noHpol
