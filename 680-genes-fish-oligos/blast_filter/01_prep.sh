
# Filter BLAST output based on homology
cat blast.out.tsv | awk '($4-$5)/30 >= 0.85' > blast.out.85percHom.tsv

# Make oligo_seq|hit_transcriptID table
cat 30mer.uniq.filter.40_70_gc.noHpol.fa | paste - - | sed 's/.//' | join -j 1 -t $'\t' - <(cat blast.out.85percHom.tsv) -o 1.2,2.2 | sed 's/\.[0-9]*$//' > blast.out.85percHom.transcripts.txt

# Identify transcripts
cut -f 2 blast.out.85percHom.transcripts.txt | sort | uniq

# Retrieve transcript gene symbol
cat transcript_list.txt | join -j 1 -t $'\t' - <(cat human_gene_symbol_tr_id_table.kf.tsv | sort -k1) > transcript_gene_list.txt

# Make oligo_seq|hit_transcriptID|Gene_Symbol table
cat blast.out.85percHom.transcripts.txt | sort -k2 | join -1 2 -2 1 - <(cat transcript_gene_list.txt) -t $'\t' -o 1.1,1.2,2.2 | sort -k1 --paralle 40 -S10G > blast.out.85percHom.transcripts.gene.txt



##########################################################################

# blastn -db /media/gire/Data/BiCro-Resources/references/Grch38_89_cdna_ncrna/Grch38.89.cdna.ncrna -query /media/gire/Data/BiCro-Analysis/Sequencing/680genes/oligo_uniq_filtered/30mer.uniq.filter.30_70_gc.noHpol.fa -outfmt 6 -out /home/gire/Desktop/30mer.uniq.filter.30_70_gc.noHpol.blast.out.tsv -word_size 10 -num_threads 35 -penalty -2 -reward 3 -gapopen 100 -gapextend 100 -evalue 1000 -strand plus

# Filter BLAST output based on homology with more relaxed threshold
cat blast.out.tsv | awk '($4-$5)/30 >= 0.7' > blast.out.70percHom.tsv
sort -k1,1 --parallel 30 -S20G blast.out.70percHom.tsv > blast.out.70percHom.sorted.tsv

# Make maxHom|oligo_seq|hit_transcriptID table
cat 30mer.uniq.filter.40_70_gc.noHpol.fa | paste - - | sed 's/^.//' | sort -k1,1 | join -j 1 -t $'\t' - <(cat 30mer.uniq.filter.40_70_gc.noHpol.blast.out.70percHom.sorted.tsv) -o 1.1,1.2,2.2 | awk '{OFS=FS="\t"; split($1, ff, ":"); split($3, gg, "."); print ff[2]"\t"$2"\t"gg[1]}' | sed 's/\.[0-9]*$//' > 30mer.uniq.filter.40_70_gc.noHpol.blast.out.70percHom.transcripts.tsv

