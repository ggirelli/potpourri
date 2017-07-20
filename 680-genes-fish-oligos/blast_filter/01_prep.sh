
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
