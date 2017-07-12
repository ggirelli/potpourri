#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 0.1.0
# Date: 20170703
# Project: COSMIC cancer gene census oligo characterization
# Description:	build database with the longest transcript of every COSMIC gene.
# 
# Note:
# 	CDS:	"Coding DNA Sequence", comprising exons but not UTRs
# 	cDNA:	"complementary DNA", gene sequence comprising exons
# 			(including UTRs) and introns
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

from biomart import BiomartServer
import numpy as np
import pandas as pd
import re

# PARAMETERS ===================================================================

# INPUT
# -----

# COSMIC tsv table
cosmic_file = "/home/gire/Desktop/BiCro-Data/Sequencing/"
cosmic_file += "Census_allMon Jul  3 09-56-47 2017.tsv"

# Gene ID conversion files
ensembl2entrez_file = "-dbOrg-Ensembl_Gene_ID__to__Gene_ID_9606"
entrez2ensembl_file = "-dbOrg-Gene_ID__to__Ensembl_Gene_ID_9606"

# OUTPUT
# ------

# Gene data table
out_gene_table_file = "gene_data.tsv"

# Gene fasta
out_gene_fasta_file = "gene_seq.fa"

# Transcript data table
out_trans_table_file = "trans_data.tsv"

# Transcript CDS fasta
out_trans_cds_fasta_file = "trans_cds_seq.fa"

# Transcript CDS+UTR fasta
out_trans_cds_utr_fasta_file = "trans_cds_utr_seq.fa"

# Exon fasta
out_exon_table_file = "exon_seq.tsv"

# Skipped genes
out_skip_gene_file = "skipped_genes.tsv"

# FUNCTIONS ====================================================================

def mk_bioDBnet_db(fname):
	# Build dictionary from bioDBnet file.
	# 
	# Args:
	# 	fname (string): path to bioDBnet table.
	# 
	# Return:
	# 	d (dict): {ensembl:entrez}
	# 	
	# Instantiate emtpy dictionary
	d = {}
	# Read file line by line
	with open(fname) as f:
		for line in f:
			# Identify key-value pairs
			(k, v) = tuple(line.strip().split('\t'))
			# If value is assigned, store it
			if not v == ' -':
				d[k] = v
	# Close file pointer
	f.close()
	# Return dictionary
	return(d)

# RUN ==========================================================================

# Gene ID conversion
# Prepare Gene ID database conversion dictionaries
# ------------------------------------------------------------------------------

ensembl2entrez = mk_bioDBnet_db(ensembl2entrez_file)
entrez2ensembl = mk_bioDBnet_db(entrez2ensembl_file)

# COSMIC
# Prepare table to be compatible with ENSEMBL
# ------------------------------------------------------------------------------

# Read COSMIC table
cosmic = pd.read_csv(cosmic_file, '\t')

# Identify genes with no synonym
no_sym_id = np.where(type('') != np.array([type(i)
	for i in np.array(cosmic['Synonyms'])]))[0]

# Remove genes with no synonym
cosmic2 = cosmic.drop(no_sym_id)
cosmic2.index = np.array(range(cosmic2.shape[0]))
print("Removed %d genes without proper coordinates." % (len(no_sym_id),))

# Identify genes without Ensembl Stable Gene ID
no_ens_id = np.where([0 == len(re.findall('ENSG', e))
	for e in np.array(cosmic2['Synonyms'])])[0].tolist()

# Remove genes without Ensembl Stable Gene ID
cosmic2 = cosmic2.drop(no_ens_id)
cosmic2.index = np.array(range(cosmic2.shape[0]))
print("Removed %d genes without Ensembl Stable Gene ID." % (len(no_ens_id),))

print("Left with %d genes (out of %d)." % (cosmic2.shape[0], cosmic.shape[0]))

# Extract ensembl gene ID
ensembl_id = []
for e in cosmic2['Synonyms']:
	k = np.nan
	for i in e.split(','):
		if i.startswith('ENSG'):
			k = i
			continue
	ensembl_id.append(k)

# BioMaRt
# Connect to the server and prepare for querying
# ------------------------------------------------------------------------------

# Connect to biomart
server = BiomartServer( "http://www.ensembl.org/biomart" )
server.verbose = False

# Check available databases
#server.show_databases()

# Select Genes database
db = server.databases['ENSEMBL_MART_ENSEMBL']

# Check available datasets (species)
#db.show_datasets()

# Select H. sapiens dataset
ds = db.datasets['hsapiens_gene_ensembl']

# show all available filters and attributes of the 'uniprot' dataset
#ds.show_filters()
#ds.show_attributes()

# Test connection and filtering
#print(ds.count())
#print(ds.count({'filters':{'chromosome_name':'1'}}))

# ENTREZGENE
# Use ENTREZ NCBI GENE ID for selection
# ------------------------------------------------------------------------------

# Point to files
out_gene_table = open(out_gene_table_file, 'w')
out_gene_fasta = open(out_gene_fasta_file, 'w')
out_trans_table = open(out_trans_table_file, 'w')
out_trans_cds_fasta = open(out_trans_cds_fasta_file, 'w')
out_trans_cds_utr_fasta = open(out_trans_cds_utr_fasta_file, 'w')
out_exon_table = open(out_exon_table_file, 'w')
out_skip_gene = open(out_skip_gene_file, 'w')

# Print headers
out_gene_table_head = 'ensembl_gene_id\tchromosome_name\tstart_position\t'
out_gene_table_head += 'end_position\tstrand\ttranscript_count\t'
out_gene_table_head += 'percentage_gene_gc_content\n'
out_gene_table.write(out_gene_table_head)

out_trans_table_head = 'ensembl_gene_id\tensembl_transcript_id\tcds_length\t'
out_trans_table_head += 'transcript_length\n'
out_trans_table.write(out_trans_table_head)

out_exon_table_head = 'ensembl_gene_id\tensembl_exon_id\texon_seq\n'
out_exon_table.write(out_exon_table_head)

# Cycle through genes
for idi in range(len(ensembl_id)):
	print("Working on %s (%d/%d)..."
		% (ensembl_id[idi], idi + 1, len(ensembl_id)))

	# Gene
	# Query for Gene characteristics
	# --------------------------------------------------------------------------

	# Retrieve gene information
	# -------------------------

	print("> Retrieving gene information...")

	# Submit query
	response = ds.search({
		'filters':{
			'ensembl_gene_id':ensembl_id[idi]
		},
		'attributes':[
			'ensembl_gene_id',				# Gene ID
			'chromosome_name',				# Chromosome
			'start_position',				# Start
			'end_position',					# End
			'strand',						# Strand
			'transcript_count',				# Number of transcripts
			'percentage_gene_gc_content'	# GC content percentage
		]
	}, header = 1)

	# Convert output
	gene_data = np.array([row.split('\t')
		for row in response.text.strip().split('\n')])

	# Check if gene was found
	if 1 == gene_data.shape[0]:
		out_skip_gene.write(ensembl_id[idi] + '\tNo gene found.\n')
		continue

	# Write
	sout = "\n".join(["\t".join(row) for row in gene_data[1:,:]]) + "\n"
	out_gene_table.write(sout)

	# Retrieve the whole gene sequence
	# --------------------------------

	print("> Retrieving gene sequence...")

	# Submit query
	response = ds.search({
		'filters':{
			'ensembl_gene_id':ensembl_id[idi]
		},
		'attributes':[
			'ensembl_gene_id',				# Gene ID
			'gene_exon_intron'				# Gene cDNA
		]
	})

	# Convert output
	gene_seq = np.array([row.split('\t')
		for row in response.text.strip().split('\n')])

	# Write
	sout = '> ' + gene_seq[0][1] + '\n' + gene_seq[0][0] + '\n'
	out_gene_fasta.write(sout)

	# Transcripts
	# Query for transcript data
	# --------------------------------------------------------------------------

	# Identify transcripts
	# --------------------

	print("> Retrieving transcript information...")

	# Submit query
	response = ds.search({
		'filters':{
			'ensembl_gene_id':ensembl_id[idi]
		},
		'attributes':[
			'ensembl_gene_id',			# Gene ID
			'ensembl_transcript_id',	# Transcript ID
			'cds_length',				# CDS length
			'transcript_length'			# CDS+UTRs length
		]
	})

	# Convert output
	trans_data = np.array([row.split('\t')
		for row in response.text.strip().split('\n')])

	# Check if transcript was found
	if 0 == trans_data.shape[0]:
		out_skip_gene.write(ensembl_id[idi] + '\tNo transcript found.')
		continue

	# Write
	sout = "\n".join(["\t".join(row) for row in trans_data]) + "\n"
	out_trans_table.write(sout)

	# Retrieve sequence of transcript with longest CDS
	# ------------------------------------------------

	print("> Retrieving transcript CDS (longest) sequence...")

	# Remove transcripts with no CDS information
	cds_data = trans_data[np.where('' != trans_data[:, 2])[0], :]

	if 0 != cds_data.shape[0]:
		# Identify transcript with longest CDS
		ltrans_cds = cds_data[np.argmax(cds_data[:, 2].astype('i')), 1]

		# Submit query to retrieve sequence
		response = ds.search({
			'filters':{
				'ensembl_transcript_id':ltrans_cds
			},
			'attributes':[
				'ensembl_gene_id',			# Gene ID
				'ensembl_transcript_id',	# Transcript ID
				'coding'					# CDS sequence
			]
		})

		# Convert output
		cds_seq = np.array([row.split('\t')
			for row in response.text.strip().split('\n')])

		# Write
		sout = '> ' + cds_seq[0][1] + '\n' + cds_seq[0][0] + '\n'
		out_trans_cds_fasta.write(sout)
	else:
		print(">> No CDS data on current gene %s ..." % (ensembl_id[idi],))
		out_skip_gene.write(ensembl_id[idi] + '\tNo CDS data found.')

	# Retrieve sequence of transcript with longest CDS+UTRs
	# -----------------------------------------------------

	print("> Retrieving transcript CDS+UTRs (longest) sequence...")

	# Remove transcripts with no CDS+UTRs information
	utr_data = trans_data[np.where('' != trans_data[:, 3])[0], :]
	
	if 0 != utr_data.shape[0]:
		# Identify transcript with longest CDS+UTRs
		ltrans_cds_utr = utr_data[np.argmax(utr_data[:, 3].astype('i')), 1]
	
		# Submit query to retrieve sequence
		response = ds.search({
			'filters':{
				'ensembl_transcript_id':ltrans_cds_utr
			},
			'attributes':[
				'ensembl_gene_id',			# Gene ID
				'ensembl_transcript_id',	# Transcript ID
				'cdna'						# CDS+UTR sequence
			]
		})
	
		# Convert output
		cds_utr_seq = np.array([row.split('\t')
			for row in response.text.strip().split('\n')])
	
		# Write
		sout = '> ' + cds_utr_seq[0][1] + '\n' + cds_utr_seq[0][0] + '\n'
		out_trans_cds_utr_fasta.write(sout)
	else:
		print(">> No CDS+UTRs data on current gene %s ..." % (ensembl_id[idi],))
		out_skip_gene.write(ensembl_id[idi] + '\tNo CDS+UTRs data found.')

	# Exons
	# Query for exon data
	# --------------------------------------------------------------------------

	print("> Retrieving exon sequences...")

	# Submit query
	response = ds.search({
		'filters':{
			'ensembl_gene_id':ensembl_id[idi]
		},
		'attributes':[
			'ensembl_gene_id',	# Gene ID
			'ensembl_exon_id',	# Exon ID
			'gene_exon'			# Exon sequence
		]
	})

	# Convert output
	exon_data = np.array([row.split('\t')
		for row in response.text.strip().split('\n')])

	# Write
	for exon in exon_data:
		sout = "\t".join([exon[1], exon[2], exon[0]]) + "\n"
		out_exon_table.write(sout)

# Close file pointers
out_gene_table.close()
out_gene_fasta.close()
out_trans_table.close()
out_trans_cds_fasta.close()
out_trans_cds_utr_fasta.close()
out_exon_table.close()
out_skip_gene.close()

# END ==========================================================================

################################################################################
