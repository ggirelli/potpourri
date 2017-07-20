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

from biomart import BiomartServer # by ggirelli, with request by POST
import numpy as np
import pandas as pd
import re
import sys

# PARAMETERS ===================================================================

# INPUT
# -----

# Gene list tsv (single-column ENSG IDs)
gene_list_file = "genelist.tsv"

# OUTPUT
# ------
out_gene_table_file = "gene_data.tsv"
out_gene_fasta_file = "gene_seq.fa"
out_trans_table_file = "trans_data.tsv"
out_trans_cds_fasta_file = "trans_cds_seq.fa"
out_trans_cds_utr_fasta_file = "trans_cds_utr_seq.fa"
out_exon_table_file = "exon_seq.tsv"
out_exon_fasta_file = "exon_seq.fa"

# FUNCTIONS ====================================================================

# RUN ==========================================================================

# GENE LIST
# ------------------------------------------------------------------------------

# Read gene list
gene_list = pd.read_csv(gene_list_file, header = None)[0].tolist()

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

# Gene
# Query for Gene characteristics
# --------------------------------------------------------------------------

# Retrieve gene information
# -------------------------

print("> Retrieving gene information...")

# Submit query
response = ds.search({
	'filters':{
		'ensembl_gene_id':gene_list
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
print(" · Reading response...")
gene_data = np.array([row.split('\t')
	for row in response.text.strip().split('\n')])

# Write
print(" · Writing output...")
sout = "\n".join(["\t".join(row) for row in gene_data]) + "\n"
out_gene_table = open(out_gene_table_file, 'w')
out_gene_table.write(sout)
out_gene_table.close()

# Retrieve the whole gene sequence
# --------------------------------

print("> Retrieving gene sequence...")

# Submit query
response = ds.search({
	'filters':{
		'ensembl_gene_id':gene_list
	},
	'attributes':[
		'ensembl_gene_id',				# Gene ID
		'gene_exon_intron'				# Gene cDNA
	]
})

# Convert output
print(" · Reading response...")
gene_seq = np.array([row.split('\t')
	for row in response.text.strip().split('\n')])

# Write
print(" · Writing output...")
sout = "\n".join(["> %s\n%s" % (row[1], row[0],) for row in gene_seq]) + "\n"
out_gene_fasta = open(out_gene_fasta_file, 'w')
out_gene_fasta.write(sout)
out_gene_fasta.close()

# Transcripts
# Query for transcript data
# --------------------------------------------------------------------------

# Identify transcripts
# --------------------

print("> Retrieving transcript information...")

# Submit query
response = ds.search({
	'filters':{
		'ensembl_gene_id':gene_list
	},
	'attributes':[
		'ensembl_gene_id',			# Gene ID
		'ensembl_transcript_id',	# Transcript ID
		'cds_length',				# CDS length
		'transcript_length'			# CDS+UTRs length
	]
}, header = 1)

# Convert output
print(" · Reading response...")
trans_data = np.array([row.split('\t')
	for row in response.text.strip().split('\n')])

# Write
print(" · Writing output...")
sout = "\n".join(["\t".join(row) for row in trans_data]) + "\n"
out_trans_table = open(out_trans_table_file, 'w')
out_trans_table.write(sout)
out_trans_table.close()

# Retrieve sequence of transcript with longest CDS
# ------------------------------------------------

print("> Retrieving transcript CDS (longest) sequence...")

ltrans_cds = []
for gene in gene_list:
	#print(" · Working on '%s'..." % (gene,))
	# Retreive CDS data for specific gene
	cds_data = trans_data[np.where(np.array(trans_data)[:,0] == gene)[0],1:3]

	# Set to 0 absent CDS data
	cds_data[np.where(cds_data[:,1] == '')[0], 1] = '0'

	# Identify transcript with longest CDS
	ltrans_cds.append(cds_data[np.argmax(cds_data[:, 1].astype('i')), 0])

# Submit query to retrieve sequence
print(" >>> Querying...")
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
print(" >>> Reading response...")
cds_seq = np.array([row.split('\t')
	for row in response.text.strip().split('\n')])

# Write
print(" · Writing output...")
sout = "\n".join(["> %s_%s\n%s" % (row[1], row[2], row[0],) for row in cds_seq]) + "\n"
out_trans_cds_fasta = open(out_trans_cds_fasta_file, 'w')
out_trans_cds_fasta.write(sout)
out_trans_cds_fasta.close()

# Retrieve sequence of transcript with longest CDS+UTRs
# -----------------------------------------------------

print("> Retrieving transcript CDS+UTRs (longest) sequence...")

ltrans_utr = []
for gene in gene_list:
	#print(" · Working on '%s'..." % (gene,))
	# Retreive CDS data for specific gene
	utr_data = trans_data[np.where(np.array(trans_data)[:,0] == gene)[0],:][:,[1, 3]]
	
	# Set to 0 absent UTR data
	utr_data[np.where(utr_data[:,1] == '')[0], 1] = '0'
	
	# Identify transcript with longest CDS
	ltrans_utr.append(utr_data[np.argmax(utr_data[:, 1].astype('i')), 0])

# Submit query to retrieve sequence
print(" >>> Querying...")
response = ds.search({
	'filters':{
		'ensembl_transcript_id':ltrans_utr
	},
	'attributes':[
		'ensembl_gene_id',			# Gene ID
		'ensembl_transcript_id',	# Transcript ID
		'cdna'					# CDS sequence
	]
})

# Convert output
print(" >>> Reading response...")
utr_seq = np.array([row.split('\t')
	for row in response.text.strip().split('\n')])

# Write
print(" · Writing output...")
sout = "\n".join(["> %s_%s\n%s" % (row[1], row[2], row[0],) for row in utr_seq]) + "\n"
out_trans_cds_utr_fasta = open(out_trans_cds_utr_fasta_file, 'w')
out_trans_cds_utr_fasta.write(sout)
out_trans_cds_utr_fasta.close()

# Exons
# Query for exon data
# --------------------------------------------------------------------------

print("> Retrieving exon sequences...")

# Submit query
response = ds.search({
	'filters':{
		'ensembl_gene_id':gene_list
	},
	'attributes':[
		'ensembl_gene_id',	# Gene ID
		'ensembl_exon_id',	# Exon ID
		'gene_exon'			# Exon sequence
	]
})

# Convert output
print(" >>> Reading response...")
exon_data = np.array([row.split('\t')
	for row in response.text.strip().split('\n')])

# Write
print(" · Writing output...")
out_exon_table = open(out_exon_table_file, 'w')
sout = "\n".join(["\t".join([exon[1], exon[2], exon[0]]) for exon in exon_data])
out_exon_table.write(sout)
out_exon_table.close()
sout = "\n".join(["> %s_%s\n%s" % (row[1], row[2], row[0],) for row in exon_data]) + "\n"
out_exon_fasta = open(out_exon_fasta_file, 'w')
out_exon_fasta.write(sout)
out_exon_fasta.close()
# END ==========================================================================

################################################################################
