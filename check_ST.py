#!/usr/bin/env python

# Parses scores from SRST(v2) and compares these with known sequence types for those samples
# 
# Author - Harriet Dashnow (h.dashnow@gmail.com)

import sys
import csv
import itertools
import os

validation_file = "/vlsci/VR0082/shared/hdashnow/testing_srst2/data/shigella-samples-1.csv"
ST_file = "/vlsci/VR0082/shared/hdashnow/testing_srst2/data/ecoli_ST.txt"
scores_file = "/vlsci/VR0082/shared/hdashnow/testing_srst2/mlst/pool10_tag2_1.fastq.gz.MLST_alleles.fasta.srst2.pileup.table.scores"

def extract_top_allele(scores_file):
	sample_name = scores_file.split(".")[0].split("/")[-1]
	hash_alleles = {}
	with open(scores_file) as scores:
		count = 0
		for line in scores:
			count += 1
			if count == 1:
				continue # Skip header
			split_line = line.split()
			allele_info = split_line[0].split("-")
			gene = allele_info[0].upper()
			allele = allele_info[1]
			score = float(split_line[1])
			coverage = float(split_line[5])

			if gene not in hash_alleles:
				hash_alleles[gene] = (allele, score, coverage)
			elif score < hash_alleles[gene][1]:
				hash_alleles[gene] = (allele, score, coverage)

	return (sample_name, hash_alleles)

def parse_ST_database(ST_filename):
	# Prepare ST database
	# Sort by gene name (alphabetical) and save alleles as space delimted string
	ST_db = {}
	gene_names = ""
	with open(ST_filename) as f:
		count = 0
		for line in f:
			count += 1
			if count == 1: # Header
				gene_names = line.split()[1:-1]
				gene_names = [ x.upper() for x in gene_names ]
				num_genes = len(gene_names)
			else:
				line_split = line.split()
				ST = line_split[0]
				alleles = " ".join(line_split[1:num_genes+1])
				if alleles not in ST_db:
					ST_db[alleles] = ST
#HD				else:
#HD					print "Warning, this combination of alleles is not unique in the database:", alleles
#HD					print "The first instance (#{0}) was used.".format(ST)
		return (ST_db, gene_names)

def get_ST(hash_alleles, ST_db, gene_names, sample_name):
	# Put sample in the correct format for comparision with the ST_db
	sample_alleles = []
	for gene in gene_names:
		if gene not in hash_alleles:
			#HD print hash_alleles
			#XXX This should produce an error
			print "Gene {0} was not found in sample {1} so sequence type cannot be identified".format(gene, sample_name)
			return None

		sample_alleles.append(hash_alleles[gene][0])
	sample_string = " ".join(sample_alleles)
	try:
		ST = ST_db[sample_string]
	except KeyError:
		#XXX Better output for this error?
		print "This combination of alleles was not found in the sequence type database:"
		print sample_name
		for key in hash_alleles:
			print key, hash_alleles[key][0]
		ST = None

	return ST
	
def parse_validation_file(validation_file):
	with open(validation_file, "rU") as f:
		count = 0
		validation_db = {}
		for line in csv.reader(f):
			count += 1
			if count == 1:
				headers = line[1:]
			else:
				sample = line[0]
				if sample not in validation_db:
					validation_db[sample] = {}
				col = 1
				for key in headers:
					validation_db[sample][key] = line[col]
					col += 1
	return validation_db

# Get databases of sequence types and validation data
ST_db, gene_names = parse_ST_database(ST_file)
validation_db = parse_validation_file(validation_file)

# Get list of all the .scores files in a given directory
directory = "/vlsci/VR0082/shared/hdashnow/testing_srst2/mlst/"
score_files = [ (directory + f) for f in os.listdir(directory) if f.endswith(".scores")]

with open("SRST2_valdiation.txt","w") as outfile:
	outfile.write("sample\tSRST2_ST\tvalidation_ST\tcorrect\n")
	
	for scores_file in score_files:
		sample_name, hash_alleles = extract_top_allele(scores_file) # Do this for all the score files
		srst2_ST = get_ST(hash_alleles, ST_db, gene_names, sample_name)
		sample = sample_name[:-2]
		validation_ST = validation_db[sample]["ST"]

		outstring = "{0}\t{1}\t{2}\t{3}\n".format(sample, srst2_ST, validation_ST, 
													srst2_ST == validation_ST)
		outfile.write(outstring)







