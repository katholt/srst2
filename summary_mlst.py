#!/usr/bin/env python

# Parses scores from SRST(v2) to produce a summary table of multple samples. 
# Assigns MLST sequence type based on a pubmlst database.
# 
# Author - Harriet Dashnow (h.dashnow@gmail.com)

import sys
import csv
import itertools
import os
from argparse import (ArgumentParser, FileType)

def parse_args():
	"Parse the input arguments, use '-h' for help"
	parser = ArgumentParser(description='SRST2 - Short Read Sequence Typer (v2): score summary')
	parser.add_argument(
		'--scores', nargs='+', type=str, required=True,
		help='One or more .scores files produced by srst2.py')
	parser.add_argument('--output', type=str, required=True, help='Output file')
	parser.add_argument(
		'--mlst', type=str, required=False,
		help='MLST definition profile e.g. from http://pubmlst.org/data/')
	return parser.parse_args() 

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
			mismatches = int(split_line[6])
			try:
				indels = int(split_line[7])
				missing = int(split_line[8])
				difference = abs(mismatches) #XXX Should be single base mismatches + indels, however these are not scored separately in srst2.py
				if missing > 0:
					difference += 1 # Score any truncation as a single difference
			except IndexError:
				difference = abs(mismatches)


			if gene not in hash_alleles:
				hash_alleles[gene] = (allele, score, coverage, difference)
			elif score < hash_alleles[gene][1]:
				hash_alleles[gene] = (allele, score, coverage, difference)

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
				#else:
					#print "Warning: this combination of alleles is not unique in the database:", alleles+".", "The first instance (#{0}) was used.".format(ST) #XXX This warning will occur for e.coli - make a supress warning option? I have commented out for now.
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

def genes_table(score_files):
	all_genes = []
	hash_samples = {}
	for scores_file in score_files:
		sample_name, hash_alleles = extract_top_allele(scores_file)
		for gene_name in hash_alleles:
			if gene_name not in all_genes:
				all_genes.append(gene_name)
		if sample_name not in hash_samples:
			hash_samples[sample_name] = hash_alleles
	return all_genes, hash_samples

def main():
	args = parse_args()
	ST_file = args.mlst
	score_files = args.scores
	
	if args.mlst != None:
		# Get databases of sequence types
		ST_db, gene_names = parse_ST_database(ST_file)

		with open(args.output, "w") as outfile:
			header = "sample\tST"
			for gene_name in gene_names:
				header = header + "\t" + gene_name
			header = header + "\n"
			outfile.write(header)
	
			for scores_file in score_files:
				sample_name, hash_alleles = extract_top_allele(scores_file)
				srst2_ST = get_ST(hash_alleles, ST_db, gene_names, sample_name)
				sample = sample_name[:-2]

				outstring = "{0}\t{1}".format(sample, srst2_ST)
				for gene_name in gene_names:
					try:
						allele_name = hash_alleles[gene_name][0]
						difference = hash_alleles[gene_name][3]
						if difference > 0:
							allele_name = allele_name + "*/" + str(difference)
						outstring = outstring + "\t" + str(allele_name)
						
					except KeyError:
						outstring = outstring + "\tNone"
				outstring = outstring + "\n"
				outfile.write(outstring)
	else:
		all_genes, hash_samples = genes_table(score_files)
		with open(args.output, "w") as outfile:
			header = "sample"
			for gene_name in all_genes:
				header = header + "\t" + gene_name
			header = header + "\n"
			outfile.write(header)
			for sample_name in hash_samples:
				outstring = sample_name
				for gene_name in all_genes:
					try:
						allele = hash_samples[sample_name][gene_name][0]
						difference = hash_samples[sample_name][gene_name][3]
						if difference > 0:
							allele = allele + "*/" + str(difference)
					except KeyError:
						allele = None
					outstring = outstring + "\t" + str(allele)
				outstring = outstring + "\n"
				outfile.write(outstring)
					


if __name__ == '__main__':
	main()

# Some testing peramters:
#validation_file = "/vlsci/VR0082/shared/hdashnow/testing_srst2/data/shigella-samples-1.csv"
#ST_file = "/vlsci/VR0082/shared/hdashnow/testing_srst2/data/ecoli_ST.txt"
#scores_file = "/vlsci/VR0082/shared/hdashnow/testing_srst2/mlst/pool10_tag2_1.fastq.gz.MLST_alleles.fasta.srst2.pileup.table.scores"





