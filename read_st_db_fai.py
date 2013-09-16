# test reading ST database after fai checking
# check if MLST delimiter we are using makes sense
# determine the number of columns to read based on the match to gene names
from argparse import (ArgumentParser, FileType)
import logging
from subprocess import call
import os, sys, re, collections, operator
from scipy.stats import binom_test, linregress
from math import log
from itertools import groupby
from operator import itemgetter

edge_a = edge_z = 2


def parse_args():
	"Parse the input arguments, use '-h' for help"

	parser = ArgumentParser(description='SRST2 - Short Read Sequence Typer (v2)')


	# Read inputs
	parser.add_argument(
		'--input_se', nargs='+',type=str, required=False,
		help='Input single end read')
	parser.add_argument(
		'--input_pe', nargs='+', type=str, required=False,
		help='Input paired end reads')
	parser.add_argument(
		'--forward', type=str, required=False, default="_1", 
			help='Designator for forward reads (e.g default is _1, expect forward reads sample_1.fastq.gz)')
	parser.add_argument(
		'--reverse', type=str, required=False, default="_2", 
			help='Designator for reverse reads (e.g default is _2, expect reverse reads sample_2.fastq.gz)')
	
	# MLST parameters
	parser.add_argument('--mlst_db', type=str, required=False, nargs=1, help='Fasta file of MLST alleles (optional)')
	parser.add_argument('--mlst_delimiter', type=str, required=False,
		help='Character(s) separating gene name from allele number in MLST database (default "-", as in arcc-1)', default="-")
		
	parser.add_argument('--mlst_definitions', type=str, required=False,
		help='ST definitions for MLST scheme (required if mlst_db supplied and you want to calculate STs)')
	parser.add_argument('--ignore_last', action="store_true", required=False,
		help='Flag to ignore last column of ST definitions table (e.g. sometimes an additional column is added to indicate clonal complex, which is not part of the ST definition).')

	# Gene database parameters
	parser.add_argument('--gene_db', type=str, required=False, nargs='+', help='Fasta file/s for gene databases (optional)')

	# Cutoffs for scoring/heuristics
	parser.add_argument('--min_coverage', type=float, required=False, help='Percent coverage cutoff for gene reporting (default 90)',default=90)
	parser.add_argument('--min_depth', type=float, required=False, help='Minimum mean depth to flag as dubious allele call (default 5)',default=5)
	parser.add_argument('--min_edge_depth', type=float, required=False, help='Minimum edge depth to flag as dubious allele call (default 2)',default=2)
	parser.add_argument('--prob_err', type=float, default=0.01, help='Probability of sequencing error (default 0.01)')

	# Mapping parameters for bowtie2
	parser.add_argument('--read_type', type=str, choices=['q', 'qseq', 'f'], default='q',
		help='Input file type (for bowtie input; default is q=fastq; other options: qseq=solexa, f=fasta).')
	parser.add_argument('--other', type=str, help='Other options for bowtie2.', required=False) 
	parser.add_argument('--mapq', type=int, default=1, help='Qual map') # ADD EXPLANATION
	parser.add_argument('--baseq', type=int, default=20, help='Qual base') # ADD EXPLANATION
	
	# Reporting options
	parser.add_argument('--output', type=str, required=True, help='Output file prefix')
	parser.add_argument('--log', metavar='FILE', type=str,
		help='log progress in FILENAME, defaults to stdout')
	parser.add_argument('--verbose', action="store_true", required=False, help='Switch on verbose reporting')

	# Run options
	parser.add_argument('--use_existing_pileup', action="store_true", required=False,
		help='Use existing pileups if available, otherwise they will be generated') # to facilitate testing of rescoring from pileups
	parser.add_argument('--use_existing_scores', action="store_true", required=False,
		help='Use existing scores files if available, otherwise they will be generated') # to facilitate testing of reporting from scores

	# Compile previous output files
	parser.add_argument('--prev_output', nargs='+', type=str, required=False,
		help='SRST2 output files to compile (any new results from this run will also be incorporated)')

	return parser.parse_args() 

def parse_ST_database(ST_filename,gene_names_from_fai):
	# Read ST definitions
	ST_db = {} # key = allele string, value = ST
	gene_names = []
	num_gene_cols_expected = len(gene_names_from_fai)
	with open(ST_filename) as f:
		count = 0
		for line in f:
			count += 1
			line_split = line.rstrip().split("\t")
			if count == 1: # Header
				gene_names = line_split[1:num_gene_cols_expected+1]
				for g in gene_names:
					if g not in gene_names_from_fai:
						print "Warning: gene " + g + " in ST definitions file isn't among those in the database " + ",".join(gene_names_from_fai)
						print " This will result in all STs being called as unknown (but allele calls will be accurate for other loci)." 
				for g in gene_names_from_fai:
					if g not in gene_names:
						print "Warning: gene " + g + " in database file isn't among those in the ST definitions: " + ",".join(gene_names)
						print " Any sequences with this gene identifer from the database will not be included in typing."
			else:
				ST = line_split[0]
				if ST not in ST_db.values():
					ST_string = " ".join(line_split[1:num_gene_cols_expected+1])
					ST_db[ST_string] = ST
				else:
					print "Warning: this ST is not unique in the ST definitions file: " + ST
		return (ST_db, gene_names)

def parse_fai(fai_file,db_type,delimiter):
	'Get sequence lengths for reference alleles - important for scoring'
	'Get gene names also, required if no MLST definitions provided'
	size = {}
	gene_clusters = [] # for gene DBs, this is cluster ID
	allele_symbols = []
	gene_cluster_symbols = {} # key = cluster ID, value = gene symbol (for gene DBs)
	unique_allele_symbols = True
	unique_gene_symbols = True
	delimiter_check = [] # list of names that may violate the MLST delimiter supplied
	with open(fai_file) as fai:
		for line in fai:
			fields = line.split('\t')
			name = fields[0] # full allele name
			size[name] = int(fields[1]) # store length
			if db_type!="mlst":
				allele_info = name.split()[0].split("__")
				gene_cluster = allele_info[0] # ID number for the cluster
				cluster_symbol = allele_info[1] # gene name for the cluster
				name = allele_info[2] # specific allele name
				if gene_cluster in gene_cluster_symbols:
					if gene_cluster_symbols[gene_cluster] != cluster_symbol:
						unique_gene_symbols = False # already seen this cluster name
				else:
					gene_cluster_symbols[gene_cluster] = cluster_symbol
			else:
				gene_cluster = name.split(delimiter)[0] # accept gene clusters raw for mlst
				# check if the delimiter makes sense
				parts = name.split(delimiter)
				if len(parts) != 2:
					delimiter_check.append(name)
				else:
					try:
						x = int(parts[1])
					except:
						delimiter_check.append(name)
			
			# check if we have seen this allele name before
			if name in allele_symbols:
				unique_allele_symbols = False # already seen this allele name
			allele_symbols.append(name)
			
			# record gene:
			if gene_cluster not in gene_clusters:
				gene_clusters.append(gene_cluster)

	if len(delimiter_check) > 0:
		print "Warning! MLST delimiter is " + delimiter + " but these genes may violate the pattern and cause problems:"
		for x in delimiter_check:
			print " " + x
			
	return size, gene_clusters, unique_gene_symbols, unique_allele_symbols

def main():
	args = parse_args()
	
	fai_file = args.mlst_db[0] + '.fai'
	print fai_file
	if not os.path.exists(fai_file):
		run_command(['samtools', 'faidx', args.mlst_db])
	
	size, gene_names, unique_gene_symbols, unique_allele_symbols = parse_fai(fai_file,"mlst",args.mlst_delimiter)
	
	print gene_names
		
	ST_db, gene_names = parse_ST_database(args.mlst_definitions, gene_names)
	
	print gene_names
	print ST_db
	
if __name__ == '__main__':
	main()
