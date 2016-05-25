#!/usr/bin/env python

# Python Version 2.7.5
#
# Basic analysis of SRST2 output
# Authors - Kathryn Holt (kholt@unimelb.edu.au)
#
# Dependencies:
#	R

from argparse import (ArgumentParser, FileType)
import logging
from subprocess import call, check_output, CalledProcessError, STDOUT
import os, sys, re, collections, operator
from scipy.stats import binom_test, linregress
from math import log
from itertools import groupby
from operator import itemgetter
from collections import OrderedDict

def parse_args():
	"Parse the input arguments, use '-h' for help"

	parser = ArgumentParser(description='Basic analysis of SRST2 output')

	# options
	parser.add_argument('--prev_output', nargs='+', type=str, required=False, help='SRST2 results files to compile')
	parser.add_argument('--sample_metadata', nargs='+', type=str, required=False, help='Sample metadata to merge into table (tab delimited)')	
	parser.add_argument('--output', type=str, required=True, help='Prefix for the output files')
	parser.add_argument('--log', action="store_true", required=False, help='Switch on logging to file (otherwise log to stdout)')
	
	return parser.parse_args() 

# take list of mlst results dictionaries and/or list of gene results dictionaries
# print table of mlst results + gene results for all samples
def compile_results(args,mlst_results,db_results,sample_metadata_hashes,compiled_output_file):

	o = file(compiled_output_file,"w")
	
	# get list of all samples and genes present in these datasets
	sample_list = [] # each entry is a sample present in at least one db
	gene_list = [] # genes read from gene db results
	variable_list = [] # metadata variables
	mlst_cols = 0
	mlst_header_string = ""
	blank_mlst_section = ""
	
	mlst_results_master = {} # compilation of all MLST results
	db_results_master = collections.defaultdict(dict) # compilation of all gene results
	metadata_master = collections.defaultdict(dict) # compilation of all metadata
	
	st_counts = {} # key = ST, value = count
	
	# store mlst data from dict (as string)
	if len(mlst_results) > 0:
	
		for mlst_result in mlst_results:
		
			# check length of the mlst string
			if "Sample" in mlst_result:
				test_string = mlst_result["Sample"]
				if mlst_cols == 0:
					mlst_header_string = test_string
			else:
				test_string = mlst_result[mlst_result.keys()[0]] # no header line?
			test_string_split = test_string.split("\t")
			this_mlst_cols = len(test_string)
			
			if (mlst_cols == 0) or (mlst_cols == this_mlst_cols):
				mlst_cols = this_mlst_cols
				blank_mlst_section = "\t" * (mlst_cols-1) # blank MLST string in case some samples missing
				# use this data
				for sample in mlst_result:
					mlst_results_master[sample] = mlst_result[sample]
					if sample not in sample_list:
						sample_list.append(sample)
			elif mlst_cols != this_mlst_cols:
				# don't process this data further
				logging.info("Problem reconciling MLST data from, first MLST results encountered had " + str(mlst_cols) + " columns, this one has " + str(this_mlst_cols) + " columns?")
				if args.mlst_db:
					logging.info("Compiled report will contain only the MLST data from this run, not previous outputs")
				else:
					logging.info("Compiled report will contain only the data from the first MLST result set provided")

	# store gene data from dict (as dict of dicts)
	if len(db_results) > 0:
		for results in db_results:
			for sample in results:
				if sample not in sample_list:
					sample_list.append(sample)
				for gene in results[sample]:
					if gene != "failed":
						db_results_master[sample][gene] = results[sample][gene]
						if gene not in gene_list:
							gene_list.append(gene)
							
	# store metadata from dict (as dict of dicts)						
	if len(sample_metadata_hashes) > 0:
		for results in sample_metadata_hashes:
			for sample in results:
				if sample not in sample_list:
					sample_list.append(sample)
				for variable in results[sample]:
#					print metadata_master
#					metadata_master[sample] = results[sample][variable]
					metadata_master[sample][variable] = results[sample][variable]
					if variable not in variable_list:
						variable_list.append(variable)
						
	print variable_list
						
	if "Sample" in sample_list:
		sample_list.remove("Sample")
	sample_list.sort()
	gene_list.sort()
	
	# print header
	header_elements = []
	if len(mlst_results) > 0:
		header_elements.append(mlst_header_string) # mlst column headings
	else:
		header_elements.append("Sample")
	if len(variable_list) > 0:
		header_elements += variable_list # metadata column headings
	if len(gene_list) > 0:
		header_elements += gene_list # gene column headings
	o.write("\t".join(header_elements)+"\n")
	
	# print results for all samples
	for sample in sample_list:
	
		sample_info = [] # first entry is mlst string OR sample name, rest are metadata or genes
		
		# print mlst if provided, otherwise just print sample name
		if len(mlst_results_master) > 0:
			if sample in mlst_results_master:
				sample_info.append(mlst_results_master[sample])
				this_st = mlst_results_master[sample].split("\t")[1]
			else:
				sample_info.append(sample+blank_mlst_section)
				this_st = "unknown"
			# record the MLST result				
			if this_st in st_counts:
				st_counts[this_st] += 1
			else:
				st_counts[this_st] = 1
		else:
			sample_info.append(sample)
		
		# get metadata if provided
		if sample in metadata_master:
			for variable in variable_list:
				if variable in metadata_master[sample]:
					sample_info.append(metadata_master[sample][variable])
				else:
					sample_info.append("-")
		else:
			for variable in variable_list:
				sample_info.append("?") # record no gene data on this strain
				
		# get gene info if provided
		if sample in db_results_master:
			for gene in gene_list:
				if gene in db_results_master[sample]:
					sample_info.append(db_results_master[sample][gene])
				else:
					sample_info.append("-")
		else:
			for gene in gene_list:
				sample_info.append("?") # record no gene data on this strain
			
		o.write("\t".join(sample_info)+"\n")
		
	o.close()
	
	logging.info("Compiled data on " + str(len(sample_list)) + " samples printed to: " + compiled_output_file)
	
	# log ST counts
	if len(mlst_results_master) > 0:
		logging.info("Detected " + str(len(st_counts.keys())) + " STs: ")
		sts = st_counts.keys()
		sts.sort()
		for st in sts:
			logging.info("ST" + st + "\t" + str(st_counts[st]))
	
	return True
	
# read data from file into dict
def read_results_from_file(infile,metadata):
	
	if metadata:
	
		# process as with a genes table
		results = collections.defaultdict(dict) # key1 = sample, key2 = gene, value = allele
		with open(infile) as f:
			header = []
			for line in f:
				line_split = line.rstrip().split("\t")
				if len(header) == 0:
					header = line_split
				else:
					sample = line_split[0]
					for i in range(1,len(line_split)):
						variable = header[i]
						results[sample][variable] = line_split[i]

		dbtype = "metadata"
		dbname = infile
		logging.info("Reading metadata from: " + infile)
				
	else:
	
		results_info = infile.split("__")
		
		if len(results_info) > 1:
	
			if results_info[-1] == "compiledResults.txt":
				dbtype = "compiled"
				dbname = results_info[0] # output identifier
			else:
				dbtype = results_info[1] # mlst or genes
				dbname = results_info[2] # database

			logging.info("Processing " + dbtype + " results from file " + infile)
	
			if dbtype == "genes":	
				results = collections.defaultdict(dict) # key1 = sample, key2 = gene, value = allele
				with open(infile) as f:
					header = []
					for line in f:
						line_split = line.rstrip().split("\t")
						if len(header) == 0:
							header = line_split
						else:
							sample = line_split[0]
							for i in range(1,len(line_split)):
								gene = header[i]
								results[sample][gene] = line_split[i]

			elif dbtype == "mlst":	
				results = {} # key = sample, value = MLST string
				with open(infile) as f:
					for line in f:
						results[line.split("\t")[0]] = line.rstrip() # store header line too (index "Sample")
				
			elif dbtype == "compiled":
				results = collections.defaultdict(dict) # key1 = sample, key2 = gene, value = allele
				with open(infile) as f:
					header = []
					mlst_cols = 0 # INDEX of the last mlst column
					n_cols = 0
					for line in f:
						line_split = line.rstrip().split("\t")
						if len(header) == 0:
							header = line_split
							n_cols = len(header)
							if header[1] == "ST":
								# there is mlst data reported
								mlst_cols = 2 # first locus column
								while header[mlst_cols] != "depth":
									mlst_cols += 1
								results["Sample"]["mlst"] = "\t".join(line_split[0:(mlst_cols+1)])
						else:
							sample = line_split[0]
							if mlst_cols > 0:
								results[sample]["mlst"] = "\t".join(line_split[0:(mlst_cols+1)])
							if n_cols > mlst_cols:
								# read genes component
								for i in range(mlst_cols+1,n_cols):
									# note i=1 if mlst_cols==0, ie we are reading all
									gene = header[i]
									if len(line_split) > i:
										results[sample][gene] = line_split[i]
									else:
										results[sample][gene] = "-"
							
		else:
			results = False
			dbtype = False
			dbname = False
			logging.info("Couldn't decide what to do with file results file provided: " + infile)
				
	return results, dbtype, dbname

	
def main():
	args = parse_args()
	if args.log is True:
		logfile = args.output + ".log"
	else:
		logfile = None
	logging.basicConfig(
		filename=logfile,
		level=logging.DEBUG,
		filemode='w',
		format='%(asctime)s %(message)s',
		datefmt='%m/%d/%Y %H:%M:%S')
	logging.info('program started')
	logging.info('command line: {0}'.format(' '.join(sys.argv)))

	# vars to store results
	mlst_results_hashes = [] # dict (sample->MLST result string) for each MLST output files created/read
	gene_result_hashes = [] # dict (sample->gene->result) for each gene typing output files created/read
	sample_metadata_hashes = [] # dict (sample->variable->value) for each variable read from metadata file

	# read in results from prior results files
	if args.prev_output:
	
		unique_results_files = list(OrderedDict.fromkeys(args.prev_output))
	
		for results_file in unique_results_files:
		
			results, dbtype, dbname = read_results_from_file(results_file, False)
			
			if dbtype == "mlst":
				mlst_results_hashes.append(results)
				
			elif dbtype == "genes":
				gene_result_hashes.append(results)
				
			elif dbtype == "compiled":
				# store mlst in its own db
				mlst_results = {}
				for sample in results:
					if "mlst" in results[sample]:
						mlst_results[sample] = results[sample]["mlst"]
						del results[sample]["mlst"]
				mlst_results_hashes.append(mlst_results)
				gene_result_hashes.append(results)
				
	# read in metadata from file
	if args.sample_metadata:
		for metadata_file in args.sample_metadata:
			results, dbtype, dbname = read_results_from_file(metadata_file, True)
			sample_metadata_hashes.append(results)

	# compile results if multiple databases or datasets provided
	if ( (len(gene_result_hashes) + len(mlst_results_hashes) + len(sample_metadata_hashes)) > 1 ):
		compiled_output_file = args.output + "__compiledResults.txt"
		compile_results(args,mlst_results_hashes,gene_result_hashes,sample_metadata_hashes,compiled_output_file)
	
	elif args.prev_output:
		logging.info('One previous output file was provided, but there is no other data to compile with.')
	
	logging.info('SRST2 has finished.')


if __name__ == '__main__':
	main()
