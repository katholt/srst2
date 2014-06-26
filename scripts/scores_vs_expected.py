#!/usr/bin/env python

# Parses scores from SRST(v2) to produce a summary table of the top scoring allele for each locus/cluster for multiple samples. 
# If a table of known/expected alleles are provided, this info is also included in the scores summary output.
# 
# Author - Harriet Dashnow (h.dashnow@gmail.com), Kathryn Holt (kholt@unimelb.edu.au)

import sys, os, collections, operator
from argparse import (ArgumentParser, FileType)

def parse_args():
	"Parse the input arguments, use '-h' for help"
	
	parser = ArgumentParser(description='SRST2 - Short Read Sequence Typer (v2): Outputs score and gene summaries. Reference fasta must have been in format cluster#__gene__allele')
	parser.add_argument(
		'--scores', nargs='+', type=str, required=True,
		help='One or more .scores files produced by srst2.py')
	parser.add_argument('--output', type=str, required=True,
		help='Prefix for output tables summarising the top scoring allele for each gene with more detailed score information, and ST table.')
	parser.add_argument('--mlst_delimiter', type=str, required=False, help='Character that separates locus symbol from allele number (e.g. "-")',default = "-")
	parser.add_argument('--known_alleles', type=str, required=False, help='File of known alleles for comparison (tab-delimited)')
	parser.add_argument('--known_STs', type=str, required=False, help='File of known STs for comparison (tab-delimited), must also provide ST definitions table (--mlst)')
	parser.add_argument('--ST_definitions', type=str, required=False, help='Table of ST definitions')
	parser.add_argument('--ignore_last', action="store_true", required=False,
		help='Ignore last column of ST profiles table (e.g. sometimes an additional column is added to indicate clonal complex, which is not part of the ST definition).')

	# Cutoffs for scoring heuristics
	parser.add_argument('--min_coverage', type=float, required=False, help='Percent coverage cutoff for gene reporting (default 90)',default=90)
	parser.add_argument('--min_depth', type=float, required=False, help='Minimum mean depth to flag as dubious allele call (default 5)',default=5)
	parser.add_argument('--min_edge_depth', type=float, required=False, help='Minimum edge depth to flag as dubious allele call (default 2)',default=2)

	return parser.parse_args() 
		
def parse_known_alleles(knowns_filename):
	# Read file of known STs and alleles and store allele numbers
	known_alleles_db = collections.defaultdict(dict) # key1 = sample, key2 = gene name (or ST), value = allele/ST number
	gene_names = []
	with open(knowns_filename) as f:
		count = 0
		for line in f:
			count += 1
			if count == 1: # Header
				gene_names = line.split()[2:]
				num_genes = len(gene_names)
			else:
				line_split = line.split()
				sample = line_split[0]
				if sample not in known_alleles_db:
					known_alleles_db[sample]["ST"] = line_split[1]
					for i in range(0,len(gene_names)):
						gene = gene_names[i]
						allele = line_split[i+2]
						known_alleles_db[sample][gene]=allele
				else:
					print "Warning: this sample is not unique in the known alleles file: " + sample
		return (known_alleles_db, gene_names)
		
def parse_ST_database(ST_filename,ignore_last):
	# Read ST definitions
	ST_db = collections.defaultdict(dict) # key1 = sample, key2 = gene name, value = allele number
	ST_dict = {} # key = allele string, value = ST
	gene_names = []
	with open(ST_filename) as f:
		count = 0
		for line in f:
			count += 1
			line_split = line.rstrip().split("\t")
			if count == 1: # Header
				if ignore_last:
					gene_names = line_split[1:-1]
				else:
					gene_names = line_split[1:]
				num_genes = len(gene_names)
			else:
				ST = line_split[0]
				if ST not in ST_db:
					for i in range(0,len(gene_names)):
						gene = gene_names[i]
						allele = line_split[i+1]
						ST_db[ST][gene]=allele
					ST_string = " ".join(line_split[1:len(gene_names)+1])
					ST_dict[ST_string] = ST
				else:
					print "Warning: this ST is not unique in the ST definitions file: " + ST

		return ST_db, ST_dict, gene_names
		
def parse_known_STs(known_STs, ST_db, gene_names):
	# Read file of known STs and store allele numbers
	known_alleles_db = collections.defaultdict(dict) # key1 = sample, key2 = gene name (or ST), value = allele/ST number
	with open(known_STs) as f:
		count = 0
		for line in f:
			count += 1
			line_split = line.split()
			if count > 1:
				sample = line_split[0]
				ST = line_split[1]
				known_alleles_db[sample]["ST"] = ST
				if ST in ST_db:
					for gene in gene_names:
						known_alleles_db[sample][gene]=ST_db[ST][gene] # assign allele numbers based on this ST value
				else:
					for gene in gene_names:
						known_alleles_db[sample][gene]="Unknown"
		return (known_alleles_db)


def read_scores_file(scores_file):
	hash_edge_depth = {}
	avg_depth_allele = {}
	coverage_allele = {}
	mismatch_allele = {}
	indel_allele = {}
	missing_allele = {}
	size_allele = {}
	next_to_del_depth_allele = {}
	mix_rates = {}
	scores = {}
	allele_info_strings = {} ### NOTE THIS LINE DOES NOT APPEAR IN THE srst2.py VERSION OF THE FUNCTION
				
	f = file(scores_file,"r")
	
	for line in f:
		line_split = line.rstrip().split("\t")
		allele = line_split[0]
		if allele != "Allele": # skip header row
			scores[allele] = float(line_split[1])
			mix_rates[allele] = float(line_split[11])
			avg_depth_allele[allele] = float(line_split[2])
			hash_edge_depth[allele] = (float(line_split[3]),float(line_split[4]))
			coverage_allele[allele] = float(line_split[5])
			size_allele[allele] = int(line_split[6])
			mismatch_allele[allele] = int(line_split[7])
			indel_allele[allele] = int(line_split[8])
			missing_allele[allele] = int(line_split[9])
			next_to_del_depth = line_split[10]
			next_to_del_depth_allele[allele] = line_split[10]
			allele_info_strings[allele] = line_split ### NOTE THIS LINE DOES NOT APPEAR IN THE srst2.py VERSION OF THE FUNCTION
			
	return hash_edge_depth, avg_depth_allele, coverage_allele, mismatch_allele, indel_allele, \
			missing_allele, size_allele, next_to_del_depth_allele, scores, mix_rates, allele_info_strings ### NOTE THE LAST OBJECT IS NOT PASSED IN srst2.py
									
def parse_scores(args, scores, hash_edge_depth, 
					avg_depth_allele, coverage_allele, mismatch_allele, indel_allele,  
					missing_allele, size_allele, next_to_del_depth_allele):
					
	# sort into hash for each gene locus
	scores_by_gene = collections.defaultdict(dict) # key1 = gene, key2 = allele, value = score
	
	for allele in scores:
		if coverage_allele[allele] > args.min_coverage:
			allele_info = allele.split(args.mlst_delimiter)
			scores_by_gene[allele_info[0]][allele] = scores[allele]
	
	# determine best allele for each gene locus/cluster
	results = {} # key = gene, value = (allele,diffs,depth)
	
	for gene in scores_by_gene:
	
		gene_hash = scores_by_gene[gene]
		scores_sorted = sorted(gene_hash.iteritems(),key=operator.itemgetter(1)) # sort by score
		(top_allele,top_score) = scores_sorted[0]
		
		# check if depth is adequate for confident call
		adequate_depth = False
		depth_problem = ""
		if hash_edge_depth[top_allele][0] > args.min_edge_depth and hash_edge_depth[top_allele][1] > args.min_edge_depth:
			if next_to_del_depth_allele[top_allele] != "NA":
				if float(next_to_del_depth_allele[top_allele]) > args.min_edge_depth:
					if avg_depth_allele[top_allele] > args.min_depth:
						adequate_depth = True
					else:
						depth_problem="depth"+str(avg_depth_allele[top_allele])
				else:
					depth_problem = "del"+str(next_to_del_depth_allele[top_allele])
			elif avg_depth_allele[top_allele] > args.min_depth:
				adequate_depth = True
			else:
				depth_problem="depth"+str(avg_depth_allele[top_allele])
		else:
			depth_problem = "edge"+str(min(hash_edge_depth[top_allele][0],hash_edge_depth[top_allele][1]))
			
		# check if there are confident differences against this allele
		differences = ""
		if mismatch_allele[top_allele] > 0:
			differences += str(mismatch_allele[top_allele])+"snp"
		if indel_allele[top_allele] > 0:
			differences += str(indel_allele[top_allele])+"indel"
		if missing_allele[top_allele] > 0:
			differences += str(missing_allele[top_allele])+"holes"
			
		if differences != "" or not adequate_depth:
			# if no SNPs or not enough depth to trust the result, no need to screen next best match
			results[gene] = (top_allele, differences, depth_problem)
		else:
			# looks good but this could be a truncated version of the real allele; check for longer versions
			truncation_override = False
			if len(scores_sorted) > 1:
				(next_best_allele,next_best_score) = scores_sorted[1]
				if size_allele[next_best_allele] > size_allele[top_allele]:
					# next best is longer, top allele could be a truncation?
					if (mismatch_allele[next_best_allele] + indel_allele[next_best_allele] + missing_allele[next_best_allele]) == 0:
						# next best also has no mismatches
						if (next_best_score - top_score)/top_score < 0.1:
							# next best has score within 10% of this one
							truncation_override = True
			if truncation_override:
				results[gene] = (next_best_allele, "trun", "") # no diffs but report this call is based on truncation test
			else:
				results[gene] = (top_allele, "", "") # no caveats to report

	return results # (allele, diffs, depth_problem)

def calculate_ST(allele_scores, ST_db, gene_names, sample_name, mlst_delimiter, avg_depth_allele):
	allele_numbers = [] # clean allele calls for determing ST. order is taken from gene names, as in ST definitions file
	alleles_with_flags = [] # flagged alleles for printing (* if mismatches, ? if depth issues)
	mismatch_flags = [] # allele/diffs
	uncertainty_flags = [] #allele/uncertainty
	depths = [] # depths for each typed loci
	
	# get allele numbers & info
	for gene in gene_names:
		if gene in allele_scores:
			(allele,diffs,depth_problem) = allele_scores[gene]
			allele_number = allele.split(mlst_delimiter)[-1]
			depths.append(avg_depth_allele[allele])
		else:
			allele_number = "-"
			diffs = ""
			depth_problem = ""
		allele_numbers.append(allele_number)

		allele_with_flags = allele_number
		if diffs != "":
			if diffs != "trun":
				allele_with_flags+="*" # trun indicates only that a truncated form had lower score, which isn't a mismatch
			mismatch_flags.append(allele+"/"+diffs)
		if depth_problem != "":
			allele_with_flags+="?"
			uncertainty_flags.append(allele+"/"+depth_problem)
		alleles_with_flags.append(allele_with_flags)
			
	# calculate ST (no flags)
	if ST_db:
		allele_string = " ".join(allele_numbers) # for determining ST
		try:
			clean_st = ST_db[allele_string]
		except KeyError:
			print "This combination of alleles was not found in the sequence type database:"
			print sample_name,
			for gene in allele_scores:
				(allele,diffs,depth_problems) = allele_scores[gene]
				print allele,
			print
			clean_st = "NF"
	else:
		clean_st = "ND"
		
	# add flags for reporting
	st = clean_st
	if len(mismatch_flags) > 0:
		if mismatch_flags!=["trun"]:
			st += "*" # trun indicates only that a truncated form had lower score, which isn't a mismatch
	else:
		mismatch_flags = ['0'] # record no mismatches
	if len(uncertainty_flags) > 0:
		st += "?"
	else:
		uncertainty_flags = ['-']
	
	# mean depth across loci
	if len(depths) > 0:
		mean_depth = float(sum(depths))/len(depths)
	else:
		mean_depth = 0
	
	return st, clean_st, alleles_with_flags, mismatch_flags, uncertainty_flags, mean_depth

def main():
	args = parse_args()
			
	known_alleles_provided = False
	
	if args.known_alleles != None:
		# Full table of known alleles provided
		known_alleles_db, gene_names = parse_known_alleles(args.known_alleles) 
		known_alleles_provided = True
		
	if args.known_STs != None:
		if args.ST_definitions != None:
			known_alleles_provided = True
			# Get databases of sequence types
			ST_db, ST_dict, gene_names = parse_ST_database(args.ST_definitions,args.ignore_last)
			# Read known STs and determine alleles
			known_alleles_db = parse_known_STs(args.known_STs, ST_db, gene_names) # full table of known alleles provided
		else:
			print "You have provided known STs but not the ST definitions. Please provide either full alleles (--known_alleles) or the ST definitions (--st_definitions)."

	# prepare output files
	score_outfile = file(args.output + ".score_summary", "w")
	st_outfile = file(args.output + ".st_summary", "w")
	score_outfile.write("Sample\tCluster\tAllele\tScore\tAvg_depth\tEdge1_depth\tEdge2_depth\tPercent_coverage\tSize\tMismatches\tIndels\tTruncated_bases\tDepthNeighbouringTruncation\tMixRate\tLeastConfident_Rate\tLeastConfident_Mismatches\tLeastConfident_Depth\tLeastConfident_Pvalue")
	st_outfile.write("\t".join(["Sample","ST"]+gene_names+["mismatches","uncertainty","depth"]))
	if known_alleles_provided:
		score_outfile.write("\tExpected\tMatchExpected\n")
		st_outfile.write("\tExpected\tMatchExpected\n")
	else:
		score_outfile.write("\n")
		st_outfile.write("\n")

	# get and print info on top scores for each sample
	for scores_file in args.scores:
	
		# read scores
		hash_edge_depth, avg_depth_allele, coverage_allele, \
			mismatch_allele, indel_allele, missing_allele, size_allele, \
			next_to_del_depth_allele, scores, mix_rates, allele_info_strings = \
				read_scores_file(scores_file)
	
		# get top scoring gene for each locus; key = locus, value = (allele, diffs, depth_problem)
		allele_scores = parse_scores(args, scores, hash_edge_depth, \
			avg_depth_allele, coverage_allele, mismatch_allele, indel_allele, \
			missing_allele, size_allele, next_to_del_depth_allele)
							
		(scores_file_path,scores_file_name) = os.path.split(scores_file)
		sample_name = scores_file_name.split(".")[0].split("__")[1]
		
		# process score info for each allele call
		for locus in gene_names:
			if locus in allele_scores:
				(allele, diffs, depth_problem) = allele_scores[locus]
				allele_number = allele.split(args.mlst_delimiter)[1]
				outstring = "\t".join([sample_name,locus] + allele_info_strings[allele])
			else:
				allele_number = "?"
				outstring = "\t".join([sample_name,locus] + ([''] * 16))
			score_outfile.write(outstring)
			if known_alleles_provided:
				if sample_name in known_alleles_db:
					expected = known_alleles_db[sample_name][locus]
					match_expected = expected == allele_number
					score_outfile.write("\t" + expected + "\t" + str(match_expected) + "\n")
				else:
					score_outfile.write("\tUnknown\tFalse\n")
			else:
				score_outfile.write("\n")
		
		# process ST information
		st,clean_st,alleles_with_flags,mismatch_flags,uncertainty_flags,mean_depth = \
			calculate_ST(allele_scores, ST_dict, gene_names, sample_name, args.mlst_delimiter, avg_depth_allele)
		st_result_string = "\t".join([sample_name,st]+alleles_with_flags+[";".join(mismatch_flags),";".join(uncertainty_flags),str(mean_depth)])
		st_outfile.write(st_result_string)
		if known_alleles_provided:
			if sample_name in known_alleles_db:
				expected = known_alleles_db[sample_name]["ST"]
				match_expected = expected == clean_st
				st_outfile.write("\t" + expected + "\t" + str(match_expected) + "\n")
			else:
				st_outfile.write("\tUnknown\tFalse\n")
		else:
			st_outfile.write("\n")

	score_outfile.close()
	st_outfile.close()

if __name__ == '__main__':
	main()

