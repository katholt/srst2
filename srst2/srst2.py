#!/usr/bin/env python

# SRST2 - Short Read Sequence Typer (v2)
# Python Version 2.7.5
#
# Authors - Michael Inouye (minouye@unimelb.edu.au), Harriet Dashnow (h.dashnow@gmail.com), 
#	Kathryn Holt (kholt@unimelb.edu.au), Bernie Pope (bjpope@unimelb.edu.au)
#
# see LICENSE.txt for the license
#
# Dependencies:
#	bowtie2	   http://bowtie-bio.sourceforge.net/bowtie2/index.shtml version 2.1.0
#	SAMtools	  http://samtools.sourceforge.net Version: 0.1.18 (Version: 0.1.19 DOES NOT WORK - loss of edge coverage)

from argparse import (ArgumentParser, FileType)
import logging
from subprocess import call, check_output, CalledProcessError
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
		
	# Gene database parameters
	parser.add_argument('--gene_db', type=str, required=False, nargs='+', help='Fasta file/s for gene databases (optional)')
	parser.add_argument('--no_gene_details', action="store_false", required=False, help='Switch OFF verbose reporting of gene typing')

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
	parser.add_argument('--save_scores', action="store_true", required=False, help='Switch ON verbose reporting of all scores')

	# Run options
	parser.add_argument('--use_existing_pileup', action="store_true", required=False,
		help='Use existing pileups if available, otherwise they will be generated') # to facilitate testing of rescoring from pileups
	parser.add_argument('--use_existing_scores', action="store_true", required=False,
		help='Use existing scores files if available, otherwise they will be generated') # to facilitate testing of reporting from scores

	# Compile previous output files
	parser.add_argument('--prev_output', nargs='+', type=str, required=False,
		help='SRST2 output files to compile (any new results from this run will also be incorporated)')

	return parser.parse_args() 


# Exception to raise if the command we try to run fails for some reason
class CommandError(Exception):
	pass

def run_command(command, **kwargs):
	'Execute a shell command and check the exit status and any O/S exceptions'
	command_str = ' '.join(command)
	logging.info('Running: {}'.format(command_str))
	try:
		exit_status = call(command, **kwargs)
	except OSError as e:
		message = "Command '{}' failed due to O/S error: {}".format(command_str, str(e))
		raise CommandError({"message": message})
	if exit_status != 0:
		message = "Command '{}' failed with non-zero exit status: {}".format(command_str, exit_status)
		raise CommandError({"message": message})


def bowtie_index(fasta_files):
	'Build a bowtie2 index from the given input fasta(s)'

	# check that both bowtie and samtools have the right versions
	check_command_version(['bowtie2', '--version'],
				'bowtie2-align version 2.1.0',
				'bowtie',
				'2.1.0')

	for fasta in fasta_files:
		built_index = fasta + '.1.bt2'
		if os.path.exists(built_index):
			logging.info('Index for {} is already built...'.format(fasta))
		else:
			logging.info('Building bowtie2 index for {}...'.format(fasta))
			run_command(['bowtie2-build', fasta, fasta])


def modify_bowtie_sam(raw_bowtie_sam):
	# fix sam flags for comprehensive pileup
	with open(raw_bowtie_sam) as sam, open(raw_bowtie_sam + '.mod', 'w') as sam_mod:
		for line in sam:
			if not line.startswith('@'):
				fields = line.split('\t')
				flag = int(fields[1])
				flag = (flag - 256) if (flag & 256) else flag
				sam_mod.write('\t'.join([fields[0], str(flag)] + fields[2:]))
			else:
				sam_mod.write(line)
	return(raw_bowtie_sam,raw_bowtie_sam + '.mod')


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
				if len(allele_info) > 2:
					gene_cluster = allele_info[0] # ID number for the cluster
					cluster_symbol = allele_info[1] # gene name for the cluster
					name = allele_info[2] # specific allele name
					if gene_cluster in gene_cluster_symbols:
						if gene_cluster_symbols[gene_cluster] != cluster_symbol:
							unique_gene_symbols = False # already seen this cluster name
					else:
						gene_cluster_symbols[gene_cluster] = cluster_symbol
				else:
					# treat as unclustered database, use whole header
					gene_cluster = cluster_symbol = name
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
		print ",".join(delimiter_check)
			
	return size, gene_clusters, unique_gene_symbols, unique_allele_symbols


def read_pileup_data(pileup_file, size, prob_err):
	with open(pileup_file) as pileup:
		prob_success = 1 - prob_err	# Set by user, default is prob_err = 0.01
		hash_alignment = {}
		hash_max_depth = {}
		hash_edge_depth = {}
		max_depth = 1
		avg_depth_allele = {}
		next_to_del_depth_allele = {}
		coverage_allele = {}
		mismatch_allele = {}
		indel_allele = {}
		missing_allele = {}
		size_allele = {}

		# Split all lines in the pileup by whitespace
		pileup_split = ( x.split() for x in pileup )
		# Group the split lines based on the first field (allele) 
		for allele, lines in groupby(pileup_split, itemgetter(0)):

			# Reset variables for new allele
			allele_line = 1 # Keep track of line for this allele
			exp_nuc_num = 0 # Expected position in ref allele
			allele_size = size[allele]
			total_depth = 0
			depth_a = depth_z = 0
			position_depths = [0] * allele_size # store depths in case required for penalties; then we don't need to track total_missing_bases
			hash_alignment[allele] = []
			total_missing_bases = 0
			total_mismatch = 0
			ins_poscount = 0
			del_poscount = 0
			next_to_del_depth = 99999

			for fields in lines:
				# Parse this line and store details required for scoring
				nuc_num = int(fields[1]) # Actual position in ref allele
				exp_nuc_num += 1
				allele_line += 1
				nuc = fields[2]
				nuc_depth = int(fields[3])
				position_depths[nuc_num-1] = nuc_depth
				if len(fields) <= 5:
					aligned_bases = ''
				else:
					aligned_bases = fields[4]

				# Missing bases (pileup skips basepairs)
				if nuc_num > exp_nuc_num:
					total_missing_bases += abs(exp_nuc_num - nuc_num)
				exp_nuc_num = nuc_num
				if nuc_depth == 0:
					total_missing_bases += 1

				# Calculate depths for this position
				if nuc_num <= edge_a:
					depth_a += nuc_depth
				if abs(nuc_num - allele_size) < edge_z:
					depth_z += nuc_depth
				if nuc_depth > max_depth:
					hash_max_depth[allele] = nuc_depth
					max_depth = nuc_depth

				total_depth = total_depth + nuc_depth

				# Parse aligned bases list for this position in the pileup
				num_match = 0
				ins_readcount = 0
				del_readcount = 0

				i = 0
				while i < len(aligned_bases):
					
					if aligned_bases[i] == "^":
						# Signifies start of a read, next char is mapping quality (skip it)
						i += 2
						continue

					if aligned_bases[i] == "+":
						i += int(aligned_bases[i+1]) + 2 # skip to next read
						ins_readcount += 1
						continue
						
					if aligned_bases[i] == "-":
						i += int(aligned_bases[i+1]) + 2 # skip to next read
						continue
						
					if aligned_bases[i] == "*":
						i += 1 # skip to next read
						del_readcount += 1
						continue

					if aligned_bases[i] == "." or aligned_bases[i] == ",":
						num_match += 1
						i += 1
						continue
					
					i += 1

				# Calculate details of this position for scoring and reporting
				
				# mismatches and indels
				num_mismatch = nuc_depth - num_match
				if num_mismatch > num_match:
					total_mismatch += 1 # record as mismatch (could be a snp or deletion)
				if del_readcount > num_match:
					del_poscount += 1
				if ins_readcount > nuc_depth / 2:
					ins_poscount += 1

				# Hash for later processing
				hash_alignment[allele].append((num_match, num_mismatch, prob_success)) # snp or deletion
				if ins_readcount > 0:
					hash_alignment[allele].append((nuc_depth - ins_readcount, ins_readcount, prob_success)) # penalize for any insertion calls at this position
				

			# Finished reading pileup for this allele
			
			# Check for missing bases at the end of the allele
			if nuc_num < allele_size:
				total_missing_bases += abs(allele_size - nuc_num)
				# determine penalty based on coverage of last 2 bases
				penalty = float(position_depths[nuc_num-1] + position_depths[nuc_num-2])/2
				m = min(position_depths[nuc_num-1],position_depths[nuc_num-2])
				hash_alignment[allele].append((0, penalty, prob_success))
				if next_to_del_depth > m:
					next_to_del_depth = m # keep track of lowest near-del depth for reporting

			# Calculate allele summary stats and save
			avg_depth = total_depth / float(allele_line)
			avg_a = depth_a / float(edge_a)   # Avg depth at 5' end, num basepairs determined by edge_a
			avg_z = depth_z / float(edge_z)	# 3'
			hash_max_depth[allele] = max_depth
			hash_edge_depth[allele] = (avg_a, avg_z)
			min_penalty = max(5, int(avg_depth))
			coverage_allele[allele] = 100*(allele_size - total_missing_bases - del_poscount)/float(allele_size) # includes in-read deletions
			mismatch_allele[allele] = total_mismatch - del_poscount # snps only
			indel_allele[allele] = del_poscount + ins_poscount # insertions or deletions
			missing_allele[allele] = total_missing_bases # truncated bases
			size_allele[allele] = allele_size
			
			# Penalize truncations or large deletions (i.e. positions not covered in pileup)
			j = 0
			while j < (len(position_depths)-2):
				# note end-of-seq truncations are dealt with above)
				if position_depths[j]==0 and position_depths[j+1]!=0:
					penalty = float(position_depths[j+1]+position_depths[j+2])/2 # mean of next 2 bases
					hash_alignment[allele].append((0, penalty, prob_success))
					m = min(position_depths[nuc_num-1],position_depths[nuc_num-2])
					if next_to_del_depth > m:
						next_to_del_depth = m # keep track of lowest near-del depth for reporting
				j += 1

			# Store depth info for reporting
			avg_depth_allele[allele] = avg_depth
			if next_to_del_depth == 99999:
				next_to_del_depth = "NA"
			next_to_del_depth_allele[allele] = next_to_del_depth

	return hash_alignment, hash_max_depth, hash_edge_depth, avg_depth_allele, coverage_allele, mismatch_allele, indel_allele, missing_allele, size_allele, next_to_del_depth_allele


def score_alleles(args,out_file_sam3, hash_alignment, hash_max_depth, hash_edge_depth, 
		avg_depth_allele, coverage_allele, mismatch_allele, indel_allele, missing_allele, 
		size_allele, next_to_del_depth_allele, run_type):
	
	if args.save_scores:
		scores_output = file(out_file_sam3 + '.table.scores', 'w')
		scores_output.write("Allele\tScore\tAvg_depth\tEdge1_depth\tEdge2_depth\tPercent_coverage\tSize\tMismatches\tIndels\tTruncated_bases\tDepthNeighbouringTruncation\tLeastConfident_Rate\tLeastConfident_Mismatches\tLeastConfident_Depth\tLeastConfident_Pvalue\n")
	
	scores = {} # key = allele, value = score
	
	for allele in hash_alignment:
		if (run_type == "mlst") or (coverage_allele[allele] > args.min_coverage):
			pvals = []
			min_pval = 1.0
			min_pval_data = (999,999) # (mismatch, depth) for position with lowest p-value
			for nuc_info in hash_alignment[allele]:
				if nuc_info is not None:
					match, mismatch, prob_success = nuc_info
					if match > 0 or mismatch > 0:
						if mismatch == 0:
							p_value = 1.0
						else:
							p_value = binom_test([match, mismatch], None, prob_success)
						# Weight pvalue by (depth/max_depth)
						max_depth = hash_max_depth[allele]
						weight = (match + mismatch) / float(max_depth)
						p_value *= weight
						if p_value == 0:
							p_value = 0.000000000000000000000000000001 ### Was getting a value error when p_value = 0.0
						if p_value < min_pval:
							min_pval = p_value
							min_pval_data = (mismatch,match + mismatch)
						p_value = -log(p_value, 10)
						pvals.append(p_value)
			# Fit linear model to observed Pval distribution vs expected Pval distribution (QQ plot)
			pvals.sort(reverse=True)
			len_obs_pvals = len(pvals)
			exp_pvals = range(1, len_obs_pvals + 1)
			exp_pvals2 = [-log(float(ep) / (len_obs_pvals + 1), 10) for ep in exp_pvals]
		
			# Slope is score
			slope, _intercept, _r_value, _p_value, _std_err = linregress(exp_pvals2, pvals)

			# Store all scores for later processing
			scores[allele] = slope
		
			# print scores for each allele, if requested
			if args.save_scores:
				if allele in hash_edge_depth:
					start_depth, end_depth = hash_edge_depth[allele]
					edge_depth_str = str(start_depth) + '\t' + str(end_depth)
				else:
					edge_depth_str = "NA\tNA"
				this_depth = avg_depth_allele.get(allele, "NA")
				this_coverage = coverage_allele.get(allele, "NA")
				this_mismatch = mismatch_allele.get(allele, "NA")
				this_indel = indel_allele.get(allele, "NA")
				this_missing = missing_allele.get(allele, "NA")
				this_size = size_allele.get(allele, "NA")
				this_next_to_del_depth = next_to_del_depth_allele.get(allele, "NA")
				scores_output.write('\t'.join([allele, str(slope), str(this_depth), edge_depth_str, 
						str(this_coverage), str(this_size), str(this_mismatch), str(this_indel), str(this_missing), str(this_next_to_del_depth), str(float(min_pval_data[0])/min_pval_data[1]),str(min_pval_data[0]),str(min_pval_data[1]),str(min_pval)]) + '\n')

	if args.save_scores:
		scores_output.close()
		
	return(scores)

# Check that an acceptable version of a command is installed
# Exits the program if it can't be found.
# - command_list is the command to run to determine the version.
# - version_identifier is the unique string we look for in the stdout of the program.
# - command_name is the name of the command to show in error messages.
# - required_version is the version number to show in error messages.
def check_command_version(command_list, version_identifier, command_name, required_version):
	try:
		command_stdout = check_output(command_list)
	except OSError as e:
		logging.error("Failed command: {}".format(' '.join(command_list)))
		logging.error(str(e))
		logging.error("Could not determine the version of {}.".format(command_name))
		logging.error("Do you have {} installed in your PATH?".format(command_name))
		exit(-1)
	except CalledProcessError as e:
		# some programs such as samtools return a non-zero exit status
		# when you ask for the version (sigh). We ignore it here.
		pass

	if version_identifier not in command_stdout:
		logging.error("Incorrect version of {} installed.".format(command_name))
		logging.error("{} version {} is required by SRST2.".format(command_name, required_version))
		exit(-1)


def run_bowtie(sample_name,fastqs,args,db_name,db_full_path):

	# check that both bowtie and samtools have the right versions
	check_command_version(['bowtie2', '--version'],
				'bowtie2-align version 2.1.0',
				'bowtie',
				'2.1.0')

	check_command_version(['samtools'],
				'Version: 0.1.8',
				'samtools',
				'0.1.8')

	command = ['bowtie2']

	if len(fastqs)==1:
		# single end
		command += ['-U', fastqs[0]]
	elif len(fastqs)==2:
		# paired end
		command += ['-1', fastqs[0], '-2', fastqs[1]]
		
	out_file = sample_name + '.' + db_name + '.srst2'
	logging.info('Output prefix set to: ' + out_file)

	command += ['-S', out_file,
				'-' + args.read_type,	# add a dash to the front of the option
				'--very-sensitive-local',
				'--no-unal',
				'-a',					 # Search for and report all alignments
				'-x', db_full_path			   # The index to be aligned to
			   ]

	if args.other:
		command += args.other.split()

	logging.info('Aligning reads to index {} using bowtie2...'.format(db_full_path))
	
	run_command(command)
	
	return(out_file)
	
def get_pileup(args,raw_bowtie_sam,bowtie_sam_mod,fasta,pileup_file):
	# Analyse output with SAMtools
	logging.info('Processing Bowtie2 output with SAMtools...')
	logging.info('Generate and sort BAM file...')
	out_file_bam = raw_bowtie_sam + ".bam"
	run_command(['samtools', 'view', '-b', '-o', out_file_bam,
				 '-q', str(args.mapq), '-S', bowtie_sam_mod])
	out_file_bam_sorted = raw_bowtie_sam + ".sorted"
	run_command(['samtools', 'sort', out_file_bam, out_file_bam_sorted])

	#XXX File deletions should be made optional. May also want to delete final sorted bam and pileup.
	logging.info('Deleting sam and bam files that are not longer needed...')
	del_filenames = [raw_bowtie_sam, bowtie_sam_mod, out_file_bam]
	for f in del_filenames:
		logging.info('Deleting ' + f)
		os.remove(f)
		
	logging.info('Generate pileup...')
	with open(pileup_file, 'w') as sam_pileup:
		run_command(['samtools', 'mpileup', '-L', '1000', '-f', fasta,
					 '-Q', str(args.baseq), '-q', str(args.mapq), out_file_bam_sorted + '.bam'],
					 stdout=sam_pileup)

def calculate_ST(allele_scores, ST_db, gene_names, sample_name, mlst_delimiter, avg_depth_allele):
	allele_numbers = [] # clean allele calls for determing ST. order is taken from gene names, as in ST definitions file
	alleles_with_flags = [] # flagged alleles for printing (* if mismatches, ? if depth issues)
	mismatch_flags = [] # allele/diffs
	uncertainty_flags = [] #allele/uncertainty
#	st_flags = [] # (* if mismatches, ? if depth issues)
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
			print "This combination of alleles was not found in the sequence type database:",
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
	if len(uncertainty_flags) > 0:
		st += "?"
	
	# mean depth across loci	
	if len(depths) > 0:
		mean_depth = float(sum(depths))/len(depths)
	else:
		mean_depth = 0
	
	return (st,clean_st,alleles_with_flags,mismatch_flags,uncertainty_flags,mean_depth)
	
def parse_ST_database(ST_filename,gene_names_from_fai):
	# Read ST definitions
	ST_db = {} # key = allele string, value = ST
	gene_names = []
	num_gene_cols_expected = len(gene_names_from_fai)
	print "Attempting to read " + str(num_gene_cols_expected) + " loci from ST database " + ST_filename
	with open(ST_filename) as f:
		count = 0
		for line in f:
			count += 1
			line_split = line.rstrip().split("\t")
			if count == 1: # Header
				gene_names = line_split[1:min(num_gene_cols_expected+1,len(line_split))]
				for g in gene_names_from_fai:
					if g not in gene_names:
						print "Warning: gene " + g + " in database file isn't among the columns in the ST definitions: " + ",".join(gene_names)
						print " Any sequences with this gene identifer from the database will not be included in typing."
						if len(line_split) == num_gene_cols_expected+1:
							gene_names.pop() # we read too many columns
							num_gene_cols_expected -= 1
				for g in gene_names:
					if g not in gene_names_from_fai:
						print "Warning: gene " + g + " in ST definitions file isn't among those in the database " + ",".join(gene_names_from_fai)
						print " This will result in all STs being called as unknown (but allele calls will be accurate for other loci)." 
			else:
				ST = line_split[0]
				if ST not in ST_db.values():
					ST_string = " ".join(line_split[1:num_gene_cols_expected+1])
					ST_db[ST_string] = ST
				else:
					print "Warning: this ST is not unique in the ST definitions file: " + ST
		print "Read ST database " + ST_filename + " successfully"
		return (ST_db, gene_names)

def get_allele_name_from_db(allele,unique_allele_symbols,unique_cluster_symbols,run_type,args):
	
	if run_type != "mlst":
		# header format: >[cluster]___[gene]___[allele]___[uniqueID] [info]
		# header format: >0__oqxB__oqxB__2 oqxB_1_EU370913; EU370913; quinolone
		allele_parts = allele.split()
		allele_info = allele_parts[0].split("__")
		if len(allele_info)>2:
			cluster_id = allele_info[0] # ID number for the cluster
			gene_name = allele_info[1] # gene name for the cluster
			allele_name = allele_info[2] # specific allele name
			seqid = allele_info[3] # unique identifier for this seq
		else:
			cluster_id = gene_name = allele_name = seqid = allele_parts[0]
		if len(allele_parts) > 1:
			annotation = " ".join(allele_parts[1:-1])
		else:
			annotation = ""
	
		if not unique_allele_symbols:	
			allele_name += "_" + seq_id
		
		if not unique_allele_symbols:	
			gene_name += "_" + cluster_id
			
	else:
		gene_name = allele.split(args.mlst_delimiter)
		allele_name = allele
		annotation = None
		seqid = None
		cluster_id = None
		
	return gene_name, allele_name, cluster_id, annotation, seqid

def parse_scores(run_type,args,scores, hash_edge_depth, 
					avg_depth_allele, coverage_allele, mismatch_allele, indel_allele,  
					missing_allele, size_allele, next_to_del_depth_allele,
					unique_cluster_symbols,unique_allele_symbols):
					
	# sort into hash for each gene locus
	scores_by_gene = collections.defaultdict(dict) # key1 = gene, key2 = allele, value = score
	
	if run_type=="mlst":	
		for allele in scores:
			allele_info = allele.split(args.mlst_delimiter)
			scores_by_gene[allele_info[0]][allele] = scores[allele]
	else:
		for allele in scores:
			if coverage_allele[allele] > args.min_coverage:
				gene_name = get_allele_name_from_db(allele,unique_allele_symbols,unique_cluster_symbols,run_type,args)[0]
				scores_by_gene[gene_name][allele] = scores[allele]
	
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
							# next best has score within 5% of this one
							truncation_override = True
			if truncation_override:
				results[gene] = (next_best_allele, "trun", "") # no diffs but report this call is based on truncation test
			else:
				results[gene] = (top_allele, "", "") # no caveats to report

	return results # (allele, diffs, depth_problem)
					

def get_readFile_components(full_file_path):
	(file_path,file_name) = os.path.split(full_file_path)
	m1 = re.match("(.*).gz",file_name)
	ext = ""
	if m1 != None:
		# gzipped
		ext = ".gz"
		file_name = m1.groups()[0]
	(file_name_before_ext,ext2) = os.path.splitext(file_name)
	full_ext = ext2+ext
	return(file_path,file_name_before_ext,full_ext)

def read_file_sets(args):	

	fileSets = {} # key = id, value = list of files for that sample
	num_single_readsets = 0
	num_paired_readsets = 0
	
	if args.input_se:
		# single end
		for fastq in args.input_se:
			(file_path,file_name_before_ext,full_ext) = get_readFile_components(fastq)
			fileSets[file_name_before_ext] = [fastq]
			num_single_readsets += 1
			
	elif args.input_pe:
		# paired end
		forward_reads = {} # key = sample, value = full path to file
		reverse_reads = {} # key = sample, value = full path to file
		num_paired_readsets = 0
		num_single_readsets = 0
		for fastq in args.input_pe:
			(file_path,file_name_before_ext,full_ext) = get_readFile_components(fastq)
			m=re.match("(.*)("+args.forward+")",file_name_before_ext)
			if m!=None:
				# store as forward read
				(baseName,read) = m.groups()
				forward_reads[baseName] = fastq
			else:
				m=re.match("(.*)("+args.reverse+")",file_name_before_ext)
				if m!=None:
				# store as reverse read
					(baseName,read) = m.groups()
					reverse_reads[baseName] = fastq
				else:
					print "Could not determine forward/reverse read status for input file " + fastq
		# store in pairs
		for sample in forward_reads:
			if sample in reverse_reads:
				fileSets[sample] = [forward_reads[sample],reverse_reads[sample]] # store pair
				num_paired_readsets += 1
			else:
				fileSets[sample] = [forward_reads[sample]] # no reverse found
				num_single_readsets += 1
				logging.info('Warning, could not find pair for read:' + forward_reads[sample])
		for sample in reverse_reads:
			if sample not in fileSets:
				fileSets[sample] = reverse_reads[sample] # no forward found
				num_single_readsets += 1
				logging.info('Warning, could not find pair for read:' + reverse_reads[sample])
				
	if num_paired_readsets > 0:
		logging.info('Total paired readsets found:' + str(num_paired_readsets))	
	if num_single_readsets > 0:
		logging.info('Total single reads found:' + str(num_single_readsets))

	return fileSets 

def read_results_from_file(infile):
	
	results_info = infile.split("__")
	if len(results_info) > 1:
	
		if results_info[-1] == "compiledResults.txt":
			dbtype = "compiled"
			dbname = results_info[0] # output identifier
		else:
			dbtype = results_info[0] # mlst or genes
			dbname = results_info[1]

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
								results[sample][gene] = line_split[i]
							
	else:
		results = False
		dbtype = False
		dbname = False
		logging.info("Couldn't decide what to do with file results file provided: " + infile)
				
	return results, dbtype, dbname
						
def read_scores_file(scores_file):
	hash_edge_depth = {}
	avg_depth_allele = {}
	coverage_allele = {}
	mismatch_allele = {}
	indel_allele = {}
	missing_allele = {}
	size_allele = {}
	next_to_del_depth_allele = {}
	scores = {}
				
	f = file(scores_file,"r")
	
	for line in f:
		line_split = line.rstrip().split("\t")
		allele = line_split[0]
		if allele != "Allele": # skip header row
			scores[allele] = float(line_split[1])
			avg_depth_allele[allele] = float(line_split[2])
			hash_edge_depth[allele] = (float(line_split[3]),float(line_split[4]))
			coverage_allele[allele] = float(line_split[5])
			size_allele[allele] = int(line_split[6])
			mismatch_allele[allele] = int(line_split[7])
			indel_allele[allele] = int(line_split[8])
			missing_allele[allele] = int(line_split[9])
			next_to_del_depth = line_split[10]
			next_to_del_depth_allele[allele] = line_split[10]

	return hash_edge_depth, avg_depth_allele, coverage_allele, mismatch_allele, indel_allele, \
			missing_allele, size_allele, next_to_del_depth_allele, scores

def run_srst2(args, fileSets, dbs, run_type):

	db_reports = [] # list of db-specific output files to return
	db_results_list = [] # list of results hashes, one per db

	for fasta in dbs:
		db_reports, db_results_list = process_fasta_db(args, fileSets, run_type, db_reports, db_results_list, fasta)

	return db_reports, db_results_list

def process_fasta_db(args, fileSets, run_type, db_reports, db_results_list, fasta):

	check_command_version(['samtools'],
				'Version: 0.1.8',
				'samtools',
				'0.1.8')

	logging.info('Processing database ' + fasta)

	db_path, db_name = os.path.split(fasta) # database
	(db_name,db_ext) = os.path.splitext(db_name)
	db_results = "__".join([run_type,db_name,args.output,"results.txt"])
	db_report = file(db_results,"w")
	db_reports.append(db_results)
	
	# Get sequence lengths and gene names
	#  lengths are needed for MLST heuristic to distinguish alleles from their truncated forms
	#  gene names read from here are needed for non-MLST dbs
	fai_file = fasta + '.fai'
	if not os.path.exists(fai_file):
		run_command(['samtools', 'faidx', fasta])
	size, gene_names, unique_gene_symbols, unique_allele_symbols = \
		parse_fai(fai_file,run_type,args.mlst_delimiter)

	# Prepare for MLST reporting
	ST_db = False
	if run_type == "mlst":
		results = {} # key = sample, value = ST string for printing
		if args.mlst_definitions:
			# store MLST profiles, replace gene names (we want the order as they appear in this file)
			ST_db, gene_names = parse_ST_database(args.mlst_definitions,gene_names)
		db_report.write("\t".join(["Sample","ST"]+gene_names+["mismatches","uncertainty","depth"]) + "\n")
		results["Sample"] = "\t".join(["Sample","ST"]+gene_names+["mismatches","uncertainty","depth"])
		
	else:
		# store final results for later tabulation
		results = collections.defaultdict(dict) #key1 = sample, key2 = gene, value = allele

	gene_list = [] # start with empty gene list; will add genes from each genedb test
	
	# Align and score each read set against this DB
	for sample_name in fileSets:
		logging.info('Processing sample ' + sample_name)
		fastq_inputs = fileSets[sample_name] # reads
		
		try:
			# try mapping and scoring this fileset against the current database
			# update the gene_list list and results dict with data from this strain
			gene_list, results = \
				map_fileSet_to_db(args,sample_name,fastq_inputs,db_name,fasta,size,gene_names,\
				unique_gene_symbols, unique_allele_symbols,run_type,ST_db,results,gene_list,db_report)
		# if we get an error from one of the commands we called
		# log the error message and continue onto the next fasta db
            	except CommandError as e:
                	logging.error(e.message)
                	# record results as unknown, so we know that we did attempt to analyse this readset
                	if run_type == "mlst":
                		st_result_string = "\t".join( [sample_name,"-"] + ["-"] * (len(gene_names)+3)) # record missing results
                		db_report.write( st_result_string + "\n")
                		logging.info(" " + st_result_string)
                		results[sample_name] = st_result_string
                	else:
                		logging.info(" failed gene detection")
                		results[sample_name]["failed"] = True # so we know that we tried this strain
                	
	if run_type != "mlst":
		# tabulate results across samples for this gene db
		logging.info('Tabulating results for database {} ...'.format(fasta))
		gene_list.sort()
		db_report.write("\t".join(["Sample"]+gene_list)+"\n") # report header row
		for sample_name in fileSets:
			db_report.write(sample_name)
			if sample_name in results:
				# print results
				if "failed" not in results[sample_name]:
					for gene in gene_list:
						if gene in results[sample_name]:
							db_report.write("\t"+results[sample_name][gene])
						else:
							db_report.write("\t-") # no hits for this gene cluster
				else:
					# no data on this, as the sample failed mapping
					for gene in gene_list:
						db_report.write("\t?") # 
						results[sample_name][gene] = "?" # record as unknown
			else:
				# no data on this because genes were not found (but no mapping errors)
				for gene in gene_list:
					db_report.write("\t?") # 
					results[sample_name][gene] = "-" # record as absent
			db_report.write("\n")

	# Finished with this database
	logging.info('Finished processing for database {} ...'.format(fasta))
	db_report.close()
	db_results_list.append(results)
	
	return db_reports, db_results_list
						
def map_fileSet_to_db(args,sample_name,fastq_inputs,db_name,fasta,size,gene_names,\
	unique_gene_symbols, unique_allele_symbols,run_type,ST_db,results,gene_list,db_report):
	
	pileup_file = sample_name + '.' + db_name + '.srst2.pileup'
	scores_file = pileup_file + '.table.scores'
	
	# Get or read scores
	
	if args.use_existing_scores and os.path.exists(scores_file):
		
		logging.info(' Using existing scores in ' + scores_file)
			
		# read in scores and info from existing scores file
		hash_edge_depth, avg_depth_allele, coverage_allele, \
				mismatch_allele, indel_allele, missing_allele, size_allele, \
				next_to_del_depth_allele, scores = read_scores_file(scores_file)
	
	else:
	
		# Get or read pileup
		
		if args.use_existing_pileup and os.path.exists(pileup_file):
			logging.info(' Using existing pileup in ' + pileup_file)

		else:
			
			# run bowtie against this db
			bowtie_sam = run_bowtie(sample_name,fastq_inputs,args,db_name,fasta)

			# Modify Bowtie's SAM formatted output so that we get secondary
			# alignments in downstream pileup
			(raw_bowtie_sam,bowtie_sam_mod) = modify_bowtie_sam(bowtie_sam)
	
			# generate pileup from sam (via sorted bam)
			get_pileup(args,raw_bowtie_sam,bowtie_sam_mod,fasta,pileup_file)

		# Get scores

		# Process the pileup and extract info for scoring and reporting on each allele
		logging.info(' Processing SAMtools pileup...')
		hash_alignment, hash_max_depth, hash_edge_depth, avg_depth_allele, coverage_allele, \
				mismatch_allele, indel_allele, missing_allele, size_allele, next_to_del_depth_allele= \
				read_pileup_data(pileup_file, size, args.prob_err)

		# Generate scores for all alleles (prints these and associated info if verbose)
		#   result = dict, with key=allele, value=score
		logging.info(' Scoring alleles...')
		scores = score_alleles(args,pileup_file, hash_alignment, hash_max_depth, hash_edge_depth, \
				avg_depth_allele, coverage_allele, mismatch_allele, indel_allele, missing_allele, \
				size_allele, next_to_del_depth_allele, run_type)
		
	# Determine best score for each gene/cluster
	
	#  result = dict, with key = gene, value = (allele,diffs,depth_problem)
	#			for MLST DBs, key = gene = locus, allele = gene-number
	#           for gene DBs, key = gene = cluster ID, allele = cluster__gene__allele__id
	#  for gene DBs, only those alleles passing the coverage cutoff are returned
	allele_scores = parse_scores(run_type,args,scores, \
			hash_edge_depth, avg_depth_allele, coverage_allele, mismatch_allele,  \
			indel_allele, missing_allele, size_allele, next_to_del_depth_allele,
			unique_gene_symbols, unique_allele_symbols)
			
	# Report results
	if run_type == "mlst" and len(allele_scores) > 0:
					
		# Calculate ST and get info for reporting
		(st,clean_st,alleles_with_flags,mismatch_flags,uncertainty_flags,mean_depth) = \
				calculate_ST(allele_scores, ST_db, gene_names, sample_name, args.mlst_delimiter, avg_depth_allele)

		# Print to database-specific report, log and save the result
		st_result_string = "\t".join([sample_name,st]+alleles_with_flags+[";".join(mismatch_flags),";".join(uncertainty_flags),str(mean_depth)])
		db_report.write( st_result_string + "\n")
		logging.info(" " + st_result_string)
		results[sample_name] = st_result_string
		
	elif run_type == "genes" and len(allele_scores) > 0:
		if args.no_gene_details:
			full_results = "__".join(["full"+run_type,db_name,args.output,"results.txt"])
			logging.info("Printing verbose gene detection results to " + full_results)
			f = file(full_results,"w")
			f.write("\t".join(["Sample","DB","gene","allele","coverage","depth","diffs","uncertainty","cluster","seqid","annotation"])+"\n")
		for gene in allele_scores:
			(allele,diffs,depth_problem) = allele_scores[gene]
			gene_name, allele_name, cluster_id, annotation, seqid = \
				get_allele_name_from_db(allele,unique_allele_symbols,unique_gene_symbols,run_type,args)
			if gene not in gene_list:
				gene_list.append(gene_name)
			results[sample_name][gene_name] = allele_name
			if diffs != "":
				results[sample_name][gene_name] += "*"
			if depth_problem != "":
				results[sample_name][gene_name] += "?"
			if args.no_gene_details:
				f.write("\t".join([sample_name,db_name,gene_name,allele_name,str(coverage_allele[allele]),str(avg_depth_allele[allele]),diffs,depth_problem,cluster_id,seqid,annotation])+"\n")
	
		# log the gene detection result
		logging.info(" " + str(len(allele_scores)) + " genes identified in " + sample_name)
	
	# Finished with this read set
	logging.info(' Finished processing for read set {} ...'.format(sample_name))
	
	return gene_list, results
	
def compile_results(args,mlst_results,db_results,compiled_output_file):

	o = file(compiled_output_file,"w")
	
	# get list of all samples and genes present in these datasets
	sample_list = [] # each entry is a sample present in at least one db
	gene_list = []
	mlst_cols = 0
	mlst_header_string = ""
	blank_mlst_section = ""
	
	mlst_results_master = {} # compilation of all MLST results
	db_results_master = collections.defaultdict(dict) # compilation of all gene results
	
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
						
	if "Sample" in sample_list:
		sample_list.remove("Sample")
	sample_list.sort()
	gene_list.sort()
	
	# print header
	header_elements = []
	if len(mlst_results) > 0:
		header_elements.append(mlst_header_string)
	else:
		header_elements.append("Sample")
	if (gene_list) > 0:
		header_elements += gene_list
	o.write("\t".join(header_elements)+"\n")
	
	# print results for all samples
	for sample in sample_list:
	
		sample_info = [] # first entry is mlst string OR sample name, rest are genes
		
		# print mlst if provided, otherwise just print sample name
		if len(mlst_results_master) > 0:
			if sample in mlst_results_master:
				sample_info.append(mlst_results_master[sample])
			else:
				sample_info.append(sample+blank_mlst_section)
		else:
			sample_info.append(sample)
		
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
	
	logging.info("Compiled data printed to: " + compiled_output_file)
	
	return True
	

def main():
	args = parse_args()
	if args.log is None:
		logfile = sys.stdout
	else:
		logfile = args.log
	logging.basicConfig(
		filename=args.log,
		level=logging.DEBUG,
		filemode='w',
		format='%(asctime)s %(message)s',
		datefmt='%m/%d/%Y %H:%M:%S')
	logging.info('program started')
	logging.info('command line: {0}'.format(' '.join(sys.argv)))

	# vars to store results
	mlst_results_hashes = [] # dict (sample->MLST result string) for each MLST output files created/read
	gene_result_hashes = [] # dict (sample->gene->result) for each gene typing output files created/read

	# parse list of file sets to analyse
	fileSets = read_file_sets(args) # get list of files to process
	
	# run MLST scoring
	if fileSets and args.mlst_db:
	
		if not args.mlst_definitions:
		
			# print warning to screen to alert user, may want to stop and restart
			print "Warning, MLST allele sequences were provided without ST definitions:"
			print " allele sequences: " + str(args.mlst_db)
			print " these will be mapped and scored, but STs can not be calculated"
			
			# log
			logging.info("Warning, MLST allele sequences were provided without ST definitions:")
			logging.info(" allele sequences: " + str(args.mlst_db))
			logging.info(" these will be mapped and scored, but STs can not be calculated")
		
		bowtie_index(args.mlst_db) # index the MLST database
		
		# score file sets against MLST database
		mlst_report, mlst_results = run_srst2(args,fileSets,args.mlst_db,"mlst")
		
		logging.info('MLST output printed to ' + mlst_report[0])
		
		#mlst_reports_files += mlst_report
		mlst_results_hashes += mlst_results
		
	# run gene detection
	if fileSets and args.gene_db:

		bowtie_index(args.gene_db) # index the gene databases
		
		db_reports, db_results = run_srst2(args,fileSets,args.gene_db,"genes")

		for outfile in db_reports:
			logging.info('Gene detection output printed to ' + outfile)
			
		gene_result_hashes += db_results
			
	# process prior results files
	if args.prev_output:
	
		for results_file in args.prev_output:
		
			results, dbtype, dbname = read_results_from_file(results_file)
			
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

	# compile results if multiple databases or datasets provided
	if ( (len(gene_result_hashes) + len(mlst_results_hashes)) > 1 ):
		compiled_output_file = args.output + "__compiledResults.txt"
		compile_results(args,mlst_results_hashes,gene_result_hashes,compiled_output_file)
	
	logging.info('SRST2 has finished.')


if __name__ == '__main__':
	main()
