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
#	bowtie2	   http://bowtie-bio.sourceforge.net/bowtie2/index.shtml version 2.1.0 or greater
#	SAMtools   http://samtools.sourceforge.net Version: 0.1.18 or greater (note optimal results are obtained with 0.1.18 rather than later versions)
#	SciPy	http://www.scipy.org/install.html
#
# Git repository: https://github.com/katholt/srst2/
# README: https://github.com/katholt/srst2/blob/master/README.md
# Questions or feature requests: https://github.com/katholt/srst2/issues
# Paper: http://genomemedicine.com/content/6/11/90


from argparse import (ArgumentParser, FileType)
import logging
from subprocess import call, check_output, CalledProcessError, STDOUT
import os, sys, re, collections, operator
from scipy.stats import binom, linregress
from math import log
from itertools import groupby
from operator import itemgetter
from collections import OrderedDict
try:
	from version import srst2_version
except:
	srst2_version = "version unknown"

edge_a = edge_z = 2


def parse_args():
	"Parse the input arguments, use '-h' for help."

	parser = ArgumentParser(description='SRST2 - Short Read Sequence Typer (v2)')

	# version number of srst2, print and then exit
	parser.add_argument('--version', action='version', version='%(prog)s ' + srst2_version)

	# Read inputs
	parser.add_argument(
		'--input_se', nargs='+',type=str, required=False,
		help='Single end read file(s) for analysing (may be gzipped)')
	parser.add_argument(
		'--input_pe', nargs='+', type=str, required=False,
		help='Paired end read files for analysing (may be gzipped)')
	parser.add_argument('--merge_paired', action="store_true", required=False, help='Switch on if all the input read sets belong to a single sample, and you want to merge their data to get a single result')
	parser.add_argument(
		'--forward', type=str, required=False, default="_1",
			help='Designator for forward reads (only used if NOT in MiSeq format sample_S1_L001_R1_001.fastq.gz; otherwise default is _1, i.e. expect forward reads as sample_1.fastq.gz)')
	parser.add_argument(
		'--reverse', type=str, required=False, default="_2",
			help='Designator for reverse reads (only used if NOT in MiSeq format sample_S1_L001_R2_001.fastq.gz; otherwise default is _2, i.e. expect forward reads as sample_2.fastq.gz')
	parser.add_argument('--read_type', type=str, choices=['q', 'qseq', 'f'], default='q',
		help='Read file type (for bowtie2; default is q=fastq; other options: qseq=solexa, f=fasta).')

	# MLST parameters
	parser.add_argument('--mlst_db', type=str, required=False, nargs=1, help='Fasta file of MLST alleles (optional)')
	parser.add_argument('--mlst_delimiter', type=str, required=False,
		help='Character(s) separating gene name from allele number in MLST database (default "-", as in arcc-1)', default="-")
	parser.add_argument('--mlst_definitions', type=str, required=False,
		help='ST definitions for MLST scheme (required if mlst_db supplied and you want to calculate STs)')
	parser.add_argument('--mlst_max_mismatch', type=str, required=False, default = "10",
		help='Maximum number of mismatches per read for MLST allele calling (default 10)')

	# Gene database parameters
	parser.add_argument('--gene_db', type=str, required=False, nargs='+', help='Fasta file/s for gene databases (optional)')
	parser.add_argument('--no_gene_details', action="store_false", required=False, help='Switch OFF verbose reporting of gene typing')
	parser.add_argument('--gene_max_mismatch', type=str, required=False, default = "10",
		help='Maximum number of mismatches per read for gene detection and allele calling (default 10)')

	# Cutoffs for scoring/heuristics
	parser.add_argument('--min_coverage', type=float, required=False, help='Minimum %%coverage cutoff for gene reporting (default 90)',default=90)
	parser.add_argument('--max_divergence', type=float, required=False, help='Maximum %%divergence cutoff for gene reporting (default 10)',default=10)
	parser.add_argument('--min_depth', type=float, required=False, help='Minimum mean depth to flag as dubious allele call (default 5)',default=5)
	parser.add_argument('--min_edge_depth', type=float, required=False, help='Minimum edge depth to flag as dubious allele call (default 2)',default=2)
	parser.add_argument('--prob_err', type=float, help='Probability of sequencing error (default 0.01)',default=0.01)
	parser.add_argument('--truncation_score_tolerance', type=float, help='%% increase in score allowed to choose non-truncated allele',default=0.15)

	# Mapping parameters for bowtie2
	parser.add_argument('--stop_after', type=str, required=False, help='Stop mapping after this number of reads have been mapped (otherwise map all)')
	parser.add_argument('--other', type=str, help='Other arguments to pass to bowtie2 (must be escaped, e.g. "\--no-mixed".', required=False)

	# Filtering parameters for initial SAM file
	parser.add_argument('--max_unaligned_overlap', type=int, default=10, help='Read discarded from alignment if either of its ends has unaligned '+\
			'overlap with the reference that is longer than this value (default 10)')

	# Samtools parameters
	parser.add_argument('--mapq', type=int, default=1, help='Samtools -q parameter (default 1)')
	parser.add_argument('--baseq', type=int, default=20, help='Samtools -Q parameter (default 20)')
	parser.add_argument('--samtools_args', type=str, help='Other arguments to pass to samtools mpileup (must be escaped, e.g. "\-A").', required=False)

	# Reporting options
	parser.add_argument('--output', type=str, required=True, help='Prefix for srst2 output files')
	parser.add_argument('--log', action="store_true", required=False, help='Switch ON logging to file (otherwise log to stdout)')
	parser.add_argument('--save_scores', action="store_true", required=False, help='Switch ON verbose reporting of all scores')
	parser.add_argument('--report_new_consensus', action="store_true", required=False, help='If a matching alleles is not found, report the consensus allele. Note, only SNP differences are considered, not indels.')
	parser.add_argument('--report_all_consensus', action="store_true", required=False, help='Report the consensus allele for the most likely allele. Note, only SNP differences are considered, not indels.')

	# Run options
	parser.add_argument('--use_existing_bowtie2_sam', action="store_true", required=False,
		help='Use existing SAM file generated by Bowtie2 if available, otherwise they will be generated') # to facilitate testing of filtering Bowtie2 output
	parser.add_argument('--use_existing_pileup', action="store_true", required=False,
		help='Use existing pileups if available, otherwise they will be generated') # to facilitate testing of rescoring from pileups
	parser.add_argument('--use_existing_scores', action="store_true", required=False,
		help='Use existing scores files if available, otherwise they will be generated') # to facilitate testing of reporting from scores
	parser.add_argument('--keep_interim_alignment', action="store_true", required=False, default=False,
		help='Keep interim files (sam & unsorted bam), otherwise they will be deleted after sorted bam is created') # to facilitate testing of sam processing
	parser.add_argument('--threads', type=int, required=False, default=1,
		help='Use multiple threads in Bowtie and Samtools')
#	parser.add_argument('--keep_final_alignment', action="store_true", required=False, default=False,
#		help='Keep interim files (sam & unsorted bam), otherwise they will be deleted after sorted bam is created') # to facilitate testing of sam processing

	# Compile previous output files
	parser.add_argument('--prev_output', nargs='+', type=str, required=False,
		help='SRST2 results files to compile (any new results from this run will also be incorporated)')

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
	check_bowtie_version()
	for fasta in fasta_files:
		built_index = fasta + '.1.bt2'
		if os.path.exists(built_index):
			logging.info('Index for {} is already built...'.format(fasta))
		else:
			logging.info('Building bowtie2 index for {}...'.format(fasta))
			run_command([get_bowtie_execs()[1], fasta, fasta])

def get_clips_cigar(cigar):
	## remove padding first if present;
	## maybe padding is never present at the edges, but it is easier to just remove
	cigar = re.sub(r'\d+P','',cigar.strip())
	x = re.search(r'^(?P<length>\d+)(?P<type>[SH])',cigar)
	if x:
		left_clip = x.groupdict()
		left_clip["length"] = int(left_clip["length"])
	else:
		left_clip = dict(length=0,type=None)

	x = re.search(r'(?P<type>[SH])(?P<length>\d+)$',cigar)
	if x:
		right_clip = x.groupdict()
		right_clip["length"] = int(right_clip["length"])
	else:
		right_clip = dict(length=0,type=None)

	return (left_clip,right_clip)

def get_end_shift_cigar(cigar):
	"""Return change in coordinate on the reference of the read end due to indels in CIGAR string"""
	shift = 0
	for edit_op in re.findall(r'(\d+)([ID])',cigar):
		shift += int(edit_op[0])*(-1 if edit_op[1] == 'I' else 1)
	return shift

def get_unaligned_read_end_lengths_sam(fields,ref_len):
	"""From SAM file line, compute clipped read length within reference"""
	left_res = 0
	right_res = 0
	if len(fields) >= 10:
		## get (clipped) start position
		ali_clipped_start = int(fields[3])
		cigar = fields[5]
		## get number and types of clipped bases on the left and right
		left_clip, right_clip = get_clips_cigar(cigar)
		left_res = min(ali_clipped_start,left_clip["length"])
		seq_start = ali_clipped_start
		if left_clip["type"] and left_clip["type"] == "S":
			seq_start -= left_clip["length"]
		## get (hard-clipped) end position as start + len(seq)
		seq_hard_clipped_end = seq_start + len(fields[9]) + get_end_shift_cigar(cigar)
		## seq end = hard end + right hard clip
		seq_end = seq_hard_clipped_end
		if right_clip["type"] and right_clip["type"] == "H":
			seq_end += right_clip["length"]
		## aligned end = hard end - right soft clip
		ali_clipped_end = seq_hard_clipped_end
		if right_clip["type"] and right_clip["type"] == "S":
			ali_clipped_end -= right_clip["length"]
		## right result = min(ref_len,right read end) - right aligned end
		right_res = min(ref_len,seq_end) - ali_clipped_end
	return (left_res,right_res)

def get_ref_length_sam(line,ref_lens):
	"""Get reference length from @ LN tag and insert into dict"""
	if line.startswith('@SQ\t'):
		ref_search = re.search(r'\tSN:(\S+)\b',line)
		if ref_search:
			ref_name = ref_search.group(1)
			assert ref_name, "Empty reference name in {}".format(line)
			len_search = re.search(r'\tLN:(\d+)\b',line)
			assert len_search,"Could not find length tag in {}".format(line)
			ref_len = int(len_search.group(1))
			if ref_name in ref_lens:
				logging.warning("Reference name is found second time in line {}".format(line))
			ref_lens[ref_name] = ref_len


def modify_bowtie_sam(raw_bowtie_sam,max_mismatch,max_unaligned_overlap):
	# fix sam flags for comprehensive pileup and filter out spurious alignments
	ref_lens = {}
	with open(raw_bowtie_sam) as sam, open(raw_bowtie_sam + '.mod', 'w') as sam_mod:
		for line in sam:
			if not line.startswith('@'):
				fields = line.split('\t')
				left_unali,right_unali = get_unaligned_read_end_lengths_sam(fields,ref_lens[fields[2].strip()])
				if left_unali > max_unaligned_overlap or right_unali > max_unaligned_overlap:
					#logging.debug("Excluding read from SAM file due to too long unaligned end overlapping the reference: {}".format(line))
					continue
				flag = int(fields[1])
				flag = (flag - 256) if (flag & 256) else flag
				m = re.search("NM:i:(\d+)\s",line)
				if m != None:
					num_mismatch = m.group(1)
					if int(num_mismatch) <= int(max_mismatch):
						sam_mod.write('\t'.join([fields[0], str(flag)] + fields[2:]))
				else:
					logging.info('Excluding read from SAM file due to missing NM (num mismatches) field: ' + fields[0])
					num_mismatch = 0
			else:
				get_ref_length_sam(line,ref_lens)
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
							unique_gene_symbols = False # already seen this cluster symbol
							logging.info( "Non-unique:" +  gene_cluster + ", " + cluster_symbol)
					else:
						gene_cluster_symbols[gene_cluster] = cluster_symbol
				else:
					# treat as unclustered database, use whole header
					gene_cluster = cluster_symbol = name.split()[0] # no spaces allowed
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

			# record gene (cluster):
			if gene_cluster not in gene_clusters:
				gene_clusters.append(gene_cluster)

	if len(delimiter_check) > 0:
		print "Warning! MLST delimiter is " + delimiter + " but these genes may violate the pattern and cause problems:"
		print ",".join(delimiter_check)

	return size, gene_clusters, unique_gene_symbols, unique_allele_symbols, gene_cluster_symbols


def read_pileup_data(pileup_file, size, prob_err, consensus_file = ""):
	with open(pileup_file) as pileup:
		prob_success = 1 - prob_err	# Set by user, default is prob_err = 0.01
		hash_alignment = {}
		hash_max_depth = {}
		hash_edge_depth = {}
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
			max_depth = 1
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
			consensus_seq = ""

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
				nuc_counts = {}

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

					elif aligned_bases[i].upper() in "ATCG":
						this_nuc = aligned_bases[i].upper()
						if this_nuc not in nuc_counts:
							nuc_counts[this_nuc] = 0
						nuc_counts[this_nuc] += 1

					i += 1

				# Save the most common nucleotide at this position
				consensus_nuc = nuc # by default use reference nucleotide
				max_freq = num_match # Number of bases matching the reference
				for nucleotide in nuc_counts:
					if nuc_counts[nucleotide] > max_freq:
						consensus_nuc = nucleotide
						max_freq = nuc_counts[nucleotide]
				consensus_seq += (consensus_nuc)

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

			# Determine the consensus sequence if required
			if consensus_file != "":
				if consensus_file.split(".")[-2] == "new_consensus_alleles":
					consensus_type = "variant"
				elif consensus_file.split(".")[-2] == "all_consensus_alleles":
					consensus_type = "consensus"
				with open(consensus_file, "a") as consensus_outfile:
					consensus_outfile.write(">{0}.{1} {2}\n".format(allele, consensus_type, pileup_file.split(".")[1].split("__")[1]))
					outstring = consensus_seq + "\n"
					consensus_outfile.write(outstring)

			# Finished reading pileup for this allele

			# Check for missing bases at the end of the allele
			if nuc_num < allele_size:
				total_missing_bases += abs(allele_size - nuc_num)
				# determine penalty based on coverage of last 2 bases
				penalty = float(position_depths[nuc_num-1] + position_depths[nuc_num-2])/2
				m = min(position_depths[nuc_num-1],position_depths[nuc_num-2])
				hash_alignment[allele].append((0, round(penalty), prob_success))
				if next_to_del_depth > m:
					next_to_del_depth = m # keep track of lowest near-del depth for reporting

			# Calculate allele summary stats and save
			avg_depth = round(total_depth / float(allele_line),3)
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
					hash_alignment[allele].append((0, round(penalty), prob_success))
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


def score_alleles(args, mapping_files_pre, hash_alignment, hash_max_depth, hash_edge_depth,
		avg_depth_allele, coverage_allele, mismatch_allele, indel_allele, missing_allele,
		size_allele, next_to_del_depth_allele, run_type,unique_gene_symbols, unique_allele_symbols):
	# sort into hash for each gene locus
	depth_by_gene = group_allele_dict_by_gene(dict( (allele,val) for (allele,val) in avg_depth_allele.items() \
			if (run_type == "mlst") or (coverage_allele[allele] > args.min_coverage) ),
			run_type,args,
			unique_gene_symbols,unique_allele_symbols)
	stat_depth_by_gene = dict(
			(gene,max(alleles.values())) for (gene,alleles) in depth_by_gene.items()
			)
	allele_to_gene = dict_of_dicts_inverted_ind(depth_by_gene)

	if args.save_scores:
		scores_output = file(mapping_files_pre + '.scores', 'w')
		scores_output.write("Allele\tScore\tAvg_depth\tEdge1_depth\tEdge2_depth\tPercent_coverage\tSize\tMismatches\tIndels\tTruncated_bases\tDepthNeighbouringTruncation\tmaxMAF\tLeastConfident_Rate\tLeastConfident_Mismatches\tLeastConfident_Depth\tLeastConfident_Pvalue\n")

	scores = {} # key = allele, value = score
	mix_rates = {} # key = allele, value = highest minor allele frequency, 0 -> 0.5

	for allele in hash_alignment:
		#stat_depth_allele = avg_depth_allele[allele]
		if (run_type == "mlst") or (coverage_allele[allele] > args.min_coverage):
			gene = allele_to_gene[allele]
			pvals = []
			min_pval = 1.0
			min_pval_data = (999,999) # (mismatch, depth) for position with lowest p-value
			mix_rate = 0 # highest minor allele frequency 0 -> 0.5
			for nuc_info in hash_alignment[allele]:
				if nuc_info is not None:
					match, mismatch, prob_success = nuc_info
					max_depth = hash_max_depth[allele]
					if match > 0 or mismatch > 0:
						# One-tailed test - prob to get that many or fewer matches
						p_value = binom.cdf(match,match+mismatch,prob_success)
						# Weight pvalue by (depth/max_depth)
						weight = (match + mismatch) / float(max_depth)
						p_value *= weight
						if p_value < min_pval:
							min_pval = p_value
							min_pval_data = (mismatch,match + mismatch)
						if p_value > 0:
							p_value = -log(p_value, 10)
						else:
							p_value = 1000
						pvals.append(p_value)
						mismatch_prop = float(match)/float(match+mismatch)
						if min(mismatch_prop, 1-mismatch_prop) > mix_rate:
							mix_rate = min(mismatch_prop, 1-mismatch_prop)
			# Fit linear model to observed Pval distribution vs expected Pval distribution (QQ plot)
			pvals.sort(reverse=True)
			len_obs_pvals = len(pvals)
			exp_pvals = range(1, len_obs_pvals + 1)
			exp_pvals2 = [-log(float(ep) / (len_obs_pvals + 1), 10) for ep in exp_pvals]

			# Slope is score
			slope, _intercept, _r_value, _p_value, _std_err = linregress(exp_pvals2, pvals)

			# Store all scores for later processing
			scores[allele] = slope
			mix_rates[allele] = mix_rate

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
						str(this_coverage), str(this_size), str(this_mismatch), str(this_indel), str(this_missing), str(this_next_to_del_depth), str(mix_rate), str(float(min_pval_data[0])/min_pval_data[1]),str(min_pval_data[0]),str(min_pval_data[1]),str(min_pval)]) + '\n')

	if args.save_scores:
		scores_output.close()

	return(scores,mix_rates)

# Check that an acceptable version of a command is installed
# Exits the program if it can't be found.
# - command_list is the command to run to determine the version.
# - version_identifier is the unique string we look for in the stdout of the program.
# - command_name is the name of the command to show in error messages.
# - required_version is the version number to show in error messages.
def check_command_version(command_list, version_identifier, command_name, required_version):
	try:
		command_stdout = check_output(command_list, stderr=STDOUT)
	except OSError as e:
		logging.error("Failed command: {}".format(' '.join(command_list)))
		logging.error(str(e))
		logging.error("Could not determine the version of {}.".format(command_name))
		logging.error("Do you have {} installed in your PATH?".format(command_name))
		exit(-1)
	except CalledProcessError as e:
		# some programs such as samtools return a non-zero exit status
		# when you ask for the version (sigh). We ignore it here.
		command_stdout = e.output

	if version_identifier not in command_stdout:
		logging.error("Incorrect version of {} installed.".format(command_name))
		logging.error("{} version {} is required by SRST2.".format(command_name, required_version))
		exit(-1)


# allow multiple specific versions that have been specifically tested
def check_bowtie_version():
	return check_command_versions([get_bowtie_execs()[0], '--version'], 'version ', 'bowtie',
								  ['2.1.0','2.2.3','2.2.4','2.2.5','2.2.6','2.2.7','2.2.8','2.2.9'])

def check_samtools_version():
	return check_command_versions([get_samtools_exec()], 'Version: ', 'samtools',
								  ['0.1.18','0.1.19','1.0','1.1','1.2','1.3','(0.1.18 is '
																			 'recommended)'])

def check_command_versions(command_list, version_prefix, command_name, required_versions):
	try:
		command_stdout = check_output(command_list, stderr=STDOUT)
	except OSError as e:
		logging.error("Failed command: {}".format(' '.join(command_list)))
		logging.error(str(e))
		logging.error("Could not determine the version of {}.".format(command_name))
		logging.error("Do you have {} installed in your PATH?".format(command_name))
		exit(-1)
	except CalledProcessError as e:
		# some programs such as samtools return a non-zero exit status
		# when you ask for the version (sigh). We ignore it here.
		command_stdout = e.output

	for v in required_versions:
		if version_prefix + v in command_stdout:
			return v

	logging.error("Incorrect version of {} installed.".format(command_name))
	logging.error("{} versions compatible with SRST2 are ".format(command_name) + ", ".join(required_versions))
	exit(-1)

def get_bowtie_execs():
	'Return the "best" bowtie2 executables'

	exec_from_environment = os.environ.get('SRST2_BOWTIE2')
	if exec_from_environment and os.path.isfile(exec_from_environment):
		bowtie2_exec = exec_from_environment
	else:
		bowtie2_exec = None

	exec_from_environment = os.environ.get('SRST2_BOWTIE2_BUILD')
	if exec_from_environment and os.path.isfile(exec_from_environment):
		bowtie2_build_exec = exec_from_environment
	elif bowtie2_exec and os.path.isfile(bowtie2_exec+'-build'):
		bowtie2_build_exec = bowtie2_exec+'-build'
	else:
		bowtie2_build_exec = 'bowtie2-build'

	if bowtie2_exec is None:
		bowtie2_exec = 'bowtie2'

	return (bowtie2_exec, bowtie2_build_exec)

def run_bowtie(mapping_files_pre,sample_name,fastqs,args,db_name,db_full_path):

	logging.info("Starting mapping with bowtie2")

	check_bowtie_version()
	check_samtools_version()

	command = [get_bowtie_execs()[0]]

	if len(fastqs)==1:
		# single end
		command += ['-U', fastqs[0]]
	elif len(fastqs)==2:
		# paired end
		command += ['-1', fastqs[0], '-2', fastqs[1]]

	sam = mapping_files_pre + ".sam"
	logging.info('Output prefix set to: ' + mapping_files_pre)

	command += ['-S', sam,
				'-' + args.read_type,	# add a dash to the front of the option
				'--very-sensitive-local',
				'--no-unal',
				'-a',					 # Search for and report all alignments
				'-x', db_full_path			   # The index to be aligned to
			   ]
	if args.threads > 1:
		command += ['--threads', str(args.threads)]

	if args.stop_after:
		try:
			command += ['-u',str(int(args.stop_after))]
		except ValueError:
			print "WARNING. You asked to stop after mapping '" + args.stop_after + "' reads. I don't understand this, and will map all reads. Please speficy an integer with --stop_after or leave this as default to map 1 million reads."

	if args.other:
		x = args.other
		x = x.replace('\\','')
		command += x.split()

	if args.use_existing_bowtie2_sam and os.path.exists(sam):
		logging.info(' Using existing Bowtie2 SAM in ' + sam)

	else:
		logging.info('Aligning reads to index {} using bowtie2...'.format(db_full_path))

		run_command(command)

	return(sam)

def get_samtools_exec():
	'Return the "best" samtools executable'

	exec_from_environment = os.environ.get('SRST2_SAMTOOLS')
	if exec_from_environment and os.path.isfile(exec_from_environment):
		return exec_from_environment
	else:
		return 'samtools'

def get_pileup(args, mapping_files_pre, raw_bowtie_sam, bowtie_sam_mod, fasta, pileup_file):
	# Analyse output with SAMtools
	samtools_exec = get_samtools_exec()
	samtools_v1 = check_samtools_version().split('.')[0] == '1'  # Usage changed in version 1.0
	logging.info('Processing Bowtie2 output with SAMtools...')
	logging.info('Generate and sort BAM file...')
	out_file_bam = mapping_files_pre + ".unsorted.bam"
	view_command = [samtools_exec, 'view']
	if args.threads > 1 and samtools_v1:
		view_command += ['-@', str(args.threads)]
	view_command += ['-b', '-o', out_file_bam, '-q', str(args.mapq), '-S', bowtie_sam_mod]
	run_command(view_command)
	out_file_bam_sorted = mapping_files_pre + ".sorted"
	sort_command = [samtools_exec, 'sort']
	if samtools_v1:
		if args.threads > 1:
			sort_command += ['-@', str(args.threads)]
		temp = mapping_files_pre + ".sort_temp"
		sort_command += ['-o', out_file_bam_sorted + '.bam', '-O', 'bam', '-T', temp, out_file_bam]
	else:  # samtools 0.x
		sort_command += [out_file_bam, out_file_bam_sorted]
	run_command(sort_command)

	# Delete interim files (sam, modified sam, unsorted bam) unless otherwise specified.
	# Note users may also want to delete final sorted bam and pileup on completion to save space.
	if not args.keep_interim_alignment:
		logging.info('Deleting sam and bam files that are not longer needed...')
		del_filenames = [raw_bowtie_sam, bowtie_sam_mod, out_file_bam]
		for f in del_filenames:
			logging.info('Deleting ' + f)
			os.remove(f)

	logging.info('Generate pileup...')
	with open(pileup_file, 'w') as sam_pileup:
		mpileup_command = [samtools_exec, 'mpileup', '-L', '1000', '-f', fasta,
					 '-Q', str(args.baseq), '-q', str(args.mapq), '-B', out_file_bam_sorted + '.bam']
		if args.samtools_args:
			x = args.samtools_args
			x = x.replace('\\','')
			mpileup_command += x.split()
		run_command(mpileup_command, stdout=sam_pileup)

def calculate_ST(allele_scores, ST_db, gene_names, sample_name, mlst_delimiter, avg_depth_allele, mix_rates):
	allele_numbers = [] # clean allele calls for determing ST. order is taken from gene names, as in ST definitions file
	alleles_with_flags = [] # flagged alleles for printing (* if mismatches, ? if depth issues)
	mismatch_flags = [] # allele/diffs
	uncertainty_flags = [] #allele/uncertainty
#	st_flags = [] # (* if mismatches, ? if depth issues)
	depths = [] # depths for each typed locus
	mafs = [] # minor allele freqencies for each typed locus

	# get allele numbers & info
	for gene in gene_names:
		if gene in allele_scores:
			(allele,diffs,depth_problem,divergence) = allele_scores[gene]
			allele_number = allele.split(mlst_delimiter)[-1]
			depths.append(avg_depth_allele[allele])
			mix_rate = mix_rates[allele]
			mafs.append(mix_rate)
		else:
			allele_number = "-"
			diffs = ""
			depth_problem = ""
			mix_rate = ""
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
				(allele,diffs,depth_problems,divergence) = allele_scores[gene]
				print allele,
			print
			clean_st = "NF"
	else:
		clean_st = "ND"

	# add flags for reporting
	st = clean_st
	if len(mismatch_flags) > 0:
		for m in mismatch_flags:
			if m.split("/")[1] != "trun":
				st = clean_st + "*" # trun indicates only that a truncated form had lower score, which isn't a mismatch
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

	# maximum maf across locus
	if len(mafs) > 0:
		max_maf = max(mafs)
	else:
		max_maf = 0

	return (st,clean_st,alleles_with_flags,mismatch_flags,uncertainty_flags,mean_depth,max_maf)

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

def get_allele_name_from_db(allele,run_type,args,unique_allele_symbols=False,unique_cluster_symbols=False):

	if run_type != "mlst":
		# header format: >[cluster]___[gene]___[allele]___[uniqueID] [info]
		allele_parts = allele.split()
		allele_detail = allele_parts.pop(0)
		allele_info = allele_detail.split("__")

		if len(allele_info)>3:
			cluster_id = allele_info[0] # ID number for the cluster
			gene_name = allele_info[1] # gene name/symbol for the cluster
			allele_name = allele_info[2] # specific allele name
			seqid = allele_info[3] # unique identifier for this seq
		else:
			cluster_id = gene_name = allele_name = seqid = allele

		if not unique_allele_symbols:
			allele_name += "_" + seqid

	else:
		gene_name = allele.split(args.mlst_delimiter)
		allele_name = gene_name[1]
		gene_name = gene_name[0]
		seqid = None
		cluster_id = None
	return gene_name, allele_name, cluster_id, seqid

def create_allele_pileup(allele_name, all_pileup_file):
	output_components = all_pileup_file.split("/")
	if len(output_components) > 1:
		all_pileup_file_name = os.path.basename(all_pileup_file)
		all_pileup_file_dir = os.path.dirname(all_pileup_file)
		outpileup = all_pileup_file_dir + '/' + allele_name + "." + all_pileup_file_name
	else:
		outpileup = allele_name + "." + all_pileup_file
	with open(outpileup, 'w') as allele_pileup:
		with open(all_pileup_file) as all_pileup:
			for line in all_pileup:
				if line.split()[0] == allele_name:
					allele_pileup.write(line)
	return outpileup


def group_allele_dict_by_gene(by_allele,run_type,args,unique_cluster_symbols=False, unique_allele_symbols=False):
	# sort into hash for each gene locus
	by_gene = collections.defaultdict(dict) # key1 = gene, key2 = allele, value = original value

	if run_type=="mlst":
		component_ind = 0 # gene_name
	else:
		component_ind = 2 # cluster_id
	for allele in by_allele:
		gene_name = get_allele_name_from_db(allele,run_type,args,unique_allele_symbols,unique_cluster_symbols)[component_ind]
		by_gene[gene_name][allele] = by_allele[allele]
	return dict(by_gene)


def dict_of_dicts_inverted_ind(dd):
	res = dict()
	for (key,val) in dd.items():
		res.update(dict((key2,key) for key2 in val))
	return res

def parse_scores(run_type,args,scores, hash_edge_depth,
					avg_depth_allele, coverage_allele, mismatch_allele, indel_allele,
					missing_allele, size_allele, next_to_del_depth_allele,
					unique_cluster_symbols,unique_allele_symbols, pileup_file):

	# sort into hash for each gene locus
	scores_by_gene = group_allele_dict_by_gene(dict( (allele,val) for (allele,val) in scores.items() \
			if coverage_allele[allele] > args.min_coverage ),
			run_type,args,
			unique_cluster_symbols,unique_allele_symbols)

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

		divergence = float(mismatch_allele[top_allele]) / float( size_allele[top_allele] - missing_allele[top_allele] )

		# check for truncated
		if differences != "" or not adequate_depth:
			# if there are SNPs or not enough depth to trust the result, no need to screen next best match
			results[gene] = (top_allele, differences, depth_problem, divergence)
		else:
			# looks good but this could be a truncated version of the real allele; check for longer versions
			truncation_override = False
			if len(scores_sorted) > 1:
				(next_best_allele,next_best_score) = scores_sorted[1]
				if size_allele[next_best_allele] > size_allele[top_allele]:
					# next best is longer, top allele could be a truncation?
					if (mismatch_allele[next_best_allele] + indel_allele[next_best_allele] + missing_allele[next_best_allele]) == 0:
						# next best also has no mismatches
						if (next_best_score - top_score)/top_score < args.truncation_score_tolerance:
							# next best has score within 10% of this one
							truncation_override = True
			if truncation_override:
				results[gene] = (next_best_allele, "trun", "", divergence) # no diffs but report this call is based on truncation test
				final_allele_reported = next_best_allele
			else:
				results[gene] = (top_allele, "", "",divergence) # no caveats to report

		# Check if there are any potential new alleles
		if depth_problem == "" and divergence > 0:
			new_allele = True
			# Get the consensus for this new allele and write it to file
			if args.report_new_consensus or args.report_all_consensus:
				new_alleles_filename = args.output + ".new_consensus_alleles.fasta"
				allele_pileup_file = create_allele_pileup(results[gene][0], pileup_file)
				read_pileup_data(allele_pileup_file, size_allele, args.prob_err, consensus_file = new_alleles_filename)
		if args.report_all_consensus:
			new_alleles_filename = args.output + ".all_consensus_alleles.fasta"
			allele_pileup_file = create_allele_pileup(results[gene][0], pileup_file)
			read_pileup_data(allele_pileup_file, size_allele, args.prob_err, consensus_file = new_alleles_filename)

	return results # (allele, diffs, depth_problem, divergence)


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
			m=re.match("(.*)(_S.*)(_L.*)(_R.*)(_.*)", file_name_before_ext)
			if m==None:
				fileSets[file_name_before_ext] = [fastq]
			else:
				fileSets[m.groups()[0]] = [fastq] # Illumina names
			num_single_readsets += 1

	elif args.input_pe:
		# paired end
		forward_reads = {} # key = sample, value = full path to file
		reverse_reads = {} # key = sample, value = full path to file
		num_paired_readsets = 0
		num_single_readsets = 0
		for fastq in args.input_pe:
			(file_path,file_name_before_ext,full_ext) = get_readFile_components(fastq)
			# try to match to MiSeq format:
			m=re.match("(.*)(_S.*)(_L.*)(_R.*)(_.*)", file_name_before_ext)
			if m==None:
				# not default Illumina file naming format, expect simple/ENA format
				m=re.match("(.*)("+args.forward+")$",file_name_before_ext)
				if m!=None:
					# store as forward read
					(baseName,read) = m.groups()
					forward_reads[baseName] = fastq
				else:
					m=re.match("(.*)("+args.reverse+")$",file_name_before_ext)
					if m!=None:
					# store as reverse read
						(baseName,read) = m.groups()
						reverse_reads[baseName] = fastq
					else:
						logging.info("Could not determine forward/reverse read status for input file " + fastq)
			else:
				# matches default Illumina file naming format, e.g. m.groups() = ('samplename', '_S1', '_L001', '_R1', '_001')
				baseName, read  = m.groups()[0], m.groups()[3]
				if read == "_R1":
					forward_reads[baseName] = fastq
				elif read == "_R2":
					reverse_reads[baseName] = fastq
				else:
					logging.info( "Could not determine forward/reverse read status for input file " + fastq )
					logging.info( "  this file appears to match the MiSeq file naming convention (samplename_S1_L001_[R1]_001), but we were expecting [R1] or [R2] to designate read as forward or reverse?" )
					fileSets[file_name_before_ext] = fastq
					num_single_readsets += 1
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

	if os.stat(infile).st_size == 0:
		logging.info("WARNING: Results file provided is empty: " + infile)
		return False, False, False

	results_info = infile.split("__")
	if len(results_info) > 1:

		if re.search("compiledResults",infile)!=None:
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
							gene = header[i] # cluster_id
							results[sample][gene] = line_split[i]

		elif dbtype == "mlst":
			results = {} # key = sample, value = MLST string
			with open(infile) as f:
				header = 0
				for line in f:
					if header > 0:
						results[line.split("\t")[0]] = line.rstrip()
						if "maxMAF" not in header:
							results[line.split("\t")[0]] += "\tNC" # empty column for maxMAF
					else:
						header = line.rstrip()
						results[line.split("\t")[0]] = line.rstrip() # store header line too (index "Sample")
						if "maxMAF" not in header:
							results[line.split("\t")[0]] += "\tmaxMAF" # add column for maxMAF

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
						if n_cols > 1:
							if header[1] == "ST":
								# there is mlst data reported
								mlst_cols = 2 # first locus column
								while header[mlst_cols] != "depth":
									mlst_cols += 1
								results["Sample"]["mlst"] = "\t".join(line_split[0:(mlst_cols+1)])
								results["Sample"]["mlst"] += "\tmaxMAF" # add to mlst header even if not encountered in this file, as it may be in others
								if header[mlst_cols+1] == "maxMAF":
									mlst_cols += 1 # record maxMAF column within MLST data, if present
							else:
								# no mlst data reported
								dbtype = "genes"
								logging.info("No MLST data in compiled results file " + infile)
						else:
							# no mlst data reported
							dbtype = "genes"
							logging.info("No MLST data in compiled results file " + infile)

					else:
						sample = line_split[0]
						if mlst_cols > 0:
							results[sample]["mlst"] = "\t".join(line_split[0:(mlst_cols+1)])
							if "maxMAF" not in header:
								results[sample]["mlst"] += "\t" # add to mlst section even if not encountered in this file, as it may be in others
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
	mix_rates = {}
	scores = {}

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

	return hash_edge_depth, avg_depth_allele, coverage_allele, mismatch_allele, indel_allele, \
			missing_allele, size_allele, next_to_del_depth_allele, scores, mix_rates

def run_srst2(args, fileSets, dbs, run_type):

	db_reports = [] # list of db-specific output files to return
	db_results_list = [] # list of results hashes, one per db

	for fasta in dbs:
		db_reports, db_results_list = process_fasta_db(args, fileSets, run_type, db_reports,
													   db_results_list, fasta)

	return db_reports, db_results_list

def samtools_index(fasta_file):
	check_samtools_version()
	fai_file = fasta_file + '.fai'
	if not os.path.exists(fai_file):
		run_command([get_samtools_exec(), 'faidx', fasta_file])
	return fai_file

def process_fasta_db(args, fileSets, run_type, db_reports, db_results_list, fasta):

	logging.info('Processing database ' + fasta)

	db_path, db_name = os.path.split(fasta) # database
	(db_name,db_ext) = os.path.splitext(db_name)
	db_results = "__".join([args.output,run_type,db_name,"results.txt"])
	db_report = file(db_results,"w")
	db_reports.append(db_results)

	# Get sequence lengths and gene names
	#  lengths are needed for MLST heuristic to distinguish alleles from their truncated forms
	#  gene names read from here are needed for non-MLST dbs
	fai_file = samtools_index(fasta)
	size, gene_names, unique_gene_symbols, unique_allele_symbols, cluster_symbols = \
		parse_fai(fai_file,run_type,args.mlst_delimiter)

	# Prepare for MLST reporting
	ST_db = False
	if run_type == "mlst":
		results = {} # key = sample, value = ST string for printing
		if args.mlst_definitions:
			# store MLST profiles, replace gene names (we want the order as they appear in this file)
			ST_db, gene_names = parse_ST_database(args.mlst_definitions,gene_names)
		db_report.write("\t".join(["Sample","ST"]+gene_names+["mismatches","uncertainty","depth","maxMAF"]) + "\n")
		results["Sample"] = "\t".join(["Sample","ST"]+gene_names+["mismatches","uncertainty","depth","maxMAF"])

	else:
		# store final results for later tabulation
		results = collections.defaultdict(dict) #key1 = sample, key2 = gene, value = allele

	gene_list = [] # start with empty gene list; will add genes from each genedb test

	# determine maximum mismatches per read to use for pileup
	if run_type == "mlst":
		max_mismatch = args.mlst_max_mismatch
	else:
		max_mismatch = args.gene_max_mismatch

	# Align and score each read set against this DB
	for sample_name in fileSets:
		logging.info('Processing sample ' + sample_name)
		fastq_inputs = fileSets[sample_name] # reads

		try:
			# try mapping and scoring this fileset against the current database
			# update the gene_list list and results dict with data from this strain
			# __mlst__ will be printed during this routine if this is a mlst run
			# __fullgenes__ will be printed during this routine if requested and this is a gene_db run
			gene_list, results = \
				map_fileSet_to_db(args, sample_name, fastq_inputs, db_name, fasta,size,gene_names,
								  unique_gene_symbols, unique_allele_symbols, run_type,ST_db,
								  results,gene_list, db_report, cluster_symbols, max_mismatch)

		# if we get an error from one of the commands we called
		# log the error message, record as failed, and continue onto the next fasta db
		except CommandError as e:
			logging.error(e.message)
			# record results as unknown, so we know that we did attempt to analyse this readset
			if run_type == "mlst":
				st_result_string = "\t".join( [sample_name,"failed"] + ["-"] * (len(gene_names) + 4)) # record missing results
				db_report.write( st_result_string + "\n")
				logging.info(" " + st_result_string)
				results[sample_name] = st_result_string
			else:
				logging.info(" failed gene detection")
				results[sample_name]["failed"] = True # so we know that we tried this strain

	if run_type != "mlst":
		# tabulate results across samples for this gene db (i.e. __genes__ file)
		logging.info('Tabulating results for database {} ...'.format(fasta))
		gene_list.sort()
		db_report.write("\t".join(["Sample"]+gene_list)+"\n") # report header row
		for sample_name in fileSets:
			db_report.write(sample_name)
			if sample_name in results:
				# print results
				if "failed" not in results[sample_name]:
					for cluster_id in gene_list:
						if cluster_id in results[sample_name]:
							db_report.write("\t"+results[sample_name][cluster_id]) # print full allele name
						else:
							db_report.write("\t-") # no hits for this gene cluster
				else:
					# no data on this, as the sample failed mapping
					for cluster_id in gene_list:
						db_report.write("\t-f") #
						results[sample_name][cluster_id] = "-f" # record as unknown as this strain failed
			else:
				# no data on this because genes were not found (but no mapping errors)
				for cluster_id in gene_list:
					db_report.write("\t-?") #
					results[sample_name][cluster_id] = "-" # record as absent
			db_report.write("\n")

	# Finished with this database
	logging.info('Finished processing for database {} ...'.format(fasta))
	db_report.close()
	db_results_list.append(results)

	return db_reports, db_results_list

def map_fileSet_to_db(args, sample_name, fastq_inputs, db_name, fasta, size, gene_names,
					  unique_gene_symbols, unique_allele_symbols, run_type, ST_db, results,
					  gene_list, db_report, cluster_symbols, max_mismatch):

	mapping_files_pre = args.output + '__' + sample_name + '.' + db_name
	pileup_file = mapping_files_pre + '.pileup'
	scores_file = mapping_files_pre + '.scores'

	# Get or read scores

	if args.use_existing_scores and os.path.exists(scores_file):

		logging.info(' Using existing scores in ' + scores_file)

		# read in scores and info from existing scores file
		hash_edge_depth, avg_depth_allele, coverage_allele, \
				mismatch_allele, indel_allele, missing_allele, size_allele, \
				next_to_del_depth_allele, scores, mix_rates = read_scores_file(scores_file)

	else:

		# Get or read pileup

		if args.use_existing_pileup and os.path.exists(pileup_file):
			logging.info(' Using existing pileup in ' + pileup_file)

		else:

			# run bowtie against this db
			bowtie_sam = run_bowtie(mapping_files_pre,sample_name,fastq_inputs,args,db_name,fasta)

			# Modify Bowtie's SAM formatted output so that we get secondary
			# alignments in downstream pileup
			(raw_bowtie_sam,bowtie_sam_mod) = modify_bowtie_sam(bowtie_sam,max_mismatch,\
					max_unaligned_overlap=args.max_unaligned_overlap)

			# generate pileup from sam (via sorted bam)
			get_pileup(args, mapping_files_pre, raw_bowtie_sam, bowtie_sam_mod, fasta, pileup_file)

		# Get scores

		# Process the pileup and extract info for scoring and reporting on each allele
		logging.info(' Processing SAMtools pileup...')
		hash_alignment, hash_max_depth, hash_edge_depth, avg_depth_allele, coverage_allele, \
				mismatch_allele, indel_allele, missing_allele, size_allele, next_to_del_depth_allele= \
				read_pileup_data(pileup_file, size, args.prob_err)

		# Generate scores for all alleles (prints these and associated info if verbose)
		#   result = dict, with key=allele, value=score
		logging.info(' Scoring alleles...')
		scores, mix_rates = score_alleles(args, mapping_files_pre, hash_alignment, hash_max_depth, hash_edge_depth, \
				avg_depth_allele, coverage_allele, mismatch_allele, indel_allele, missing_allele, \
				size_allele, next_to_del_depth_allele, run_type,unique_gene_symbols, unique_allele_symbols)

	# GET BEST SCORE for each gene/cluster
	#  result = dict, with key = gene, value = (allele,diffs,depth_problem)
	#			for MLST DBs, key = gene = locus, allele = gene-number
	#			for gene DBs, key = gene = cluster ID, allele = cluster__gene__allele__id
	#  for gene DBs, only those alleles passing the coverage cutoff are returned

	allele_scores = parse_scores(run_type, args, scores, \
			hash_edge_depth, avg_depth_allele, coverage_allele, mismatch_allele, \
			indel_allele, missing_allele, size_allele, next_to_del_depth_allele,
			unique_gene_symbols, unique_allele_symbols, pileup_file)

	# REPORT/RECORD RESULTS

	# Report MLST results to __mlst__ file
	if run_type == "mlst" and len(allele_scores) > 0:

		# Calculate ST and get info for reporting
		(st,clean_st,alleles_with_flags,mismatch_flags,uncertainty_flags,mean_depth,max_maf) = \
				calculate_ST(allele_scores, ST_db, gene_names, sample_name, args.mlst_delimiter, avg_depth_allele, mix_rates)

		# Print to MLST report, log and save the result
		st_result_string = "\t".join([sample_name,st]+alleles_with_flags+[";".join(mismatch_flags),";".join(uncertainty_flags),str(mean_depth),str(max_maf)])
		db_report.write( st_result_string + "\n")
		logging.info(" " + st_result_string)
		results[sample_name] = st_result_string

		# Make sure scores are printed if there was uncertainty in the call
		scores_output_file = mapping_files_pre + '.scores'
		if uncertainty_flags != ["-"] and not args.save_scores and not os.path.exists(scores_output_file):
			# print full score set
			logging.info("Printing all MLST scores to " + scores_output_file)
			scores_output = file(scores_output_file, 'w')
			scores_output.write("Allele\tScore\tAvg_depth\tEdge1_depth\tEdge2_depth\tPercent_coverage\tSize\tMismatches\tIndels\tTruncated_bases\tDepthNeighbouringTruncation\tMmaxMAF\n")
			for allele in scores.keys():
				score = scores[allele]
				scores_output.write('\t'.join([allele, str(score), str(avg_depth_allele[allele]), \
					str(hash_edge_depth[allele][0]), str(hash_edge_depth[allele][1]), \
					str(coverage_allele[allele]), str(size_allele[allele]), str(mismatch_allele[allele]), \
					str(indel_allele[allele]), str(missing_allele[allele]), str(next_to_del_depth_allele[allele]), str(round(mix_rates[allele],3))]) + '\n')
			scores_output.close()

	# Record gene results for later processing and optionally print detailed gene results to __fullgenes__ file
	elif run_type == "genes" and len(allele_scores) > 0:
		if args.no_gene_details:
			full_results = "__".join([args.output,"fullgenes",db_name,"results.txt"])
			logging.info("Printing verbose gene detection results to " + full_results)
			if os.path.exists(full_results):
				f = file(full_results,"a")
			else:
				f = file(full_results,"w") # create and write header
				f.write("\t".join(["Sample","DB","gene","allele","coverage","depth","diffs","uncertainty","divergence","length", "maxMAF","clusterid","seqid","annotation"])+"\n")
		for gene in allele_scores:
			(allele,diffs,depth_problem,divergence) = allele_scores[gene] # gene = top scoring alleles for each cluster
			gene_name, allele_name, cluster_id, seqid = \
				get_allele_name_from_db(allele,run_type,args,unique_allele_symbols,unique_gene_symbols)

			# store for gene result table only if divergence passes minimum threshold:
			if divergence*100 <= float(args.max_divergence):
				column_header = cluster_symbols[cluster_id]
				results[sample_name][column_header] = allele_name
				if diffs != "":
					results[sample_name][column_header] += "*"
				if depth_problem != "":
					results[sample_name][column_header] += "?"
				if column_header not in gene_list:
					gene_list.append(column_header)

			# write details to full genes report
			if args.no_gene_details:

				# get annotation info
				header_string = os.popen(" ".join(["grep",allele,fasta]))
				try:
					header = header_string.read().rstrip().split()
					header.pop(0) # remove allele name
					if len(header) > 0:
						annotation = " ".join(header) # put back the spaces
					else:
						annotation = ""

				except:
					annotation = ""

				f.write("\t".join([sample_name,db_name,gene_name,allele_name,str(round(coverage_allele[allele],3)),str(avg_depth_allele[allele]),diffs,depth_problem,str(round(divergence*100,3)),str(size_allele[allele]),str(round(mix_rates[allele],3)),cluster_id,seqid,annotation])+"\n")

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
	st_counts = {} # key = ST, value = count

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
			this_mlst_cols = len(test_string_split)
			if (mlst_cols == 0) or (mlst_cols == this_mlst_cols):
				mlst_cols = this_mlst_cols
				blank_mlst_section = "\t?" * (mlst_cols-1) # blank MLST string in case some samples missing
				# use this data
				for sample in mlst_result:
					mlst_results_master[sample] = mlst_result[sample]
					if sample not in sample_list:
						sample_list.append(sample)
			elif mlst_cols != this_mlst_cols:
				# don't process this data further
				logging.info("Problem reconciling MLST data from two files, first MLST results encountered had " + str(mlst_cols) + " columns, this one has " + str(this_mlst_cols) + " columns?")
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
				st_data_split = mlst_results_master[sample].split("\t")
				if len(st_data_split) > 1:
					this_st = st_data_split[1]
					sample_info.append(mlst_results_master[sample])
				else:
					sample_info.append(sample+blank_mlst_section)
					this_st = "unknown" # something wrong with the string
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


def main():
	args = parse_args()

	# Check output directory
	output_components = args.output.split("/")
	if len(output_components) > 1:
		output_dir = "/".join(output_components[:-1])
		if not os.path.exists(output_dir):
			try:
				os.makedirs(output_dir)
				print "Created directory " + output_dir + " for output"
			except:
				print "Error. Specified output as " + args.output + " however the directory " + output_dir + " does not exist and our attempt to create one failed."

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

	# Delete consensus file if it already exists (so can use append file in functions)
	if args.report_new_consensus or args.report_all_consensus:
		new_alleles_filename = args.output + ".consensus_alleles.fasta"
		if os.path.exists(new_alleles_filename):
			os.remove(new_alleles_filename)

	# vars to store results
	mlst_results_hashes = [] # dict (sample->MLST result string) for each MLST output files created/read
	gene_result_hashes = [] # dict (sample->gene->result) for each gene typing output files created/read

	# parse list of file sets to analyse
	fileSets = read_file_sets(args) # get list of files to process

	if args.merge_paired:
		mate1 = [] # list of forward read files
		mate2 = [] # list of reverse read files
		for prefix in fileSets:
			reads = fileSets[prefix] # forward, reverse as list
			mate1.append(reads[0])
			mate2.append(reads[1])
		fileSets.clear() # remove all individual read sets
		fileSets["combined"] = [",".join(mate1),",".join(mate2)] # all input reads belong to same strain, ie single file set
		logging.info('Assuming all reads belong to single strain. A single combined result will be returned.')

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
		mlst_report, mlst_results = run_srst2(args, fileSets, args.mlst_db, "mlst")

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

		unique_results_files = list(OrderedDict.fromkeys(args.prev_output))

		for results_file in unique_results_files:

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

	elif args.prev_output:
		logging.info('One previous output file was provided, but there is no other data to compile with.')

	logging.info('SRST2 has finished.')


if __name__ == '__main__':
	main()
