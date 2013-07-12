#!/usr/bin/env python

# Parses scores from SRST(v2) to produce a summary table of multple samples. 
# 
# Author - Harriet Dashnow (h.dashnow@gmail.com)

import sys
import csv
import itertools
import os
from argparse import (ArgumentParser, FileType)

def parse_args():
	"Parse the input arguments, use '-h' for help"
	parser = ArgumentParser(description='SRST2 - Short Read Sequence Typer (v2): Outputs score and gene summaries. Reference fasta must have been in format cluster#__gene__allele')
	parser.add_argument(
		'--scores', nargs='+', type=str, required=True,
		help='One or more .scores files produced by srst2.py')
	parser.add_argument('--output', type=str, required=True, help='Base name for output files')
	parser.add_argument('--gene_summary', action="store_true", required=False,
		help='Output is a table summarising alleles for each gene.')
	parser.add_argument('--score_summary', action="store_true", required=False,
		help='Output is a table summarising the top scoring allele for each gene with more detailed score information.')
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
			allele_info = split_line[0].split("__")
			cluster = allele_info[0]
			gene = allele_info[1]
			allele = allele_info[2]
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

			full_cluster = cluster + "_" + gene
			#Only consider alleles with 90%+ coverage
			if coverage < 90:
				continue


			info_string = "\t".join([sample_name, full_cluster, gene, allele]) + "\t" + "\t".join(split_line[1:]) + "\n"
			if full_cluster not in hash_alleles:
				hash_alleles[full_cluster] = (score, gene, allele, info_string, difference)
			elif score < hash_alleles[full_cluster][0]:
				hash_alleles[full_cluster] = (score, gene, allele, info_string, difference)

	return (sample_name, hash_alleles)

def genes_table(score_files):
	all_clusters = []
	hash_samples = {}
	for scores_file in score_files:
		sample_name, hash_alleles = extract_top_allele(scores_file)
		for cluster_name in hash_alleles:
			if cluster_name not in all_clusters:
				all_clusters.append(cluster_name)
		if sample_name not in hash_samples:
			hash_samples[sample_name] = hash_alleles
	return all_clusters, hash_samples

def main():
	args = parse_args()
	if args.gene_summary == False and args.score_summary == False:
		print "You did not specify an output type. Select --gene_summary or --score_summary or both."
	score_files = args.scores

	if args.score_summary:
		with open(args.output + ".score_summary", "w") as score_outfile:
			header = "Sample\tCluster\tGene\tAllele\tScore\tAvg_depth\tEdge1_depth\tEdge2_depth\tPercent_coverage\tMismatches\tIndels\tMissing_bases\n"
			score_outfile.write(header)
			for scores_file in score_files:
				sample_name, hash_allele = extract_top_allele(scores_file)
				for cluster in hash_allele:
					outstring = hash_allele[cluster][3]
					score_outfile.write(outstring)

	if args.gene_summary:
		with open(args.output + ".gene_summary", "w") as gene_outfile:
			all_clusters, hash_samples = genes_table(score_files)
			header = "Sample"
			for cluster in all_clusters:
				header = header + "\t" + cluster
			gene_outfile.write(header + "\n")
			for sample in hash_samples:
				outstring = sample
				for cluster in all_clusters:
					try:
						allele = hash_samples[sample][cluster][2]
						difference = hash_samples[sample][cluster][4]
						if difference > 0:
							allele = allele + "*/" + str(difference)
					except KeyError:
						allele = "None"
					outstring = outstring + "\t" + allele
				gene_outfile.write(outstring + "\n")


if __name__ == '__main__':
	main()

# Some testing peramters:
#scores_file = "/vlsci/VR0082/shared/hdashnow/testing_srst2/mlst/pool10_tag2_1.fastq.gz.MLST_alleles.fasta.srst2.pileup.table.scores"

