import re as re
import sys as sys
import argparse as ap
import maptools as mt
import scoretools as st
import Bio.SeqIO

def printScores(scores, fastqs, regex):
	""""Print the scores for variants with significant matches in tab
	delimited format; sample_name st locus_1_variant similarity 
	least_significant_match corrected_significance_threshold #SNPs
	minimum_read_depth
	"""

	#Get the fastq identifiers
	sample = re.search("(.*)\.fastq", fastqs[0].split("/")[-1]).group(1)
	sample = re.sub("_\d+$", "", sample)

	try:
		#Determine the set of unique locus names based on the given regex
		loci = {re.search(regex, variant).group(1) for variant in scores.keys()}
		
		#Split the scores by locus
		scoresByLocus = splitScoresByLocus(scores, loci)
		
		for locus in scoresByLocus:
			for variant in st.sortKeysByScore(scoresByLocus[locus]):
				print "\t".join([sample, locus] + [str(s) for s in scores[variant]])
	except ValueError:
		print >> sys.stderr, "Failed to match locus for variant, scores not sorted by locus"
		for variant in st.sortKeysByScore(scores):
			match = re.search(regex, variant)
			if match:
				locus = match.group(1)
			else:
				locus = "unknown"
			
			print "\t".join([sample, "unmatched", variant] + [str(s) for s in scores[variant]])


def splitScoresByLocus(scores, loci):
	"""Split a dictionary of scores, keyed by variant name,
	into a dictionary of dictionaries keyed by loci name,
	then allele name."""

	scoresByLocus = dict()

	for locus in loci:
		scoresByLocus[locus] = dict()

	for variant in scores.keys():
		for locus in loci:
			if variant.startswith(locus):
				scoresByLocus[locus][variant] = scores[variant]

	return scoresByLocus


def findMaxLengths(variantFilename, regex):
	"""Find the maximum length of the variants for each locus"""
	maxLengths = dict()
	with open(args.var) as fin:
		for rec in Bio.SeqIO.parse(fin, "fasta"):
			locus = re.search(regex, rec.id).group(1)
			if locus in maxLengths:
				maxLengths[locus] = max(maxLengths[locus], len(rec.seq))
			else:
				maxLengths[locus] = len(rec.seq)
	
	return maxLengths


def main(args):
	"""Map a set of fastq files to a set of reference sequences and score each mapping"""
	maxLengths = findMaxLengths(args.var, args.regex)

	#Paired fastq files are given one after another since pattern matching is used 
	#to pass the set of file names and the shell sorts the filenames. Therefore we
	#define the step through the fastq list based on whether the reads are paired or not
	step = 1 if args.unpaired else 2
	for i in xrange(0, len(args.fastq), step):
		fastqs = args.fastq[i:i+step]
	
		#Get the sample name by first splitting on / to separate the file
		#name from the path if given, then take everything before the fastq
		#extensiom. Finally, remove the identifier for the pair number.
		sample = re.search("(.*)\.fastq", fastqs[0].split("/")[-1]).group(1)
		sample = re.sub("_\d+$", "", sample)

		#Map the reads to the reference sequences
		mt.mapReads(args.mapper, args.var, fastqs, args.base)
		
		#Print the header
		print ",".join(("sample", "locus", "variant", "pid", "entropy", "binomial", "sig", "mixed", "mismatches", "length", "min. depth", "mean depth", "max depth")) 
		
		#Calculate and print the scores, ignore unscored reads (_)
		if args.procs > 1:
			for score in st.scoreVariantsMultiProc(args.base + ".pileup", args.ranges, 0.01, maxLengths, args.mismatches, args.depth, args.pid, args.sig, args.procs):
				#Determine the locus name based on the given regex
				locus = re.search(args.regex, score.variant).group(1) 
				print ",".join([sample, locus, str(score)])
		else:
			for score in st.scoreVariants(args.base + ".pileup", args.ranges, args.regex, 0.01, maxLengths, args.mismatches, args.depth, args.pid, args.sig):
				#Determine the locus name based on the given regex
				locus = re.search(args.regex, score.variant).group(1) 
				print ",".join([sample, locus, str(score)])


if __name__ == "__main__":
	parser = ap.ArgumentParser()
	parser.add_argument("--mapper", default="bowtie")
	parser.add_argument("--unpaired", action="store_true", help="Fastq files are unpaired") 
	parser.add_argument("--ranges", default=None, help="position of alleles with reference sequences")
	parser.add_argument("--sig", type=float, default=0.05, help="minimum statistical significance threshold")
	parser.add_argument("--depth", type=int, default=0, help="minimum read depth to score variant")
	parser.add_argument("--pid", type=float, default=0, help="minimum % identity to score variant")
	parser.add_argument("--mismatches", type=int, default=float('inf'), help="maximum number of mismatches to score variant")
	parser.add_argument("--base", required=True, help="base name for output files")
	parser.add_argument("--var", required=True, help="name of variant fasta file")
	parser.add_argument("--fastq", required=True, nargs="+", help="list of (paired end) fastq files (fastq.gz)")
	parser.add_argument("--regex", default="^(.*)-(\d+)", help="Regex to extract the locus from the variant names")
	parser.add_argument("--procs", type=int, default=1, help="Number of processers to spawn for multiprocessing")
	args = parser.parse_args()

	main(args)
