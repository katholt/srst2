# input = consensus alleles fasta files output by srst2 (*.all_consensus_alleles.fasta), containing >1 locus per file
# output = one fasta file of seqs per locus
# Author: Kat Holt (kholt@unimelb.edu.au)

# modules
import string, re, collections
import os, sys, subprocess
from argparse import (ArgumentParser, FileType)

# BioPython modules for reading and writing sequences
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

def parse_args():
	"Parse the input arguments, use '-h' for help."

	parser = ArgumentParser(description='SRST2 - Short Read Sequence Typer (v2)')

	parser.add_argument(
		'--input', nargs='+', type=str, required=False,
		help='Input files (should be *.all_consensus_alleles.fasta)')

	parser.add_argument('--pre', type=str, required=True, help='Prefix for output files')
	parser.add_argument('--type', type=str, required=True, help='either mlst or gene')
	parser.add_argument('--mlst_delimiter', type=str, required=False, default="-",
		help='For MLST seqs, specifcy the character(s) separating gene name from allele number in MLST database (default "-", as in arcc-1)')

	return parser.parse_args() 


def main():
	args = parse_args()
	
	# seqs = collections.defaultdict(dict) # key1 = locus, key2 = strain, value = sequence
	
	# read in seqs and print each new locus to a file
	for fasta_file in args.input:
		for record in SeqIO.parse(open(fasta_file, "r"), "fasta"):
			x = record.description # locus strain
			strain = x.split()[1]
			locus = x.split()[0]
			if args.type == "mlst":
				locus_id = locus.split(args.mlst_delimiter)[0]
				allele = locus
			elif args.type == "gene":
				locus_id = locus.split("__")[0] + "__" + locus.split("__")[1]
				allele = locus.split("__")[2] + "__" + locus.split("__")[3]
			out_file = args.pre + "__" + locus_id + ".fasta"
			o = file(out_file,"a") # open seq file for appending
			o.write(">" + strain + " " + locus + "\n")
			o.write(str(record.seq) + "\n")
			o.close()

if __name__ == '__main__':
	main()
