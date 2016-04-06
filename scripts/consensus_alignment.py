'''
This script combines allele sequences in every sample into a FASTA file per locus, organising consensus sequences according to loci.

Input:
	FASTA files of consensus allele sequences produced by SRST2 (*.all_consensus_alleles.fasta or *.new_consensus_alleles.fasta), containing >1 locus per file
	Note that it is recommended to process *.all_consensus_alleles.fasta files because variants of MLST alleles are not recorded in *.new_consensus_alleles.fasta files.
Output: one FASTA file of seqs per locus

Notes about the format of FASTA definition lines in FASTA files produced by SRST2:
	* for an MLST allele: >[gene name]_[allele number].consensus [sample name]
	* for an allele of other genes: >[cluster number]__[gene name]_[product/phenotype]__[allele name]__[internal sequence number].consensus [sample name]
	* alleles of both MLST genes and the other kind of genes will be present in the same FASTA file if SRST2 is used for MLST and genetic profiling simultaneously.

Notes about arguments:
	The argument "type" specifies types of alleles present in input FASTA files. The program will try to parse descriptions of alleles of both MLST genes
	and the other kind of genes if this argument is left blank. You may use "mlst" if all FASTA files only consist of allele sequences of MLST genes, and use
	"gene" if these files only consist of the other kind of allele sequences. The program does work successfully if you use a wrong value for this argument.
	
Author: Kathryn Holt (kholt@unimelb.edu.au), Yu Wan (wanyuac@gmail.com)
'''

# modules
import string, os, sys
from argparse import ArgumentParser

# BioPython modules for reading and writing sequences
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_args(): # Parse the input arguments, use "-h" for help.
	parser = ArgumentParser(description = "SRST2 - Short Read Sequence Typer (v2)")
	parser.add_argument("--input", nargs = "+", type = str, required = False, help = "Input files (recommended to be *.all_consensus_alleles.fasta)")
	parser.add_argument("--pre", type = str, required = True, help = "Prefix for output files")
	parser.add_argument("--type", type = str, default = "", required = False, help = "mlst, gene, or leave it as a blank (default)")
	parser.add_argument("--mlst_delimiter", type = str, required = False, default = "-",
		help = "For MLST seqs, specifcy the character(s) separating gene name from allele number in MLST database (default \"-\", as in arcc-1)")
	parser.add_argument("--ext", type = str, default = "fasta", required = False, help = "The file name extension. Default: fasta")
	return parser.parse_args()

def parse_descr_mlst(locus, delimiter):
	# Now, locus = [gene name]_[allele number]
	locus_id = locus.split(delimiter)[0]  # For example, gapA_1 becomes gapA and 1.
	allele = locus
	return locus_id, allele

def parse_descr_gene(locus):
	# locus = [cluster number]__[gene name]_[product/phenotype]__[allele name]__[internal sequence number]
	fields = locus.split("__")
	locus_id = fields[0] + "__" + fields[1]  # locus_id = [cluster number]__[gene name]_[product/phenotype]
	allele = fields[2] + "__" + fields[3]  # allele = [allele name]__[internal sequence number]
	return locus_id, allele

def main():
	args = parse_args()
	'''
	read in every sequence and print each new locus to a file
	Otherwise, allele sequences are appended to corresponding extant files.
	'''
	for fasta_file in args.input:
		for record in SeqIO.parse(open(fasta_file, "r"), "fasta"):
			locus, sample = record.description.split(" ") # locus sample
			locus = locus.split(".")[0]  # remove the suffix ".consensus"
			
			if args.type == "mlst":
				locus_id, allele = parse_descr_mlst(locus, args.mlst_delimiter)
			elif args.type == "gene":
				locus_id, allele = parse_descr_gene(locus)
			else:
				if locus[0].isdigit():  # The description of an allele other than MLST alleles always starts with a digit, which is the cluster number.
					locus_id, allele = parse_descr_gene(locus)
				else:
					locus_id, allele = parse_descr_mlst(locus, args.mlst_delimiter)
			
			out_file = args.pre + "__" + locus_id + "." + args.ext
			o = file(out_file,"a") # open a seq file for appending
			o.write(">" + sample + " " + locus_id + " " + allele + "\n")
			o.write(str(record.seq) + "\n")
			o.close()

if __name__ == "__main__":
	main()
