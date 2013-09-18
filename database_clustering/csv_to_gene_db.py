# take csv table detailing clustering etc and sequences for gene DB, write as fasta
# expected csv file format:
# seqID,ResFinderDBname,accession,clusterID,gene,allele,antibiotic_class,DNA,protein
# headers in output will be srst2 compatible, ie [clusterID]__[gene]__[allele]__[seqID] [accession]; [antibiotic_class]

# modules
import string, re, collections
import os, sys, subprocess
from optparse import OptionParser

# BioPython modules for reading and writing sequences
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)
	
	parser.add_option("-t", "--table_file", action="store", dest="table_file", help="table to read (csv)", default="")

	return parser.parse_args()

if __name__ == "__main__":

	(options, args) = main()
	
	if options.table_file == "":
		DoError("Please specify input table using -t")

	# read contents of a table and print as fasta
	f = file(options.table_file,"r")
	header = []
	for line in f:
		fields = line.rstrip().split(",")
		if len(header) > 0:
			unique_ID = fields[0]
			cluster = fields[3]
			gene = fields[4]
			allele = fields[5]
			id = "__".join([cluster,gene,allele,unique_ID]) ## this is the format for SRST2 detection and typing
			description = fields[1]
			if fields[2] != description:
				description += "; " + fields[2]
			description += "; " +  fields[6]
			seq = fields[7]
			record = SeqRecord(Seq(seq,
                   IUPAC.unambiguous_dna),
                   id=id, 
                   description=description)
			print record.format("fasta"),
		else:
			header = fields
	f.close()
	