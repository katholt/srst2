# take csv table detailing clustering etc and sequences for gene DB, write as fasta
# expected csv file format:
# seqID,clusterid,gene,allele,(DNAseq),other....
# headers in output will be srst2 compatible, ie [clusterID]__[gene]__[allele]__[seqID] [other stuff]
# sequence can be read from a specified column or from a fasta file (specify which column contains fasta header to match in seqs file)

# Author: Kat Holt (kholt@unimelb.edu.au)

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
	
	parser.add_option("-t", "--table", action="store", dest="table_file", help="table to read (csv)", default="")
	parser.add_option("-o", "--out", action="store", dest="output_file", help="output file (fasta)", default="")
	parser.add_option("-s", "--seq_col", action="store", dest="seq_col", help="column number containing sequences", default="")
	parser.add_option("-f", "--fasta", action="store", dest="fasta_file", help="fasta file to read sequences from (must specify which column in the table contains the sequence names that match the fasta file headers)", default="")	
	parser.add_option("-c", "--headers_col", action="store", dest="headers_col", help="column number that contains the sequence names that match the fasta file headers", default="")	

	return parser.parse_args()

if __name__ == "__main__":

	(options, args) = main()
	
	seqid_col = False
	seqs_file_col = False
	
	input_seqs = {}
		
	if options.table_file == "":
		DoError("Please specify input table using -t")
	if options.output_file == "":
		DoError("Please specify output fasta file using -o")
	if options.seq_col != "":
		print "Reading DNA sequences from table, column" + options.seq_col
		seqid_col = int(options.seq_col)
	elif options.fasta_file != "":
		if options.headers_col == "":
			DoError("Please specify which column of the table contains identifiers that match the headers in the fasta file")
		seqs_file_col = int(options.headers_col)
		print "Reading DNA sequences from fasta file: " + options.fasta_file
		for record in SeqIO.parse(open(options.fasta_file, "r"), "fasta"):
			input_seqs[record.id] = record.seq
			
	else:
		print DoError("Where are the sequences? If they are in the table, specify which column using -s. Otherwise provide a fasta file of sequence using -f and specify which column contains sequence identifiers that match the fasta headers, using -h")

	# read contents of a table and print as fasta
	f = file(options.table_file,"r")
	o = open(options.output_file, "w")
	header = []
	for line in f:
		fields = line.rstrip().split(",")
		if len(header) > 0:
			seqID = fields[0]
			cluster = fields[1]
			gene = fields[2]
			allele = fields[3]
			db_id = "__".join([cluster,gene,allele,seqID]) ## this is the format for SRST2 detection and typing

			if seqid_col:
				seq = fields.pop(seqid_col-1)
				record = SeqRecord(Seq(seq,
					   IUPAC.unambiguous_dna),
					   id=db_id, 
					   description=db_id)
			elif seqs_file_col:
				seqs_file_id = fields.pop(seqs_file_col-1)
				if seqs_file_id in input_seqs:
					record = SeqRecord(input_seqs[seqs_file_id],id=db_id, description=db_id)
				else:
					print "Warning, couldn't find a sequence in the fasta file matching this id: " + seqs_file_id
				
			else:
				"??"
				
			# add annotation from other columns
			if len(fields) > 4:
				description = ";".join(fields[4:len(fields)])
				record.description = description

			count = SeqIO.write(record, o, "fasta")
			
		else:
			header = fields
			
	f.close()
	o.close()
	
	