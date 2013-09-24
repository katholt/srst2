'''
Extract virulence genes by genus from the VFDB database at http://www.mgc.ac.cn/VFs/Down/CP_VFs.ffn.gz
'''

import sys, re
from argparse import ArgumentParser

# BioPython modules for reading and writing sequences
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

def parse_args():
	parser = ArgumentParser(description='Extract virulence genes by genus from the VFDB database available at http://www.mgc.ac.cn/VFs/Down/CP_VFs.ffn.gz')

	parser.add_argument('--infile',
						required = True,
						help = 'Raw VFDB sequences file (fasta, e.g. download from http://www.mgc.ac.cn/VFs/Down/CP_VFs.ffn.gz)')
	parser.add_argument('--genus',
						required = False,
						help = 'Genus to extract (if not specified, all genera will be extracted to individual files)')
	return parser.parse_args()  

def main():
	args = parse_args()
	
	db = {} # key = genus, value = list of sequences
	
	for record in SeqIO.parse(open(args.infile, "r"), "fasta"):
		full_name = record.description
		genus = full_name.split("[")[-1].split()[0]
		if (not args.genus) or (genus == args.genus):
			if genus in db:
				db[genus].append(record)
			else:
				db[genus] = [record]

	# Save all alleles from the same genus to separate fasta files
	for genus in db:
		records = db[genus] # list of records
		SeqIO.write(records, (genus + ".fsa"), "fasta")

if __name__ == '__main__':
	sys.exit(main())
