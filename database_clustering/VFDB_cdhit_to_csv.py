'''
Created on 19/06/2013
@author: Harriet Dashnow, Kat Holt
Output a comma-separated table of alleles, their inferred gene name (from annotation) 
and their 90% cluster from CD-HIT analysis
'''

import sys, re
from Bio import SeqIO
from argparse import ArgumentParser


def parse_args():
	parser = ArgumentParser(description='Parse cdhit results')

	parser.add_argument('--cluster_file',
						required = True,
						help = 'cd hit output file (.clstr)')
	parser.add_argument('--infile',
						required = True,
						help = 'raw sequences file that was input to cdhit (fasta)')
	parser.add_argument('--outfile',
						required = True,
						help = 'output file (csv)')
	return parser.parse_args()  

def main():
	args = parse_args()
	outfile = file(args.outfile,"w")
	outfile.write("seqID,clusterid,gene,allele,DNA,annotation\n")
	 
	database = {} # key = clusterid, value = list of seqIDs
	seq2cluster = {} # key = seqID, value = clusterid
	
	for line in open(args.cluster_file):
		if line.startswith(">"):
			ClusterNr = line.split()[1]
			continue
 
		line_split =  line.split()
		seqID = line_split[2].strip("*.>")
		
		if ClusterNr not in database:
			database[ClusterNr] = []
		if seqID not in database[ClusterNr]:
			database[ClusterNr].append(seqID) # for virulence gene DB, this is the unique ID R0xxx
		seq2cluster[seqID] = ClusterNr
	
	for record in SeqIO.parse(open(args.infile, "r"), "fasta"):
		full_name = record.description
		genus = full_name.split("[")[1].split()[0]
		id_bits = re.sub("[()]","",full_name.split("-")[0]).split() # 'R004852 fliL VP2243 '
		seqID = id_bits[0] # R004852
		gene = id_bits[1] # fliL
		if len(id_bits) > 2:
			allele = id_bits[1]+"_"+id_bits[2] # fliL_VP2243
		else:
			allele = id_bits[1]
		clusterid = seq2cluster[seqID]
		outstring = ",".join([seqID, clusterid, gene, allele, str(record.seq), re.sub(",","",record.description)]) + "\n"
		outfile.write(outstring)
	outfile.close()


if __name__ == '__main__':
	sys.exit(main())
