'''
Created on 19/06/2013
@author: Harriet Dashnow, Kat Holt
Output a tab-separated table of alleles, their inferred gene name (from annotation) 
and their 90% cluster using CD-HIT
Also output fasta files containing all genes sharing the same basic name (prefix before '-' in gene id)
'''

import sys
from Bio import SeqIO
from argparse import ArgumentParser


def parse_args():
	parser = ArgumentParser(description='Parse cdhit results')

	parser.add_argument('--cluster_file',
						required = True,
						help = 'cd hit output file (.clstr)')
	parser.add_argument('--infasta_file',
						required = True,
						help = 'raw sequences file (fasta)')
	parser.add_argument('--outfile',
						required = True,
						help = 'output file (csv)')
	return parser.parse_args()

def remove_trailing_numbers(string):
	index = len(string) - 1	   
	while index >= 0:
		if string[index].isdigit():
			string = string[:index]
		else:
			break
		index -= 1
	return string

def extract_allele_name(gene_allele):
	gene_name = ""
	i = 0
	for part in gene_allele.split("-"):
		
		if not part.isdigit():
			if i !=0:
				gene_name += "-"
			gene_name += part
		i += 1
	
	return remove_trailing_numbers(gene_name)	  

def main():
	args = parse_args()
	outfile = file(args.outfile,"w")
	outfile.write("seqID,clusterid,gene,allele,cluster_contains_multiple_genes,gene_found_in_multiple_clusters\n")
	 
	database = {}
	full_database = {}
	
	for line in open(args.cluster_file):
		if line.startswith(">"):
			ClusterNr = line.split()[1]
			continue
 
		line_split =  line.split()
		QueryLabel = line_split[2].strip("*.>")
		
		if ClusterNr not in full_database:
			full_database[ClusterNr] = []
		if QueryLabel not in full_database[ClusterNr]:
			full_database[ClusterNr].append(QueryLabel)
		
		cluster_label = QueryLabel.split("-")[0]
		
		if ClusterNr not in database:
			database[ClusterNr] = []
		if cluster_label not in database[ClusterNr]:
			gene_dict = {'gene':cluster_label, 'allele':QueryLabel}
			database[ClusterNr].append(gene_dict)
	
	max_cluster_size = 0
	largest_cluster = ""
	for cluster in full_database:
		cluster_size = len(full_database[cluster])
		if cluster_size > max_cluster_size:
			max_cluster_size = cluster_size
			largest_cluster = cluster
	
	all_genes = []
	multiple_genes = {}
	multiples_count = 0
	total_count = 0

	span_genes = []
	clusters_with_multiple_genes = []
		
	for cluster in database:
		total_count += 1
		cluster_genes = []
		
		for gene in database[cluster]:
			gene_name = gene['gene']
			
			if gene_name not in multiple_genes:
				multiple_genes[gene_name] = []
			multiple_genes[gene_name].append(cluster)
			
			if gene_name not in cluster_genes:
				cluster_genes.append(gene_name)
			
		if len(cluster_genes) > 1:
			multiples_count += 1
			clusters_with_multiple_genes.append(cluster)
		
	
	for gene in multiple_genes:
		unique_clusters = list(set(multiple_genes[gene]))
		if len(unique_clusters) > 1:
			span_genes.append(gene)
	
	genes_printed_to_file = []
	uniqueID = 0	
	for cluster in database:
		
		for gene in database[cluster]:
			
			if gene not in genes_printed_to_file:
				uniqueID += 1
				gene_name = gene['gene']
				allele = gene['allele']
				if cluster in clusters_with_multiple_genes:
					cluster_flag = "yes"
				else:
					cluster_flag = "no"
				if gene_name in span_genes:
					gene_flag = "yes"
				else:
					gene_flag = "no"
				
				outstring = "{0:05d},{1},{2},{3},{4},{5}\n".format(uniqueID, cluster, gene_name, allele, cluster_flag, gene_flag)
				outfile.write(outstring)
				genes_printed_to_file.append(gene)

	outfile.close()
	
	# Save all alleles with the same gene name to separate fasta files
	records_by_gene = {}
	for seq_record in SeqIO.parse(args.infasta_file, "fasta"):
		gene = seq_record.id.split("-")[0]
		if gene not in records_by_gene:
			records_by_gene[gene] = []
		records_by_gene[gene].append(seq_record)

	for gene in records_by_gene:
		SeqIO.write(records_by_gene[gene], (gene + ".fsa"), "fasta")
		


if __name__ == '__main__':
	sys.exit(main())
