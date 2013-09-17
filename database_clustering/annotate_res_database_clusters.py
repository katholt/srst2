'''
Created on 28/06/2013
@author: Harriet Dashnow
Output a tab-separated table of alleles, their inferred gene name (from annotation) 
and their 90% cluster using CD-HIT
Also output fasta files containing genes with a given gene name
'''

import sys
from Bio import SeqIO

def main():

    cluster_info_file = "resistance_gene_DB_annotated_final.csv"
    database_fasta_file = "/home/UNIMELB/hdashnow/resistance_database/all_fixed.fsa"
    out_fasta_file = "all_fixed_annotated.fsa"

    # Read in the cluster file info for each allele: cluster # and cluster annotation
    cluster_info = {}
    with open(cluster_info_file) as clust_info:
        for f in clust_info:
            for line in f.split()[1:]:
                fields = line.split(",")
                allele = fields[1]
                clust_number = fields[5]
                clust_annotation = fields[8]
                unique_id = fields[0]
                if allele not in cluster_info:
                    cluster_info[allele] = (clust_number, clust_annotation)
                else:
                    pass #report duplicates?
                    
    # Erase output file if it already exists, then open it for appending
    with open(out_fasta_file, "w") as out_fasta:
        pass
    with open(out_fasta_file, "a") as out_fasta:

        # Read in fasta database, annotate then write out to a output file
        with open(database_fasta_file, "r") as database_fasta:

            # For each entry in fasta database, write it out to a new file with the addition of cluster#__cluster_annotation__...
            for record in SeqIO.parse(database_fasta, "fasta"):

                # Look up relevent cluster info
                allele = record.id
                clust_number = cluster_info[allele][0]
                clust_annotation = cluster_info[allele][1]
                
                # Edit the ID with cluster info
                record.id = record.name = record.description =  "{0}__{1}__{2}".format(clust_number, clust_annotation, allele)
                print record.id
                sys.exit()

                # Write to a new file 
                SeqIO.write(record, out_fasta, "fasta")

if __name__ == '__main__':
    sys.exit(main())
