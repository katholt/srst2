'''
Created on 19/06/2013
@author: Harriet Dashnow
Output a tab-separated table of alleles, their inferred gene name (from annotation) 
and their 90% cluster using CD-HIT
Also output fasta files containing genes with a given gene name
'''

import sys
from Bio import SeqIO

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
    cluster_file = "/home/UNIMELB/hdashnow/resistance_database/cdhit/all_fixed_cdhit90.clstr"
    infasta_file = "/home/UNIMELB/hdashnow/resistance_database/all_fixed.fsa"
    outfile = open("resistance_gene_DB.txt", "w")
    outfile.write("uniqueID\tResFinderDBname\tgene\tallele\taccession\tcluster90%\tcluster_contains_multiple_genes\tgene_found_in_multiple_clusters\n")
     
    database = {}
    full_database = {}
    
    for line in open(cluster_file):
        if line.startswith(">"):
            ClusterNr = line.split()[1]
            continue
 
        line_split =  line.split()
        #ClusterNr = line_split[0]
        QueryLabel = line_split[2].strip("*.>")
        
        if ClusterNr not in full_database:
            full_database[ClusterNr] = []
        if QueryLabel not in full_database[ClusterNr]:
            full_database[ClusterNr].append(QueryLabel)
        
        gene_allele = QueryLabel.split("_")[0]
        
# get allele             
        
        raw_gene_name = extract_allele_name(gene_allele)
        #print QueryLabel
        raw_allele = QueryLabel.split(raw_gene_name)[1].split("_")[1]
        accession = QueryLabel.split("_")[-1]
        
        if ClusterNr not in database:
            database[ClusterNr] = []
        if raw_gene_name not in database[ClusterNr]:
            gene_dict = {'ResFinderDB': QueryLabel,'gene':raw_gene_name, 'allele':raw_allele, 'accession':accession}
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

#        if len(database[cluster]) != 1:
#            multiple_genes[cluster] = database[cluster]
#            multiples_count += 1
#            clusters_with_multiple_genes.append(cluster)
        
        for gene in database[cluster]:
            #gene_ID = gene['ResFinderDB']
            gene_name = gene['gene']
            
#            if gene_name in all_genes:
#                span_genes.append(gene_ID)
#                       
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
            
            ResFinderDB = gene['ResFinderDB']
            if ResFinderDB not in genes_printed_to_file:
                uniqueID += 1
                gene_name = gene['gene']
                allele = gene['allele']
                accession = gene['accession']
                if cluster in clusters_with_multiple_genes:
                    cluster_flag = "yes"
                else:
                    cluster_flag = "no"
                if gene_name in span_genes:
                    gene_flag = "yes"
                else:
                    gene_flag = "no"
                
                outstring = "{0:05d}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(uniqueID, ResFinderDB, 
                            gene_name, allele, accession, cluster, cluster_flag, gene_flag)
                outfile.write(outstring)
                genes_printed_to_file.append(ResFinderDB)

    # Save all alleles with the same gene name to separate fasta files
    records_by_gene = {}
    for seq_record in SeqIO.parse(infasta_file, "fasta"):
        gene = extract_allele_name(seq_record.id.split("_")[0])
        if gene not in records_by_gene:
            records_by_gene[gene] = []
        records_by_gene[gene].append(seq_record)

    for gene in records_by_gene:
        SeqIO.write(records_by_gene[gene], (gene + ".fsa"), "fasta")

if __name__ == '__main__':
    sys.exit(main())
