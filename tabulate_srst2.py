#!/usr/bin/env python

import sys
import os.path

# e.g. shigella-samples.csv
st_true = {}
with open(sys.argv[1]) as true_st:
    for line in true_st:
        fields = line.split(',')
        samplename = fields[0]
        this_st = fields[1]
        if this_st in ["", "?", "unknown"] or 'SLV' in this_st or 'DLV' in this_st:
            this_st = "NA"
        st_true[samplename] = this_st

# e.g. ../data/pub_mlst/Shigella-sonnei/ecoli_ST_table.txt
# XXX we should really use the CSV library for this
st_map2 = {}
with open(sys.argv[2]) as st_file:
    for count, line in enumerate(st_file):
        if count == 0:
            gene_order = line.upper().split('\t')
            # drop the first column of the heading
            gene_order = gene_order[1:]
        else:
            fields = line.split('\t')
            st_map2[fields[0]] = fields[1:] 


class AlleleInfo(object):
    def __init__(self, avg_depth, depth_a, depth_z, score):
        self.avg_depth = avg_depth
        self.depth_a = depth_a
        self.depth_z = depth_z
        self.score = score

hash = {}
with open(sys.argv[3]) as list:
    for line in list:
        line = line.strip()
        with open(line) as file:
            # chop off any path prefix in the filename
            _path, filename = os.path.split(line) 
            # chop of any extension after the dot
            filenameparts = filename.split('.')
            sample = filenameparts[0]
            allele = filenameparts[1].upper()

            #skip the first line
            next(file)
            for line in file:
                score = avg_depth = depth_a = depth_z = "NA"
                allelenum, score, avg_depth, depth_a, depth_z = line.strip().split('\t')
                score = float(score)
                if depth_a == '': depth_a = "NA"
                if depth_z == '': depth_z = "NA"

                info = AlleleInfo(avg_depth, depth_a, depth_z, score)

                if sample in hash:
                    if allele in hash[sample]:
                        hash[sample][allele][allelenum] = info 
                    else:
                        hash[sample][allele] = { allelenum : info }
                else:
                    hash[sample] = { allele : { allelenum : info }}


print("Sample\tGene\tAllele_SRST2\tAllele_True\tScore\tAvg_depth\tEdge1_depth\tEdge2_depth\tConcordant_alleles")

for sample in hash:
    this_true_hash = {}
    for count, gene in enumerate(gene_order):
        if sample in st_true:
            this_true_st = st_true[sample]
            this_true_alleles = st_map2[this_true_st]
            this_true_hash[gene] = this_true_alleles[count]
        else:
            sys.stderr.write("No 'true' ST for $sample!")
            this_true_hash[gene] = "NA"

    for allele in hash[sample]:
        for count, allelenum in enumerate(sorted(hash[sample][allele])): 
            info = hash[sample][allele][allelenum]
            if count == 0:
                best_score = info.score

            if count > 0 and best_score < info.score:
                continue 

            if allele in this_true_hash:
                a, b = allelenum.split('-')
                if b == this_true_hash[allele]:
                    concord = "1"
                elif this_true_hash[allele] == "NA":
                    concord = "NA"
                else:
                    concord = "0"
                print('\t'.join([sample, allele, allelenum, this_true_hash[allele],
                                str(info.score), info.avg_depth, info.depth_a, info.depth_z, concord]))
            else:
                print('\t'.join([sample, allele, allelenum, "NA",
                                str(info.score), info.avg_depth, info.depth_a, info.depth_z, "NA"]))
