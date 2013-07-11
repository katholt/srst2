# Replacement function for: 
# pileup_binomial_scoring(sam_file, size)
# Harriet Dashnow 28/6/13

pileup_file = "/vlsci/VR0082/hdashnow/SRST_github/srst2_etc_output/6581_3#8.all_fixed_annotated.fsa.srst2.pileup"
fasta = "/vlsci/VR0082/hdashnow/database_preprocessing/all_fixed_annotated.fsa"
prob_err = 0.01
edge_a = edge_z = 2

from itertools import groupby
from operator import itemgetter
import sys

def sequence_lengths_for_ref_alleles(fai_file):
    'Get sequence lengths for reference alleles - important for scoring'
    size = {}
    with open(fai_file) as fai:
        for line in fai:
            fields = line.split('\t')
            size[fields[0]] = int(fields[1])
    return size

def pileup_binomial_scoring(pileup_file, size):
    with open(pileup_file) as pileup:

        allele_line = 1
        exp_nuc_num = 1
        # XXX not actually set by user even though the comment in perl says it is
        prob_success = 1 - prob_err    # Set by user, default is prob_err = 0.01
        hash_alignment = {}
        hash_max_depth = {}
        hash_edge_depth = {}
        total_depth = 0
        depth_a = 0   # Gene 5' depth
        depth_z = 0   # Gene 3' depth
        max_depth = 1
        avg_depth_allele = {}
        coverage_allele = {}

        # Split all lines in the pileup by whitespace
        pileup_split = [ x.split() for x in pileup ]
        # Group the split lines based on the first field (allele) 
        for allele, lines in groupby(pileup_split, itemgetter(0)):

            # Reset variables for new allele
            allele_line = 1 # Keep track of line for this allele
            exp_nuc_num = 1 # Expected position in ref allele
            total_depth = 0
            depth_a = depth_z = 0
            #num_indel = 0
            #this_nuc_num = nuc_num
            #max_depth = 0
            #max_depth = nuc_depth
            hash_alignment[allele] = []
            allele_size = size[allele]
            total_indels = 0

            for fields in lines:
                nuc_num = int(fields[1]) # Actual position in ref allele
                exp_nuc_num += 1
                allele_line += 1
                nuc = fields[2]
                nuc_depth = int(fields[3])
                if len(fields) <= 5:
                    aligned_bases = ''
                else:
                    aligned_bases = fields[4]
                
                # Get info from this line:

                # Missing bases (alignment skips basepairs)
                if nuc_num > exp_nuc_num:
                    total_indels += abs(exp_nuc_num - nuc_num)
                exp_nuc_num = nuc_num

                # Calculate allele info: average depth -> penalise missing
                if exp_nuc_num <= edge_a:
                    depth_a += nuc_depth
                if abs(exp_nuc_num - allele_size) < edge_z:
                    depth_z += nuc_depth
                if nuc_depth > max_depth:
                    hash_max_depth[allele] = nuc_depth
                    max_depth = nuc_depth

                total_depth = total_depth + nuc_depth

                num_match = 0

                i = 0
                while i < len(aligned_bases):
                    if aligned_bases[i] == "^":
                        # Signifies start of a read, next char is mapping quality (skip it)
                        i += 2
                        continue

                    if aligned_bases[i] == "." or aligned_bases[i] == ",":
                        num_match += 1

                    i += 1

                num_mismatch = nuc_depth - num_match

                # Hash for later processing in R
                hash_alignment[allele].append((num_match, num_mismatch, prob_success))
 
            # Check for missing bases at the end of the allele
            if nuc_num < allele_size:
                total_indels += abs(allele_size - nuc_num)

            # Calculate allele summary stats and save
            avg_depth = total_depth / allele_line
            avg_a = depth_a / edge_a   # Avg depth at 5' end, num basepairs determined by edge_a
            avg_z = depth_z /edge_z    # 3'
            hash_max_depth[allele] = max_depth
            hash_edge_depth[allele] = (avg_a, avg_z)
            min_penalty = max(5, avg_depth)
            coverage_allele[allele] = 100*(allele_size - total_indels)/float(allele_size)

            # Penalize insertions/deletions and truncations 
            for j in range(total_indels):
                # Maintain a penalty of at least 5 mismatches if there's a gap
                # Save in hash for later processing in R
                hash_alignment[allele].append((0, min_penalty, prob_success))

            if allele in avg_depth_allele:
                avg_depth_allele[allele] += avg_depth
            else:
                avg_depth_allele[allele] = avg_depth

    return hash_alignment, hash_max_depth, hash_edge_depth, avg_depth_allele, coverage_allele

# Get sequence lengths for reference alleles - important for scoring
size = sequence_lengths_for_ref_alleles(fasta + '.fai')

coverage = pileup_binomial_scoring(pileup_file, size)[4]
for allele in coverage:
    if coverage[allele] != 100:
        print allele, coverage[allele]

