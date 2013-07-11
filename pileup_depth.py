#!/usr/bin/env python

# SRST2 - Short Read Sequence Typer (v2)
# 
# Author - Michael Inouye (minouye@unimelb.edu.au)
# Translated to Python by Bernie Pope (bjpope@unimelb.edu.au)
#
# Dependencies:
#    bowtie2       http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
#    SAMtools      http://samtools.sourceforge.net
#    R             http://www.r-project.org 

from argparse import (ArgumentParser, FileType)
import logging
from subprocess import call
import os
import sys
from scipy.stats import binom_test, linregress
from math import log

edge_a = edge_z = 2
prob_err = 0.01

#pileup file (not sure why it's called out_file_sam3)
out_file_sam3 = "/vlsci/VR0082/hdashnow/SRST_github/srst2_etc_output/6581_3#8.all_fixed_annotated.fsa.srst2.pileup"
size = 972

def pileup_binomial_scoring(sam_file, size):
    # Use pileup for binomial-based scoring
    with open(sam_file) as pileup:
        count = 1
	line_count = 1
        # XXX not actually set by user even though the comment in perl says it is
        prob_success = 1 - prob_err    # Set by user, default is prob_err = 0.01
        hash_alignment = {}
        hash_max_depth = {}
        hash_edge_depth = {}
        total_depth = 0
        depth_a = 0   # Gene 5' depth
        depth_z = 0   # Gene 3' depth
        num_indel = 0
        this_allele = None
        max_depth = 1
        avg_depth_allele = {}

        for line in pileup:
	    #print line
            fields = line.split()
            allele = fields[0]
            nuc_num = int(fields[1])
            nuc = fields[2]
            nuc_depth = int(fields[3])
            allele_size = size
            if len(fields) <= 5:
                aligned_bases = ''
            else:
                aligned_bases = fields[4]
                
            if count > 1 and allele != this_allele:
                # Not the first line and allele alignment has changed in pileup
                ### Check depth here, print total_depth, count
		#print total_depth, count
		avg_depth = total_depth / float(line_count)
		#print avg_depth
                avg_a = depth_a / edge_a   # Avg depth at 5' end, num basepairs determined by edge_a
                avg_z = depth_z /edge_z    # 3'
                hash_max_depth[this_allele] = max_depth
                hash_edge_depth[this_allele] = (avg_a, avg_z)
                        
                min_penalty = max(5, avg_depth)

		#HD#print line
		#print total_depth, count, avg_depth
                # Penalize insertions/deletions and truncations 
                num_missing = abs(this_allele_size - (this_nuc_num - 1)) + num_indel

                for j in range(num_missing):
                    # Maintain a penalty of at least 5 mismatches if there's a gap
                    # Save in hash for later processing in R
                    hash_alignment[this_allele].append((0, min_penalty, prob_success))

                if allele in avg_depth_allele:
                    avg_depth_allele[this_allele] += avg_depth
                else:
                    avg_depth_allele[this_allele] = avg_depth
		### Check above that this is saved correctly (allele is unique?)                        
		#all_alleles_depth = []
		#for allele in avg_depth_allele:
		#	all_alleles_depth.append(avg_depth_allele[allele])
		#print all_alleles_depth
                # Reset counters and indicators
                count = 1
                line_count = 1
		total_depth = 0
                depth_a = depth_z = 0
                num_indel = 0
                this_allele = allele
                this_nuc_num = nuc_num
                max_depth = 0
                
            if count == 1:
                # First line of file or of new allele
                this_allele = allele
                this_nuc_num = nuc_num
                this_allele_size = allele_size
                max_depth = nuc_depth
                hash_alignment[this_allele] = []

            # XXX assumes nuc_num > count, what if count > nuc_num?
            # could that even happen? Check the pileup format.
            if allele == this_allele and nuc_num != count:
                # Same allele but alignment skips basepairs
                num_indel += abs(count - nuc_num)
                count = nuc_num     ### Check that this is the right way to index count
                this_nuc_num = nuc_num
                # add dummy entries for the missing nucleotide positions
                hash_alignment[allele] += [None] * num_indel

            if allele == this_allele:
                # Same allele, calculate depths and update counters
                if count > 0 and count < edge_a + 1:
                    depth_a = depth_a + nuc_depth
                if abs(count - allele_size) < edge_z:
                    depth_z = depth_z + nuc_depth

                if nuc_depth > max_depth:
                    hash_max_depth[allele] = nuc_depth
                    max_depth = nuc_depth

                total_depth = total_depth + nuc_depth       ### Check here
                this_nuc_num += 1
                count += 1 
		line_count += 1

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
	print avg_depth_allele["162__aadA__aadA1_1_X02340"]
        return hash_alignment, hash_max_depth, hash_edge_depth, avg_depth_allele





pileup_binomial_scoring(out_file_sam3, size)

