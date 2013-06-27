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


def parse_args():
    "Parse the input arguments, use '-h' for help"

    parser = ArgumentParser(description='SRST2 - Short Read Sequence Typer (v2)')
    parser.add_argument(
        '--input_se', type=str, required=False,
        help='Input single end read')
    parser.add_argument(
        '--input_pe', nargs=2, type=str, required=False,
        help='Input paired end reads')
    parser.add_argument('--output', type=str, required=True, help='Output file')
    parser.add_argument('--build', type=str, required=True, nargs='+',
        help='Fasta file(s) to build bowtie2 indices')
    parser.add_argument(
        '--verbose', type=int, default=0, required=False, help='Verbosity')
    parser.add_argument(
        '--input_type', type=str, choices=['q', 'qseq', 'f'], default='q',
        help='Input file type.')
    parser.add_argument( '--other', type=str, help='Other options for bowtie2.') 
    parser.add_argument( '--mapq', type=int, default=1, help='Qual map') 
    parser.add_argument( '--baseq', type=int, default=20, help='Qual base') 
    parser.add_argument( '--log', metavar='FILE', type=str,
        help='log progress in FILENAME, defaults to stdout')
    return parser.parse_args() 


def run_command(command, **kwargs):
    'Execute a shell command and check the exit status and any O/S exceptions'
    command_str = ' '.join(command)
    logging.info('Running: {}'.format(command_str))
    try:
        exit_status = call(command, **kwargs)
    except OSError as e:
        exit("Command '{}' failed due to O/S error: {}".format(command_str, str(e)))
    if exit_status != 0:
        exit("Command '{}' failed with non-zero exit status: {}".format(command_str, exit_status))


def bowtie_index(fasta_files):
    'Build a bowtie2 index from the given input(s)'

    for fasta in fasta_files:
        built_index = fasta + '.1.bt2'
        if os.path.exists(built_index):
            logging.info('Index for {} is already built...'.format(fasta))
        else:
            logging.info('Building bowtie2 index for {}...'.format(fasta))
            run_command(['bowtie2-build', fasta, fasta])


def modify_bowtie_sam(output_file):
    with open(output_file) as sam, open(output_file + '.mod', 'w') as sam_mod:
        # XXX test parsing better
        for line in sam:
            if not line.startswith('@'):
                fields = line.split('\t')
                flag = int(fields[1])
                flag = (flag - 256) if (flag & 256) else flag
                sam_mod.write('\t'.join([fields[0], str(flag)] + fields[2:]))
            else:
                sam_mod.write(line)


def sequence_lengths_for_ref_alleles(fai_file):
    'Get sequence lengths for reference alleles - important for scoring'
    size = {}
    with open(fai_file) as fai:
        for line in fai:
            fields = line.split('\t')
            size[fields[0]] = int(fields[1])
    return size


def pileup_binomial_scoring(sam_file, size):
    # Use pileup for binomial-based scoring
    with open(sam_file) as pileup:
        count = 1
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
            fields = line.split()
            allele = fields[0]
            nuc_num = int(fields[1])
            nuc = fields[2]
            nuc_depth = int(fields[3])
            allele_size = size[fields[0]]
            if len(fields) <= 5:
                aligned_bases = ''
            else:
                aligned_bases = fields[4]
                
            if count > 1 and allele != this_allele:
                # Not the first line and allele alignment has changed in pileup
                ### Check depth here, print total_depth, count
		avg_depth = total_depth / count
                if this_allele == "blaCTX-M_blaCTX-M-71_1_FJ815436":
                        print "line 132:", total_depth, count, avg_depth
                avg_a = depth_a / edge_a   # Avg depth at 5' end, num basepairs determined by edge_a
                avg_z = depth_z /edge_z    # 3'
                hash_max_depth[this_allele] = max_depth
                hash_edge_depth[this_allele] = (avg_a, avg_z)
                        
                min_penalty = max(5, avg_depth)

		#HD#print line
		#HD#print total_depth, count, avg_depth
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
		if allele == "blaCTX-M_blaCTX-M-71_1_FJ815436":
			print avg_depth
                # Reset counters and indicators
                count = 1
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

                if allele == "blaCTX-M_blaCTX-M-71_1_FJ815436":
                        print total_depth, count, total_depth/count

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

        print avg_depth_allele["blaCTX-M_blaCTX-M-71_1_FJ815436"]
        return hash_alignment, hash_max_depth, hash_edge_depth, avg_depth_allele


def score_alleles(out_file_sam3, hash_alignment, hash_max_depth, hash_edge_depth, avg_depth_allele):
    # Table for R processing
    with open(out_file_sam3 + '.table.scores', 'w') as scores:
        scores.write("Allele\tScore\tAvg_depth\tEdge1_depth\tEdge2_depth\n")
        for allele in hash_alignment:
            pvals = []
            for nuc_info in hash_alignment[allele]:
                if nuc_info is not None:
                    match, mismatch, prob_success = nuc_info
                    if match > 0 or mismatch > 0:
                        if mismatch == 0:
                            p_value = 1.0
                        else:
                            p_value = binom_test([match, mismatch], None, prob_success)
                        # Weight pvalue by (depth/max_depth)
                        max_depth = hash_max_depth[allele]
                        weight = (match + mismatch) / float(max_depth)
                        p_value *= weight
			if p_value == 0:
			    p_value = 0.000000000000000000000000000001 ### Was getting a value error when p_value = 0.0
                        p_value = -log(p_value, 10)
                        pvals.append(p_value)
            # Fit linear model to observed Pval distribution vs expected Pval distribution (QQ plot)
            pvals.sort(reverse=True)
            len_obs_pvals = len(pvals)
            exp_pvals = range(1, len_obs_pvals + 1)
            exp_pvals2 = [-log(float(ep) / (len_obs_pvals + 1), 10) for ep in exp_pvals]
            # Slope is score
            slope, _intercept, _r_value, _p_value, _std_err = linregress(exp_pvals2, pvals)
            if allele in hash_edge_depth:
                start_depth, end_depth = hash_edge_depth[allele]
                edge_depth_str = str(start_depth) + '\t' + str(end_depth)
            else:
                edge_depth_str = "NA\tNA"
            this_depth = avg_depth_allele.get(allele, "NA")
            scores.write('\t'.join([allele, str(slope), str(this_depth), edge_depth_str]) + '\n')
            if allele == "blaCTX-M_blaCTX-M-71_1_FJ815436":
            	print allele, this_depth
		print avg_depth_allele[allele]
		sys.exit()

def run_bowtie_on_indices(args):
    'Run bowtie2 on the newly built index/indices'

    for fasta in args.build:
        path, fasta_name = os.path.split(fasta) 
        out_file = args.output + '.' + fasta_name + '.srst2'
        logging.info('Output prefix set to: ' + out_file)
        command = ['bowtie2']

        if args.input_se:
            # single end
            command += ['-U', args.input_se]
        elif args.input_pe:
            # paired end
            command += ['-1', args.input_pe[0], '-2', args.input_pe[1]]
        else:
            # XXX should fix this error message and possibly test it earlier
            exit('One of input_se or input_pe must be defined')

        command += ['-S', out_file,
                    '-' + args.input_type,    # add a dash to the front of the option
                    '--very-sensitive-local',
                    '--no-unal',
                    '-a',                     # Search for and report all alignments
                    '-x', fasta               # The index to be aligned to
                   ]

        if args.other:
            command += args.other.split()

        logging.info('Aligning reads to index {} using bowtie2...'.format(fasta))
#        run_command(command)

        # Modify Bowtie's SAM formatted output so that we get secondary
        # alignments in downstream pileup

#        modify_bowtie_sam(out_file)

        # Analyse output with SAMtools
        logging.info('Processing Bowtie2 output with SAMtools...')
        logging.info('Generate and sort BAM file...')
        out_file_sam1 = out_file + ".bam"
#        run_command(['samtools', 'view', '-b', '-o', out_file_sam1,
#                     '-q', str(args.mapq), '-S', out_file + '.mod'])
        out_file_sam2 = out_file + ".sorted"
#        run_command(['samtools', 'sort', out_file_sam1, out_file_sam2])
#        run_command(['samtools', 'faidx', fasta])

        logging.info('Generate pileup...')
        out_file_sam3 = out_file + '.pileup'
#        with open(out_file_sam3, 'w') as sam_pileup:
#            run_command(['samtools', 'mpileup', '-L', '250', '-f', fasta,
#                         '-Q', str(args.baseq), '-q', str(args.mapq), out_file_sam2 + '.bam'],
#                         stdout=sam_pileup)

        # Process SAMtools output
        logging.info('Processing SAMtools output...')

        # Get sequence lengths for reference alleles - important for scoring
        size = sequence_lengths_for_ref_alleles(fasta + '.fai')

        hash_alignment, hash_max_depth, hash_edge_depth, avg_depth_allele = \
            pileup_binomial_scoring(out_file_sam3, size)

        logging.info('Scoring alleles...')
        score_alleles(out_file_sam3, hash_alignment, hash_max_depth, hash_edge_depth, avg_depth_allele)
        logging.info('Finished processing for index file {} ...'.format(fasta))


def main():
    args = parse_args()
    if args.log is None:
        logfile = sys.stdout
    else:
        logfile = args.log
    logging.basicConfig(
        filename=args.log,
        level=logging.DEBUG,
        filemode='w',
        format='%(asctime)s %(message)s',
        datefmt='%m/%d/%Y %H:%M:%S')
    logging.info('program started')
    logging.info('command line: {0}'.format(' '.join(sys.argv)))
    bowtie_index(args.build)
    run_bowtie_on_indices(args)
    logging.info('SRST2 has finished.')


if __name__ == '__main__':
    main()
