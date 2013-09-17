#!/usr/bin/env python
import string, re, collections
import os, sys, subprocess		
from subprocess import call
from argparse import (ArgumentParser, FileType)

def parse_args():
	"Parse the input arguments, use '-h' for help"

	parser = ArgumentParser(description='Submit SRST2 jobs through SLURM')

	# Job details
	parser.add_argument(
		'--walltime', type=str, required=False, help='wall time (default 0-1:0 = 1 h)', default="0-1:0")
	parser.add_argument(
		'--memory', type=str, required=False, help='mem (default 4096 = 4gb)', default="4096")
	parser.add_argument(
		'--rundir', type=str, required=False, help='directory to run in (default current dir)')
	
	# SRST2 inputs
	parser.add_argument(
		'--script', type=str, required=True, help='SRST2 script (/vlsci/VR0082/shared/srst2_sep/srst2_1509_reporting2.py)', 
		default="/vlsci/VR0082/shared/srst2_sep/srst2_1509_reporting2.py")
	parser.add_argument(
		'--output', type=str, required=True, help='identifier for outputs (will be combined with read set identifiers)')
	parser.add_argument(
		'--input_se', nargs='+', type=str, required=False, help='Input single end reads')
	parser.add_argument(
		'--input_pe', nargs='+', type=str, required=False, help='Input paired end reads')
	parser.add_argument(
		'--forward', type=str, required=False, default="_1", 
			help='Designator for forward reads (e.g default is _1, expect forward reads sample_1.fastq.gz)')
	parser.add_argument(
		'--reverse', type=str, required=False, default="_2", 
			help='Designator for reverse reads (e.g default is _2, expect reverse reads sample_2.fastq.gz)')
	parser.add_argument('--mlst_db', type=str, required=False, nargs=1, help='Fasta file of MLST alleles (optional)')
	parser.add_argument('--gene_db', type=str, required=False, nargs='+', help='Fasta file/s for gene databases (optional)')
	parser.add_argument(
		'--other_args', type=str, required=False, help='string containing all other arguments to pass to function')
		
	return parser.parse_args() 

def get_readFile_components(full_file_path):
	(file_path,file_name) = os.path.split(full_file_path)
	m1 = re.match("(.*).gz",file_name)
	ext = ""
	if m1 != None:
		# gzipped
		ext = ".gz"
		file_name = m1.groups()[0]
	(file_name_before_ext,ext2) = os.path.splitext(file_name)
	full_ext = ext2+ext
	return(file_path,file_name_before_ext,full_ext)

def read_file_sets(args):	

	fileSets = {} # key = id, value = list of files for that sample
	num_single_readsets = 0
	num_paired_readsets = 0
	
	if args.input_se:
		# single end
		for fastq in args.input_se:
			(file_path,file_name_before_ext,full_ext) = get_readFile_components(fastq)
			fileSets[file_name_before_ext] = [fastq]
			num_single_readsets += 1
			
	if args.input_pe:
		# paired end
		forward_reads = {} # key = sample, value = full path to file
		reverse_reads = {} # key = sample, value = full path to file
		num_paired_readsets = 0
		num_single_readsets = 0
		for fastq in args.input_pe:
			(file_path,file_name_before_ext,full_ext) = get_readFile_components(fastq)
			m=re.match("(.*)("+args.forward+")",file_name_before_ext)
			if m!=None:
				# store as forward read
				(baseName,read) = m.groups()
				forward_reads[baseName] = fastq
			else:
				m=re.match("(.*)("+args.reverse+")",file_name_before_ext)
				if m!=None:
				# store as reverse read
					(baseName,read) = m.groups()
					reverse_reads[baseName] = fastq
				else:
					print "Could not determine forward/reverse read status for input file " + fastq
		# store in pairs
		for sample in forward_reads:
			if sample in reverse_reads:
				fileSets[sample] = [forward_reads[sample],reverse_reads[sample]] # store pair
				num_paired_readsets += 1
			else:
				fileSets[sample] = [forward_reads[sample]] # no reverse found
				num_single_readsets += 1
				print 'Warning, could not find pair for read:' + forward_reads[sample]
		for sample in reverse_reads:
			if sample not in fileSets:
				fileSets[sample] = reverse_reads[sample] # no forward found
				num_single_readsets += 1
				print 'Warning, could not find pair for read:' + reverse_reads[sample]

	if num_paired_readsets > 0:
		print 'Total paired readsets found:' + str(num_paired_readsets)

	if num_single_readsets > 0:
		print 'Total single reads found:' + str(num_single_readsets)

	return fileSets 

def run_command(command, **kwargs):
	'Execute a shell command and check the exit status and any O/S exceptions'
	command_str = ' '.join(command)
	print 'Running: {}'.format(command_str)
	try:
		exit_status = call(command, **kwargs)
	except OSError as e:
		exit("Command '{}' failed due to O/S error: {}".format(command_str, str(e)))
	if exit_status != 0:
		exit("Command '{}' failed with non-zero exit status: {}".format(command_str, exit_status))

def bowtie_index(fasta_files):
	'Build a bowtie2 index from the given input fasta(s)'
	for fasta in fasta_files:
		built_index = fasta + '.1.bt2'
		if os.path.exists(built_index):
			print 'Bowtie index for {} is already built...'.format(fasta)
		else:
			print 'Building bowtie2 index for {}...'.format(fasta)
			run_command(['bowtie2-build', fasta, fasta])

def samtools_index(fasta_files):
	for fasta in fasta_files:
		fai_file = fasta + '.fai'
		if os.path.exists(fai_file):
			print 'Samtools index for {} is already built...'.format(fasta)
		else:
			print 'Building samtools index for {}...'.format(fasta)
			run_command(['samtools', 'faidx', fasta])
def main():

	args = parse_args()

	if not args.rundir:
		args.rundir = os.getcwd()

	# parse list of file sets to analyse
	fileSets = read_file_sets(args) # get list of files to process
	
	# build indexes to avoid issues later
	if args.mlst_db:
		bowtie_index(args.mlst_db)
		samtools_index(args.mlst_db)
		
	if args.gene_db:
		bowtie_index(args.gene_db)
		samtools_index(args.gene_db)
	
	# build and submit commands
	for sample in fileSets:		
		cmd = "#!/bin/bash"
		cmd += "\n#SBATCH -p main"
		cmd += "\n#SBATCH --job-name=srst2" + sample + args.output
		cmd += "\n#SBATCH --ntasks=1"
		cmd += "\n#SBATCH --mem-per-cpu=" + args.memory
		cmd += "\n#SBATCH --time=" + args.walltime
		cmd += "\ncd " + args.rundir
		cmd += "\nmodule load bowtie2-intel/2.1.0"
		cmd += "\nmodule load samtools-intel/0.1.18"
		cmd += "\nmodule load python-gcc"
		cmd += "\n" + args.script
		fastq = fileSets[sample]
		if len(fastq) > 1:
			cmd += " --input_pe " + fastq[0] + " " + fastq[1]
			cmd += " --forward " + args.forward
			cmd += " --reverse " + args.reverse
		else:
			cmd += " --input_se " + fastq[0]
		if args.mlst_db:
			cmd += " --mlst_db " + args.mlst_db[0]
		if args.gene_db:
			cmd += " --gene_db " + " ".join(args.gene_db)
		cmd += " --output " + sample + "_" + args.output
		cmd += " --log log_" + sample + "_" + args.output + ".log"
		cmd += " " + args.other_args
		
		# print and run command
		print cmd
		os.system('echo "' + cmd + '" | sbatch')

if __name__ == '__main__':
	main()