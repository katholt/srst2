#!/usr/bin/env python
'''
This script generates SRST2 jobs for the SLURM scheduling system (http://slurm.schedmd.com/). It
allows many samples to be processed in parallel. After they all complete, the results can be
merged together using SRST2's --prev_output argument.

Some of the specifics are set up for the cluster used by Kat Holt's lab, so modifications may be
necessary to make it run properly on a different cluster using SLURM.
'''

import string, re, collections
import os, sys, subprocess
from subprocess import call, check_output, CalledProcessError, STDOUT
from argparse import (ArgumentParser, FileType)
import logging

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
	parser.add_argument(
		'--threads', type=int, required=False, help='number of CPUs per job', default=1)

	# SRST2 inputs
	parser.add_argument(
		'--script', type=str, required=True, help='path to srst2.py')
	parser.add_argument(
		'--output', type=str, required=True, help='identifier for outputs (will be combined with read set identifiers)')
	parser.add_argument(
		'--input_se', nargs='+', type=str, required=False, help='Single end read file(s) for analysing (may be gzipped)')
	parser.add_argument(
		'--input_pe', nargs='+', type=str, required=False, help='Paired end read files for analysing (may be gzipped)')
	parser.add_argument(
		'--forward', type=str, required=False, default="_1",
			help='Designator for forward reads (only used if NOT in MiSeq format sample_S1_L001_R1_001.fastq.gz; otherwise default is _1, i.e. expect forward reads as sample_1.fastq.gz)')
	parser.add_argument(
		'--reverse', type=str, required=False, default="_2",
			help='Designator for reverse reads (only used if NOT in MiSeq format sample_S1_L001_R2_001.fastq.gz; otherwise default is _2, i.e. expect forward reads as sample_2.fastq.gz)')
	parser.add_argument(
		'--other_args', type=str, required=True, help='string containing all other arguments to pass to srst2')

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
			m=re.match("(.*)(_S.*)(_L.*)(_R.*)(_.*)", file_name_before_ext)
			if m==None:
				fileSets[file_name_before_ext] = [fastq]
			else:
				fileSets[m.groups()[0]] = [fastq] # Illumina names
			num_single_readsets += 1

	elif args.input_pe:
		# paired end
		forward_reads = {} # key = sample, value = full path to file
		reverse_reads = {} # key = sample, value = full path to file
		num_paired_readsets = 0
		num_single_readsets = 0
		for fastq in args.input_pe:
			(file_path,file_name_before_ext,full_ext) = get_readFile_components(fastq)
			# try to match to MiSeq format:
			m=re.match("(.*)(_S.*)(_L.*)(_R.*)(_.*)", file_name_before_ext)
			if m==None:
				# not default Illumina file naming format, expect simple/ENA format
				m=re.match("(.*)("+args.forward+")$",file_name_before_ext)
				if m!=None:
					# store as forward read
					(baseName,read) = m.groups()
					forward_reads[baseName] = fastq
				else:
					m=re.match("(.*)("+args.reverse+")$",file_name_before_ext)
					if m!=None:
					# store as reverse read
						(baseName,read) = m.groups()
						reverse_reads[baseName] = fastq
					else:
						print "Could not determine forward/reverse read status for input file " + fastq
			else:
				# matches default Illumina file naming format, e.g. m.groups() = ('samplename', '_S1', '_L001', '_R1', '_001')
				baseName, read = m.groups()[0], m.groups()[3]
				if read == "_R1":
					forward_reads[baseName] = fastq
				elif read == "_R2":
					reverse_reads[baseName] = fastq
				else:
					print "Could not determine forward/reverse read status for input file " + fastq
					print "  this file appears to match the MiSeq file naming convention (samplename_S1_L001_[R1]_001), but we were expecting [R1] or [R2] to designate read as forward or reverse?"
					fileSets[file_name_before_ext] = fastq
					num_single_readsets += 1
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

class CommandError(Exception):
	pass

def run_command(command, **kwargs):
	'Execute a shell command and check the exit status and any O/S exceptions'
	command_str = ' '.join(command)
	print 'Running: {}'.format(command_str)
	try:
		exit_status = call(command, **kwargs)
	except OSError as e:
		message = "Command '{}' failed due to O/S error: {}".format(command_str, str(e))
		raise CommandError({"message": message})
	if exit_status != 0:
		message = "Command '{}' failed with non-zero exit status: {}".format(command_str, exit_status)
		raise CommandError({"message": message})

def check_bowtie_version():
	check_command_versions([get_bowtie_execs()[0], '--version'], 'version ', 'bowtie',
						   ['2.1.0','2.2.3','2.2.4','2.2.5','2.2.6','2.2.7','2.2.8','2.2.9'])

def check_samtools_version():
	check_command_versions([get_samtools_exec()], 'Version: ', 'samtools',
						   ['0.1.18','0.1.19','1.0','1.1','1.2','1.3','(0.1.18 is recommended)'])

def check_command_versions(command_list, version_prefix, command_name, required_versions):
	try:
		command_stdout = check_output(command_list, stderr=STDOUT)
	except OSError as e:
		logging.error("Failed command: {}".format(' '.join(command_list)))
		logging.error(str(e))
		logging.error("Could not determine the version of {}.".format(command_name))
		logging.error("Do you have {} installed in your PATH?".format(command_name))
		exit(-1)
	except CalledProcessError as e:
		# some programs such as samtools return a non-zero exit status
		# when you ask for the version (sigh). We ignore it here.
		command_stdout = e.output

	version_ok = False
	for v in required_versions:
		if version_prefix + v in command_stdout:
			version_ok = True

	if not version_ok:
		logging.error("Incorrect version of {} installed.".format(command_name))
		logging.error("{} versions compatible with SRST2 are ".format(command_name) + ", ".join(required_versions))
		exit(-1)

def get_bowtie_execs():
	'Return the "best" bowtie2 executables'

	exec_from_environment = os.environ.get('SRST2_BOWTIE2')
	if exec_from_environment and os.path.isfile(exec_from_environment):
		bowtie2_exec = exec_from_environment
	else:
		bowtie2_exec = None

	exec_from_environment = os.environ.get('SRST2_BOWTIE2_BUILD')
	if exec_from_environment and os.path.isfile(exec_from_environment):
		bowtie2_build_exec = exec_from_environment
	elif bowtie2_exec and os.path.isfile(bowtie2_exec+'-build'):
		bowtie2_build_exec = bowtie2_exec+'-build'
	else:
		bowtie2_build_exec = 'bowtie2-build'

	if bowtie2_exec is None:
		bowtie2_exec = 'bowtie2'

	return (bowtie2_exec, bowtie2_build_exec)

def bowtie_index(fasta_files):
	'Build a bowtie2 index from the given input fasta(s)'
	check_bowtie_version()
	for fasta in fasta_files:
		built_index = fasta + '.1.bt2'
		if os.path.exists(built_index):
			print 'Bowtie 2 index for {} is already built...'.format(fasta)
		else:
			print 'Building bowtie2 index for {}...'.format(fasta)
			run_command([get_bowtie_execs()[1], fasta, fasta])

def get_samtools_exec():
	'Return the "best" samtools executable'

	exec_from_environment = os.environ.get('SRST2_SAMTOOLS')
	if exec_from_environment and os.path.isfile(exec_from_environment):
		return exec_from_environment
	else:
		return 'samtools'

def samtools_index(fasta_files):
	'Build a samtools faidx index from the given input fasta(s)'
	check_samtools_version()
	for fasta in fasta_files:
		built_index = fasta + '.fai'
		if os.path.exists(built_index):
			print 'Samtools index for {} is already built...'.format(fasta)
		else:
			print 'Building samtools faidx index for {}...'.format(fasta)
			run_command([get_samtools_exec(), 'faidx', fasta])

def main():

	args = parse_args()

	if args.rundir:
		if not os.path.exists(args.rundir):
			os.makedirs(args.rundir)
	else:
		args.rundir = os.getcwd()

	# parse list of file sets to analyse
	fileSets = read_file_sets(args) # get list of files to process

	# make sure the databases are formated for bowtie2 and samtools before running the jobs
	db = []
	m = re.search( r'(--mlst_db) (.*?) .*', args.other_args)
	if m != None:
		db.append(m.group(2))
	g = re.search( r'(--gene_db) (.*?) --', args.other_args)
	if g != None:
		db += g.group(2).split()
	else:
		g = re.search( r'(--gene_db) (.*?)$', args.other_args)
		if g != None:
			db += g.group(2).split()
	bowtie_index(db)
	samtools_index(db)

	# build and submit commands
	for sample in fileSets:
		cmd = "#!/bin/bash"
		cmd += "\n#SBATCH -p sysgen"
		cmd += "\n#SBATCH --job-name=srst2_" + sample + "_" + args.output
		cmd += "\n#SBATCH --nodes=1"
		cmd += "\n#SBATCH --ntasks=1"
		cmd += "\n#SBATCH --cpus-per-task=" + str(args.threads)
		cmd += "\n#SBATCH --mem-per-cpu=" + args.memory
		cmd += "\n#SBATCH --time=" + args.walltime
		cmd += "\ncd " + args.rundir
		cmd += "\nmodule load srst2"
		cmd += "\n" + args.script
		fastq = fileSets[sample]
		if len(fastq) > 1:
			cmd += " --input_pe " + fastq[0] + " " + fastq[1]
			cmd += " --forward " + args.forward
			cmd += " --reverse " + args.reverse
		else:
			cmd += " --input_se " + fastq[0]
		cmd += " --output " + sample + "_" + args.output
		cmd += " --log"
		if args.threads > 1:
			cmd += " --threads " + str(args.threads)
		cmd += " " + args.other_args

		# print and run command
		print cmd
		print ''
		os.system('echo "' + cmd + '" | sbatch')
		print ''

if __name__ == '__main__':
	main()
