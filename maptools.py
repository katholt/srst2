"""Copyright David Savage, 2013.

This module is a wrapper around bwa, bowtie2, and samtools,
and provides functionality to quickly perform mapping of short
read sequence data to a set of reference sequences, to
create a sorted, indexed bam file describing this mapping, and
to provide a means for reading this bam file into an in-memory pileup
"""

import subprocess as sp
import argparse as ap
import sys as sys
import re as re
import os as os

class Pileup(object):
	def __init__(self):
		self.pid = None #Reference id
		self.pos = [] #Position in reference
		self.ref = [] #Reference sequence
		self.cons = [] #Consensus for reads
		self.depth = [] #Read depth
		self.cross = [] #Bases across reads
		self.mis = [] #Number of mismatched reads
		self.q = [] #Quality scores


def _callConsensus(pileup):
	pileup.cons = []
	pileup.mis = []

	for i, cross in enumerate(pileup.cross):
		#Count the occurence of each base in the cross section
		counts = dict(A=0, G=0, C=0, T=0, N=0)
		for c in cross:
			if c in counts:
				counts[c] += 1

		#Sort and call the consensus
		sortedCounts = sorted(counts, key=counts.get, reverse=True)
		pileup.cons.append(sortedCounts[0])
		pileup.mis.append(sum((counts[key] for key in sortedCounts[1:])))


def parseCrossSection(cs, ref):
	"""Parse the cross section string from a pileup
	"""
	parsed = []
	indel = False
	indelLength = 0

	cs = cs.upper()

	#Remove the read start and end characters
	cs = re.sub("\^.|\$", "", cs)
	
	#Find the lengths of any indels
	indelLengths = re.findall("[-+]([0-9]+)", cs)
	
	#Remove the indel strings
	for length in indelLengths:
		cs = re.sub("[-+][0-9]+[ACGTN]{%d}" % int(length) , "", cs)

	#Parse the remaining characters
	for c in cs:
		#base matches reference
		if c == "." or c == ",":
			parsed.append(ref)
		#base doesn't match reference (SNP or read extension)
		else:
			parsed.append(c)

	return parsed


#Parse a pileup file, creating a generator
#that yields the pileup for each reference 
def parsePileup(pileupf):
	#with open(pileupFn) as pileupf:
	pileup = Pileup()

	for line in pileupf:
		fields = line.split()
		
		#If we hit the next reference in the pileup	
		if pileup.pid == None:
			pileup.pid = fields[0]
		elif pileup.pid != fields[0]:
			#Call the consensus and return the pileup
			_callConsensus(pileup)
			yield pileup

			#Reset for the next reference
			pileup = Pileup()
			pileup.pid = fields[0]

		#Parse each of the fields and append a dict
		#representing the cross-section to the list
		pileup.pos.append(int(fields[1]))
		pileup.ref.append(fields[2])
		pileup.depth.append(int(fields[3]))
		pileup.cross.append(parseCrossSection(fields[4], fields[2]))
		pileup.q.append(fields[5])

	_callConsensus(pileup)
	yield pileup


def splitAlternativeHits(line):
	"""Process the given read line (from a bwa generated SAM)
	by detecting alternative locations that the read maps to,
	returning a list of lines giving these alternative locations
	"""

	newLines = []

	#Split the line into individual fields
	fields = line.split()
	
	if len(fields) > 11:
		#Look for alternative hits denoted by bwa using the xa tag
		for field in fields[11:]:
			if field.startswith("XA"):
				alts = field.split(":")[2].split(";")

				#Split gives an empty set at the end, because
				#bwa puts a ; after the last entry
				for alt in alts[:-1]:
					variant, pos, cigar, nm = alt.split(',')
					pos = int(pos)

					if pos > 0:
						r = (fields[0], '0', variant, str(pos), fields[4], cigar, '*', '0', '0', fields[9], fields[10])
					else:
						r = (fields[0], '16', variant, str(-pos), fields[4], cigar, '*', '0', '0', fields[9], fields[10])

					newLines.append("\t".join(r) + "\n")

	return newLines


def fixDuplicateBWAHits(samFn, outFn, mapper="bwa"):
	with open(outFn, "w") as fout:
		with open(samFn) as fin:
			for line in fin:
				fout.write(line)
				for nl in splitAlternativeHits(line):
					fout.write(nl)


def fixDuplicateBowtieHits(samFn, outFn, mapper="bwa"):
	with open(outFn, "w") as fout:
		with open(samFn) as fin:
			seqLengths = dict()
			for line in fin:
				if line.startswith("@"):
					match = re.search("SN:(.*?)\tLN:(\d+)", line) 
					if match:
						seqLengths[match.group(1)] = int(match.group(2)) 
					
					fout.write(line)
				else:
					fields = line.split()
					flag, cigar, start = (int(fields[1]), fields[5], int(fields[3]))

					if not flag & 4 and flag & 256:
						fields[1] = str(flag - 256)
					
					fields[3], fields[5] = undoClipping(start, cigar, len(fields[9]), seqLengths[fields[2]])

					fout.write("\t".join(fields) + "\n")


def undoClipping(start, cigar, seqLength, refLength):
	match = re.search("^(\d)S", cigar)
	if match:
		startClip = int(match.group(1))
		if start - startClip > 1:
			start = start - startClip
			cigar = re.sub("S", "M", cigar, count=1)
	
	match = re.search("(\d)S$", cigar)
	if match and start + seqLength < refLength:
		endClip = match.group(1)
		cigar = re.sub(endClip + "S$", endClip + "M", cigar, count=1)

	return str(start), cigar


def mapReadsBowtie(fastqs, ref, samFilename, tmpDir):
	#Perform indexing of the reference sequence
	#if it hasn't previously been performed
	if not os.path.exists(ref + ".index.1.b2"):
		sp.call(["bowtie2-build", "--quiet", ref, ref + ".index"])

	ref += ".index"
	
	if len(fastqs) == 2:
		sp.call(["bowtie2", "--quiet", "--very-sensitive-local", "-a", "--no-unal", "-S", tmpDir + "/map.sam", "-x", ref, "-1", fastqs[0], "-2", fastqs[1]])
	else:
		sp.call(["bowtie2", "--quiet", "--very-sensitive-local", "-a", "--no-unal", "-S", tmpDir + "/map.sam", "-x", ref, "-U", fastqs[0]], stdout=None, stderr=None)

	#Fill out each alternative hit so that each variant has a full map
	fixDuplicateBowtieHits(tmpDir + "/map.sam", samFilename)


def mapReadsBWA(fastqs, ref, samFilename, tmpDir):
	#Index the reference
	sp.call(["bwa", "index", ref])
	
	if len(fastqs) == 2:
		#Perform the alignment of each paired end file
		sp.call(["bwa", "aln", "-f", tmpDir + "/map1.sai", ref, fastqs[0]])
		sp.call(["bwa", "aln", "-f", tmpDir + "/map2.sai", ref, fastqs[1]])
			
		#Combine alignments into a single SAM file
		sp.call(["bwa", "sampe", "-n", "99999", "-f", tmpDir + "/map.sam", ref, tmpDir + "/map1.sai", tmpDir + "/map2.sai"] + fastqs)
		#tmpDir + "/map1.sai", tmpDir + "/map2.sai", fastq1, fastq2])
	else:
		sp.call(["bwa", "aln", "-f", tmpDir + "/map.sai", ref, fastqs[0]])
		sp.call(["bwa", "samse", "-n", "99999", "-f", tmpDir + "/map.sam", ref] + sais + fastqs)

	#Duplicate reads and fill out each alternative hit
	fixDuplicateBWAHits(tmpDir + "/map.sam", samFilename)


def convertToBam(samFilename, bamFilename, deleteSam=False):
	#Convert to bam
	with open(bamFilename, "w") as fout:
		sp.call(["samtools", "view", "-Sb", samFilename], stdout=fout)
	
	if deleteSam:
		os.remove(samFilename)


def mapReads(mapper, ref, fastqs, base):
	#Make a temporary directory to hold the intermediate files
	tmpDir = os.path.abspath(base + "tmp")
	try:
		os.makedirs(tmpDir)
	except OSError:
		print >> sys.stderr, "Couldn't make tmp directory" + tmpDir
	
	#Generate the sam files from each mapping
	if mapper == "bwa":
		mapReadsBWA(fastqs, ref, base + ".sam", tmpDir)	
	elif mapper == "bowtie":
		mapReadsBowtie(fastqs, ref, base + ".sam", tmpDir)

	#Convert the sam to bam
	convertToBam(base + ".sam", base + ".bam", False)

	#Sort and index the bam file
	sp.call(["samtools", "sort", base + ".bam", base + "_sorted"])
	sp.call(["samtools", "index", base + "_sorted.bam", base + "_sorted.bam.index"])

	#Create the pileup file
	with open(base + ".pileup", "w") as f:
		sp.call(["samtools", "mpileup", "-f", ref, base + "_sorted.bam"], stdout=f)

	#Remove the temporary files
	for f in os.listdir(tmpDir):
		os.remove(tmpDir + "/" + f)
	os.rmdir(tmpDir)


if __name__ == "__main__":
	parser = ap.ArgumentParser()
	parser.add_argument("--mapper", default="bwa")
	parser.add_argument("--ref", help="name of reference fasta file")
	parser.add_argument("--fastq1", help="first set of paired ends (fastq.gz)")
	parser.add_argument("--fastq2", help="second set of paired ends (fastq.gz)")
	parser.add_argument("--base", help="base name for output files")
	args = parser.parse_args()

	mapReads(args.mapper, args.ref, args.fastq1, args.fastq2, args.base)
	
	
