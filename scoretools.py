import multiprocessing
import itertools
import pickle as pickle
import argparse as ap
import scipy as sp
import scipy.stats
import sys as sys
import re as re

class Score(object):
	def __init__(self, variant, pid, escore, sscore, sig, mixed, mismatches, length, minDepth, meanDepth, maxDepth):
		self.variant = variant
		self.pid = pid
		self.escore = escore
		self.sscore = sscore
		self.sig = sig
		self.mixed = mixed
		self.mismatches = mismatches
		self.length = length
		self.minDepth = minDepth
		self.meanDepth = meanDepth
		self.maxDepth = maxDepth


	def __str__(self):
		ordered = [self.variant, self.pid, self.escore, self.sscore, self.sig,
			self.mixed, self.mismatches, self.length, self.minDepth, self.meanDepth, self.maxDepth]
		return ",".join([str(v) for v in ordered])


	def __cmp__(self, other):
		#Changed: swapped mismatche for pid and
		#switched the order of length and mixed
		
		#if self.pid != other.pid:
		#	return cmp(other.pid, self.pid)
		if self.mismatches != other.mismatches:
			return cmp(self.mismatches, other.mismatches)
		elif self.length != other.length:
			return cmp(other.length, self.length)
		elif self.mixed != other.mixed:
			return cmp(self.mixed, other.mixed)
		else:
			return cmp(other.sscore, self.sscore)


	def isSignificant(self):
		return self.sscore > self.sig


def sortKeysByScore(scores):
	return sorted(scores.keys(), key=scores.get)

#Everything to do with parsing the pileup file
#could go into a separate module
class Pileup(object):
	def __init__(self, id):
		self.id = id #Reference id
		self.pos = [] #Position in reference
		self.ref = [] #Reference sequence
		self.cons = [] #Consensus for reads
		self.depth = [] #Read depth
		self.cross = [] #Bases across reads
		self.mis = [] #Number of mismatched reads
		self.counts = [] #Count for each base
		self.q = [] #Quality scores

	
	def addCrossSection(self, pos, ref, depth, cross, q):
		ref = "N" if ref == "R" else ref
		#Parse the string representing the cross section
		self.cross.append(_parseCrossSection(cross, ref))
		
		#Call the consensus and add the called nucleotide and
		#the count of all nucleotides to the appropriate lists
		cons, counts = _callConsensus(self.cross[-1])
		self.cons.append(cons)
		self.counts.append(counts)

		#Add the remaining elements of the cross section
		self.pos.append(pos)
		self.ref.append(ref)
		self.depth.append(depth)
		self.q.append(q)
		

def _callConsensus(cross):
	#Count the occurence of each base in the cross section
	counts = dict(A=0, G=0, C=0, T=0, N=0, X=0)
	for c in cross:
		if c in counts:
			counts[c] += 1

	#Sort and call the consensus
	sortedNucs = sorted(counts, key=counts.get, reverse=True)
	return sortedNucs[0], counts


def _parseCrossSection(cs, ref):
	"""Parse the cross section string from a pileup"""
	parsed = []
	indel = False
	indelLength = 0

	cs = cs.upper()

	#Remove the read start and end characters
	cs = re.sub("\^.|\$", "", cs)
	
	#Find the lengths of any indels
	indelLengths = re.findall("[-+]([0-9]+)", cs)
	indels = re.findall("[-+]([0-9]+)", cs)
	
	#Remove the indel strings
	for length in indelLengths:
		cs = re.sub("[-+][0-9]+[ACGTN]{%d}" % int(length) , "", cs)

	#Parse the remaining characters
	for c in cs:
		#base matches reference
		if c == "." or c == ",":
			parsed.append(ref)
		elif c == "*":
			parsed.append("X")
		else:
			parsed.append(c)

	return parsed


#Parse a pileup file, creating a generator
#that yields the pileup for each reference 
def parsePileup(pileupFn):
	with open(pileupFn) as pileupf:
		pileup = None#Pileup()

		for line in pileupf:
			fields = line.split()

			#If we hit the next reference in the pileup	
			if pileup == None:
				pileup = Pileup(fields[0])
			elif pileup.id != fields[0]:
				#Return the pileup and then reset for the next one
				yield pileup
				pileup = Pileup(fields[0])

			#Parse each of the fields and append the cross-section to the pileup
			pileup.addCrossSection(int(fields[1]), fields[2], int(fields[3]), fields[4], fields[5])
	
	#Yield the final pileup
	yield pileup


def isMixed(ref, counts, depth, pErrorAndBase, sig):
	return counts[ref] > 0 and any(sp.stats.binom_test(counts[nuc], depth, pErrorAndBase) < sig for nuc in counts if nuc != "N" and nuc != ref)


def scoreRef(ref, counts, depth, pError):
	"""Calculate the p-value for a binomial test. This is the 
	probability of observing k reads with base r if in fact the true base was
	some q != r. If p < 0.05, we can reject the hypothesis that q != r
	"""
	return sp.stats.binom_test(counts[ref], depth, pError / 3.)


def scoreRefH(ref, counts, depth, maxLength):
	m = sum(d == 0 for d in depth)
	sumLogs = sum(sp.log2((c[r] + 1) / float(t + 4))**2 for r, c, t in zip(ref, counts, depth)) 
	return len(ref)**2 / (maxLength * sumLogs * (m + 1))
	

def scorePileup(pileup, pmin, pmax, pError, maxLength, maxMismatches=float('inf'), minDepth=0, minPID=0, minSignificance=0.05):
	length = float(pmax - pmin)
	ref = pileup.ref[pmin:pmax]
	cons = pileup.cons[pmin:pmax]
	counts = pileup.counts[pmin:pmax]
	depth = pileup.depth[pmin:pmax]

	matches = [r == c or r == "N" for r, c in zip(ref, cons)]
	pid = sum(matches) / float(len(ref))
	numMismatches = len(ref) - sum(c[r] > 0 for r, c in zip(ref, counts))
	
	#If the pileup meets the requirements for the percent identity, minimum depth and minimum coverage
	if pid < minPID or min(depth) < minDepth or numMismatches > maxMismatches:
		logBFSig, numMixed, meanSig, hscore = (-1, -1, -1, -1)
	else:
		#Calculate the bonferroni corrected significance threshold
		bfsig = minSignificance / float(length)
		logBFSig = -sp.log10(bfsig)
		
		#Calculate the p-values for each base where there are reads matching the reference and then take the mean
		pvalues = [scoreRef(r, c, d, pError) for r, c, d in zip(ref, counts, depth) if c[r] > 0]
		meanSig = -sp.log10(sp.mean(pvalues))
		
		#Calculate the number of bases where there is significant evidence for at least two alternate
		#variants (i.e. there is significant number of reads matching some base q != ref)
		numMixed = sum(isMixed(r, c, d, pError / 3., bfsig) for r, c, d in zip(ref, counts, depth))

		hscore = scoreRefH(ref, counts, depth, maxLength) 
	
	return Score(pileup.id, pid, hscore, meanSig, logBFSig, numMixed, numMismatches, len(ref), min(depth), sp.mean(depth), max(depth))


def scoreVariants(pileupFn, rangesFn, regex, pError, varLengths, maxMismatches, minDepth, minSimilarity, minSignificance):
	#Unpickle the ranges
	if rangesFn:
		ranges = pickle.load(open(rangesFn))
	else:
		ranges = None
	
	for p in parsePileup(pileupFn):
		if ranges != None:
			pmin, pmax = ranges[p.id]
		else:
			pmin, pmax = 0, len(p.ref)
		
		#locus = re.search("(.*)-\d", p.id).group(1)
		locus = re.search(regex, p.id).group(1)
		maxLength = varLengths[locus]
		
		#Score the pileup
		yield scorePileup(p, pmin, pmax, pError, maxLength, maxMismatches, minDepth, minSimilarity, minSignificance)


def mpHelperFunction((p, pError, varLengths, maxMismatches, minDepth, minSimilarity, minSignificance)):
	pmin, pmax = 0, len(p.ref)
	locus = re.search("(.*)-\d", p.id).group(1)
	maxLength = varLengths[locus]
	return scorePileup(p, pmin, pmax, pError, maxLength, maxMismatches, minDepth, minSimilarity, minSignificance)


def scoreVariantsMultiProc(pileupFn, rangesFn, pError, varLengths, maxMismatches, minDepth, minSimilarity, minSignificance, procs):
	pool = multiprocessing.Pool(processes=procs)
	ps = parsePileup(pileupFn)
	scores = []
	
	chunk = 100
	perrs, vls, mgs, mds, msims, msigs = ([pError] * chunk,
			[varLengths] * chunk, [maxMismatches] * chunk, [minDepth] * chunk, [minSimilarity] * chunk, [minSignificance] * chunk)

	while True:
		res = pool.map(mpHelperFunction, itertools.izip(itertools.islice(ps, chunk), perrs, vls, mgs, mds, msims, msigs))
		if res:
			scores.extend(res)
		else:
			break
		break
	return scores


if __name__ == "__main__":
	parser = ap.ArgumentParser()
	parser.add_argument("--pileup", help="The pileup file to process")
	parser.add_argument("--ranges", default="")
	args = parser.parse_args()

	#loci = ["aroe", "ddl", "gdh", "gki", "recP", "spi", "xpt"]
	#loci = ["ADK", "FUMC", "GYRB", "ICD", "MDH", "PURA", "RECA"]

	#Calculate the scores
	#scores = scoreVariants(args.pileup, args.ranges)
	
	#Split the scores by locus and then write to file, sorted by score
	#scoresByLocus = splitScoresByLocus(scores, loci)

	#for locus in scoresByLocus.keys():
	#	writeScores(scoresByLocus[locus], sys.stdout, sort=True)

