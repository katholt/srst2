import re as re
import sys as sys
import argparse as ap
import scoretools as st

def scoreFromStringList(fields):
	typedFields = [fields[0]] + [float(v) if re.search("\.", v) else int(v) for v in fields[1:]]
	return st.Score(*typedFields)


def printSummaryHeader(loci, fout):
	fout.write("Sample,ST,mixed," + ",".join(loci) + ",exact match\n")


def printScoresHeader(loci, fout):
	locusHeaders = ["percent identity", "entropy", "binomial", "significance", "mixed bases", "mismatches", "length", "min. read depth", "mean read depth", "max read depth"]
	fout.write("Sample,ST,all significant,all exact")
	for locus in loci:
		fout.write("," +  ",".join([locus] + locusHeaders))
	fout.write("\n")


def rankVariantCombinations(variants, scoresByLocus, loci):
	maxVars = max(len(vs) for vs in variants)
	sortedSTs = []
	
	for i in xrange(maxVars):
		if i == 0:
			base = [vs[0] for vs in variants]
			score = (0, 1, 1)
			base.append(score)
			sortedSTs.append(base)
		else:
			alternatives = [vs[i] if len(vs) > i else None for vs in variants]
			for j, alt in enumerate(alternatives):
				if alt != None:
					#Substitute the next alternative
					altsubs = base[:]
					altsubs[j] = alt
					
					#Claculate the change to the score from substituting the alternative
					score, bestScore = scoresByLocus[loci[j]][alt], scoresByLocus[loci[j]][base[j]]
					change = (score.mismatches - bestScore.mismatches, score.length / bestScore.length, score.sscore / bestScore.sscore)

					if 0 < change[2] <= 1.0:
						altsubs[-1] = change
						sortedSTs.append(altsubs)
	
	return sorted(sortedSTs, cmp=vcmp)


def vcmp(a, b):
	s1, s2 = a[-1], b[-1]
	if s1[0] != s2[0]:#Change in num. mismatches
		return cmp(s1[0], s2[0])
	elif s1[1] != s2[1]:#Change in length
		return cmp(s2[1], s1[1])
	else:#Change in significance
		return cmp(s2[2], s1[2])


def entropyCmp(a, b):
	return cmp(b.escore, a.escore)


def callSTsForSample(scoresByLocus, loci, sts, sample, scoresOut):
	"""Call the STs for the given sample. For each sample,
	the ST for top ranked variants across all loci is determined.
	Remaining alternatives are then given in order of the likelihood
	that they are present in the sample."""
	
	#Generate a list of lists, with each sub-list giving the sorted alternative variants for the corresponding locus
	variantsRankedBySScore = [sorted(scoresByLocus[locus].keys(), key=scoresByLocus[locus].get) if locus in scoresByLocus else ["none"] for locus in loci]  
	variantsRankedByEScore = [sorted(scoresByLocus[locus].keys(), key=scoresByLocus[locus].get, cmp=entropyCmp) if locus in scoresByLocus else ["none"] for locus in loci]  
	
	#Rank the possible combinations
	stsSortedBySScore = rankVariantCombinations(variantsRankedBySScore, scoresByLocus, loci)
	stsSortedByEScore = rankVariantCombinations(variantsRankedByEScore, scoresByLocus, loci)
	sortedSTs = stsSortedBySScore
	
	return ([resolveST(vs[:len(loci)], sts) for vs in sortedSTs], [vs[:len(loci)] for vs in sortedSTs])


def resolveST(variants, sts):
	if "none" in variants:
		st = "none"
	elif "mixed" in variants:
		st = "mixed"
	else:
		try:
			st = str(sts[tuple([int(re.search("-(\d+)", v).group(1)) for v in variants])])
		except KeyError:
			st = "novel-ST"

	return st


def writeSummary(sample, scores, calledSTs, calledVariants, loci, fout):
	mixedStr = "mixed" if len(calledSTs) > 1 else "unmixed"
	novelVariants = [v for locus, v in zip(loci, calledVariants[0]) if v in scores[locus] and scores[locus][v].mismatches > 0]
	if len(novelVariants) == 1:
		novelStr = "contains mismatches in " + novelVariants[0] + " possibly indicating a novel variant"
	elif len(novelVariants) > 1:
		novelStr = "contains mismatches in " + ",".join(novelVariants) + " possibly indicating novel variants"
	else:
		novelStr = "exact match"

	fout.write(",".join([sample, calledSTs[0], mixedStr] + calledVariants[0] + [novelStr]) + "\n")


def writeScores(sample, scores, calledSTs, variants, loci, fout):
	#Go through each called ST
	for st, vs in zip(calledSTs, variants):
		#Test if all variants were present and had significant scores, and if all varinats were exactly matched
		allSignificant = all(v in scores[locus] and scores[locus][v].isSignificant() for locus, v in zip(loci, vs)) 
		allExact = all(v in scores[locus] and scores[locus][v].mismatches == 0 for locus, v in zip(loci, vs)) 
		
		fout.write(",".join([sample, st, str(allSignificant), str(allExact)]))
		
		for locus, v in zip(loci, vs):
			try:
				score = scores[locus][v]
				fout.write("," + str(score))
			except KeyError:
				fout.write(",none")
		
		fout.write("\n")


def parseSTDatabase(databaseFilename):
	"""Parse a tab delimited ST database file, creating
	a dictionary of STs keyed by a tuple giving the
	allele at each loci, and a list of locus names"""
	sts = dict()

	with open(databaseFilename) as fin:
		loci = fin.readline().split()[1:]
		
		for line in fin:
			fields = [int(f) for f in line.split()]
			if len(fields) > 0:
				sts[tuple(fields[1:])] = fields[0]

	return sts, loci


def readVariantScores(fin, filterNonSignificant, filterThreshold):
	scores = dict()
	
	for line in fin:
		if not line.startswith("sample"):
			#Get the sample, locus, variant and associated score
			fields = line.strip().split(",")
			sample, locus, variant = fields[:3]
			score = scoreFromStringList(fields[2:])
			
			#If the score for the variant is significant
			if (not filterNonSignificant) or score.isSignificant():#(score.escore > 0 and score.escore >= filterThreshold):
				#Make sure the sample and locus are in the dictionaries
				if not sample in scores:
					scores[sample] = dict()
				if not locus in scores[sample]:
					scores[sample][locus] = dict()
				
				#Add the score
				scores[sample][locus][variant] = score 

	return scores
	

def main(args):
	#Read in the ST database, this gives the names of
	#each loci, and the variants associated with each st
	sts, loci = parseSTDatabase(args.st)
	scores = readVariantScores(sys.stdin, args.filter, args.filter_threshold)

	with open(args.base + args.summary, "w") as summaryOut:
		with open(args.base + args.scores, "w") as scoresOut:
			#Print the file headers
			printSummaryHeader(loci, summaryOut)
			printScoresHeader(loci, scoresOut)
			
			for sample in sorted(scores.keys()):
				#Call the STs for the current sample
				calledSTs, variants = callSTsForSample(scores[sample], loci, sts, sample, scoresOut)

				#Output the top ranked ST to the summary, indicating a possible mixture
				#if there are multiple called STs. Output all STs to the scores file 
				writeSummary(sample, scores[sample], calledSTs, variants, loci, summaryOut)
				writeScores(sample, scores[sample], calledSTs, variants, loci, scoresOut)


if __name__ == "__main__":
	parser = ap.ArgumentParser()
	parser.add_argument("--st", help="sequence type database file")
	parser.add_argument("--base", required=True, help="base name for output files")
	parser.add_argument("--regex", default="-(\d+)", help="Regex used to extract the variant number from the variant name")
	parser.add_argument("--summary", default="summary.csv", help="Name of summary file, appended to base")
	parser.add_argument("--scores", default="scores.csv", help="Name of score file, appended to base")
	parser.add_argument("--filter", action="store_true", help="filter out non-significant scores")
	parser.add_argument("--filter_threshold", type=int, default=10, help="threshold for significance filter") 
	args = parser.parse_args()

	main(args)
