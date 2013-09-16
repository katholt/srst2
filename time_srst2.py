# detect time taken for srst

# modules
import string, re, collections, datetime
import os, sys, subprocess
from optparse import OptionParser

# BioPython modules for reading and writing sequences
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)
	
	return parser.parse_args()

if __name__ == "__main__":

	(options, args) = main()

	print "log sample start stop seconds"
	for log in args:
		print log,
		sample = os.popen("grep 'command line' " + log + " | awk '{print $10}'").read().rstrip()
		if sample == "":
			print "error",
		else:
			print sample,
		start = os.popen("grep 'program started' " + log + " | awk '{print $2}'").read().rstrip()
		if start == "":
			print "error",
		else:
			print start,
		stop = os.popen("grep 'SRST2 has finished' " + log + " | awk '{print $2}'").read().rstrip()
		if stop == "":
			print "error",
			print "NA"
		else:
			print stop,
			if start != "":
				seconds = (datetime.datetime.strptime(stop, "%H:%M:%S") - datetime.datetime.strptime(start, "%H:%M:%S")).total_seconds()
				print str(seconds)
			else:
				print "NA"