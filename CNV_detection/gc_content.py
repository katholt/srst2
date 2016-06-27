"""
This script calculates the GC content of every consensus sequence in collated FASTA files produced by the script consensus_alignment.py
in the SRST2 package. The script produces a table of GC contents as an output. Note that the format of sequence definition lines in the
FAST files must be ">[sample name] [locus] [allele]" (space-delimited).

Example command: python gc_content.py -s *.fasta -p Kp

Author: Yu Wan (wanyuac@gmail.com)
Last update: 8 April 2016
"""

from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC  # for calculating GC contents

def parse_args():
    parser = ArgumentParser(description = "Calculate GC contents of consensus sequences")
    parser.add_argument("--seq", "-s", dest = "seq", nargs = "+", type = str, required = True, help = "Input FASTA files containing collated consensus sequences")
    parser.add_argument("--prefix", "-p", dest = "prefix", type = str, required = False, default = "", help = "The prefix of the output file")
    return parser.parse_args()

def main():
    args = parse_args()
    output_filename = args.prefix + "__" + "gc.txt"
    open(output_filename, "w").close()  # erase the content of an existing output file or create a new file
    with open(output_filename, "a") as output: # open the file again for appending lines
        print >> output, "\t".join(["Sample", "Locus", "Allele", "GC"])  # the header line
        for fasta_file in args.seq:
            for record in SeqIO.parse(open(fasta_file, "r"), "fasta"):
                sample, locus, allele = record.description.split(" ")
                gc = round(GC(str(record.seq)), 4)
                print >> output, "\t".join([sample, locus, allele, str(gc)])

if __name__ == "__main__":
    main()