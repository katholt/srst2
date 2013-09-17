
import rpy2.robjects as robjects
r = robjects.r

from Bio.Align.Applications import MuscleCommandline
import os
import subprocess as sub

fasta_directory = "/home/UNIMELB/hdashnow/resistance_database/by_gene"

# Get the filenames of all files in the input directory "fasta_directory"
files = []
for f in os.listdir(fasta_directory):
	if not f.startswith("."):
            if f.endswith(".fsa") or f.endswith(".fasta"):
		files.append('"'+os.path.join(fasta_directory,f)+'"') # need to surround with " " due to () in filenames
print "Number of input files:", len(files)

# Run muscle on each fasta file to produce an alignment for each
# Save the alignment file names so that they can be used in R ape
for f in files:
	outfilename = (f+".aln").replace("(","").replace(")","")
	#print outfilename
        alignment = MuscleCommandline(input=f, out=outfilename)
	#alignment()


r('''
	require(ape, quietly=TRUE)

	all_files = list.files("/home/UNIMELB/hdashnow/resistance_database/by_gene",
                       pattern="aln$", full.names=TRUE)

	pdf(width=9,height=12,file="trees.pdf")

	for (filename in all_files) {
		#print(filename)
		aln = read.dna(filename, format="fasta") 
		if (dim(aln)[1] >= 3) {
                    d = dist(aln)
		    #print(aln)
		    tree=nj(d)
		    plot(tree, cex=0.6)
                    #print(filename)
		    min_filename = tail(strsplit(filename,"/")[[1]],1) # This needs to be made robust
		    #print(min_filename)
		    gene_name = strsplit(min_filename,"\\\.")[[1]][1]
		    #print(cluster)
		    title(paste("Gene: ", gene_name))
    		    #dev.off()
                }
	}
	dev.off()
''')
