# Replacement function for: 
# pileup_binomial_scoring(sam_file, size)
# Harriet Dashnow 28/6/13

pileup_file = "/vlsci/VR0082/hdashnow/SRST_github/srst2_etc_output/6581_3#8.all_fixed_annotated.fsa.srst2.pileup"
size = 537
prob_err = 0.01

def pileup_binomial_scoring(pileup_file, size):
	with open(pileup_file) as pileup:
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

	#return hash_alignment, hash_max_depth, hash_edge_depth, avg_depth_allele

pileup_binomial_scoring(pileup_file, size)
