import sys
test_allele = "234__blaCTX-M__blaCTX-M-28_1_EU531513"
filename = "/vlsci/VR0082/hdashnow/SRST_github/srst2_etc_output/6581_3#8.all_fixed_annotated.fsa.srst2.pileup"
f = open(filename, "r")

count = 0
total_depth = 0
for line in f:
	if test_allele in line:
		depth = int(line.split()[3])
		pos = int(line.split()[1])
		count += 1
		total_depth += depth
avg_depth = total_depth/float(count)
print test_allele, avg_depth
#print count
#print pos
