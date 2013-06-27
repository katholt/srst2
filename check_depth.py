import sys

filename = "6581_3#8.all.renamed.nr.fsa.srst2.pileup"
f = open(filename, "r")

count = 0
total_depth = 0
for line in f:
	if line.startswith("aadA_aadA15_1_DQ393783"):
		depth = int(line.split()[3])
		pos = int(line.split()[1])
		count += 1
		total_depth += depth
avg_depth = total_depth/float(count)
print avg_depth
print count
print pos
