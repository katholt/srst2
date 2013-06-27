#!/bin/bash
# Created by the VLSCI job script generator
# Fri Jun 14 2013 10:12:38 GMT+1000 (EST)

# Queue for the job:
#PBS -q batch

# The name of the job:
#PBS -N test_entropy 

# The total amount of memory in gigabytes used by the whole job:
#PBS -l mem=4gb

# Send yourself an email when the job (a)borts abnormally (b)egins (e)nds successfully
#PBS -m abe

# The maximum running time of the job in hours:mins:secs
#PBS -l walltime=0:40:00
# Run the job from the directory where it was launched:
cd $PBS_O_WORKDIR

# The modules to load:

# The job command(s):
bamtools coverage -in 6581_3#8.all.renamed.nr.fsa.srst2.sorted.bam -out 6581_3#8.all.renamed.nr.fsa.srst2.sorted.bam.coverage 2> /dev/null

python /vlsci/VR0082/hdashnow/python_code/getCoverByChrLength.py /vlsci/VR0082/hdashnow/raw_data/reference/all.renamed.nr.fsa 6581_3#8.all.renamed.nr.fsa.srst2.sorted.bam.coverage 6581_3#8.all.renamed.nr.fsa.srst2.sorted.bam.coverage.bychr

python /vlsci/VR0082/hdashnow/python_code/parse_pileup_CLI_2.py -i 6581_3#8.all.renamed.nr.fsa.srst2.pileup -o 6581_3#8.all.renamed.nr.fsa.srst2.pileup.entropyscore -c 6581_3#8.all.renamed.nr.fsa.srst2.sorted.bam.coverage.bychr
