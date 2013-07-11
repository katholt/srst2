srst
====

Short Read Sequence Typing for Bacterial Pathogens


SRST Version 2
==============

srst2.py

Authors - Michael Inouye (minouye@unimelb.edu.au), Harriet Dashnow (h.dashnow@gmail.com), Bernie Pope (bjpope@unimelb.edu.au)

Dependencies:
    Python Version 2.7.3
    bowtie2    http://bowtie-bio.sourceforge.net/bowtie2/index.shtml version 2.1.0
    SAMtools          http://samtools.sourceforge.net Version: 0.1.18 (Version: 0.1.19 DOES NOT WORK - loss of edge coverage)


Running srst2.py for MLST with paired end reads:
/vlsci/VR0082/hdashnow/SRST_github/srst/srst2.py --input_pe reads_1.fastq.gz reads_2.fastq.gz --output output_base_name --build MLST_alleles.fa --log log_file.log --input_type q

Some examples of running srst2.py in a Bpipe pipeline:
validate_srst2_MLST_Shigella.pipe
validate_srst2_resistance_ecoli.pipe	
validate_srst2_resistance_staph_single.pipe
Note: if using a Torque queueing system (e.g. on VLSCI Merri) put bpipe.config in the folder you are rinning the .pipe file from.


