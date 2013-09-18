SRST2
====

Short Read Sequence Typing (for Bacterial Pathogens)

This program is designed to take Illumina sequence data, a MLST database and/or a database
of gene sequences (e.g. resistance genes, virulence genes, etc) and report the presence of
STs and/or reference genes.

Dependencies:

python (v2.7.5), scipy

bowtie2 v2.1.0     http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

SAMtools v0.1.18   http://samtools.sourceforge.net (NOTE 0.1.19 DOES NOT WORK)



Authors - Michael Inouye (minouye@unimelb.edu.au)
		- Harriet Dashnow (h.dashnow@gmail.com)
		- Kathryn Holt (kholt@unimelb.edu.au)
		- Bernie Pope (bjpope@unimelb.edu.au)
		
================

Basic usage - MLST
====

1 - Gather your input files:

(i) sequence reads

(ii) a fasta sequence database to match to. For MLST, this means a fasta file of all allele sequences. If you want to assign STs, you also need a tab-delim file which defines the ST profiles as a combination of alleles. You can retrieve these files automatically from pubmlst.org/data/ using the script provided:

getmlst.py --species "Escherichia coli"

2 - Run MLST calling:

srst2.py --input_pe strainA_1.fastq.gz strainA_2.fastq.gz --output test 
	--mlst_db Escherichia_coli.fasta --mlst_definitions ecoli.txt
	--gene_db resistance.fasta

3 - Check the outputs:

(i) MLST results are output in: "mlst__Escherichia_coli__strainA_test__results.txt"

Sample  ST      adk     fumC    gyrB    icd     mdh     purA    recA    mismatches      uncertainty     depth

strainA     152     11      63      7       1       14      7       7                       25.8319955826


All usage options
====
SRST2 - Short Read Sequence Typer (v2)

optional arguments:
  -h, --help            show this help message and exit
  
  --input_se INPUT_SE [INPUT_SE ...]
                        Input single end read
                        
  --input_pe INPUT_PE [INPUT_PE ...]
                        Input paired end reads
                        
  --forward FORWARD     Designator for forward reads (e.g default is _1,
                        expect forward reads sample_1.fastq.gz)
                        
  --reverse REVERSE     Designator for reverse reads (e.g default is _2,
                        expect reverse reads sample_2.fastq.gz)
                        
  --mlst_db MLST_DB     Fasta file of MLST alleles (optional)
  
  --mlst_delimiter MLST_DELIMITER
                        Character(s) separating gene name from allele number
                        in MLST database (default "-", as in arcc-1)
                        
  --mlst_definitions MLST_DEFINITIONS
                        ST definitions for MLST scheme (required if mlst_db
                        supplied and you want to calculate STs)
                        
  --ignore_last         Flag to ignore last column of ST definitions table
                        (e.g. sometimes an additional column is added to
                        indicate clonal complex, which is not part of the ST
                        definition).
                        
  --gene_db GENE_DB [GENE_DB ...]
                        Fasta file/s for gene databases (optional)
                        
  --min_coverage MIN_COVERAGE
                        Percent coverage cutoff for gene reporting (default
                        90)
                        
  --min_depth MIN_DEPTH
                        Minimum mean depth to flag as dubious allele call
                        (default 5)
                        
  --min_edge_depth MIN_EDGE_DEPTH
                        Minimum edge depth to flag as dubious allele call
                        (default 2)
                        
  --prob_err PROB_ERR   Probability of sequencing error (default 0.01)
  --read_type {q,qseq,f}
                        Input file type (for bowtie input; default is q=fastq;
                        other options: qseq=solexa, f=fasta)
                        
  --other OTHER         Other options for bowtie2
  
  --mapq MAPQ           Qual map
  
  --baseq BASEQ         Qual base
  
  --output OUTPUT       Output file prefix
  
  --log FILE            log progress in FILENAME, defaults to stdout
  
  --verbose             Switch on verbose reporting
  
  --use_existing_pileup
                        Use existing pileups if available, otherwise they will
                        be generated
                        
  --use_existing_scores
                        Use existing scores files if available, otherwise they
                        will be generated
                        
  --prev_output PREV_OUTPUT [PREV_OUTPUT ...]
                        SRST2 output files to compile (any new results from
                        this run will also be incorporated)

More basic usage examples
====

Run single read sets against MLST database and resistance database

srst2.py --input_pe pool11_tag2_1.fastq.gz pool11_tag2_2.fastq.gz 
	--output pool11_tag2_Shigella 
	--log log_pool11_tag2_Shigella.log 
	--gene_db /vlsci/VR0082/shared/srst2_sep/resistance.fasta 
	--mlst_db Escherichia_coli.fasta 
	--mlst_definitions ecoli.txt 
	--ignore_last 
	--verbose 

------------

Run against multiple read sets in serial

srst2.py --input_pe *.fastq.gz
	--output Shigella 
	--log Shigella.log 
	--gene_db /vlsci/VR0082/shared/srst2_sep/resistance.fasta 
	--mlst_db Escherichia_coli.fasta 
	--mlst_definitions ecoli.txt 
	--ignore_last 
	--verbose 

------------

Run against new read sets, merge with previous reports (individual or compiled)

srst2.py --input_pe strainsY-Z*.fastq.gz
	--output strainsA-Z
	--log shigella_new.log 
	--gene_db /vlsci/VR0082/shared/srst2_sep/resistance.fasta 
	--mlst_db Escherichia_coli.fasta 
	--mlst_definitions ecoli.txt 
	--ignore_last 
	--verbose
	--prev_output genes__resistance__strainA_results.txt
		  mlst__Escherichia_coli__strainA_results.txt
		  genes__resistance__strainB_results.txt
		  mlst__Escherichia_coli__strainB_results.txt
		  strainsC-X__compiledResults.txt

------------

Run against Enterococcus reads, where read names are different from the usual _1.fastq and _2.fastq

python srst2.py --input_pe strain_R1.fastq.gz strain_R2.fastq.gz 
	--forward _R1 --reverse _R2 
	--output strainA --log strainA.log 
	--gene_db /vlsci/VR0082/shared/srst2_sep/resistance.fasta 
	--mlst_db Enterococcus_faecium.fasta 
	--mlst_definitions efaecium.txt 
	--verbose
	
Compile results from completed runs
====

python srst2.py --prev_output *compiledResults.txt --output Shigella_report
 
Running lots of jobs and compiling results
====

Run against multiple read sets: submitting 1 job per readset to SLURM queueing system

python slurm_srst2.py --script srst2.py 
	--output Shigella1609 
	--input_pe *.fastq.gz 
	--other_args '--gene_db resistance.fasta 
	--mlst_db Escherichia_coli.fasta 
	--mlst_definitions ecoli.txt 
	--verbose 
	--ignore_last' 
	--walltime 0-1:0 
		> job_sub_list.txt

------------

SLURM job script usage options

optional arguments:
  -h, --help            show this help message and exit
  
  --walltime WALLTIME   wall time (default 0-0:30 = 30 minutes)
  
  --memory MEMORY       mem (default 4096 = 4gb)
  
  --rundir RUNDIR       directory to run in (default current dir)
  
  --script SCRIPT       SRST2 script (/vlsci/VR0082/shared/srst2_sep/srst2_150
                        9_reporting2.py)
                        
  --output OUTPUT       identifier for outputs (will be combined with read set
                        identifiers)
                        
  --input_se INPUT_SE [INPUT_SE ...]
                        Input single end reads
                        
  --input_pe INPUT_PE [INPUT_PE ...]
                        Input paired end reads
                        
  --forward FORWARD     Designator for forward reads (e.g default is _1,
                        expect forward reads sample_1.fastq.gz)
                        
  --reverse REVERSE     Designator for reverse reads (e.g default is _2,
                        expect reverse reads sample_2.fastq.gz)
                        
  --other_args OTHER_ARGS
                        string containing all other arguments to pass to
                        function
