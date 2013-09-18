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

(i) sequence reads (this example uses paired reads in gzipped fastq format, see below for options)

(ii) a fasta sequence database to match to. For MLST, this means a fasta file of all allele sequences. If you want to assign STs, you also need a tab-delim file which defines the ST profiles as a combination of alleles. You can retrieve these files automatically from pubmlst.org/data/ using the script provided:

getmlst.py --species "Escherichia coli"

2 - Run MLST:

srst2.py --input_pe strainA_1.fastq.gz strainA_2.fastq.gz --output test 
	--mlst_db Escherichia_coli.fasta --mlst_definitions ecoli.txt

3 - Check the outputs:

(i) MLST results are output in: "mlst__Escherichia_coli__strainA_test__results.txt"

Sample  ST      adk     fumC    gyrB    icd     mdh     purA    recA    mismatches      uncertainty     depth

strainA     152     11      63      7       1       14      7       7                       25.8319955826

Basic usage - Resistance genes
====

1 - Gather your input files:

(i) sequence reads (this example uses paired reads in gzipped fastq format, see below for options)

(ii) a fasta sequence database to match to. For resistance genes, this means a fasta file of all the resistance genes/alleles that you want to screen for, clustered into gene groups. A suitable database, which combines sequences from ResFinder and CARD, is distributed with SRST2 (resistance.fasta).

2 - Run gene detection:

srst2.py --input_pe strainA_1.fastq.gz strainA_2.fastq.gz --output test --gene_db resistance.fasta

3 - Check the outputs:

(i) Gene detection results are output in: "genes__resistance__strainA_test__results.txt"

Sample  aadA    dfrA    sul2    tet(B)

strainA     aadA1-5 dfrA1_1 sul2_2  tet(B)_4

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

Input read formats and options
====

Any number of readsets can be provided using --input_se (for single end reads) and/or --input_pe (for paired end reads). You can provide both types of reads at once. Note however that if you do this, each readset will be typed one by one (in serial). So if it takes 2 minutes to type each read set, and you provide 10 read sets, it will take 20 minutes to get your results. The better way to proces lots of samples quickly is to give each one its own srst2.py job (e.g. submitted simultaneously to your job scheduler or server); then compile the results into a single report using "srst2.py --prev_output *results.txt --output all". That way each readset's 2 minutes of analysis is occurring in parallel on different nodes, and you'll get your results for all 10 samples in 2 minutes rather than 20.

Read formats - reads can be in any format readable by bowtie2. The format is passed on to the bowtie2 command via the --read_type flag in srst2. The default format is fastq (passed to bowtie 2 as q); other options are qseq=solexa, f=fasta. So to use fasta reads, you would need to tell srst2 this via '--read_type f'.

Reads may be gzipped.

Paired reads - bowtie2 requires forward and reverse reads to be supplied in separate files, e.g strainA_1.fastq.gz and strainA_2.fastq.gz. srst2 attempts to sort reads supplied via --input_pe into read pairs, based on the suffix that occurs before the file extension (which in this example would be .fastq.gz). So if you supplied --input_pe strainA_1.fastq.gz strainB_1.fastq.gz strainA_2.fastq.gz strainB_2.fastq.gz, srst2 would sort these into two pairs (strainA_1.fastq.gz, strainA_2.fastq.gz) and (strainB_1.fastq.gz, strainB_2.fastq.gz) and pass each pair to bowtie2 for mapping. By default, the suffixes are assumed to "_1" for forward reads and "_2" for reverse reads, but you can tell srst2 if you have other conventions, via --forward and --reverse. E.g. if your files were named strainA_read1.fastq.gz and strainA_read2.fastq.gz, you would use these commands: --input_pe strainA_read1.fastq.gz strainA_read2.fastq.gz --forward _read1 --reverse _read2.

Sample names are taken from the first part of the read file name (before the suffix if you have paired reads). E.g. 'strainA_1.fastq.gz' is assumed to belong to a sample called "strainA"; 'strainB_C_1.fastq.gz" would be assumed to belong to a sample called "strainB_C". These sample names will be used to name all output files, and will appear in the results files.

MLST Database format
====
MLST databases are specified by allele sequences, and a profiles table. These can be downloaded from the public databases, ready to use with srst2, using the provided script getmlst.py (see above).

1 - Allele sequences file, fasta format. 

--mlst_db alleles.fasta

This should contain ALL allele sequences for the MLST scheme; i.e. if there are 7 loci in the scheme, then the sequences of all alleles from all 7 loci should appear together in this file. If you have one file per locus, just cat these together first. If you use getmlst.py, this is done for you.

--mlst_delimiter '-'

The names of the alleles (i.e. the fasta headers) are critical for a functioning MLST scheme and therefore for correct calling of STs. There are two key components to every fasta header: the name of the locus (e.g. in E. coli these are adk, fumC, gyrB, icd, mdh, purA, recA) and the number assigned to the allele (1, 2, 3 etc). These are usually separated by a delimiter like '-' or '_'; e.g. in the E. coli scheme, alleles are named adk-1, fumC-10, etc. For ST calling to work properly, srst2 needs to know what the delimiter is. By default, we assume it is '-' as this is the most common; however some schemes use '_' (e.g. in the C. difficile scheme, the first allele is 'adk_1', so you would need to set --mlst_delimiter '_' on the srst2 command line). If you use getmlst.py, it will remind you of this and try to guess for you what the most likely delimiter is.

2 - MLST definitions file, tab delimited format.

--mlst_definitions

This is the file that tells you the ST number that is assigned to known combinations of alleles. Column 1 is the ST, and subsequent columns are the loci that make up the scheme. The names of these loci must match the allele names in the sequences database, e.g. adk, fumC, gyrB, icd, mdh, purA, recA in the E. coli scheme. If you download data from pubmlst this should not be a problem. Sometimes there are additional columns in this file, e.g. a column assigning STs to clonal complexes. srst2 will ignore any columns that don't match gene names found in the allele sequences file.

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
                        
Known issues
====

Reference indexing - srst2 uses bowtie2 for mapping reads to reference sequences. To do this, srst2 must first check the index exists, call bowtie2-build to generate the index if it doesn't already exist, and then call bowtie2 to map the reads to this indexed reference. Occasionallly bowtie2 will return an Error message saying that it doesn't like the index. This seems to be due to the fact that if you submit multiple srst2 jobs to a cluster at the same time, they will all test for the presence of the index and, if index files are present, will proceed with mapping... but this doesn't mean the indexing process is actually finished, and so errors will arise. The simple way out of this is, if you are running lots of srst2 jobs, FIRST index your reference(s) for bowtie2 and samtools (using 'bowtie2-build ref.fasta ref.fasta' and 'samtools faidx ref.fasta'), then submit your srst2 jobs.
