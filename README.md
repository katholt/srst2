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

Gene database format
====

In addition to MLST, srst2 can do gene/allele detection. This works by mapping reads to each of the reference sequences in a fasta file(s) (provided through --gene_db) and reporting details of all genes that are covered above a certain level (--min_coverage, 90% by default). 

If the input database contains different alelles of the same gene, srst2 can report just the best matching allele for that gene (much like with MLST we report the best matching allele for each locus in the scheme). This makes the output manageable, as you will get one column per gene/locus (e.g. blaCTX-M) which reports the specific allele that was called in each sample (e.g. blaCTX-M-15 in sample A, blaCTX-M-13 in sample B).

To do this properly, srst2 needs to know which of the reference sequences are alleles of the same gene. This is done by adhering to the following format in the naming of the sequences (i.e. the headers in the fasta sequence for the database):

>[clusterUniqueIdentifier]__[clusterSymbol]__[alleleSymbol]__[alleleUniqueIdentifier]

e.g. in the resistance gene database provided, the first entry is:

>344__blaOXA__blaOXA-181__1

Note these are separated by two underscores. The individual components are:

clusterUniqueIdentifier = 344;  unique identifier for this cluster (uniquely identifes the cluster)
clusterSymbol = blaOXA;  gene symbol for this cluster (may be shared by multiple clusters)
alleleSymbol = blaOXA-181;  full name of this allele
alleleUniqueIdentifier = 1;  uniquely identifies the sequence

Ideally the alleleSymbol would be unique (as it is in the reference.fasta file provided). However it doesn't have to be: if allele symbols are not unique, then srst2 will use the combination '[alleleSymbol]__[alleleUniqueIdentifier]' to  uniquely identify the sequence in the resulting reports, so that you can trace exactly which sequence was present in each sample.

Additional gene annotation can appear on the header line, after a space. This additional info will be printed in the full genes report, but not in the compiled results files.

e.g. for the blaOXA sequence above, the full header is actually:

>344__blaOXA__blaOXA-181__1 blaOXA-181_1_HM992946; HM992946; betalactamase

------------
If you want to use your own database of allele sequences, with the reporting behaviour described, you will need to assign your sequences to clusters and use the header format specified above.

To facilitate this, use the scripts provided in the database_clustering directory provided with srst2.

------------
You can also use unclustered sequences. This is perfectly fine for gene detection applications, where you have one representative allele sequence for each gene, and you simply want to know which samples contain a sequence that is similar to this one (e.g. detecting plasmid replicons, where there is one target sequence per replicon).

However, this won't work well for allele typing. If the sequence database contains multiple allele sequences for the same gene, then all of these that are covered above the length threshold (default 90%) will be reported in the output, which makes for messy reporting. If you do this, you would probably find it most useful to look at the full gene results table rather than looking at the compiled results output.

Output files
====

MLST results

If MLST sequences and profiles were provided, STs will be printed in tab-delim format to a file called "mlst__[sequenceFileName]__[sample]__results.txt", e.g.: "mlst__Escherichia_coli__strainA__results.txt".

The format looks like this:

Sample  ST      adk     fumC    gyrB    icd     mdh     purA    recA    mismatches      uncertainty     depth

strainA     1502    6       63      7       1       14      7       7                       12.3771855735

- Each locus has a column in which the best scoring allele number is printed. * indicates the best scoring allele has >=1 mismatche (SNP or indel, according to majority rules consensus of the aligned reads vs the allele sequence). Details of the mismatches are given in the mismatches column. This often means you have a novel allele. ? indicates uncertainty in the result because the best scoring allele has some low-depth bases; either the the first or last 2 bases of the allele had <N reads mapped, or a truncation was called in which neigbhbouring bases were coverd with <N reads, or the average depth across the whole allele was <X. N is set by the parameter --min_edge_depth (default 2), X is set by --min_depth (default 5). The source of the uncertainty is printed to the uncertainty column. If no allele could be assigned, the entry is set to ‘-’. 

- If the combination of these alleles appears in the ST definitions file provided by --mlst_definitions, this ST will be printed in the ST column. "NF" indicates the allele combination was not found; "ND" indicates ST calculations were not done (because no ST definitions were provided). Here, * next to the ST indicates that there were mismatches against at least one of the alleles. This suggests that you have a novel variant of this ST rather than a precise match to this ST. ? indicates that there was uncertainty in at least one of the alleles. In all cases, the ST is calculated using the best scoring alleles, whether or not there are mismatches or uncertainty in those calls.

- The 'mismatches' column gives details of any mismatches (defined by majority rules consensus of the aligned reads vs the allele sequence) against the top scoring allele that is reported in the corresponding locus column. Possibilities are: (i) snps, adk-1/1snp indicates there was 1 SNP against the adk-1 allele; (ii) indels, adk-1/2indel indicates there were 1 indels (insertion or deletion calls in the alignment); (iii) holes, adk-1/5holes indicates there were 5 sections of the allele sequence that were not covered in the alignment and pileup (e.g. due to truncation at the start or end of the gene, or large deletions within the gene).

- The 'uncertainty' column gives details of parts of the top scoring alleles for which the depth of coverage was too low to give confidence in the result, this may be zero or any number up to the specified cutoff levels set via --min_edge_depth and --min_depth. Possibilities are considered in this order: (i) edge depth, adk-1/edge1.0 indicates that the mean read depth across either the first 2 or last 2 bases of the assigned allele was 1.0; this is monitored particularly because coverage at the ends of the allele sequences is dependent on bowtie2 to properly map reads that overhang the ends of the allele sequences, which is not as confident as when the whole length of a read maps within the gene (reported if this value is below the cutoff specified (default --min_edge_depth 2), low values can be interpreted as indicating uncertainty in the result as we can’t confidently distinguish alleles that differ at these low-covered bases); (ii) truncations, adk-1/del1.0 indicates that a truncation or large deletion was called for which the neighbouring 2 bases were covered to depth 1.0, this can be interpreted as indicating there is only very weak evidence for the deletion, as it is likely just due to random decline in coverage at this point (reported if this value is below the cutoff specified (default --min_edge_depth 2); (iii) average depth, adk-1/depth3.5 indicates that the mean read depth across the length of the assigned allele was 3.5 (reported if this value is below the cutoff specified, which by default is --min_depth 5)

- The 'depth' column indicates the mean read depth across the length of all alleles which were assigned a top scoring allele number (i.e. excluding any which are recorded as '-'). So if there are 7 loci with alleles called, this number represents the mean read depth across those 7 loci. If say, 2 of the 7 alleles were not called (recorded as ‘-’), the mean depth is that of the 5 loci that were called.

------------
Gene typing results files report the details of sequences provided in fasta files via --genes_db that are detected above the minimum %coverage theshold set by --min_coverage (default 90).

Two output files are produced:

1. A detailed report, fullgenes__[db]__[sample]__results.txt, with one row per gene per sample:

Sample  DB      gene    allele  coverage        depth   diffs   uncertainty     cluster seqid   annotation

strainA     resistance      dfrA    dfrA1_1 100.0   6.79368421053           edge1.5 590     137     

strainA     resistance      aadA    aadA1-5 100.0   10.6303797468                   162     1631    

strainA     resistance      sul2    sul2_9  100.0   79.01992966                     265     1763    

strainA     resistance      blaTEM  blaTEM-1_5      100.0   70.8955916473                   258     1396    

strainA     resistance      tet(A)  tet(A)_4        97.6666666667   83.5831202046   28holes edge0.0 76      1208    

strainB     resistance      strB    strB1   100.0   90.0883054893                   282     1720    

strainB     resistance      strA    strA4   100.0   99.0832298137                   325     1142 

- coverage indicates the % of the gene length that was covered (if clustered DB, then this is the highest coverage of any members of the cluster)
- uncertainty is as above

2. A tabulated summary report of samples x genes, genes__[db]__[sample]__results.txt:

Sample  aadA    blaTEM  dfrA    strA    strB    sul2    tet(A)

strainA     aadA1-5 blaTEM-1_5      dfrA1_1?        -   -   sul2_9  tet(A)_4*?

strainB     -     -      -        strA4   strB1   -  -

The first column indicates the sample name, all other columns report the genes/alleles that were detected in the sample set. If multiple samples were input, or if previous outputs were provided for compiling results, then all the genes detected in ANY of the samples will have their own column in this table.

If you were using a clustered gene database (such as the resistance.fasta database provided), the name of each cluster (i.e. the basic gene symbol) will be printed in the column headers, while specific alleles will be printed in the sample rows.

- * indicates mismatches
- ? indicates uncertainty due to low depth in some parts of the gene
- - indicates the gene was not detected (> %coverage threshold, --min_coverage 90)

------------
Combined results

If more then one database is provided for typing (via --mlst_db and/or --gene_db), or if previous results are provided for merging with the current run which contain data from >1 database (via --prev_output), then an additional table summarizing all the database results is produced. This is named "[sample]__compiledResults.txt" and is a combination of the MLST style table plus the tabulated gene summary (file 2 above).

Sample  ST      adk     fumC    gyrB    icd     mdh     purA    recA    mismatches      uncertainty     depth   aadA    blaTEM  dfrA    strA    strB    sul2    tet(A)
sampleA     152*     11      63*      7       1       14      7       7                       21.3139900892   aadA1-5 blaTEM-1_5      dfrA1_1?        strA4   strB1   sul2_9  tet(A)_4*?

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
