# SRST2

Short Read Sequence Typing for Bacterial Pathogens

This program is designed to take Illumina sequence data, a MLST database and/or a database
of gene sequences (e.g. resistance genes, virulence genes, etc) and report the presence of
STs and/or reference genes.

Authors - Michael Inouye, Harriet Dashnow, Bernie Pope, Ryan Wick, Kathryn Holt (University of Melbourne)
		
How to cite - The peer-reviewed open-access paper is available in Genome Medicine: http://genomemedicine.com/content/6/11/90

Story-behind-the-paper is [here](http://holtlab.net/2014/12/27/behind-the-paper-srst2-for-short-read-sequence-typing-of-bacterial-pathogens/)
		
Problems? Please post an issue here in github: https://github.com/katholt/srst2/issues.

To be notifed of updates, join the SRST2 google group at https://groups.google.com/forum/#!forum/srst2.

Contents
----
[Current release](https://github.com/katholt/srst2#current-release)

[Installation](https://github.com/katholt/srst2#installation)

[Basic usage - MLST](https://github.com/katholt/srst2#basic-usage---mlst)

[Basic usage - Resistance genes](https://github.com/katholt/srst2#basic-usage---resistance-genes)

[All usage options](https://github.com/katholt/srst2#all-usage-options)

[Input read formats and options](https://github.com/katholt/srst2#input-read-formats-and-options)

[MLST Database format](https://github.com/katholt/srst2#mlst-database-format)

[Gene databases](https://github.com/katholt/srst2#gene-databases)

[Output files](https://github.com/katholt/srst2#output-files)

* [MLST results](https://github.com/katholt/srst2#mlst-results)

* [Gene typing](https://github.com/katholt/srst2#gene-typing)

* [Combined results](https://github.com/katholt/srst2#combined-results)

* [Mapping results](https://github.com/katholt/srst2#mapping-results)

[Printing consensus sequences](https://github.com/katholt/srst2#printing-consensus-sequences)

[More basic usage examples](https://github.com/katholt/srst2#more-basic-usage-examples)

[Compile results from completed runs](https://github.com/katholt/srst2#compile-results-from-completed-runs)

[Running lots of jobs and compiling results](https://github.com/katholt/srst2#running-lots-of-jobs-and-compiling-results)

[Known issues](https://github.com/katholt/srst2#known-issues)

[Generating SRST2-compatible clustered database from raw sequences](https://github.com/katholt/srst2#generating-srst2-compatible-clustered-database-from-raw-sequences)

* [Using the VFBD Virulence Factor Database with SRST2](https://github.com/katholt/srst2#using-the-vfdb-virulence-factor-database-with-srst2)

* [Using the EcOH database for serotyping E. coli with SRST2](https://github.com/katholt/srst2#using-the-ecoh-database-for-serotyping-e-coli-with-srst2)

[Typing the LEE pathogenicity island of E. coli](https://github.com/katholt/srst2#typing-the-lee-pathogenicity-island-of-e-coli)

[Plotting output in R](https://github.com/katholt/srst2#plotting-output-in-r)

[Example - Shigella sonnei public data](example.txt)


Current release - v0.2.0 - July 28, 2016
-----

Dependencies:
* python (v2.7.5 or later)
* scipy, numpy   http://www.scipy.org/install.html
* bowtie2 (v2.1.0 or later)   http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
* SAMtools v0.1.18   https://sourceforge.net/projects/samtools/files/samtools/0.1.18/ (NOTE: later versions can be used, but better results are obtained with v0.1.18, especially at low read depths (<20x))

**NOTE AMR gene database has been updated since this release to include all the new mcr genes and correct allele nomenclature according to this [Partridge et al, JAC 2018](https://academic.oup.com/jac/article/73/10/2625/5055843), you can download the current one (ARGannot_r3.fasta) from the data directory

-----------

Updates in v0.2.0

1. Some improvements to allele calling, particularly for Klebsiella MLST locus mdh, kindly contributed by [andreyto](https://github.com/andreyto). Includes rejection of read alignments that are clipped on both ends (likely to be spurious) and minor bug fixes associated with depth calculations.
2. Updated E. coli serotype database to remove duplicate sequences.
3. Added mcr-2 colistin resistance gene to `ARGannot.r1.fasta` resistance gene database.
4. A `--threads` option was added, which makes SRST2 call Bowtie and Samtools with their threading options. The resulting speed up is mostly due to the Bowtie mapping step which parallelises very well.
5. The `VFDB_cdhit_to_csv.py` script was updated to work with the new VFDB FASTA format.
6. Versions of Bowtie2 up to 2.2.9 are now supported. Samtools v1.3 can now be used as well, however v0.1.18 is still the recommended version (for reasons discussed below).
7. Added `scripts/qsub_srst2.py` to generate SRST2 jobs for the Grid Engine (qsub) scheduling system (http://gridscheduler.sourceforge.net/). Thanks to Ramon Fallon from the University of St Andrews for putting this together. Some of the specifics are set up for his cluster, so modifications may be necessary to make it run properly on a different cluster using Grid Engine.
8. Various other small bug fixes!

-----------

Updates in v0.1.8

1. /data directory includes files for subtyping of the LEE pathogenicity island of E. coli, as per [Ingle et al, 2016, Nature Microbiology](http://www.nature.com/articles/nmicrobiol201510). [Instructions below](https://github.com/katholt/srst2#typing-the-lee-pathogenicity-island-of-e-coli)
2. Resistance gene database updates:
  * Fixed `ARGannot.r1.fasta` to include proper mcr1 DNA sequence.
  * Added columns to the `ARGannot_clustered80.csv` table, to indicate classes of beta-lactamases included in the `ARGannot.r1.fasta` database according to the [NCBI beta-lactamase resource](http://www.ncbi.nlm.nih.gov/pathogens/beta-lactamase-data-resources/) (new location for the Lahey list).
3. Fixed some issues with handling of missing data (i.e. where there were no hits to MLST and/or no hits to genes) when compiling results into a table via `--prev_output`. This could result in misalignment of gene columns in previous versions.

-----------

Updates in v0.1.7

1. Use the following environment variables to specify your prefered samtools and bowtie2 executables (thanks to Ben Taylor for this):
  * SRST2_SAMTOOLS
  * SRST2_BOWTIE2
  * SRST2_BOWTIE2_BUILD
2. Added mcr1, the plasmid-borne colisting resistance gene to the included ARG-Annot-based resistance gene DB (`ARGannot.r1.fasta`)
3. Fixed a problem with writing consensus files that occurred when a directory structure was specified using `--output` (bug introduced in v0.1.6)

-----------

Updates in v0.1.6

1. The original validation of SRST2 (see [paper](http://genomemedicine.com/content/6/11/90)) was performed with bowtie2 version 2.1.0 and samtools v0.1.18.
  * bowtie2: SRST2 has now been tested on the tutorial example and other test data sets using the latest versions of bowtie2, 2.2.3 and 2.2.4, which gave identical results to those obtained with bowtie2 v2.1.0. Therefore, the SRST2 code will now run if any of these versions of bowtie2 are available: 2.1.0, 2.2.3 or 2.2.4. 
  * samtools: SRST2 has now been tested on the Staph & Salmonella test data sets used in the paper, and will work with newer samtools versions (tested up to v1.1). Note however that SRST2 still works best with [samtools v0.1.18](https://sourceforge.net/projects/samtools/files/samtools/0.1.18/), due to small changes in the mapping algorithms in later versions that result in some loss of reads at the ends of alleles. This has most impact at low read depths, however we do recommend using v0.1.18 for optimum results.
2. Minor fixes to the ARG-Annot database of resistance genes, including removal of duplicate sequences and fixes to gene names (thanks to Wan Yu for this). Old version remains unchanged for backwards compatibility, but we recommend using the revised version (located in `data/ARGannot.r1.fasta`).
3. Added EcOH database for serotyping E. coli (thanks to Danielle Ingle for this). See [Using the EcOH database for serotyping E. coli with SRST2](https://github.com/katholt/srst2#using-the-ecoh-database-for-serotyping-e-coli-with-srst2) and [this MGen paper](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000064).
4. Fixed a problem where, when analysing multiple read sets in one SRST2 call against a gene database in which cluster ids don't match gene symbols, individual gene clusters appear multiple times in the output. The compile function was unaffected and remains unchanged.
5. Fixed behaviour so that including directory paths in `--output` parameter works (thanks to nyunyun for contributing most of this fix). E.g. `--output test_dir/test` will create output files prefixed with `test`, located in `test_dir/`, and all SRST2 functions should work correctly including consensus allele calling. If `test_dir/` doesn't exist, we attempt to create it; if this is not possible the user is alerted and SRST2 stops.
6. Fixed problem when using a gene database with a simple fasta header (ie not clustered for SRST2; note best results are achieved by pre-clusering your sequence database beforehand) (thanks to cglambert for this one).
7. Fixes contributed by ppcherng (thanks!): 
  * Fixed KeyErrors that occured when a given seqID was not found in the seq2cluster dictionary, which tended to happen if the FASTA file (gene database) contained empty entries that only have a header and no sequence.
  * Note v0.1.5 included addition of ppcherng's utility scripts to help automate creation of SRST2-compatible gene databases from VFDB.
8. Added new parameter `--samtools_args` to pass additional options to samtools mpileup (e.g. SionBayliss requested this in order to use `-A` option in samtools mpileup to include anomalous reads).
9. Fixed problem with consensus sequence reporting of truncated alleles (issue #39).
10. Added basic instructions for the R scripts provided for plotting output. See [Plotting output in R](https://github.com/katholt/srst2#plotting-output-in-r)

-----------

Updates in v0.1.5

1. Optionally switch on reporting of pileups and consensus sequences (fasta) for novel alleles (`--report_new_consensus`) or for all alleles (`--report_all_consensus`). See [Printing consensus sequences](https://github.com/katholt/srst2#printing-consensus-sequences)
2. Post-process consensus sequences from a set of strains, to generate one file per locus containing all/new consensus sequences. See [Collate consensus sequences](https://github.com/katholt/srst2/blob/master/README.md#collate-consensus-sequences-output-by-srst2-run-on-multiple-strains--loci-into-one-file-per-locus)
3. Some enhancements to getmlst.py script to handle some more unusual scheme names (force download of specific schemes that have non-unique  names, handle forward slashes in names).
4. Fixed an issue where, if multiple readsets analysed in serial in a SRST2 run, the fullgenes report would only contain the results for the last readset. Fullgenes report now contains gene output for all readsets.
5. Added option (`--merge_paired`) to accommodate cases where users have multiple read sets for the same sample. If this flag is used, SRST2 will assume that all the input reads belong to the same sample, and outputs will be named as `[prefix]__combined.xxx`, where SRST2 was run using `--output [prefix]`. If the flag is not used, SRST2 will operate as usual and assume that each read pair is a new sample, with output files named as `[prefix]__[sample].xxx`, where [sample] is taken from the base name of the reads fastq files. Note that if you have lots of multi-run read sets to analyse, the ease of job submission will depend heavily on how your files are named and you will need to figure out your own approach to manage this (ie there is no way to submit multiple sets of multiple reads).
6. The original validation of SRST2 (see [paper](http://genomemedicine.com/content/6/11/90)) was performed with bowtie2 version 2.1.0. SRST2 has now been tested on the tutorial example using the latest versions of bowtie2, 2.2.3 and 2.2.4, which gave identical results to those obtained with bowtie2 v2.1.0. Therefore, the SRST2 code will now run if any of these versions of bowtie2 are available: 2.1.0, 2.2.3 or 2.2.4. (Note however that there are still incompatibilities with the recent release of samtools, so you will need to stick to [samtools v0.1.18](https://sourceforge.net/projects/samtools/files/samtools/0.1.18/) unless you want to modify the SRST2 code to allow later versions, and are happy with a dramatic loss in accuracy!)
7. Thanks to ppcherng for adding utility scripts to help automate creation of SRST2-compatible gene databases from VFDB.

-----------

Updates in v0.1.4

1. No longer store sam and unsorted bam (can be retained via the `--keep_interim_alignment` flag)
2. Added options to specify a maximum number of mismatches to allow during mapping; this is specified separately for mlst and genes, so that it is possible to relax the stringency of gene detection in the same run as a high-accuracy MLST test.
  * `--mlst_max_mismatch`
  * `--gene_max_mismatch`
  * Default value for both is 10 mismatches.
3. The highest minor allele frequency (MAF) of variants encountered in the alignment is now calculated and reported for each allele (in the scores file) and also at the gene level and ST level, to facilitate checking for mixed/contaminated read sets. 
This value is in the range 0 -> 0.5; with e.g. 0 indicating no variation between reads at any aligned base (i.e. at all positions in the alignment, all aligned reads agree on the same base call; although this agreed base may be different from the reference); and 0.25 indicating there is at least one position in the alignment at which all reads do not agree, and the least common variant (either match or mismatch to the reference) is present in 25% of reads. This value is printed, for all alleles, to the scores file. Note this is different to the ‘LeastConfident’ information printed to scores, which presents the strongest evidence for mismatch *compared to the reference*, i.e. between 0 -> 1.
The highest such value for each gene/cluster/locus is reported in the fullgenes output table.
The highest such value across all MLST loci is reported in the mlst output table.
Note that all compiled reports will now include a maxMAF column; if you provide MLST or compiled reports from previous versions without this columns, the value “NC” will be inserted in the maxMAF column to indicate “not calculated”. This ensures the updated SRST2 (v0.1.4+) is backwards compatible with previous SRST2 outputs; do be aware though that the older versions of SRST2 (<v0.1.4) will not be forwards-compatible with output generated by more recent versions (v0.14 onwards).
4. Added R code for plotting SRST2 output in R (plotSRST2data.R). Instructions will be added to the read me.
5. Added formatted versions of the ARG-Annot resistance gene database, PlasmidFinder database and 18 plasmid replicon sequences to the /data directory. See /data/README.md for details and citations. It is recommended to use `ARGannot.r1.fasta` for detection of acquired resistance genes.

-----------
 
Updates in v0.1.3

1. Fixed a bug that occurred while trying to type genes from a user-supplied database (see issue #5, thanks to Scott Long)
2. Fixed a bug in gene detection reporting - genes are now correctly reported by cluster, rather than by gene symbol (see issue #7)
3. Added maximum divergence option for reporting (`--max_divergence`), default is now to report only hits with <10% divergence from the database (see issue #8)
4. added parameter to pass to bowtie2 parameter `-u N` to stop mapping after the first N reads. Default behaviour remains to map all reads. However, for large read sets (e.g. >100x), extra reads do not help and merely increase the time taken for mapping and scoring, and you may want to limit to the first million reads (100x of a 2 Mbp genome) using `--stop_after 1000000`.


# Installation

### 1 - Install dependencies

* python (v2.7.5)
* scipy http://www.scipy.org/install.html
* bowtie2 v2.1.0 http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
* SAMtools v0.1.18 https://sourceforge.net/projects/samtools/files/samtools/0.1.18/ (NOTE 0.1.19 DOES NOT WORK)

N.B. If you have multiple versions of samtools or bowtie2 installed, you can pick which one `srst2` or `slurm_srst2` should use by setting the following environment variables.

* `SRST2_SAMTOOLS`
* `SRST2_BOWTIE2`
* `SRST2_BOWTIE2_BUILD`

If these aren't set or are missing, they will default to looking in your `PATH` for `samtools`, `bowtie2` and `bowtie2-build`. The exception is `SRST2_BOWTIE2_BUILD` which, if it is not set or missing, will try adding `-build` to `SRST2_BOWTIE2` if it exists, otherwise it defaults to looking in your `PATH`


### 2 - Get and install the code

Make sure you have installed [git](https://help.github.com/articles/set-up-git) and [pip](http://www.pip-installer.org/).

Clone the git repository: `git clone https://github.com/katholt/srst2`

and then install with pip: `pip install srst2/`

OR do both at once: `sudo pip install git+https://github.com/katholt/srst2`

### 3 - Test that the programs are installed properly

```
srst2 --version
getmlst.py -h
slurm_srst2.py -h
```

The downloaded directory also contains things that might be useful for SRST2 users:

* `data/` contains a databases for resistance genes and plasmids
* `database_clustering/` contains scripts and instructions for formatting of other gene databases for use with SRST2
* [example.txt](example.txt) contains a tutorial/example on running SRST2 using public data 



# Basic usage - MLST

### 1 - Gather your input files

(i) sequence reads (this example uses paired reads in gzipped fastq format, see below for options)

(ii) a fasta sequence database to match to. For MLST, this means a fasta file of all allele sequences. If you want to assign STs, you also need a tab-delim file which defines the ST profiles as a combination of alleles. You can retrieve these files automatically from pubmlst.org/data/ using the script provided:

```
getmlst.py --species 'Escherichia coli#1'
```

### 2 - Run MLST

```
srst2 --input_pe strainA_1.fastq.gz strainA_2.fastq.gz --output strainA_test --log --mlst_db Escherichia_coli#1.fasta --mlst_definitions profiles_csv --mlst_delimiter _
```

### 3 - Check the outputs

(i) MLST results are output in: `strainA_test__mlst__Escherichia_coli#1__results.txt`

Sample | ST | adk | fumC | gyrB | icd | mdh | purA | recA | mismatches | uncertainty | depth | maxMAF
:---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---:
strainA | 152 | 11 | 63 | 7 | 1 | 14 | 7 | 7 | 0 | - | 25.8319955826 | 0.125


# Basic usage - Resistance genes 

### 1 - Gather your input files

(i) sequence reads (this example uses paired reads in gzipped fastq format, see below for options)

(ii) a fasta sequence database to match to. For resistance genes, this means a fasta file of all the resistance genes/alleles that you want to screen for, clustered into gene groups. Some suitable databases are distributed with SRST2 (in the /data directory); we recommend using `/data/ARGannot_r3.fasta` for acquired resistance genes (note this latest version has been added since the last SRST2 release so you may need to download it directly from the /data directory in this github repository).

### 2 - Run gene detection

```
srst2 --input_pe strainA_1.fastq.gz strainA_2.fastq.gz --output strainA_test --log --gene_db resistance.fasta
```

### 3 - Check the outputs

(i) Gene detection results are output in: `strainA_test__genes__resistance__results.txt`

Sample | aadA | dfrA | sul2 | tet(B)
:---: | :---: | :---: | :---: | :---:
strainA | aadA1-5 | dfrA1_1 | sul2_2 | tet(B)_4


# All usage options

```
srst2 -h

SRST2 - Short Read Sequence Typer (v2)

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  --input_se INPUT_SE [INPUT_SE ...]
                        Single end read file(s) for analysing (may be gzipped)
  --input_pe INPUT_PE [INPUT_PE ...]
                        Paired end read files for analysing (may be gzipped)
  --merge_paired        Switch on if all the input read sets belong to a
                        single sample, and you want to merge their data to get
                        a single result
  --forward FORWARD     Designator for forward reads (only used if NOT in
                        MiSeq format sample_S1_L001_R1_001.fastq.gz; otherwise
                        default is _1, i.e. expect forward reads as
                        sample_1.fastq.gz)
  --reverse REVERSE     Designator for reverse reads (only used if NOT in
                        MiSeq format sample_S1_L001_R2_001.fastq.gz; otherwise
                        default is _2, i.e. expect forward reads as
                        sample_2.fastq.gz
  --read_type {q,qseq,f}
                        Read file type (for bowtie2; default is q=fastq; other
                        options: qseq=solexa, f=fasta).
  --mlst_db MLST_DB     Fasta file of MLST alleles (optional)
  --mlst_delimiter MLST_DELIMITER
                        Character(s) separating gene name from allele number
                        in MLST database (default "-", as in arcc-1)
  --mlst_definitions MLST_DEFINITIONS
                        ST definitions for MLST scheme (required if mlst_db
                        supplied and you want to calculate STs)
  --mlst_max_mismatch MLST_MAX_MISMATCH
                        Maximum number of mismatches per read for MLST allele
                        calling (default 10)
  --gene_db GENE_DB [GENE_DB ...]
                        Fasta file/s for gene databases (optional)
  --no_gene_details     Switch OFF verbose reporting of gene typing
  --gene_max_mismatch GENE_MAX_MISMATCH
                        Maximum number of mismatches per read for gene
                        detection and allele calling (default 10)
  --min_coverage MIN_COVERAGE
                        Minimum %coverage cutoff for gene reporting (default
                        90)
  --max_divergence MAX_DIVERGENCE
                        Maximum %divergence cutoff for gene reporting (default
                        10)
  --min_depth MIN_DEPTH
                        Minimum mean depth to flag as dubious allele call
                        (default 5)
  --min_edge_depth MIN_EDGE_DEPTH
                        Minimum edge depth to flag as dubious allele call
                        (default 2)
  --prob_err PROB_ERR   Probability of sequencing error (default 0.01)
  --truncation_score_tolerance TRUNCATION_SCORE_TOLERANCE
                        % increase in score allowed to choose non-truncated
                        allele
  --stop_after STOP_AFTER
                        Stop mapping after this number of reads have been
                        mapped (otherwise map all)
  --other OTHER         Other arguments to pass to bowtie2 (must be escaped,
                        e.g. "\--no-mixed".
  --max_unaligned_overlap MAX_UNALIGNED_OVERLAP
                        Read discarded from alignment if either of its ends
                        has unaligned overlap with the reference that is
                        longer than this value (default 10)
  --mapq MAPQ           Samtools -q parameter (default 1)
  --baseq BASEQ         Samtools -Q parameter (default 20)
  --samtools_args SAMTOOLS_ARGS
                        Other arguments to pass to samtools mpileup (must be
                        escaped, e.g. "\-A").
  --output OUTPUT       Prefix for srst2 output files
  --log                 Switch ON logging to file (otherwise log to stdout)
  --save_scores         Switch ON verbose reporting of all scores
  --report_new_consensus
                        If a matching alleles is not found, report the
                        consensus allele. Note, only SNP differences are
                        considered, not indels.
  --report_all_consensus
                        Report the consensus allele for the most likely
                        allele. Note, only SNP differences are considered, not
                        indels.
  --use_existing_bowtie2_sam
                        Use existing SAM file generated by Bowtie2 if
                        available, otherwise they will be generated
  --use_existing_pileup
                        Use existing pileups if available, otherwise they will
                        be generated
  --use_existing_scores
                        Use existing scores files if available, otherwise they
                        will be generated
  --keep_interim_alignment
                        Keep interim files (sam & unsorted bam), otherwise
                        they will be deleted after sorted bam is created
  --threads THREADS     Use multiple threads in Bowtie and Samtools
  --prev_output PREV_OUTPUT [PREV_OUTPUT ...]
                        SRST2 results files to compile (any new results from
                        this run will also be incorporated)
```

# Input read formats and options

Any number of readsets can be provided using `--input_se` (for single end reads) and/or `--input_pe` (for paired end reads). You can provide both types of reads at once. Note however that if you do this, each readset will be typed one by one (in serial). So if it takes 2 minutes to type each read set, and you provide 10 read sets, it will take 20 minutes to get your results. The better way to proces lots of samples quickly is to give each one its own SRST2 job (e.g. submitted simultaneously to your job scheduler or server); then compile the results into a single report using `srst2 --prev_output *results.txt --output all`. That way each readset's 2 minutes of analysis is occurring in parallel on different nodes, and you'll get your results for all 10 samples in 2 minutes rather than 20.

### Read formats 
Reads can be in any format readable by bowtie2. The format is passed on to the bowtie2 command via the `--read_type` flag in SRST2. The default format is fastq (passed to bowtie 2 as q); other options are qseq=solexa, f=fasta. So to use fasta reads, you would need to tell SRST2 this via `--read_type f`.

Reads may be gzipped.

### Read names 
SRST2 can parse Illumina MiSeq reads files; we assume that files with names in the format `XXX_S1_L001_R1_001.fastq.gz` and `XXX_S1_L001_R2_001.fastq.gz` are the forward and reverse reads from a sample named "XXX". So, you can simply use `srst2 --input_pe XXX_S1_L001_R1_001.fastq.gz XXX_S1_L001_R2_001.fastq.gz` and SRST2 will recognise these as forward and reverse reads of a sample named "XXX". If you have single rather than paired MiSeq reads, you would use `srst2 --input_se XXX_S1_L001_R1_001.fastq.gz`.

### Paired reads
If you have paired reads that are named in some way other than the Illumina MiSeq format, e.g. from the SRA or ENA public databases, you need to tell SRST2 how to pass these to bowtie2.
bowtie2 requires forward and reverse reads to be supplied in separate files, e.g `strainA_1.fastq.gz` and `strainA_2.fastq.gz`. SRST2 attempts to sort reads supplied via `--input_pe` into read pairs, based on the suffix (_1, _2 in this example) that occurs before the file extension (.fastq.gz in this example). So if you supplied `--input_pe strainA_1.fastq.gz strainB_1.fastq.gz strainA_2.fastq.gz strainB_2.fastq.gz`, SRST2 would sort these into two pairs (`strainA_1.fastq.gz`, `strainA_2.fastq.gz`) and (`strainB_1.fastq.gz`, `strainB_2.fastq.gz`) and pass each pair on to bowtie2 for mapping. By default, the suffixes are assumed to be `_1` for forward reads and `_2` for reverse reads, but you can tell SRST2 if you have other conventions, via `--forward` and `--reverse`. E.g. if your files were named `strainA_read1.fastq.gz` and `strainA_read2.fastq.gz`, you would use these commands: `--input_pe strainA_read1.fastq.gz strainA_read2.fastq.gz --forward _read1 --reverse _read2`. 

### Sample names
Sample names are taken from the first part of the read file name (before the suffix if you have paired reads). E.g. `strainA_1.fastq.gz` is assumed to belong to a sample called "strainA"; `strainB_C_1.fastq.gz` would be assumed to belong to a sample called "strainB_C". These sample names will be used to name all output files, and will appear in the results files.

### If you have multiple read sets for the same sample
The flag `--merge_paired` tells SRST2 to assume that all the input reads belong to the same sample. Outputs will be named as `[prefix]__combined.xxx`, where SRST2 was run using `--output [prefix]`. If this flag is not used, SRST2 will operate as usual and assume that each read pair is a new sample, with output files named as `[prefix]__[sample].xxx`, where [sample] is taken from the base name of the reads fastq files. Note that if you have lots of multi-run read sets to analyse, the ease of job submission will depend heavily on how your files are named and you will need to figure out your own approach to manage this (ie there is no way to submit multiple sets of multiple reads).


# MLST Database format

MLST databases are specified by allele sequences, and a profiles table. These can be downloaded from the public databases, ready to use with SRST2, using the provided script getmlst.py (see above).

### Allele sequences file, fasta format. 

`--mlst_db alleles.fasta`

This should contain ALL allele sequences for the MLST scheme; i.e. if there are 7 loci in the scheme, then the sequences of all alleles from all 7 loci should appear together in this file. If you have one file per locus, just cat these together first. If you use getmlst.py, this is done for you.

`--mlst_delimiter '-'`

The names of the alleles (i.e. the fasta headers) are critical for a functioning MLST scheme and therefore for correct calling of STs. There are two key components to every fasta header: the name of the locus (e.g. in E. coli these are adk, fumC, gyrB, icd, mdh, purA, recA) and the number assigned to the allele (1, 2, 3 etc). These are usually separated by a delimiter like '-' or '_'; e.g. in the E. coli scheme, alleles are named adk-1, fumC-10, etc. For ST calling to work properly, SRST2 needs to know what the delimiter is. By default, we assume it is '-' as this is the most common; however some schemes use '_' (e.g. in the C. difficile scheme, the first allele is 'adk_1', so you would need to set `--mlst_delimiter '_'` on the SRST2 command line). If you use getmlst.py, it will remind you of this and try to guess for you what the most likely delimiter is.

### MLST definitions file, tab delimited format.

`--mlst_definitions`

This is the file that tells you the ST number that is assigned to known combinations of alleles. Column 1 is the ST, and subsequent columns are the loci that make up the scheme. The names of these loci must match the allele names in the sequences database, e.g. adk, fumC, gyrB, icd, mdh, purA, recA in the E. coli scheme. If you download data from pubmlst this should not be a problem. Sometimes there are additional columns in this file, e.g. a column assigning STs to clonal complexes. SRST2 will ignore any columns that don't match gene names found in the allele sequences file.

# Gene databases

In addition to MLST, SRST2 can do gene/allele detection. This works by mapping reads to each of the reference sequences in a fasta file(s) (provided through `--gene_db`) and reporting details of all genes that are covered above a certain level (`--min_coverage`, 90% by default). 

If the input database contains different alelles of the same gene, SRST2 can report just the best matching allele for that gene (much like with MLST we report the best matching allele for each locus in the scheme). This makes the output manageable, as you will get one column per gene/locus (e.g. blaCTX-M) which reports the specific allele that was called in each sample (e.g. blaCTX-M-15 in sample A, blaCTX-M-13 in sample B).

We have provided some databases of resistance genes, plasmid genes and [E. coli serotyping loci](https://github.com/katholt/srst2#using-the-ecoh-database-for-serotyping-e-coli-with-srst2) in /data, ready for use with SRST2. We recommend using `data/ARGannot_r3.fasta` for detecting resistance genes, but you can also use `data/ResFinder.fasta` (this is the same as `data/resistance.fasta` in earlier distributions of SRST2).

You can however format any sequence set for screening with SRST2. [See instructions below](https://github.com/katholt/srst2#generating-srst2-compatible-clustered-database-from-raw-sequences).

# Output files

### MLST results

If MLST sequences and profiles were provided, STs will be printed in tab-delim format to a file called `[outputprefix]__mlst__[db]__results.txt`, e.g.: `strainArun1__mlst__Escherichia_coli__results.txt`.

The format looks like this:

Sample | ST | adk | fumC | gyrB | icd | mdh | purA | recA | mismatches | uncertainty | depth | maxMAF
:---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---:
strainA | 1502 | 6 | 63 | 7 | 1 | 14 | 7 | 7 | 0 | - | 12.3771855735 | 0.275

Each locus has a column in which the best scoring allele number is printed. 

\* indicates the best scoring allele has >=1 mismatch (SNP or indel, according to majority rules consensus of the aligned reads vs the allele sequence). Details of the mismatches are given in the mismatches column. This often means you have a novel allele. 

? indicates uncertainty in the result because the best scoring allele has some low-depth bases; either the the first or last 2 bases of the allele had <N reads mapped, or a truncation was called in which neigbhbouring bases were coverd with <N reads, or the average depth across the whole allele was <X. N is set by the parameter `--min_edge_depth` (default 2), X is set by `--min_depth` (default 5). The source of the uncertainty is printed to the uncertainty column. 

\- indicates that no allele could be assigned (generally means there were no alleles that were >90% covered by reads) 

If the combination of these alleles appears in the ST definitions file provided by `--mlst_definitions`, this ST will be printed in the ST column. "NF" indicates the allele combination was not found; "ND" indicates ST calculations were not done (because no ST definitions were provided). Here, * next to the ST indicates that there were mismatches against at least one of the alleles. This suggests that you have a novel variant of this ST rather than a precise match to this ST. ? indicates that there was uncertainty in at least one of the alleles. In all cases, the ST is calculated using the best scoring alleles, whether or not there are mismatches or uncertainty in those calls.

The *mismatches* column gives details of any mismatches (defined by majority rules consensus of the aligned reads vs the allele sequence) against the top scoring allele that is reported in the corresponding locus column. Possibilities are: (i) snps, adk-1/1snp indicates there was 1 SNP against the adk-1 allele; (ii) indels, adk-1/2indel indicates there were 1 indels (insertion or deletion calls in the alignment); (iii) holes, adk-1/5holes indicates there were 5 sections of the allele sequence that were not covered in the alignment and pileup (e.g. due to truncation at the start or end of the gene, or large deletions within the gene).

The *uncertainty* column gives details of parts of the top scoring alleles for which the depth of coverage was too low to give confidence in the result, this may be zero or any number up to the specified cutoff levels set via `--min_edge_depth` and `--min_depth`. Possibilities are considered in this order: (i) edge depth, adk-1/edge1.0 indicates that the mean read depth across either the first 2 or last 2 bases of the assigned allele was 1.0; this is monitored particularly because coverage at the ends of the allele sequences is dependent on bowtie2 to properly map reads that overhang the ends of the allele sequences, which is not as confident as when the whole length of a read maps within the gene (reported if this value is below the cutoff specified (default `--min_edge_depth 2`), low values can be interpreted as indicating uncertainty in the result as we can’t confidently distinguish alleles that differ at these low-covered bases); (ii) truncations, adk-1/del1.0 indicates that a truncation or large deletion was called for which the neighbouring 2 bases were covered to depth 1.0, this can be interpreted as indicating there is only very weak evidence for the deletion, as it is likely just due to random decline in coverage at this point (reported if this value is below the cutoff specified (default `--min_edge_depth 2`); (iii) average depth, adk-1/depth3.5 indicates that the mean read depth across the length of the assigned allele was 3.5 (reported if this value is below the cutoff specified, which by default is `--min_depth 5`)

The *depth* column indicates the mean read depth across the length of all alleles which were assigned a top scoring allele number (i.e. excluding any which are recorded as '-'). So if there are 7 loci with alleles called, this number represents the mean read depth across those 7 loci. If say, 2 of the 7 alleles were not called (recorded as ‘-’), the mean depth is that of the 5 loci that were called.

The *maxMAF* column reports the highest minor allele frequency (MAF) of variants encountered across the MLST locus alignments. This value is in the range 0 -> 0.5; with e.g. 0 indicating no variation between reads at any aligned base (i.e. at all positions in the alignment, all aligned reads agree on the same base call; although this agreed base may be different from the reference); and 0.25 indicating there is at least one position in the alignment at which all reads do not agree, and the least common variant (either match or mismatch to the reference) is present in 25% of reads. This value is printed, for all alleles, to the scores file; in the MLST file, the value reported is the highest MAF encountered across the MLST loci.

------------
### Gene typing

Gene typing results files report the details of sequences provided in fasta files via `--genes_db` that are detected above the minimum %coverage theshold set by `--min_coverage` (default 90).

Two output files are produced:

1 - A detailed report, `[outputprefix]__fullgenes__[db]__results.txt`, with one row per gene per sample:

Sample | DB | gene | allele | coverage | depth | diffs | uncertainty | cluster | seqid | annotation
:---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---:
strainA | resistance | dfrA | dfrA1_1 | 100.0 | 6.79368421053 | | edge1.5 | 590 | 137 | 
strainA | resistance | aadA | aadA1-5 | 100.0 | 10.6303797468 | | | 162 | 1631 | 
strainA | resistance | sul2 | sul2_9 | 100.0 | 79.01992966 | | | 265 | 1763 | 
strainA | resistance | blaTEM | blaTEM-1_5 | 100.0 | 70.8955916473 | | | 258 | 1396 | 
strainA | resistance | tet(A) | tet(A)_4 | 97.6666666667 | 83.5831202046 | 28holes | edge0.0 | 76 | 1208 | 
strainB | resistance | strB | strB1 | 100.0 | 90.0883054893 | | | 282 | 1720 | 
strainB | resistance | strA | strA4 | 100.0 | 99.0832298137 | | | 325 | 1142 |

*coverage* indicates the % of the gene length that was covered (if clustered DB, then this is the highest coverage of any members of the cluster)

*uncertainty* is as above

2 - A tabulated summary report of samples x genes, `[outputprefix]__genes__[db]__results.txt`:

Sample | aadA | blaTEM | dfrA | strA | strB | sul2 | tet(A)
:---: | :---: | :---: | :---: | :---: | :---: | :---: | :---:
strainA | aadA1-5 | blaTEM-1_5 | dfrA1_1? | - | - | sul2_9 | tet(A)_4*?
strainB | - | - | - | strA4 | strB1 | - | -

The first column indicates the sample name, all other columns report the genes/alleles that were detected in the sample set. If multiple samples were input, or if previous outputs were provided for compiling results, then all the genes detected in ANY of the samples will have their own column in this table.

If you were using a clustered gene database (such as the `ARGannot_r3.fasta` database provided), the name of each cluster (i.e. the basic gene symbol) will be printed in the column headers, while specific alleles will be printed in the sample rows.

\* indicates mismatches

? indicates uncertainty due to low depth in some parts of the gene

\- indicates the gene was not detected (> %coverage threshold, `--min_coverage 90`)

------------
### Combined results

If more then one database is provided for typing (via `--mlst_db` and/or `--gene_db`), or if previous results are provided for merging with the current run which contain data from >1 database (via `--prev_output`), then an additional table summarizing all the database results is produced. This is named `[outputprefix]__compiledResults.txt` and is a combination of the MLST style table plus the tabulated gene summary (file 2 above).

Sample | ST | adk | fumC | gyrB | icd | mdh | purA | recA | mismatches | uncertainty | depth | maxMAF | aadA | blaTEM | dfrA | strA | strB | sul2 | tet(A)
:---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---:
sampleA | 152* | 11 | 63* | 7 | 1 | 14 | 7 | 7 | 0 | - | 21.3 | 0.05 | aadA1-5 | blaTEM-1_5 | dfrA1_1? | strA4 | strB1 | sul2_9 | tet(A)_4*?

------------
### Mapping results

The bowtie2 alignment of reads to each input database is stored in `[outputprefix]__[sample].[db].sorted.bam` and the samtools pileup of the alignment is stored in `[outputprefix]__[sample].[db].pileup`. 

If you used `--save_scores`, a table of scores for each allele in the database is printed to `[outputprefix]__[sample].[db].scores`.


# Printing consensus sequences

SRST2 can optionally report consensus sequences & pileups for novel alleles, or for all alleles. Note that only SNPs are included in the consensus fasta sequence as these can be reliably extracted from the read alignments; if there are indels (reported in the output tables) these will not be included in the consensus fasta file and we suggest you assembly the mapped reads (these are in the bam file).

By default, no allele sequences are generated, the results are simply tabulated.

### Report consensus sequences for novel alleles 

```
--report_new_consensus
```

For all samples and loci where the top scoring allele contains SNPs:

* a pileup file will be generated for the top scoring allele, with the name `[allele].[output]__[readset].[database].pileup`
* the consensus sequence will be printed to a fasta file with the name `[output].new_consensus_alleles.fasta`
* fasta headers will be in the format `>[allele].variant [sample]`


IN ADDITION TO THE NOVEL ALLELES FILE OUTLINED ABOVE, the following will ALSO occur:

* a pileup file will be generated for the top scoring allele for each sample at each locus, with the name `[allele].[output]__[readset].[database].pileup`
* the consensus sequence for the top scoring allele for each sample at each locus will be printed to a fasta file with the name `[output].all_consensus_alleles.fasta`
* fasta headers will be in the format `>[allele].variant [sample]`


### Report consensus sequences for all alleles

```
--report_all_consensus
```

### Collate consensus sequences output by SRST2 run on multiple strains & loci, into one file per locus

For genes:

```
python srst2/scripts/consensus_alignment.py --in *.all_consensus_alleles.fasta --pre test --type gene
```

For mlst:

```
python srst2/scripts/consensus_alignment.py --in *.all_consensus_alleles.fasta --pre test --type mlst --mlst_delimiter _
```

# More basic usage examples

Run single read sets against MLST database and resistance database

```
srst2 --input_pe pool11_tag2_1.fastq.gz pool11_tag2_2.fastq.gz 
	--output pool11_tag2_Shigella --log 
	--gene_db /vlsci/VR0082/shared/srst2_sep/resistance.fasta 
	--mlst_db Escherichia_coli.fasta 
	--mlst_definitions ecoli.txt
```

------------

Run against multiple read sets in serial

```
srst2 --input_pe *.fastq.gz
	--output Shigella --log
	--gene_db /vlsci/VR0082/shared/srst2_sep/resistance.fasta 
	--mlst_db Escherichia_coli.fasta 
	--mlst_definitions ecoli.txt
```

------------

Run against new read sets, merge with previous reports (individual or compiled)

```
	srst2 --input_pe strainsY-Z*.fastq.gz
		--output strainsA-Z --log
		--gene_db /vlsci/VR0082/shared/srst2_sep/resistance.fasta 
		--mlst_db Escherichia_coli.fasta 
		--mlst_definitions ecoli.txt
		--prev_output ShigellaA__genes__resistance__results.txt
        	 ShigellaA__mlst__Escherichia_coli__results.txt
        	 ShigellaB__genes__resistance__results.txt
        	 ShigellaB__mlst__Escherichia_coli__results.txt
        	 ShigellaC-X__compiledResults.txt
```
------------

Run against Enterococcus reads, where read names are different from the usual `_1.fastq` and `_2.fastq`

```
srst2 --input_pe strain_R1.fastq.gz strain_R2.fastq.gz 
	--forward _R1 --reverse _R2 
	--output strainA --log 
	--gene_db /vlsci/VR0082/shared/srst2_sep/resistance.fasta 
	--mlst_db Enterococcus_faecium.fasta 
	--mlst_definitions efaecium.txt
```
	
# Compile results from completed runs

```
srst2 --prev_output *compiledResults.txt --output Shigella_report
```


# Running lots of jobs and compiling results

Run against multiple read sets: submitting 1 job per readset to SLURM queueing system.

```
slurm_srst2.py --script srst2 
	--output test 
	--input_pe *.fastq.gz 
	--other_args '--gene_db resistance.fasta 
	--mlst_db Escherichia_coli.fasta 
	--mlst_definitions ecoli.txt 
	--save_scores' 
	--walltime 0-1:0 
		> job_sub_list.txt
```
		
The results from all the separate runs can then be compiled together using:		

```
srst2 --prev_output *compiledResults.txt --output Shigella_report
```

------------

SLURM job script usage options

```
slurm_srst2.py -h

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
                        Single end read file(s) for analysing (may be gzipped)
                        
  --input_pe INPUT_PE [INPUT_PE ...]
                        Paired end read files for analysing (may be gzipped)
                        
  --forward FORWARD     Designator for forward reads (only used if NOT in
                        MiSeq format sample_S1_L001_R1_001.fastq.gz; otherwise
                        default is _1, i.e. expect forward reads as
                        sample_1.fastq.gz)
                        
  --reverse REVERSE     Designator for reverse reads (only used if NOT in
                        MiSeq format sample_S1_L001_R2_001.fastq.gz; otherwise
                        default is _2, i.e. expect forward reads as
                        sample_2.fastq.gz)
                        
  --other_args OTHER_ARGS
                        string containing all other arguments to pass to srst2
```

# Known issues

Reference indexing - SRST2 uses bowtie2 for mapping reads to reference sequences. To do this, SRST2 must first check the index exists, call `bowtie2-build` to generate the index if it doesn't already exist, and then call bowtie2 to map the reads to this indexed reference. Occasionallly bowtie2 will return an Error message saying that it doesn't like the index. This seems to be due to the fact that if you submit multiple SRST2 jobs to a cluster at the same time, they will all test for the presence of the index and, if index files are present, will proceed with mapping... but this doesn't mean the indexing process is actually finished, and so errors will arise. 

The simple way out of this is, if you are running lots of SRST2 jobs, FIRST index your reference(s) for bowtie2 and samtools (using `bowtie2-build ref.fasta ref.fasta` and `samtools faidx ref.fasta`), then submit your SRST2 jobs. The slurm_srst2.py script takes care of this for you by formatting the databases before submitting any SRST2 jobs.


# Generating SRST2-compatible clustered database from raw sequences

### Gene database format
In addition to MLST, SRST2 can do gene/allele detection. This works by mapping reads to each of the reference sequences in a fasta file(s) (provided through `--gene_db`) and reporting details of all genes that are covered above a certain level (`--min_coverage`, 90% by default). 

If the input database contains different alelles of the same gene, SRST2 can report just the best matching allele for that gene (much like with MLST we report the best matching allele for each locus in the scheme). This makes the output manageable, as you will get one column per gene/locus (e.g. blaCTX-M) which reports the specific allele that was called in each sample (e.g. blaCTX-M-15 in sample A, blaCTX-M-13 in sample B).

To do this properly, SRST2 needs to know which of the reference sequences are alleles of the same gene. This is done by adhering to the following format in the naming of the sequences (i.e. the headers in the fasta sequence for the database):

```
>[clusterUniqueIdentifier]__[clusterSymbol]__[alleleSymbol]__[alleleUniqueIdentifier]
```

e.g. in the resistance gene database provided, the first entry is:

```
>344__blaOXA__blaOXA-181__1
```

Note these are separated by two underscores. The individual components are:

clusterUniqueIdentifier = 344;  unique identifier for this cluster (uniquely identifes the cluster)
clusterSymbol = blaOXA;  gene symbol for this cluster (may be shared by multiple clusters)
alleleSymbol = blaOXA-181;  full name of this allele
alleleUniqueIdentifier = 1;  uniquely identifies the sequence

Ideally the alleleSymbol would be unique (as it is in the `reference.fasta` file provided). However it doesn't have to be: if allele symbols are not unique, then SRST2 will use the combination `[alleleSymbol]__[alleleUniqueIdentifier]` to  uniquely identify the sequence in the resulting reports, so that you can trace exactly which sequence was present in each sample.

Additional gene annotation can appear on the header line, after a space. This additional info will be printed in the full genes report, but not in the compiled results files.

e.g. for the blaOXA sequence above, the full header is actually:

	>344__blaOXA__blaOXA-181__1 blaOXA-181_1_HM992946; HM992946; betalactamase


### Sourcing suitable gene databases

To get started, we have provided a resistance gene database (`data/ARGannot_r3.fasta`) and code (`database_clustering/`) to extract virulence factors for a genus of interest from the Virulence Factor DB (detailed instructions below).

If you want to use your own database of allele sequences, with the reporting behaviour described, you will need to assign your sequences to clusters and use this header format. To facilitate this, use the scripts provided in the database_clustering directory provided with SRST2, and follow the instructions below.

You can also use unclustered sequences. This is perfectly fine for gene detection applications, where you have one representative allele sequence for each gene, and you simply want to know which samples contain a sequence that is similar to this one (e.g. detecting plasmid replicons, where there is one target sequence per replicon). However, this won't work well for allele typing. If the sequence database contains multiple allele sequences for the same gene, then all of these that are covered above the length threshold (default 90%) will be reported in the output, which makes for messy reporting. If you do this, you would probably find it most useful to look at the full gene results table rather than looking at the compiled results output.

### If you already know which alleles belong to the same gene family and should be clustered for reporting:

Use this info to generate the appropriate headers for SRST2 to read. The provided script `csv_to_gene_db.py` can help:

```
csv_to_gene_db.py -h

Usage: csv_to_gene_db.py [options]

Options:

  -h, --help            show this help message and exit
  
  -t TABLE_FILE, --table=TABLE_FILE
                        table to read (csv)
                        
  -o OUTPUT_FILE, --out=OUTPUT_FILE
                        output file (fasta)
                        
  -s SEQ_COL, --seq_col=SEQ_COL
                        column number containing sequences
                        
  -f FASTA_FILE, --fasta=FASTA_FILE
                        fasta file to read sequences from (must specify which
                        column in the table contains the sequence names that
                        match the fasta file headers)
                        
  -c HEADERS_COL, --headers_col=HEADERS_COL
                        column number that contains the sequence names that
                        match the fasta file headers
```

The input table should be comma-separated (csv, unix newline characters) and have these columns:

seqID,clusterid,gene,allele,(DNAseq),other....

which will be used to make headers of the required form `[clusterID]__[gene]__[allele]__[seqID] [other stuff]`

If you have the sequences as a column in the table, specify which column they are in using `-s`:

```
csv_to_gene_db.py -t genes.csv -o genes.fasta -s 5
```

Alternatively, if you have sequences in a separate fasta file, you can provide this file via `-f`. You will also need to have a column in the table that links the rows to unique sequences, specify which column this is using `-c`:

```
csv_to_gene_db.py -t genes.csv -o genes.fasta -f rawseqs.fasta -c 5
```

### Clustering sequences

If your sequences are not already assigned to gene clusters, you can do this automatically using CD-HIT (http://weizhong-lab.ucsd.edu/cd-hit/).

1 - Run CD-HIT to cluster the sequences at 90% nucleotide identity:

```
cdhit-est -i rawseqs.fasta -o rawseqs_cdhit90 -d 0 > rawseqs_cdhit90.stdout
```

2 - Parse the cluster output and tabulate the results, check for inconsistencies between gene names and the sequence clusters, and generate individual fasta files for each cluster to facilitate further checking:

```
python cdhit_to_csv.py --cluster_file rawseqs_cdhit90.clstr --infasta raw_sequences.fasta --outfile rawseqs_clustered.csv
```

For comparing gene names to cluster assignments, this script assumes very basic nomenclature of the form gene-allele, ie a gene symbol followed by '-' followed by some more specific allele designation. E.g. adk-1, blaCTX-M-15. The full name of the gene (adk-1, blaCTX-M-15) will be stored as the allele, and the bit before the '-' will be stored as the name of the gene cluster (adk, blaCTX). This won't always give you exactly what you want, because there really are no standards for gene nomenclature! But it will work for many cases, and you can always modify the script if you need to parse names in a different way. Note though that this only affects how sensible the gene cluster nomenclature is going to be in your srst2 results, and will not affect the behaviour of the clustering (which is purely sequence based using cdhit) or srst2 (which will assign a top scoring allele per cluster, the cluster name may not be perfect but the full allele name will always be reported anyway).

3 - Convert the resulting csv table to a sequence database using:

```
csv_to_gene_db.py -t rawseqs_clustered.csv -o seqs_clustered.fasta -f rawseqs.fasta -c 4
```

The output file, `seqs_clustered.fasta`, should now be ready to use with srst2 (`--gene_db seqs_clustered.fasta`).

If there are potential inconsistencies detected at step 2 above (e.g. multiple clusters for the same gene, or different gene names within the same cluster), you may like to investigate further and change some of the cluster assignments or cluster names. You may find it useful to generate neighbour joining trees for each cluster that contains >2 genes, using align_plot_tree_min3.py

### Screening for resistance genes with SRST2
A preliminary set of resistance genes based on the ResFinder database and CARD is included with SRST2: `data/ARGannot_r3.fasta`. The fasta file is ready for use with SRST2. The CSV table contains the same sequence information, but in tabular format for easier parsing/editing.

An easy way to add sequences to this database would be to add new rows to the table, and then generate an updated fasta file using:

```
csv_to_gene_db.py -t rawseqs_clustered.csv -o seqs_clustered.fasta -s rawseqs.fasta -c 5
```

### Using the VFDB Virulence Factor Database with SRST2

The VFDB houses sets of virulence genes for a range of bacterial genera, see http://www.mgc.ac.cn/VFs/.

To type these virulence genes using SRST2, download the full set of sequences from the VFDB website (http://www.mgc.ac.cn/VFs/Down/VFDB_setB_nt.fas.gz) and follow these steps to generate SRST2-compatible files for your genus of interest.

1 - Extract virulence genes by genus from the main VFDB file, `CP_VFs.ffn`:

```
python VFDBgenus.py --infile CP_VFs.ffn --genus Clostridium
```

or, to get all availabel genera in separate files:

```
python VFDBgenus.py --infile CP_VFs.ffn
```

2 - Run CD-HIT to cluster the sequences for this genus, at 90% nucleotide identity:

```
cd-hit -i Clostridium.fsa -o Clostridium_cdhit90 -c 0.9 > Clostridium_cdhit90.stdout
```

3 - Parse the cluster output and tabulate the results using the specific Virulence gene DB compatible script:

```
python VFDB_cdhit_to_csv.py --cluster_file Clostridium_cdhit90.clstr --infile Clostridium.fsa --outfile Clostridium_cdhit90.csv
```

4 - Convert the resulting csv table to a SRST2-compatible sequence database using:

```
python csv_to_gene_db.py -t Clostridium_cdhit90.csv -o Clostridium_VF_clustered.fasta -s 5
```
    
The output file, `Clostridium_VF_clustered.fasta`, should now be ready to use with srst2 (`--gene_db Clostridium_VF_clustered.fasta`).

### Using the EcOH database for serotyping E. coli with SRST2

Details can be found in [this MGen paper](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000064).

The EcOH database includes genes for identifying O and H types in E. coli, see `/data/EcOH.fasta`

O types are represented by the presence of two loci (either wzy and wzy, or wzm and wzt). Note that allelic variation is possible but does not impact serotype in a predictable way, so typing calls should be made based on the presence of genes rather than allele assignments (i.e. it is generally safe to ignore *? characters)

H types are represented by alleles of fliC or flnA flagellin genes (one allele per H type). Note that it is possible to carry both a fliC allele and a flnA allele, such strains should be considered phase variable for flagellin).

#### Basic Usage – Serotyping E. coli

```
srst2 --input_pe strainA_1.fastq.gz strainA_2.fastq.gz --output strainA_serotypes --log --gene_db EcOH.fasta
```

#### Example

Note these reads sets are each ~100 MB each, so you need to download 400 MB of test data to run this example

```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR178/ERR178148/ERR178148_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR178/ERR178148/ERR178148_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR178/ERR178156/ERR178156_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR178/ERR178156/ERR178156_2.fastq.gz

srst2 --input_pe ERR178148*.fastq.gz ERR178156*.fastq.gz --output serotypes --log --gene_db EcOH.fasta
```

Results will be output in: `[prefix]__genes__EcOH__results.txt`

Output from the above example would appear in: `serotypes__genes__EcOH__results.txt`

Sample | fliC | wzm | wzt | wzx | wzy
:---: | :---: | :---: | :---: | :---: | :---:
ERR178148 | fliC-H31_32* | - | - | wzx-O131_165* | wzy-O131_386
ERR178156 | fliC-H33_34 | wzm-O9_93* | wzt-O9_107* | - | -

#### Interpretation

ERR178148 has matching wzx and wzy hits for O131, and fliC allele H31, thus the predicted serotype is O131:H31. 

ERR178156 has matching wzm and wzt hits for O9 and fliC allele H33, thus the predicted serotype is O9:H33. 

Note that each O antigen type is associated with loci containing EITHER wzx and wzy, OR wzm and wzt genes.

No alleles for flnA were detected in these strains, indicating they are not phase variable for flagellin.


### Typing the LEE pathogenicity island of E. coli

Details can be found in [this Nature Microbiology paper](http://www.nature.com/articles/nmicrobiol201510).

The LEE typing database is based on analysis of >250 LEE-containing E. coli genomes and includes 7 loci (eae (i.e. intimin), tir, espA, espB, espD, espH, espZ). The data is provided as a MLST-style database, in which combinations of alleles are assigned to a LEE subtype, to facilitate a common nomenclature for LEE subtypes. However please note the database does not contain every known allele and is not intended to function as a typical MLST scheme. Rather, each sequence in the database represents a cluster of closely related alleles that have been assigned to the same locus type. So you will generally not get precise matches to a known allele, and this is not something to worry about - the goal is to identify the locus type and the overall LEE subtype. 

Also note that the LEE scheme includes three distinct lineages: Lineage 1 consists of LEE subtypes 1-2; Lineage 2 consists of LEE subtypes 3-8; Lineage 3 consists of LEE subtypes 9-30. See the paper for more details. 

The database files (sequences + MLST-style profile definitions) are in the SRST2 distribution under /data.

#### Basic Usage – Typing E. coli LEE variants

```
srst2 --input_pe strainA_1.fastq.gz strainA_2.fastq.gz --output strainA_LEE --log --mlst_db LEE_mlst.fasta --mlst_definitions LEE_STscheme.txt
```

#### Example

Note these reads sets are each ~100 MB each, so you need to download 400 MB of test data to run this example

```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR178/ERR178148/ERR178148_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR178/ERR178148/ERR178148_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR178/ERR178156/ERR178156_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR178/ERR178156/ERR178156_2.fastq.gz

srst2 --input_pe ERR178156*.fastq.gz ERR124656*.fastq.gz --output LEE --log --mlst_db LEE_mlst.fasta --mlst_definitions LEE_STscheme.txt
```

Results will be output in: `[prefix]__mlst__LEE_mlst__results.txt`

Output from the above example would appear in: `LEE__mlst__LEE_mlst__results.txt`

Sample | ST | eae | tir | espA | espB | espD | espH | espZ | mismatches | uncertainty | depth | maxMAF
:---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---:
ERR178148 | 30* | 19 | 10* | 9* | 5* | 6* | 4* | 13* | tir-10/33snp3indel;espA-9/3snp;espB-5/4snp;espD-6/9snp;espH-4/1snp;espZ-13/3snp | - | 59.8217142857 | 0.285714285714
ERR178156 | 9* | 6* | 5* | 5* | 3* | 3* | 4* | 8 | eae-6/4snp;tir-5/18snp;espA-5/6snp;espB-3/4snp;espD-3/8snp;espH-4/3snp | - | 33.5062857143 | 0.444444444444

#### Interpretation

ERR178148 and ERR178156 carry LEE subtypes 30 and 9 respectively, which both belong to LEE lineage 3. The imprecise matches are to be expected and can be taken as a confident assignment of LEE type.

### Plotting output in R

Some R functions are provided in scripts/plotSRST2data.R for plotting SRST2 output to produce images like those in the paper (e.g. Figure 8: http://www.genomemedicine.com/content/6/11/90/figure/F8)

These functions require the `ape` package to be installed.

Example usage:

```
# load the functions in R
source("srst2/scripts/plotSRST2data.R")

1. EXAMPLE FROM FIGURE 8A (http://www.genomemedicine.com/content/6/11/90/figure/F8): viewing resistance genes in individual strains that have been analysed for MLST + resistance genes

# read in a compiled MLST + genes table output by SRST2
Ef_JAMA<-read.delim("srst2/data/EfaeciumJAMA__compiledResults.txt", stringsAsFactors = F)

# Check column names. Sample names are in column 1 (strain_names=1 in the function call below), MLST data is in columns 2 to 9 (mlst_columns = 2:9), while gene presence/absence information is recorded in columns 13 to 31 (gene_columns = 13:31). 

colnames(Ef_JAMA)
 [1] "Sample"         "ST"             "AtpA"           "Ddl"            "Gdh"            "PurK"          
 [7] "Gyd"            "PstS"           "Adk"            "mismatches"     "uncertainty"    "depth"         
[13] "Aac6.Aph2_AGly" "Aac6.Ii_AGly"   "Ant6.Ia_AGly"   "Aph3.III_AGly"  "Dfr_Tmt"        "ErmB_MLS"      
[19] "ErmC_MLS"       "MsrC_MLS"       "Sat4A_AGly"     "TetL_Tet"       "TetM_Tet"       "TetU_Tet"      
[25] "VanA_Gly"       "VanH.Pt_Gly"    "VanR.A_Gly"     "VanS.A_Gly"     "VanX.M_Gly"     "VanY.A_Gly"    
[31] "VanZ.A_Gly"    

# Make a tree based on the MLST loci, and plot gene content as a matrix. 
# cluster=T turns on hierarchical clustering of the columns (=genes); labelHeight, infoWidth and treeWidth control the relative dimensions of the plotting areas available for the tree, printed MLST information, and gene labels.

geneContentPlot(m=Ef_JAMA, mlst_columns = 2:9, gene_columns = 13:31, strain_names=1, cluster=T, labelHeight=40, infoWidth=15, treeWidth=5)

2. EXAMPLE FROM FIGURE 7C (http://www.genomemedicine.com/content/6/11/90/figure/F7): viewing summaries of resistance genes by ST, in a set of isolates that have been analysed for MLST + resistance genes

# read in a compiled MLST + genes table output by SRST2
Ef_Howden<-read.delim("srst2/data/EfaeciumHowden__compiledResults.txt", stringsAsFactors = F)

# Make a tree of STs based on the MLST alleles and plot this; Group isolates by ST and calculate the number of strains in each ST that contain each resistance gene;  plot these counts as a heatmap. Note that barplots of the number of strains in each ST are also plotted to the right.
geneSTplot(d,mlst_columns=8:15,gene_columns=17:59,plot_type="count",cluster=T)

# Same as above but plot the rate of detection of each gene within each ST, not the raw counts
geneSTplot(d,mlst_columns=8:15,gene_columns=17:59,plot_type="rate",cluster=T)

# To suppress SNPs, i.e. collapse ST1 and ST1* into a single group for summarisation at the clonal complex level, set suppressSNPs=T.

# To suppress uncertainty due to low depth, i.e. collapse ST1 and ST1? into a single group for summarisation at the clonal complex level, set suppressUncertainty=T.
```

Note, heatmap colours can be set via the `matrix.colours` parameter in both of these functions. The default value is matrix.colours=colorRampPalette(c("white","yellow","blue"),space="rgb")(100), i.e. white=0% gene frequency, yellow = 50% and blue = 100%.
