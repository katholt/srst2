Databases formatted for use with SRST2
====

See documentation at http://katholt.github.io/srst2/ for details on formatting other sequences databases for use.

The latest MLST databases can be retrieved from pubmlst.org using getmlst.py (distributed with SRST2 in the /scripts directory).

# 1. Antimicrobial resistance

## ARG-ANNOT (ARGannot_r2.fasta) - recommended

Citation: Gupta et al, Antimicrob Agents Chemother 2014; 58(1):212. DOI: [10.1128/AAC.01310-13](http://aac.asm.org/content/58/1/212.abstract)

### Features
* This is the second revision of the SRST2-formatted ARG-ANNOT database. See [update summary](#summary_r2) for details.
* Quality: ARG-ANNOT has been manually curated and appears to be the highest quality sequence set in terms of its informativeness (accession & position of original sequence; antibiotic class) and internal consistency. It contains almost all sequences found in CARD and ResFinder, but has been manually curated and had many redundancies removed.
* Annotations (following the ARG-ANNOT documentation): The nucleotide sequences included in this database from different antibiotics classes are abbreviated as **AGly**: aminoglycosides, **Bla**: beta-lactamases, **Fos**: fosfomycin, **Flq**: fluoroquinolones, **Gly**: glycopeptides, **MLS**: macrolide-lincosamide- streptogramin, **Phe**: phenicols, **Rif**: rifampicin, **Sul**: sulfonamides, **Tet**: tetracyclines and **Tmt**: trimethoprim. A unified nomenclature system was followed in which the name contains all of the information regarding gene class, gene name, accession number, gene location in the sequence and gene size. For example, *(AGly) AadA1:M95287:3320-4111:792* tells the researcher that the class of antibiotics is AGly: aminoglycosides, the gene name is *AadA1*, the accession number is M95287, the gene location is 3320-4111, and the gene size is 792 bp.

### Update summary<a name="summary_r2"></a>
* This update aims to incorporate revisions and new sequences from the latest version of the ARG-ANNOT database ([ARG-ANNOT Nt V3](http://en.mediterranee-infection.com/arkotheque/client/ihumed/_depot_arko/articles/1424/arg-annot-nt-v3-march2017_doc.fasta), accessed in July 2017) into ARGannot.r1.fasta.
* ARGannot_r2.fasta is a non-redundant (in terms of identical nucleotide sequences) database designed for screening antimicrobial resistance genes.
* For colistin resistance determinants, six novel alleles (*Mcr1-2* - *Mcr1-7*) of the *Mcr1* gene and another two genes (*Mcr3* and *Mcr4*) were added to this database.
* Checked and corrected gene/allele names, accession numbers, genome coordinates and sequences. The file ARGannot\_clustered80\_r2.csv in this repository provides a table of sequence annotations.
* Removed *TetR* and *TetRG* alleles from the previous database because they are regulatory genes for *tetA* and *tetG* genes but themselves do not confer resistance.
* Allele sequences were clustered under a 80% nucleotide identity using CD-HIT-EST. None of sequences from the ARGannot.r1.fasta database dropped out from its original cluster. As such, cluster IDs remain the same for sequences from the previous version and were assigned to new sequences accordingly. The same rule applies to sequence IDs as well.
* See ARGannot_r2_changes.log in this repository for more details.

**Notice** Since this database is used for gene screen at the nucleotide level, different alleles encoding the same amino acid sequences were retained, which is different from the ARG-ANNOT Nt V3 database. 

### Previous versions
**ARGannot.r1.fasta**
To make this version, redundant sequences (those of 100% nucleotide identity) were removed from ARGannot.fasta. Furthermore, minor manual corrections of gene or allele names, sequences, accession numbers, coordinates and gene sizes were applied. This version was Released in July 2015.

**ARGannot.fasta**
This is the first release of the SRST2-formatted ARG-ANNOT database, which was published with the [SRST2 paper](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-014-0090-6) in 2014.

* Sequences were downloaded (May 2014) from: [IHU Méditerranée Infection](http://www.mediterranee-infection.com/article.php?laref=282&titer=arg-annot).  
* Clustered at 80% nucleotide identity using [CD-HIT-EST](https://github.com/weizhongli/cdhit/releases).  
* Manually curated annotations.

## ResFinder (ResFinder.fasta)

Citation: Zankari et al, J Antimicrob Chemother 2012; 67:2640–2644 DOI: [10.1093/jac/dks261](https://academic.oup.com/jac/article/67/11/2640/707208/Identification-of-acquired-antimicrobial)

Processing history:

* Sequences were downloaded (2013) from the [Center for Genomic Epidemiology (CGE)](http://cge.cbs.dtu.dk/services/data.php).
* Clustered at 90% nucleotide identity using CD-HIT-EST.
* Manual correction of gene names and removal of redundant sequences; some nomenclature corrected to common usage.

# 2. Plasmid replicons

### PlasmidFinder (PlasmidFinder.fasta)

Citation: Carattoli et al, Antimicrob Agents Chemother 2014; 58(7):3895-3903 DOI: [10.1128/AAC.02412-14](http://aac.asm.org/content/58/7/3895)

Processing history:
 
* Sequences were downloaded in June 2014 from [CGE](http://cge.cbs.dtu.dk/services/data.php).
* Clustered at 80% nucleotide identity using CD-HIT-EST.
* Minor manual correction of gene names and removal of duplicate sequences.

### Eighteen replicons targeted by Caratolli in 2005  (Plasmid18Replicons.fasta)
Citation: Carattoli et al, Identification of plasmids by PCR-based replicon typing, Journal of Microbiological Methods 63 (2005) 219–228 DOI: [10.1016/j.mimet.2005.03.018](http://www.sciencedirect.com/science/article/pii/S0167701205001132?via%3Dihub)

Processing history: Sequences for target amplicons were retrieved from EMBL following the primers and sequence accessions listed in the Table 1 of this paper.

# 3. Virulence factors

See [Using the VFDB Virulence Factor Database with SRST2](https://github.com/katholt/srst2#using-the-vfdb-virulence-factor-database-with-srst2) for details on formatting the VFDB and other databases.
