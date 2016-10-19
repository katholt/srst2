Databases formatted for use with SRST2
====

See documentation at http://katholt.github.io/srst2/ for details on formatting other sequences databases for use.

Current MLST databases can be retrieved from pubmlst.org using getmlst.py (distributed with srst2 in the /scripts directory).

====
RESISTANCE GENES

1. ARG-Annot (ARGannot.r1.fasta) - *RECOMMENDED

Citation: Gupta et al, Antimicrob Agents Chemother 2014; 58(1):212. DOI: 10.1128/AAC.01310-13

Processing history: 
- Sequences downloaded (May 2014) from: http://www.mediterranee-infection.com/article.php?laref=282&titer=arg-annot
- Clustered at 80% nucleotide identity using CDhit-est
- Minor manual correction of gene names and accession formats

Notes:
- Quality: ARG-Annot has been manually curated and appears to be the highest quality sequence set in terms of its informativeness (accession & position of original sequence; antibiotic class) and internal consistency. It contains almost all sequences found in CARD and ResFinder, but has been manually curated and had many redundancies removed.
- Annotation (from ARG-Annot documentation): The nucleotide sequences included in this database from different antibiotics classes are abbreviated as AGly: aminoglycosides, Bla: beta-lactamases, Fos: fosfomycin, Flq: fluoroquinolones, Gly: glycopeptides, MLS: macrolide-lincosamide- streptogramin, Phe: phenicols, Rif: rifampicin, Sul: sulfonamides, Tet: tetracyclines and Tmt: trimethoprim. A unified nomenclature system was followed in which the name contains all of the information regarding gene class, gene name, accession number, gene location in the sequence and gene size. For example, (AGly) AadA1:M95287:3320-4111:792 tells the researcher that the class of antibiotics is AGly: aminoglycosides, the gene name is AadA1, the accession number is M95287, the gene location is 3320-4111, and the gene size is 792.

----------

2. ResFinder (ResFinder.fasta)

Citation: Zankari et al, J Antimicrob Chemother 2012; 67:2640–2644 DOI: 10.1093/jac/dks261

Processing history: 
- Sequences downloaded (2013) from: http://cge.cbs.dtu.dk/services/data.php
- Clustered at 90% nucleotide identity using CDhit-est
- Manual correction of gene names and removal of redundant sequences; some nomenclature corrected to common usage

====

PLASMID REPLICONS

1. 18 replicons targeted by Caratolli 2005, "Identification of plasmids by PCR-based replicon typing" (Plasmid18Replicons.fasta)

Citation: Carattoli et al, Journal of Microbiological Methods 63 (2005) 219–228 DOI: 10.1016/j.mimet.2005.03.018

Processing history:
- Sequences for target amplicons retrieved from EMBL using the primers and sequence accessions listed in Table 1 of the paper

----------

2. PlasmidFinder (PlasmidFinder.fasta)

Citation: Carattoli et al, Antimicrob Agents Chemother 2014; 58(7):3895-3903 DOI: 10.1128/AAC.02412-14

Processing history: 
- Sequences downloaded (June 2014) from: http://cge.cbs.dtu.dk/services/data.php
- Clustered at 80% nucleotide identity using CDhit-est
- Minor manual correction of gene names and removal of duplicate sequences

====

VIRULENCE FACTORS

See documentation at http://katholt.github.io/srst2/ for details on formatting VFDB and other databases.