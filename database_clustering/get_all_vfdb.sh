#!/bin/bash
#this is a utility bash script that automates generation of all the VFDB gene databases for use with srst2.py
#script assumes you already have python, and cd-hit installed somewhere on the $PATH
#this script MUST be in the same folder as all the other database_clustering python scripts
#example usage:
#/srst2/database_clustering/get_all_vfdb.sh ./CP_VFs.ffn ./VFDB

VFDBFILE=$(readlink -e $1)
OUTPUTFOLDER=$2
#get the srst2/database_clustering folder where all the other python scripts live side-by-side with this one
DBCLUSTERINGSCRIPTFOLDER=$(dirname $(readlink -e $0))

#if the specified output folder doesn't exist, then create it
if [ ! -d ${OUTPUTFOLDER} ]; then
  mkdir ${OUTPUTFOLDER}
fi
cd ${OUTPUTFOLDER}

#extract virulence genes from all available genera into separate files
python ${DBCLUSTERINGSCRIPTFOLDER}/VFDBgenus.py --infile ${VFDBFILE}

#loop over each genus' *.fsa file and generate the gene database fasta file
for FSAFILE in *.fsa; do
  FILENAME=${FSAFILE##*/}
  GENUS=${FILENAME%.*}
  mkdir ${GENUS}

  echo Generating virulence gene database for ${GENUS}

  #Run CD-HIT to cluster the sequences for this genus, at 90% nucleotide identity
  cd-hit -i ${FILENAME} -o ${GENUS}/${GENUS}_cdhit90 -c 0.9 > ${GENUS}/${GENUS}_cdhit90.stdout

  #Parse the cluster output and tabulate the results using the specific Virulence gene DB compatible script:
  python ${DBCLUSTERINGSCRIPTFOLDER}/VFDB_cdhit_to_csv.py --cluster_file ${GENUS}/${GENUS}_cdhit90.clstr --infile ${FILENAME} --outfile ${GENUS}/${GENUS}_cdhit90.csv

  #Convert the resulting csv table to a SRST2-compatible sequence
  python ${DBCLUSTERINGSCRIPTFOLDER}/csv_to_gene_db.py -t ${GENUS}/${GENUS}_cdhit90.csv -o ${GENUS}/${GENUS}_VF_clustered.fasta -s 5

  #move the original *.fsa file to the created genus subfolder
  mv ${FILENAME} ${GENUS}/${FILENAME}
done
