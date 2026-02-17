# SPARTA Data analysis
Code used for data analysis for the original SPARTA paper

## Overview
This pipeline builds a custom database for taxonomic classification of long-read 16S rRNA gene sequencing data.\
The database is based on NCBI RefSeq 16S TargetedLoci sequences and supports TaxID-aware classification and downstream lineage (taxonomy tree) annotation.\
(https://www.ncbi.nlm.nih.gov/nuccore?term=33175%5BBioProject%5D+OR+33317%5BBioProject%5D)\
Important: Steps 1–3 must be run from the same working directory to ensure all intermediate files are found correctly.

## 1. ete3_build
Create an ETE3 NCBI taxonomy SQLite database at taxdump/ncbi.sqlite.\
Run example:
```console
python ete3_build.py
```

## 2. Database construction
Build a custom Kraken2 16S database from NCBI RefSeq TargetedLoci 16S FASTAs (Bacteria + Archaea).\
Run example:
```console
python build_kraken2_refseq_16s_db
```
If needed, standard database can be constructed with kraken buitin-database\
(https://github.com/DerrickWood/kraken2)


## 3. Long-read 16S rRNA gene sequencing taxonomic classification
This script takes a folder of FASTA/FASTQ files, runs Kraken2 on each one to do taxonomic classification, and saves the Kraken report. 
It then runs Bracken on those Kraken reports to estimate abundances (at the level you set, like species),
and converts the Bracken results into a clean CSV that includes the full taxonomy lineup (Kingdom → Species) using an ETE3 NCBI taxonomy database. 
Run example:
```console
python taxa_classification.py input_directory_fastq
```

## System Requirements
No additional hardware is required.\
The pipeline has been tested on:\
- Linux Ubuntu 16.04
- macOS Ventura 13.2.1\
Database construction requires approximately 4-5 Gb of storage. 


## Dependencies:

python 3.9.22\
ete3 3.1.3\
kraken2 2.1.5\
bracken 3.1\
pandas 2.2.3


## Run time
Installing packages should only take a few minutes.

On a normal computer, the demos in the repository should also complete within a few minutes.
