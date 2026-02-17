# SPARTA Data analysis
Code used for data analysis for the original SPARTA paper

## Overview
This pipeline builds a custom database for taxonomic classification of long-read 16S rRNA gene sequencing data.\
The database is based on NCBI RefSeq 16S TargetedLoci sequences and supports TaxID-aware classification and downstream lineage (taxonomy tree) annotation.\
(https://www.ncbi.nlm.nih.gov/nuccore?term=33175%5BBioProject%5D+OR+33317%5BBioProject%5D)\
Important: Steps 1–3 must be run from the same working directory to ensure all intermediate files are found correctly.

## 1. ete3_build
Example:
'''bash
python ete3_build.py
'''
`python ete3_build.py`

## Dependencies:

python 3.9.22\
ete3 3.1.3\
kraken2 2.1.5\
bracken 3.1\
pandas 2.2.3


## Run time
Installing packages should only take a few minutes.

On a typical computer, the demos in the repository should also complete within a few minutes.
