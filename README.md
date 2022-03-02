# COMP 383 Mini-Project

## Introduction
This python wrapper will take the accession number for a set of Illumina reads and assemble the genome, predict the coding regions, and annotate the functionality of these regions. 

## Installation
Download the python file and place it in your home directory. It can be run from a command line terminal. 

## Required Dependencies
* [SRA Toolkit](https://github.com/ncbi/sra-tools)
* [SPAdes](https://cab.spbu.ru/software/spades/)
* [GeneMarkS-2](http://exon.gatech.edu/GeneMark/license_download.cgi)

## Command Line Options
* --a: Accession number of Illumina reads (default: SRR8185130)
* --o: Directory for file outputs (default: /home/mkill/results) Note: Ensure directory exists before running command
* --h: Help

## Relevant Outputs
* miniproject.log: Summary of assembly and annotation
* predicted_functionality.csv: Best hit functionality matches for each predicted CDS

# Test Data
The wrapper was tested using the accession 'SRR8185130'.


