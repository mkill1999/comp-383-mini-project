# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 17:21:29 2022

@author: mck1999
"""
import os
import argparse

parser = argparse.ArgumentParser(description = 'A script to automate the assembly and annotation of illumina reads.')
parser.add_argument('--a', default = 'SRR8185310', help = 'Illumina accession to be retrieved')
parser.add_argument('--o', default = '/home/mkill/results', help ='Directory for output')
args = parser.parse_args()
accession = args.a
output_dir = args.o

#accession = 'SRR8185310'
#output_dir = '/home/mkill/results'
log = open(output_dir + '/miniproject.log', 'w')

#1: retrieve Illumina reads
sra_toolkit_path='/home/mkill/sra/sratoolkit.2.11.2-ubuntu64/bin/' #path of sratoolkit
os.system(sra_toolkit_path + 'prefetch ' + accession + ' -O ' + output_dir) #fetch accession from NCBI
os.system(sra_toolkit_path + 'fasterq-dump ' + output_dir + '/' + accession + ' -O ' + output_dir) #get fastq file from accession
print(sra_toolkit_path + 'fasterq-dump ' + output_dir + '/' + accession + ' -outdir ' + output_dir)
#writes out to results directory, name is the SRR#.fastq

#2: assemble genome using SPAdes
spades_path='/home/mkill/SPAdes/SPAdes-3.15.4-Linux/bin/' #path of spades
spades_command='python3 ' + spades_path + 'spades.py -k 55,77,99,127 -t 2 --only-assembler -s ' + output_dir + '/' + accession + '.fastq -o ' + output_dir #command for spades
os.system(spades_command) #run spades command
log.write(spades_command + '\n') #write spade command to log
#writes out to results directory, name is contigs.fasta

#3: calculate number of contigs > 1000
f = open(output_dir + '/contigs.fasta', 'r') #open file of all scaffolds
inp = f.readlines()  #read in input
fasta=dict() #creat dictionary for sequences
i=0 #initialize for while loop
while i < len(inp): #for every element in inp
    if inp[i][0]=='>': #if line is header
        name = inp[i] #set header as 'name'
        fasta[name]=inp[i+1] #create dictionary entry of header: beginning of sequence
        i+=2 #go to next line after first line of sequence
    else: #if line is sequence
        fasta[name]+=inp[i] #append sequence to existing sequence
        i += 1 #go to next line
outp = open(output_dir + '/contigs_1000.fasta', 'w') #open output file for filtered contigs

fasta_1000 = dict()
counter=0 #intiatlize counter for large sequences
for i in fasta: #for each dicttionary entry
    if len(fasta[i])>1000: #if sequence is > 1000
        counter+=1 #increase count of long contigs
        outp.write(i+fasta[i]) #write header and sequence to output file
        fasta_1000[i] = fasta[i] #write to dict of sequences > 1000
outp.close() #close output file     
log.write("There are " + str(counter) + " contigs > 1000 in the assembly. \n") #write result to log
#writes out to results directory, name is contigs_1000.fasta 

#4: Calculate the length of the assembly (> 1000 bp)
total = 0 #initilize total to 0
for z in fasta_1000.values(): #for each contig > 1000
    total += len(z) #add length of sequence to total
log.write("There are " + str(total) + " bp in the assembly. \n") #write result to log

#5: predict protein sequences with GeneMarkS-2
genemark_path = '/home/mkill/gms2_linux_64/gms2.pl' #path for gene-mark
fasta_path = output_dir + '/contigs_1000.fasta' #path for contigs > 1000
os.system('perl ' + genemark_path + ' --seq ' + fasta_path + ' --genome-type bacteria --output ' + output_dir + '/gms_output.lst --faa ' + output_dir + '/gms_protein.fasta') #run genemark command
#writes out to results directory, name is gms_protein.fasta

#6: predict functionality from amino acid sequences
ref_path = '/home/mkill/mini-project/ecoli-ref' #path of reference database
predicted_path = output_dir + '/gms_protein.fasta' #path of predicted AA sequences
blast_output = output_dir + '/predicted_functionality.csv' #path for outputsv'
blast_cmd = 'blastp -query ' + predicted_path + ' -db ' + ref_path + ' -out ' + blast_output + ' -outfmt "10 qseqid sseqid pident qcovs" -max_target_seqs 1' #blast command
os.system(blast_cmd) #run blast command

#7:identify discrepency in CDS
ref_cds = 4140
genemark_fasta = open(predicted_path, 'r')
genemark_seq = genemark_fasta.read()
seq_cds = genemark_seq.count('>')
log.write('GeneMarkS found ' + str(seq_cds-ref_cds) + ' additional CDS than the RefSeq. \n')
log.close()