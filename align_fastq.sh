#!/bin/bash
REFPATH=../ncbi-genomes-2022-06-01/GCF_000021165.1_ASM2116v1_cds_from_genomic.fna

head $REFPATH

for file in *.fastq;
do
       	echo $file	
	bwa mem -t 8 -T -15 -V ${REFPATH} ${file} | samtools view -bS > ~/vannini_data_aligned_scorefilt/${file}.bam;
done

