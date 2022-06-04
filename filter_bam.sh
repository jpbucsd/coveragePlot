#!/bin/bash

MINSCORE=60

for file in *.bam
do
	samtools view -bq $MINSCORE $file > filtered_ex/ex_filtered_$file
done

