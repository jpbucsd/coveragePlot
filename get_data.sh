#!/bin/bash

while read line; \
do
	echo "Now downloading ${line}..."
	fastq-dump -O vannini_data $line	
done < SraAccList.txt
