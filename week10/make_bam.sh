#!/bin/bash
for FILENAME in *.sam 
do 
	samtools sort -o ${FILENAME%.sam}.bam -O bam $FILENAME
done