#!/bin/bash

#make index file
bwa index sacCer3.fa 
#makes SAM file for each fastq
for SAMPLE in *.fastq 
do 
	bwa mem -R "@RG\tID:${SAMPLE%.fastq}\tSM:${SAMPLE%.fastq}" sacCer3.fa $SAMPLE >${SAMPLE%.fastq}.sam
done 
#bam files from all sam files
for FILENAME in *.sam 
do 
	samtools sort -o ${FILENAME%.sam}.bam -O bam $FILENAME
done
# find all variants
freebayes -f sacCer3.fa  -p 1 --genotype-qualities *.bam > yeast_vars.vcf
# filter variants based on genotype 
vcffilter -f "QUAL >20" yeast_vars.vcf >yeast_filtered.vcf 
#decompose complex haplotypes
vcfallelicprimitives -kg yeast_filtered.vcf > yeast_decomposed.vcf
#variant annotation
snpeff ann -lof R64-1-1.86 yeast_decomposed.vcf -o vcf > yeast_annotated.vcf
#1st 1000 lines of the file 
head -n 1000 yeast_annotated.vcf > yeast_annoted_short.vcf