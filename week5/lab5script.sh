#!/bin/bash 
makes index file 
bowtie2-build chr19.fa chr19 
for sample in g1e/*.fastq
do 
	#maps reads
	bowtie2 -x chr19 -U ${sample} -S ${sample%.fastq}.sam -p 6
	#makes a bam file
	samtools view -bSo ${sample%.fastq}.bam ${sample%.fastq}.sam
	#makes sorted file
	samtools sort ${sample%.fastq}.bam -o ${sample%.fastq}.sorted.bam
	#indexes file
	samtools index ${sample%.fastq}.sorted.bam
done

#finds macs peak
for ct in G1E ER4
do 
	macs2 callpeak -t g1e/CTCF_${ct}.bam -c g1e/input_${ct}.bam --format=BAM --name CTCF_${ct} -g 61420004

done 
#finds differentially expressed regions 
bedtools subtract -a CTCF_ER4_peaks.narrowPeak -b CTCF_G1E_peaks.narrowPeak >lost.bed 
bedtools subtract -a CTCF_G1E_peaks.narrowPeak -b CTCF_ER4_peaks.narrowPeak >gained.bed

#splits up 
for feature in promoter intron exon
do 
	grep "${feature}" Mus_musculus.GRCm38.94_features.bed>${feature}.bed 
done 
for ct in G1E ER4 
do
	for feature in promoter intron exon
	do 
		bedtools intersect -a CTCF_${ct}_peaks.narrowPeak -b ${feature}.bed  >${ct}_${feature}_overlap.bed
	done
done 
