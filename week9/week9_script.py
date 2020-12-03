#!/usr/bin/env python2
from __future__ import division
import hifive
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math
import hifive
import pyBigWig


##Heatmap

hic=hifive.HiC("hic_project.txt",'r')

chr13=hic.cis_heatmap('chr13',1000000,datatype='fend',arraytype='full',diagonalincluded=True)
#Enrichment value with pseudocounts
enrichment=(chr13[:,:,0]+1)/(chr13[:,:,1]+1)
#log the encrichment soces
log_enrichment=np.log(enrichment)

fig,ax=plt.subplots(figsize=(10,10))
ax.set_title("Chr 13 Heatmap")
ax=sns.heatmap(log_enrichment)
ax.set_xlabel("X-Coord (Mb)")
ax.set_ylabel("Y-Coord (Mb)")
plt.savefig("ch13_heatmap.png")

##Compartment Analysis 

Comp = hifive.hic_domains.Compartment(hic, 100000, chroms=['chr13'], out_fname='tmp.hdf5')
Comp.write_eigen_scores('hic_comp.bed')
X = Comp.positions['chr13']
#extracts the start coordinate of each composition
X_start=[x[0] for x in X]
Y = Comp.eigenv['chr13']

fig,ax=plt.subplots(figsize=(10,10))
sns.scatterplot(X_start,Y)
ax.set_title("Comparment Scores")
ax.set_xlabel("X-Coord (Mb)")
ax.set_ylabel("Eigen Scores")
plt.savefig("Comparment_scores.png")

genes=open('data/WT_fpkm.bed')
gene_dict={}
for line in genes: 
	a=line.split('\t')
	#gene name: length, expression, start,end 
	gene_dict[a[3]]=[int(a[2])-int(a[1]),float(a[4][:-1]),int(a[1]),int(a[2])]

compa=open('a_compartment_genes.bed')
compartment_a_genes=[]
a_dict={}
for line in compa:
	a=line.split('\t')
	gene=a[3]
	start=int(a[1])
	end=int(a[2])
	if gene in a_dict:

		a_dict[gene]+=end-start
	else:
		a_dict[gene]=end-start
for gene in a_dict:
	if gene in gene_dict:
		comp=a_dict[gene]	
		size=comp/(gene_dict[gene][0])
		if size>0.5:
			compartment_a_genes.append([gene,gene_dict[gene][1]])

compb=open('b_compartment_genes.bed')
compartment_b_genes=[]
b_dict={}
for line in compb:
	b=line.split('\t')
	gene=b[3]
	start=int(b[1])
	end=int(b[2])
	if gene in b_dict:
		b_dict[gene]+=end-start
	else:
		b_dict[gene]=end-start
for gene in b_dict:
	if gene in gene_dict:
		comp=b_dict[gene]	
		size=comp/(gene_dict[gene][0])
		if size>0.5:
			compartment_b_genes.append([gene,gene_dict[gene][1]])
a_ge=[]
b_ge=[]
#Removing outliers: massaging data for our nefarious purposes 
for a in compartment_a_genes:
	if a[1]<10:
		a_ge.append(a[1])
for b in compartment_b_genes:
	if b[1]<10:
		b_ge.append(b[1])

fig, ax = plt.subplots()
ax.violinplot(dataset=[a_ge, b_ge],showmedians=True)
ax.set_xticks([1,2])
ax.set_xticklabels(["Comp A","Comp B"])
ax.set_ylabel("Gene Expression")
ax.set_title("Gene Expression across compartments")
plt.savefig("violinplot.png")


bw = pyBigWig.open('data/WT_H3K27me3.bw')
a_h327=[]
a_exp=[]
#Expression vs. Repression
for a in compartment_a_genes:
	gene=a[0]
	start=gene_dict[gene][2]
	end=gene_dict[gene][3]
	gene_expression=gene_dict[gene][1]
	#Caps gene expression at 10x enriched
	if gene_expression<10:
		h3_val=bw.stats('chr13', start, end, type='sum')
		#corrects for nonetype
		if h3_val is None:
			h3_val=0
		a_h327+=h3_val
		a_exp.append(gene_expression)
b_h327=[]
b_exp=[]
#Expression vs. Repression
for b in compartment_b_genes:
	gene=b[0]
	start=gene_dict[gene][2]
	end=gene_dict[gene][3]
	gene_expression=gene_dict[gene][1]
	#Caps gene expression at 10x enriched
	if gene_expression<10:
		h3_val=bw.stats('chr13', start, end, type='sum')
		#corrects for nonetype
		if h3_val is None:
			h3_val=0
		b_h327+=h3_val
		b_exp.append(gene_expression)
#subplot for Compartment A

fig,ax=plt.subplots(1,1)
ax.scatter(a_exp,a_h327,alpha=0.5)
ax.set_title("Expression vs. H327me3 in Compartment A")
ax.set_xlabel('Gene Expression')
ax.set_ylabel("H327me3 Score")
plt.savefig("A_Compartment_H3K27me3.png")

#subplot for Compartment B
fig,ax=plt.subplots(1,1)
ax.scatter(b_exp,b_h327,alpha=0.5)
ax.set_title("Expression vs. H327me3 in Compartment B")
ax.set_xlabel('Gene Expression')
ax.set_ylabel("H327me3 Score")
plt.savefig("B_Compartment_H3K27me3.png")
