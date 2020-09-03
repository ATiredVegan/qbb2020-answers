#!/usr/bin/env python3

from fasta_iterator_class import FASTAReader

import sys
#gest user import
target_name=sys.argv[1]
query_name=sys.argv[2]
k = int(sys.argv[3])

#imports files from user command 
target=FASTAReader(open(target_name))
query=FASTAReader(open(query_name))

target_kmers={}
#iterates through target and establishes a dictionary with each unique location in each of the target sequences
for seq_id, sequence in target:
    for i in range(0, len(sequence) - k+1):
        kmer = sequence[i:i + k].upper()
        target_kmers.setdefault(kmer, [])
        target_kmers[kmer].append((i,seq_id))
#format output string
output=["target_sequence_name \t target_start \t query_start \t k-mer"]
#iterates through query and tries to find matches. Each match=newline
for seq_id,sequence in query:
    for start in range(0,len(sequence)-k+1):
        kmer=sequence[start:start+k].upper()
        if kmer in target_kmers:
            for entry in target_kmers[kmer]:
                output.append("\n"+entry[1]+"\t"+str(entry[0])+"\t"+str(start)+"\t"+kmer)

#returns first 100 lines.                
with open('alignment.dat','w') as f: 
    for line in output[:1001]:
        f.write("".join(line))