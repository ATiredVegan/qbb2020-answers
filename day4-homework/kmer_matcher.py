#!/usr/bin/env python3

from fasta_iterator_class import FASTAReader

import sys
target_name=sys.argv[1]
query_name=sys.argv[2]
k = int(sys.argv[3])

target=FASTAReader(open("subset.fa"))


target_kmers={}
for seq_id, sequence in target:
    for i in range(0, len(sequence) - k+1):
        kmer = sequence[i:i + k].upper()
        target_kmers.setdefault(kmer, [])
        target_kmers[kmer].append((i,seq_id))
query=FASTAReader(open("droYak2_seq.fa"))
output=["target_sequence_name \t target_start \t query_start \t k-mer"]
for seq_id,sequence in query:
    for start in range(0,len(sequence)-k+1):
        kmer=sequence[start:start+k].upper()
        if kmer in target_kmers:
            for entry in target_kmers[kmer]:
                output.append("\n"+entry[1]+"\t"+str(entry[0])+"\t"+str(start)+"\t"+kmer)

                
with open('alignment.dat','w') as f: 
    for line in output[:1001]:
        f.write("".join(line))