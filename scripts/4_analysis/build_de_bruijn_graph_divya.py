import random
import os
import csv


def read_fasta(fasta_content):
    name, seq = None, []
    for line in fasta_content:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name.replace('>',''), ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name.replace('>',''), ''.join(seq))


def rev_compl(seq):
    reverse_code = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(reverse_code[n] for n in reversed(seq))


     

ribotype = '999255'
kmer_size = '21'


colors = {'L1 Logistic Regression': 'red', 'L2 Logistic Regression': 'orange', 'Random Forest': 'blue'}

ribotype_colors = {'010': 'pink', '027': 'red', '078': 'orange', '106':'blue', '125':'green'}

# get sequences
all_seq = list()
with open('ribotype_' + ribotype + '.fas') as fasta_content:
    for name, seq in read_fasta(fasta_content):
        all_seq.append(seq)

print(len(all_seq))

ribotype_kmers = dict()


# get reference kmers
d_good_kmers = dict()
file_content = csv.reader(open('with_truth_collated_important_kmers_k' + kmer_size + '.csv'), delimiter=',')
for line in file_content: 
    if line[0] != '':
        #if line[1] == ribotype:
        if 1 == 1:
            method = line[2]
            kmer = line[3]
            size = 50 + int(round(100 * float(line[4])))
            truth_ribotype = line[5]
            ribotype_kmers[kmer] = truth_ribotype

            if truth_ribotype in ribotype_colors:
                truth_ribotype_color = ribotype_colors[truth_ribotype]
            else:
                truth_ribotype_color = 'black'
            # save it
            if not kmer in d_good_kmers:
                d_good_kmers[kmer] = [colors[method], size, truth_ribotype_color]
            else:
                d_good_kmers[kmer][0] = 'black'
            # save reverse-complement
            revcompl = rev_compl(kmer)
            if not revcompl in d_good_kmers:
                d_good_kmers[revcompl] = [colors[method], size, truth_ribotype_color]
            else:
                d_good_kmers[revcompl][0] = 'black'
            
            ribotype_kmers[revcompl] = truth_ribotype

print(len(d_good_kmers))


# build graph
out = open('graph_k' + kmer_size + '_ribotype_' + ribotype + '.txt', 'w+')
out.write('node1	node2	color	size	truth_ribotype\n')

edges_done = set()
nodes_done = set()

int_kmer_size = int(kmer_size)
for seq in all_seq:

    for kmer in d_good_kmers:
        if kmer in seq:
            print(kmer)
            print(ribotype_kmers[kmer])

    previous_kmer = ''
        
    # extract kmers
    for n in range(len(seq) - int_kmer_size +1):
        kmer = seq[n:n + int_kmer_size]
        # save it in graph
        if previous_kmer == '':
            previous_kmer = kmer
            # save node color if a good kmer
            if kmer in d_good_kmers and not kmer in nodes_done:
                out.write(kmer + '	' + kmer + '	' + d_good_kmers[kmer][0] + '	' + str(d_good_kmers[kmer][1]) +  '	' + d_good_kmers[kmer][2] +'\n')
                nodes_done.add(kmer)
        else:
            edge = previous_kmer + kmer
            if not edge in edges_done:
                out.write(previous_kmer + '	' + kmer + '	\n')
                edges_done.add(edge)
            previous_kmer = kmer
            # save node color if a good kmer
            if kmer in d_good_kmers and not kmer in nodes_done:
                out.write(kmer + '	' + kmer + '	' + d_good_kmers[kmer][0] + '	' + str(d_good_kmers[kmer][1]) + '	' + d_good_kmers[kmer][2] +'\n')
                nodes_done.add(kmer)
                

out.close()






