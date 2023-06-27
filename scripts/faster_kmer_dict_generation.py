import sys
import subprocess
import pandas as pd
import os
from multiprocessing import Pool as ThreadPool  
import collections
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map  # or thread_map


kmer_size = 21 # ACTIONABLE: to be declared

os.chdir('/rds/general/user/dds122/ephemeral/kmer_'+str(kmer_size)+'_txt/')

# import data containing accession | ribotype of all top 10 ribotype samples
df = pd.read_csv('/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/noNA_downsized.csv')
n_threads = 1

# load library kmers
global kmer_2_id
kmer_2_id = pd.read_csv('/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/kmer_id21.csv')

def create_csv(file):
    id_count = collections.defaultdict(int)

    accession = file.replace('_1.txt', '')

    output = open(file, "r")
    Lines = output.readlines()
    for line in Lines:
        insert = line.strip().split('\t')
        if insert and len(insert) == 2:
            kmer = insert[0]
            if kmer in kmer_2_id['kmer'].tolist():
                id_ = int(kmer_2_id[kmer_2_id['kmer']==kmer]['id'])
                id_count[id_] += int(insert[1])

    output = open(file.replace('_1.txt','_2.txt'), "r")
    Lines = output.readlines()
    for line in Lines:
        insert = line.strip().split('\t')
        if insert and len(insert) == 2:
            kmer = insert[0]
            if kmer in kmer_2_id['kmer'].tolist():
                id_ = int(kmer_2_id[kmer_2_id['kmer']==kmer]['id'])
                id_count[id_] += int(insert[1])
    
    
    out = open('/rds/general/user/dds122/ephemeral/kmer_21_df/' + accession + '.txt', 'w+')
    for id_, count_ in id_count.items():
        out.write(str(id_) + '	' + str(count_) + '\n')
    out.close()
    
    pbar.update(1)
    
   
file_list = []
for file in os.listdir('/rds/general/user/dds122/ephemeral/kmer_21_txt/'):
    if file.endswith('_1.txt'):
        file_list.append(file)

pbar = tqdm(total=len(file_list))

## process samples        
pool = ThreadPool(n_threads) 
tmp_res = pool.map_async(create_csv, file_list, chunksize=None)
output = tmp_res.get()
pool.close() 
pool.join() 
