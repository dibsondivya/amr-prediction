import sys
import subprocess
import pandas as pd
import os
from multiprocessing import Pool as ThreadPool  

kmer_size = 21 # ACTIONABLE: to be declared

os.chdir('/rds/general/user/dds122/ephemeral/kmer_'+str(kmer_size)+'_txt/')

# import data containing accession | ribotype of all top 10 ribotype samples
df = pd.read_csv('/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/noNA_downsized.csv')
n_threads = 1

def create_csv(file):
    sample_id = 0
    kmer_id = {}
    id_count = {}
    id_sample = {}
    id_ribotype = {}

    if "fastq.gz" in file:
        filename = file.replace('.fastq.gz.txt', '')
    else:
        filename = file.replace('.txt', '')

    output = open(file, "r")
    Lines = output.readlines()
    for line in Lines:
        if len(line.strip().split('\t')) == 2:
            (kmer, counts) = line.strip().split('\t')
            kmer_id[kmer] = sample_id
            id_count[sample_id] = int(counts)
            id_sample[sample_id] = filename
            if "_1" in filename:
                accession = filename.replace('_1', '')
            else:
                accession = filename.replace('_2', '')
            id_ribotype[sample_id] = int(df[df['accession'] == accession]['ribotype'])
            sample_id += 1

    kmer_id_df = pd.DataFrame(list(kmer_id.items()),columns = ['kmer','id']) 
    kmer_id_df.to_csv('/rds/general/user/dds122/ephemeral/kmer_21_df/kmer_id/' + str(filename) +'.csv', index=False)
    id_count_df = pd.DataFrame(list(id_count.items()),columns = ['id','count']) 
    id_count_df.to_csv('/rds/general/user/dds122/ephemeral/kmer_21_df/id_count/' + str(filename) +'.csv', index=False)
    id_sample_df = pd.DataFrame(list(id_sample.items()),columns = ['id','sample']) 
    id_sample_df.to_csv('/rds/general/user/dds122/ephemeral/kmer_21_df/id_sample/' + str(filename) +'.csv', index=False)
    id_ribotype_df = pd.DataFrame(list(id_ribotype.items()),columns = ['id','ribotype']) 
    id_ribotype_df.to_csv('/rds/general/user/dds122/ephemeral/kmer_21_df/id_ribotype/' + str(filename) +'.csv', index=False)
    
file_list = []
for file in os.listdir('/rds/general/user/dds122/ephemeral/kmer_21_txt/'):
    file_list.append(file)
    
## process samples        
pool = ThreadPool(n_threads) 
tmp_res = pool.map_async(create_csv, file_list, chunksize=None)
output = tmp_res.get()
pool.close() 
pool.join() 
