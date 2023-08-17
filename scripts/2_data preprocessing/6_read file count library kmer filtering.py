import sys
import subprocess
import os
import csv
from multiprocessing import Pool as ThreadPool  
import collections
from tqdm import tqdm


kmer_size = 31 # ACTIONABLE: to be declared
os.chdir('/rds/general/user/dds122/ephemeral/kmer_'+str(kmer_size)+'_txt/')
n_threads = 4

# load library kmers
global kmer_2_id
kmer_2_id = dict()
file_content = csv.reader(open('/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/kmer_id'+str(kmer_size)+'.csv'), delimiter=',')
for line in file_content:
    kmer = line[0]
    id_ = line[1]
    kmer_2_id[kmer] = id_
    


def create_csv(file):
    id_count = collections.defaultdict(list)

    accession = file.replace('_1.txt', '')
    
    output = open(file, "r")
    Lines = output.readlines()
    for line in Lines:
        insert = line.strip().split('\t')
        if insert and len(insert) == 2:
            kmer = insert[0]
            if kmer in kmer_2_id:
                id_ = kmer_2_id[kmer]
                id_count[id_] += int(insert[1])
    
    try:
        output = open(file.replace('_1.txt','_2.txt'), "r")
        Lines = output.readlines()
        for line in Lines:
            insert = line.strip().split('\t')
            if insert and len(insert) == 2:
                kmer = insert[0]
                if kmer in kmer_2_id:
                    id_ = kmer_2_id[kmer]
                    id_count[id_].extend(insert[1])
    except:
        print('Error with file ' + file.replace('_1.txt','_2.txt'))

    out = open('/rds/general/user/dds122/ephemeral/kmer_'+str(kmer_size)+'_df/' + accession + '.txt', 'w+')
    for id_, count_ in id_count.items():
        out.write(str(id_) + '	' + count_ + '\n')
    out.close()

    
    pbar.update(1)
    
   
file_list = []
for file in os.listdir('/rds/general/user/dds122/ephemeral/kmer_'+str(kmer_size)+'_txt/'):
    if file.endswith('_1.txt'):
        file_list.append(file)

pbar = tqdm(total=len(file_list))


## process samples        
pool = ThreadPool(n_threads) 
tmp_res = pool.map_async(create_csv, file_list)
output = tmp_res.get()
pool.close() 
pool.join() 
