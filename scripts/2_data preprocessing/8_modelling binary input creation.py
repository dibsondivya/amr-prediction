################################################## get unique kmers in library ##############################################################
#### get kmers that are unique to ribotype in library
import sys
import subprocess
import pandas as pd
import os
import csv
from multiprocessing import Pool as ThreadPool  
import collections
from collections import Counter
from tqdm import tqdm

# initialize
kmer_size = 21 # ACTIONABLE: to be declared
os.chdir('/rds/general/user/dds122/ephemeral/kmer_'+str(kmer_size)+'_txt/')
n_threads = 4

# load library kmers
global library_kmer_2_id
library_kmer_2_id = dict()
file_content = csv.reader(open('/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/kmer_id'+str(kmer_size)+'.csv'), delimiter=',')
for line in file_content:
    kmer = line[0]
    id_ = line[1]
    library_kmer_2_id[kmer] = id_

global library_id_2_ribotype
library_id_2_ribotype = dict()
global library_kmer_2_ribotype
library_kmer_2_ribotype = dict()
file_content = csv.reader(open('/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/id_ribotype'+str(kmer_size)+'.csv'), delimiter=',')
for line in file_content:
    id_ = line[0]
    ribotype_list = line[1]
    #print(ribotype_list)
    #print(ribotype_list.count(','))
    if ribotype_list.count(',') == 0: # UNIQUE KMERS; remove this to get all library kmers
        library_id_2_ribotype[id_] = ribotype_list.replace("['", '').replace("']", '') # for k=21 only
        #library_id_2_ribotype[id_] = ribotype_list.replace("{'", '').replace("'}", '') # for k=17 and k=31
        kmer_list = list(library_kmer_2_id.keys())
        id_list = list(library_kmer_2_id.values())
        kmer_position = id_list.index(id_)
        kmer_ = kmer_list[kmer_position]
        library_kmer_2_ribotype[kmer_] = ribotype_list.replace("['", '').replace("']", '') # for k=21 only
        #library_kmer_2_ribotype[kmer_] = ribotype_list.replace("{'", '').replace("'}", '') # for k=17 and k=31

# get distinct library kmers per ribotype
print(Counter(library_kmer_2_ribotype.values()))

################################################## get kmers in read files ##############################################################
# combine all txt one huge dictionary
import json
import os
import math
from tqdm import tqdm

# Define output filename
OutputFilename = 'new10_kmer'+str(kmer_size)+'_dict.txt'

# Define path to input and output files
InputPath  = '/rds/general/user/dds122/ephemeral/kmer_'+str(kmer_size)+'_df/'
OutputPath = '/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/'

# Define output file
out_file = os.path.join(OutputPath,OutputFilename)

readfile_dict = dict()

def append_record(fn, record):
    with open(fn, 'a') as f:
        json.dump(record, f)

# Loop through each file in input directory
for fn in tqdm(os.listdir(InputPath)):
    # Define full filename
    in_file = os.path.join(InputPath,fn)
    if os.path.isfile(in_file) and fn.endswith('.txt'):
        with open(in_file, 'r') as file_in:
            content = file_in.readlines()
            id_present = collections.defaultdict(int)
            for line in content:
                insert = line.strip().split('\t')
                id_ = insert[0]
                if id_ in library_id_2_ribotype: # if kmer exists in library; replacing it with shrunk_library_id_2_ribotype does not give sane result!
                      id_present[id_] = 1
            if len(id_present) != 0:
                readfile_dict[fn.replace('.txt', '')] = id_present

# save to txt file
out = open(out_file, 'w+')
for id_, presence_ in id_present.items():
    out.write(str(id_) + '	' + str(presence_) + '\n')
out.close()     

################################################## get per kmer counts of unique k-mers per readfile ##############################################################
readfile_by_kmer = dict()

for k, v in readfile_dict.items(): # for every readfile
    kmer_count = collections.defaultdict(int)
    for id_, presence_ in v.items():
        kmer_list = list(library_kmer_2_id.keys())
        id_list = list(library_kmer_2_id.values())
        kmer_position = id_list.index(id_)
        kmer_ = kmer_list[kmer_position]
        extracted_count = presence_
        #extracted_count = int(extracted_count)
        kmer_count[kmer_] = extracted_count

    missing_kmers = list(set(library_kmer_2_ribotype.keys()) - set(kmer_count.keys()))
    for missing_kmer in missing_kmers:
        kmer_count[missing_kmer] = 0
    
    readfile_by_kmer[k] = kmer_count
    
################################################## account for files with no kmers belonging to library ##############################################################
for file in os.listdir('/rds/general/user/dds122/ephemeral/kmer_'+str(kmer_size)+'_df/'):
    filename = file.replace('.txt', '')
    if filename not in readfile_by_kmer.keys() and file.endswith('.txt'):
        kmer_count = collections.defaultdict(int)
        for kmer_ in library_kmer_2_ribotype.keys():
            kmer_count[kmer_] = 0 # tried w 0
        readfile_by_kmer[filename] = kmer_count
################################################## export data in x and y, ready for prediction ##############################################################

# prep and export x data 
## x_train index
training = pd.read_csv('/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/X_train_readfile_by_absence_kmer_df_k17.csv', index_col=0)
testing = pd.read_csv('/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/X_test_readfile_by_absence_kmer_df_k17.csv', index_col=0)

readfile_by_kmer_df = pd.DataFrame.from_dict(readfile_by_kmer,orient='index')
X_train = readfile_by_kmer_df.loc[training.index]
X_test = readfile_by_kmer_df.loc[testing.index]
X_train.to_csv('/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/10ribo_X_train_readfile_by_kmer_absence_df_k'+str(kmer_size)+'.csv', index = True, header = True)
X_test.to_csv('/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/10ribo_X_test_readfile_by_kmer_absence_df_k'+str(kmer_size)+'.csv', index = True, header = True)

# prep and export y data
df = pd.read_csv('/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/noNA_downsized.csv')
df = df.loc[df['accession'].isin(readfile_by_kmer_df.index.tolist())] 
df = df.drop(columns='source_code')
df = df.set_index('accession')
df.index.name = None
df = df.loc[readfile_by_kmer_df.index]
df.to_csv('/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/readfile_by_absence_kmer_trueribotype_df_k'+str(kmer_size), index = True, header = True)

