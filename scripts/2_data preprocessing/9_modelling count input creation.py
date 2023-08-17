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
    if ribotype_list.count(',') == 0: # if only one ribotype exists for this kmer id. # UNIQUE KMERS; remove this to get all library kmers
        library_id_2_ribotype[id_] = ribotype_list.replace("{'", '').replace("'}", '') 
        kmer_list = list(library_kmer_2_id.keys())
        id_list = list(library_kmer_2_id.values())
        kmer_position = id_list.index(id_)
        kmer_ = kmer_list[kmer_position]
        library_kmer_2_ribotype[kmer_] = ribotype_list.replace("{'", '').replace("'}", '')
        
################################################## get kmers counts list in read files ##############################################################
# combine all txt one huge dictionary
import json
import os
import math
from tqdm import tqdm

# Define output filename
OutputFilename = 'count_new10_kmer'+str(kmer_size)+'_dict.txt'

# Define path to input and output files
InputPath  = '/rds/general/user/dds122/ephemeral/kmer_list_'+str(kmer_size)+'_df/'
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
            id_count = collections.defaultdict(list)
            for line in content:
                insert = line.strip().split('\t')
                id_ = insert[0]
                if id_ in library_id_2_ribotype: # if kmer exists in library; replacing it with shrunk_library_id_2_ribotype does not give sane result!
                    if id_ in id_count:
                        id_count[id_].extend(insert[1])
                    else:
                        id_count[id_] = insert[1]
            if len(id_count) != 0:
                readfile_dict[fn.replace('.txt', '')] = id_count


# save to txt file
out = open(out_file, 'w+')
for id_, count_ in readfile_dict.items():
    out.write(str(id_) + '	' + str(count_) + '\n')
out.close()     
################################################## get kmers counts per readfile as dict ##############################################################

readfile_by_kmer = dict()

for k, v in tqdm(readfile_dict.items()): # for every readfile
    kmer_count = collections.defaultdict(list)
    for id_, count_ in v.items():
        kmer_list = list(library_kmer_2_id.keys())
        id_list = list(library_kmer_2_id.values())
        kmer_position = id_list.index(id_)
        kmer_ = kmer_list[kmer_position]
        kmer_count[kmer_] = count_.replace("'", '')

    missing_kmers = list(set(library_kmer_2_ribotype.keys()) - set(kmer_count.keys()))
    for missing_kmer in missing_kmers:
        kmer_count[missing_kmer] = '[0]'
                                       
    readfile_by_kmer[k] = kmer_count

################################################## account for files with no kmers belonging to library ##############################################################
for file in os.listdir('/rds/general/user/dds122/ephemeral/kmer_'+str(kmer_size)+'_df/'):
    filename = file.replace('.txt', '')
    if filename not in readfile_by_kmer.keys() and file.endswith('.txt'):
        kmer_count = collections.defaultdict(int)
        for kmer_ in library_kmer_2_ribotype.keys():
            kmer_count[kmer_] = '[0]' # tried w 0
        readfile_by_kmer[filename] = kmer_count
        
################################################## get sum/mean/median counts of unique k-mers from prev step ##############################################################
import statistics
import pandas as pd

readfile_by_kmer_sumcount = dict()
readfile_by_kmer_meancount = dict()
readfile_by_kmer_mediancount = dict()


for k, v in tqdm(readfile_by_kmer.items()): # for every readfile
    # initialie dicts  
    kmer_sumcount = collections.defaultdict(int)
    kmer_meancount = collections.defaultdict(int)
    kmer_mediancount = collections.defaultdict(int)
    
    for kmer_, count_ in v.items():
        count_ = count_.strip('][').split(', ')
        count_ = [ int(x) for x in count_ ] # convert to int
        
        sum_count = sum(count_)
        mean_count = statistics.mean(count_)
        median_count = statistics.median(count_)
              
        kmer_sumcount[kmer_] = sum_count
        kmer_meancount[kmer_] = mean_count
        kmer_mediancount[kmer_] = median_count

    readfile_by_kmer_sumcount[k] = kmer_sumcount
    readfile_by_kmer_meancount[k] = kmer_meancount    
    readfile_by_kmer_mediancount[k] = kmer_mediancount  

################################################## export data in x and y, ready for prediction ##############################################################

## x_train index
training = pd.read_csv('/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/X_train_readfile_by_absence_kmer_df_k17.csv', index_col=0)
testing = pd.read_csv('/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/X_test_readfile_by_absence_kmer_df_k17.csv', index_col=0)

# prep and export x data 
readfile_by_kmer_sumcount_df = pd.DataFrame.from_dict(readfile_by_kmer_sumcount,orient='index')
X_train = readfile_by_kmer_sumcount_df.loc[training.index]
X_test = readfile_by_kmer_sumcount_df.loc[testing.index]
X_train.to_csv('/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/10ribo_X_train_readfile_by_kmer_sumcount_df_k'+str(kmer_size)+'.csv', index = True, header = True)
X_test.to_csv('/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/10ribo_X_test_readfile_by_kmer_sumcount_df_k'+str(kmer_size)+'.csv', index = True, header = True)

readfile_by_kmer_meancount_df = pd.DataFrame.from_dict(readfile_by_kmer_meancount,orient='index')
X_train = readfile_by_kmer_meancount_df.loc[training.index]
X_test = readfile_by_kmer_meancount_df.loc[testing.index]
X_train.to_csv('/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/10ribo_X_train_readfile_by_kmer_meancount_df_k'+str(kmer_size)+'.csv', index = True, header = True)
X_test.to_csv('/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/10ribo_X_test_readfile_by_kmer_meancount_df_k'+str(kmer_size)+'.csv', index = True, header = True)

readfile_by_kmer_mediancount_df = pd.DataFrame.from_dict(readfile_by_kmer_mediancount,orient='index')
X_train = readfile_by_kmer_mediancount_df.loc[training.index]
X_test = readfile_by_kmer_mediancount_df.loc[testing.index]
X_train.to_csv('/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/10ribo_X_train_readfile_by_kmer_mediancount_df_k'+str(kmer_size)+'.csv', index = True, header = True)
X_test.to_csv('/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/10ribo_X_test_readfile_by_kmer_mediancount_df_k'+str(kmer_size)+'.csv', index = True, header = True)

# prep and export y data
df = pd.read_csv('/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/noNA_downsized.csv')
df = df.loc[df['accession'].isin(readfile_by_kmer_sumcount.keys())] 
df = df.drop(columns='source_code')
df = df.set_index('accession')
df.index.name = None
y_train = df.loc[training.index]
y_test = df.loc[testing.index]
y_train.to_csv('/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/10ribo_y_train_readfile_by_all_kmer_trueribotype_df_k'+str(kmer_size)+'.csv', index = True, header = True)
y_test.to_csv('/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/10ribo_y_test_readfile_by_all_kmer_trueribotype_df_k'+str(kmer_size)+'.csv', index = True, header = True)
