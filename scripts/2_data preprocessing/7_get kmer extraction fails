## determine successful n failed k-mer extraction counts 

# get counts of successful file extractions

import pandas as pd
import csv
import collections 
import os
from tqdm import tqdm
from collections import Counter


readfile_accession_2_ribotype = collections.defaultdict(int)
file_content = csv.reader(open('/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/noNA_downsized.csv'), delimiter=',')
for line in file_content:
    accession = line[0]
    ribotype = line[2]
    if (ribotype == '27') | (ribotype == '78') | (ribotype == '14'):
        readfile_accession_2_ribotype[accession] = int(ribotype)
    
kmer_size = 17
# Define path to input and output files
InputPath  = '/rds/general/user/dds122/ephemeral/kmer_'+str(kmer_size)+'_df/'

present_accession_2_ribotype = collections.defaultdict(int)
absent_accession_2_ribotype = collections.defaultdict(int)

for fn in os.listdir(InputPath):
    # Define full filename
    in_file = os.path.join(InputPath,fn)
    if os.path.isfile(in_file) and fn.endswith('.txt'):
        #print("  Adding: " + fn)
        accession_ = fn.replace('.txt', '')
        if accession_ in readfile_accession_2_ribotype:
            present_accession_2_ribotype[accession_] = readfile_accession_2_ribotype[accession_]
            
print(Counter(present_accession_2_ribotype.values()))
missing_list = list(set(readfile_accession_2_ribotype.keys()) - set(present_accession_2_ribotype.keys()))
for accession_ in missing_list:
    if accession_ in readfile_accession_2_ribotype:
        absent_accession_2_ribotype[accession_] = readfile_accession_2_ribotype[accession_]
Counter(absent_accession_2_ribotype.values())
