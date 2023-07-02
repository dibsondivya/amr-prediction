# combine all txt one huge dictionary
import json
import os
import math
from tqdm import tqdm

kmer_size = 21

# Define output filename
OutputFilename = 'top3_kmer'+str(kmer_size)+'_dict.txt'

# Define path to input and output files
InputPath  = '/rds/general/user/dds122/ephemeral/kmer_'+str(kmer_size)+'_df/'
OutputPath = '/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/'

# Define output file
out_file = os.path.join(OutputPath,OutputFilename)

id_count = collections.defaultdict(int)

def append_record(fn, record):
    with open(fn, 'a') as f:
        json.dump(record, f)

# Loop through each file in input directory
for fn in tqdm(os.listdir(InputPath)):
    # Define full filename
    in_file = os.path.join(InputPath,fn)
    if os.path.isfile(in_file) and fn.endswith('.txt'):
        print("  Adding: " + fn)
        with open(in_file, 'r') as file_in:
            content = file_in.readlines()
            for line in content:
                insert = line.strip().split('\t')
                id_ = insert[0]
                if id_ in id_count: # add to dictionary
                    id_count[id_] += int(insert[1])
                else:
                    id_count[id_] = int(insert[1])

# save to txt file
out = open(out_file, 'w+')
for id_, count_ in id_count.items():
    out.write(str(id_) + '	' + str(count_) + '\n')
out.close()  
