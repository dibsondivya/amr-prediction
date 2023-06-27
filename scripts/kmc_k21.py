import sys
import subprocess
import os
import pandas as pd
from multiprocessing import Pool as ThreadPool 

df = pd.read_csv('/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/noNA_downsized.csv')
# subset to only top three ribotypes
df = df[(df['ribotype'] == 27) | (df['ribotype'] == 78) | (df['ribotype'] == 14)]
SRA_ID_list = df['accession'].tolist()

sample_id = 0
kmer_id = {}
id_count = {}
id_sample = {}
id_ribotype = {}

kmer_size = 21 # ACTIONABLE: to be changed as required; 11, 17, 21 and 31
n_threads = 1

def mapping(accession):
    filename = accession+'_1.fastq.gz'
    try:
        subprocess.check_output('kmc -k' + str(kmer_size) + ' ' + 'reads/' + filename + ' /rds/general/user/dds122/ephemeral/kmer_'+str(kmer_size)+'/' + filename + ' /rds/general/user/dds122/ephemeral/tmp', shell=True)
        subprocess.check_output('kmc_tools transform /rds/general/user/dds122/ephemeral/kmer_'+str(kmer_size)+'/'+filename+ ' dump /rds/general/user/dds122/ephemeral/kmer_'+str(kmer_size)+'_txt/' + accession+'_1.txt', shell=True)
    except subprocess.CalledProcessError as error:
        sys.stdout.write("%s\n" % (str(error)))
    filename = accession+'_2.fastq.gz'
    try:
        subprocess.check_output('kmc -k' + str(kmer_size) + ' ' + 'reads/' + filename + ' /rds/general/user/dds122/ephemeral/kmer_'+str(kmer_size)+'/' + filename + ' /rds/general/user/dds122/ephemeral/tmp', shell=True)
        subprocess.check_output('kmc_tools transform /rds/general/user/dds122/ephemeral/kmer_'+str(kmer_size)+'/'+filename+ ' dump /rds/general/user/dds122/ephemeral/kmer_'+str(kmer_size)+'_txt/' + accession+'_2.txt', shell=True)
    except subprocess.CalledProcessError as error:
        sys.stdout.write("%s\n" % (str(error)))

## process samples        
pool = ThreadPool(n_threads) 
tmp_res = pool.map_async(mapping, SRA_ID_list, chunksize=None)
output = tmp_res.get()
pool.close() 
pool.join() 
