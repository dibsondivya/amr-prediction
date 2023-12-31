import sys
import subprocess
import os
import pandas as pd
from multiprocessing import Pool as ThreadPool 
from tqdm import tqdm

os.chdir('/rds/general/user/dds122/')
kmer_size = 21
n_threads = 4


# get existing files
file_list = []
for file in os.listdir('/rds/general/user/dds122/ephemeral/kmer_'+str(kmer_size)+'_txt/'):
    if file.endswith('_1.txt'):
        file_list.append(file.replace('_1.txt', ''))
# get total wanted files
df = pd.read_csv('/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/noNA_downsized.csv')
df = df[(df['ribotype'] == 27) | (df['ribotype'] == 78) | (df['ribotype'] == 14)]
SRA_ID_list = df['accession'].tolist()


def mapping(accession):
    filename = accession+'_1.fastq.gz'
    try:
        subprocess.check_output('kmc -k' + str(kmer_size) + ' ' + 'projects/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/reads/' + filename + ' ephemeral/kmer_'+str(kmer_size)+'/' + filename + ' ephemeral/tmp', shell=True)
        subprocess.check_output('kmc_tools transform ephemeral/kmer_'+str(kmer_size)+'/'+filename+ ' dump ephemeral/kmer_'+str(kmer_size)+'_txt/' + accession+'_1.txt', shell=True)
    except subprocess.CalledProcessError as error:
        sys.stdout.write("%s\n" % (str(error)))
    pbar.update(1)

input_list = list(set(SRA_ID_list) - set(file_list)) # find missing _1 files
pbar = tqdm(total=len(input_list))

## process samples        
pool = ThreadPool(n_threads) 
tmp_res = pool.map_async(mapping, input_list, chunksize=None)
output = tmp_res.get()
pool.close() 
pool.join() 
