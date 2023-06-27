import sys
import collections
import csv
import gzip
import os
import json
import subprocess
import statistics
import operator
from multiprocessing import Pool as ThreadPool 
import os.path
import time
import random


input_file = 'list_SRA_ID.txt'
n_threads = 1

global directory
directory = './reads/'


def mapping(SRA_ID):
    ## get FTP link(s)

    info = True
    try:
        subprocess.check_output('ffq ' + SRA_ID + ' -o ' + SRA_ID + '.json 2>&1', shell=True)
    except:
        info = False  
    
    ## extract links from json file
    if info:
        f = open(SRA_ID + '.json')
        data = json.load(f)
        l_links = list()
        
        if SRA_ID[0:3] == "ERS":
            for k,v in data[SRA_ID]['experiments'].items():
                name = k
            for k,v in data[SRA_ID]['experiments'][name]['runs'].items():
                err_name = k
            file_name = data[SRA_ID]['experiments'][name]['runs'][err_name]
        else:
            file_name = data[SRA_ID]
        
        #for d in data[SRA_ID]['files']['ftp']:
        for d in file_name['files']['ftp']:
            if d['filetype'] == 'fastq':
                link = d['url']
                l_links.append(link)
    
        ## case no link or more than 2 links
        if not l_links or len(l_links) > 2:
            print('error: found ' + str(len(l_links)) + ' FTP links for run ' + SRA_ID)
    
        else:
            ## download reads
            for link in l_links:
                if '_1.fastq.gz' in link:
                    subprocess.check_output('wget ' + link + ' -O ' + directory + SRA_ID + '_1.fastq.gz 2>&1', shell=True)
                elif '_2.fastq.gz' in link:
                    subprocess.check_output('wget ' + link + ' -O ' + directory + SRA_ID + '_2.fastq.gz 2>&1', shell=True)
                else:
                    subprocess.check_output('wget ' + link + ' -O ' + directory + SRA_ID + '.fastq.gz 2>&1', shell=True)

        
        os.system('rm ' + SRA_ID + '.json')
        time.sleep(2)
    
    else:
        print('error: the SRA ID ' + SRA_ID + ' seems invalid')




l_final = list()
file_content = csv.reader(open(input_file), delimiter='	')
for line in file_content:
    for SRA_ID in line:
        SRA_ID = line[0]
        l_final.append(SRA_ID)


## process samples        
pool = ThreadPool(n_threads) 
tmp_res = pool.map_async(mapping, l_final, chunksize=1)
output = tmp_res.get()
pool.close() 
pool.join()



