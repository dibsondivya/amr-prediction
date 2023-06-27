import pandas as pd
import os

kmer_size = 17

sample_id = 0
kmer_id = {}
id_count = {}
id_sample = {}
id_ribotype = {}

df = pd.read_csv('/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/noNA_downsized.csv')

os.chdir('/rds/general/user/dds122/ephemeral/kmer_'+str(kmer_size)+'_txt/')

for file in os.listdir('/rds/general/user/dds122/ephemeral/kmer_'+str(kmer_size)+'_txt/'):
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
kmer_id_df.to_csv('/rds/general/project/hda-22-23/live/Summer_projects/dds122/kmer_id_k' + str(kmer_size) +'.csv', index=False)
id_count_df = pd.DataFrame(list(id_count.items()),columns = ['id','count']) 
id_count_df.to_csv('/rds/general/project/hda-22-23/live/Summer_projects/dds122/id_count_k' + str(kmer_size) +'.csv', index=False)
id_sample_df = pd.DataFrame(list(id_sample.items()),columns = ['id','sample']) 
id_sample_df.to_csv('/rds/general/project/hda-22-23/live/Summer_projects/dds122/id_sample_k' + str(kmer_size) +'.csv', index=False)
id_ribotype_df = pd.DataFrame(list(id_ribotype.items()),columns = ['id','ribotype']) 
id_ribotype_df.to_csv('/rds/general/project/hda-22-23/live/Summer_projects/dds122/id_ribotype_k' + str(kmer_size) +'.csv', index=False)
