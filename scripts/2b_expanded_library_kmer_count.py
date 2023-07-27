import os
import sys
import subprocess
import pandas as pd
from sklearn.preprocessing import MultiLabelBinarizer


# convert multi fasta to single line fasta
#Cmd = 'seqtk seq pone.0106545.s006.fasta > single_pone.fasta'
## seqtk seq <multiline fasta> > <nameofnew>

# run kmer counter
os.chdir('/rds/general/projectdf = df[(df['ribotype'] == 106)]
/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/')
k = 31
fastaFile = 'single_pone_added.fasta'
#kmerCmd = 'kmer-counter --fasta --k=%d %s > kmer_count_output%d.txt' % (k, fastaFile, k)
#subprocess.check_output(kmerCmd, shell=True)

# open output of kmer count
with open('/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/kmer_count_output'+str(k)+'.txt', 'r') as file:
#with open('kmer_count_output'+str(k)+'.txt', 'r') as file:
    data = file.read().replace('\n\n', '\n')

dict_output = {}
count = 0
for line in data.split('\n'):
    if len(line.strip().split('\t')) == 2:
        (sequence_id, count_dict) = line.strip().split('\t')
        sequence_id = sequence_id[1:]
        kmers = dict((k,int(v)) for (k,v) in [d.split(':') for d in count_dict.split(' ')])
        dict_output[sequence_id] = kmers
    else: # if count is on a new line
        if line.startswith('>'):
            start_of_line = line
            count += 1
        else:
            end_of_line = line
            count += 1
        
        if count%2 == 0:
            new_line = start_of_line + end_of_line
            (sequence_id, count_dict) = new_line.strip().split('\t')
            sequence_id = sequence_id[1:]
            kmers = dict((k,int(v)) for (k,v) in [d.split(':') for d in count_dict.split(' ')])
            dict_output[sequence_id] = kmers
    
# get ribotype information
new_result_df = pd.DataFrame(list(dict_output.items()),columns = ['sequence_id','kmers_counts'])
#new_col = ['012', '
', '029', '029', '033', '063/1', '126', '247', '434', '434', '441/3', '444', '484', '519', 
#           '519', '524', '542', '542', 'Strain630', 'StrainCD196', '012', '012', '012', '012', '012', '012', '029', '029', '029', '029', '033', '033', '033', '033', '045/1', '045', '063', '063', '063', 
#           '063/1', '063/1', '066', '066/2', '066/2', '066', '078/4', '078ecdc', '126', '247', '247', '247', '247', '251', '251', '434', '434', '441', '444', '444', '444', '484', '519', 
#           '519', '519', '524', '524', '524', '542', '542', '542', '542', 'AI-56', 'StrainCD196', 
#           'StrainCD196', '045/1', 'Strain027', 'AI5', 'CloneB-14', '176', 'CloneA-6', '078', 'CloneAT-16', 'CloneAT-3', 'CloneB-12', 'CloneAT-17', 'CloneAT-11', '001', 'CloneAT-7', 'CloneAT-9', 
#           'CloneAT-8', 'CloneAT-14', 'CloneAT-12', 'CloneAT-1', 'CloneA-15', 'CloneB-11']
new_col = ['012', '027', '029', '029', '033', '063/1', '126', '247', '434', '434', '441/3', '444', '484', '519', 
           '519', '524', '542', '542', 'Strain630', 'StrainCD196', '012', '012', '012', '012', '012', '012', '029', '029', '029', '029', '033', '033', '033', '033', '045/1', '045', '063', '063', '063', 
           '063/1', '063/1', '066', '066/2', '066/2', '066', '078/4', '078ecdc', '126', '247', '247', '247', '247', '251', '251', '434', '434', '441', '444', '444', '444', '484', '519', 
           '519', '519', '524', '524', '524', '542', '542', '542', '542', 'AI-56', 'StrainCD196', 
           'StrainCD196', '045/1', 'Strain027', 'AI5', 'CloneB-14', '176', 'CloneA-6', '078', 'CloneAT-16', 'CloneAT-3', 'CloneB-12', 'CloneAT-17', 'CloneAT-11', '001', 'CloneAT-7', 'CloneAT-9', 
           'CloneAT-8', 'CloneAT-14', 'CloneAT-12', 'CloneAT-1', 'CloneA-15', 'CloneB-11', 
           '002', '002', '002', '002', '002', '002', '002', '002', '002', '002', '002', '002', 
           '010', '010', '010', '010', '010', '010', '010', '010', '010', '010', '010', '010', 
           '012', '012', '012', '012', '012', '012', '012', '012', '012', '012', '012', '012', 
           '017', '017', '017', '017', '017', '017', '017', '017', '017', '017', '017', 
           '027', '027', '027', '027', '027', '027', '027', '027', '027', '027', '027', '027', 
           '056', '056', '056', '056', '056', '056', '056', '056', '056', '056', '056', '056', 
           '078', '078', '078', '078', '078', '078', '078', '078', '078', '078', '078', 
           '106', '106', '106', '106', '106', '106', '106', '106', '106', '106', '106', '106', 
           '125', '125', '125', '125', '125', '125', '125', '125', '125', '125', '125', '125', 
           '255', '255', '255', '255', '255', '255', '255', '255', '255', '255', '255']

new_result_df['ribotype'] = new_col
new_result_df

    # since all kmers occur only once per sequence id
new_list_kmers = list()
for i in dict_output.values():
   new_list_kmers.append(list(i.keys()))

new_result_df['list_kmers'] = new_list_kmers
new_result_df

mlb = MultiLabelBinarizer()
new_OHE_df = pd.DataFrame(mlb.fit_transform(new_result_df['list_kmers']),columns=mlb.classes_, index=new_result_df.index)

# create dictionaries
uh_kmer_id_dict = dict()
uh_id_count_dict = dict()

id = 0
for kmer in new_OHE_df.columns.tolist():
    uh_kmer_id_dict[kmer] = id
    count = new_OHE_df[kmer].sum()
    uh_id_count_dict[id] = count
    
    id += 1
    
#len(uh_kmer_id_dict) # 1684

# include reverse complements
def reverse_complement(st):
    nn = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
    return "".join(nn[n] for n in reversed(st))

for kmer in new_OHE_df.columns.tolist():
    rc = reverse_complement(kmer)
    
    if rc in uh_kmer_id_dict:
        print('yes')
        print(rc)
        id_count = uh_kmer_id_dict.get(rc)
        og_kmer_id = uh_kmer_id_dict.get(kmer)
        kmer_count = uh_id_count_dict.get(og_kmer_id)
        uh_id_count_dict[id_count] += kmer_count
    else:
        current_largest_id = uh_kmer_id_dict.get(max(uh_kmer_id_dict, key=uh_kmer_id_dict.get))
        id_count = current_largest_id+1
        uh_kmer_id_dict[rc] = id_count
        og_kmer_id = uh_kmer_id_dict.get(kmer)
        kmer_count = uh_id_count_dict.get(og_kmer_id)
        uh_id_count_dict[id_count] = kmer_count

uh_kmer_id_df = pd.DataFrame(list(uh_kmer_id_dict.items()),columns = ['kmer','id']) 
uh_id_count_df = pd.DataFrame(list(uh_id_count_dict.items()),columns = ['id','count']) 

new_col_kmer = list()
for listofkmers in new_result_df['list_kmers']:
    include_rc_kmer = list()
    for kmer in listofkmers:
        include_rc_kmer.append(kmer)
        include_rc_kmer.append(reverse_complement(kmer))
    new_col_kmer.append(include_rc_kmer)
    
new_result_df['list_kmer_rc'] = new_col_kmer

mlb = MultiLabelBinarizer()
OHE_withrc_df = pd.DataFrame(mlb.fit_transform(new_result_df['list_kmer_rc']),columns=mlb.classes_, index=new_result_df.index)


id_ribotype_dict = dict()
for kmer in uh_kmer_id_df['kmer'].tolist():
    id = int(uh_kmer_id_df[uh_kmer_id_df['kmer'] == kmer]['id'])
    if kmer in new_OHE_df:
        index_list = new_OHE_df.index[new_OHE_df[kmer] == 1].tolist() # get index
        ribotype_list = new_result_df[new_result_df.index.isin(index_list)]['ribotype'].tolist()
        ribotype_list = set(ribotype_list) # get distinct
        id_ribotype_dict[id] = ribotype_list
        
multiple_id_ribotype_dict = dict()
for kmer in uh_kmer_id_df['kmer'].tolist():
    id = int(uh_kmer_id_df[uh_kmer_id_df['kmer'] == kmer]['id'])
    if kmer in OHE_withrc_df:
        index_list = OHE_withrc_df.index[OHE_withrc_df[kmer] == 1].tolist() # get index
        ribotype_list = new_result_df[new_result_df.index.isin(index_list)]['ribotype'].tolist()
        multiple_id_ribotype_dict[id] = ribotype_list
        
id_ribotype_df = pd.DataFrame(list(id_ribotype_dict.items()),columns = ['id','ribotype']) 
id_ribotype_df.to_csv('id_ribotype'+str(k)+'.csv', index=False)
multiple_id_ribotype_df = pd.DataFrame(list(multiple_id_ribotype_dict.items()),columns = ['id','ribotype']) 
multiple_id_ribotype_df.to_csv('multiple_id_ribotype'+str(k)+'.csv', index=False)
uh_kmer_id_df.to_csv('kmer_id'+str(k)+'.csv', index=False)
uh_id_count_df.to_csv('id_count'+str(k)+'.csv', index=False)
new_result_df.to_csv('result_df'+str(k)+'.csv', index=True)
OHE_withrc_df.to_csv('OHE_withrc'+str(k)+'.csv', index=True)
