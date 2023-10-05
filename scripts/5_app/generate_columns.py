import pandas as pd
import csv

k_sizes = ['17', '21']
for k in k_sizes:
    training = pd.read_csv('library/training_k'+k+'.csv', index_col = 0)
    out = open('library/k'+k+'_col_list.txt', 'w+')
    for kmer_ in training.columns.tolist():
        out.write(kmer_ + '\n')
    out.close()