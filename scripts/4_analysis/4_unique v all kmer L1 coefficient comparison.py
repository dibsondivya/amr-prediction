import shap
from matplotlib_venn import venn2
import pandas as pd
import numpy as np
import joblib
import matplotlib.pyplot as plt
import os
from sklearn.linear_model import LogisticRegression

label_list = [1, 2, 5, 14, 15, 17, 20, 27, 78, 106]
filepath = '/rds/general/user/dds122/projects/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/'
kmer = 'kmer'
predicting = 'binary'

k = '17'
method = 'l1'
unique_X_train = pd.read_csv(filepath+'10ribo_X_train_readfile_by_'+kmer+'_absence_df_k'+k+'.csv', index_col=0)
unique_X_test= pd.read_csv(filepath+'10ribo_X_test_readfile_by_'+kmer+'_absence_df_k'+k+'.csv', index_col=0)
y_train = pd.read_csv(filepath+'10ribo_y_train_readfile_by_all_kmer_trueribotype_df_k'+k+'.csv', index_col=0)
y_test = pd.read_csv(filepath+'10ribo_y_test_readfile_by_all_kmer_trueribotype_df_k'+k+'.csv', index_col=0)

kmer = 'all'
all_X_train = pd.read_csv(filepath+'10ribo_X_train_readfile_by_'+kmer+'_absence_df_k'+k+'.csv', index_col=0)
all_X_test= pd.read_csv(filepath+'10ribo_X_test_readfile_by_'+kmer+'_absence_df_k'+k+'.csv', index_col=0)

for i in range(len(label_list)):
    kmer = 'kmer'
    model = joblib.load(filepath+'models/k'+k+'_'+kmer+'_'+predicting+'_'+method+'.joblib') # load ovo model
    regular_l1 = LogisticRegression(**model.estimator.get_params()) # feed optimized parameters into logreg 
    regular_l1.fit(unique_X_train, y_train)
    unique_kmer_list = unique_X_train.iloc[:,regular_l1.coef_[i]!=0].columns.tolist()

    kmer = 'all'
    model = joblib.load(filepath+'models/k'+k+'_'+kmer+'_'+predicting+'_'+method+'.joblib') # load ovo model
    all_l1 = LogisticRegression(**model.estimator.get_params()) # feed optimized parameters into logreg 
    all_l1.fit(all_X_train, y_train)
    all_kmer_list = all_X_train.iloc[:,all_l1.coef_[0]!=0].columns.tolist()
    
    venn2([set(unique_kmer_list), set(all_kmer_list)], set_labels=('Unique', 'All'))
    plt.title('Comparing Selected K-mers from L1 Logistic Regression for K='+k+" for Ribotype "+str(label_list[i]))
    plt.savefig(filepath+'venn/k'+k+'_'+predicting+'_'+method+'_ribo'+str(label_list[i])+'.png')
    plt.show()

