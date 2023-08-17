import shap
import venn
import pandas as pd
import numpy as np
import joblib
import matplotlib.pyplot as plt
import os
from sklearn.linear_model import LogisticRegression


## note: analysing rf, l1 and l2 only
kmer = '21'
predicting = 'all' # or 'kmer'

readfile_by_kmer_df = pd.read_csv('/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/10ribo_X_test_readfile_by_all_absence_df_k'+str(kmer_size)+'.csv', index_col = 0)
X_train = pd.read_csv('/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/10ribo_X_train_readfile_by_all_absence_df_k'+str(kmer_size)+'.csv', index_col = 0)
y_train = pd.read_csv('/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/10ribo_y_train_readfile_by_all_kmer_trueribotype_df_k'+str(kmer_size)+'.csv', index_col = 0)

feature_imp = pd.DataFrame(columns=['class', 'model', 'features', 'importance']) # initialize

########################## for l1 ##########################
method = 'l1'

model = joblib.load(filepath+'models/k'+k+'_'+kmer+'_'+predicting+'_'+method+'.joblib') # load model
log_reg = LogisticRegression(**model.estimator.get_params()) # feed optimized parameters into logreg
log_reg.fit(X_train, y_train)
explainer = shap.LinearExplainer(log_reg, readfile_by_kmer_df)
shap_values = explainer.shap_values(readfile_by_kmer_df)

def shapley_feature_ranking(j, X):
    feature_order = np.argsort(np.mean(np.abs(shap_values[j]), axis=0))
    return pd.DataFrame(
        {"class": log_reg.classes_[j],
         "model": "L1 Logistic Regression",
            "features": [X.columns[i] for i in feature_order[::-1][:3]],
            "importance": [
                np.mean(np.abs(shap_values[j]), axis=0)[i] for i in feature_order[::-1][:3]
            ],
        }
    )
for i in range(len(log_reg.classes_)):
    output = shapley_feature_ranking(i, readfile_by_kmer_df)
    feature_imp = feature_imp.append(output)
    
########################## for l2 ##########################
feature_imp2 = feature_imp
method = 'l2'

model2 = joblib.load(filepath+'models/k'+k+'_'+kmer+'_'+predicting+'_'+method+'.joblib') # load model
log_reg2 = LogisticRegression(**model2.estimator.get_params()) # feed optimized parameters into logreg
log_reg2.fit(X_train, y_train)
explainer = shap.LinearExplainer(log_reg2, readfile_by_kmer_df)
shap_values = explainer.shap_values(readfile_by_kmer_df)

def shapley_feature_ranking(j, X):
    feature_order = np.argsort(np.mean(np.abs(shap_values[j]), axis=0))
    return pd.DataFrame(
        {"class": log_reg2.classes_[j],
         "model": "L2 Logistic Regression",
            "features": [X.columns[i] for i in feature_order[::-1][:3]],
            "importance": [
                np.mean(np.abs(shap_values[j]), axis=0)[i] for i in feature_order[::-1][:3]
            ],
        }
    )
for i in range(len(log_reg2.classes_)):
    output = shapley_feature_ranking(i, readfile_by_kmer_df)
    feature_imp2 = feature_imp2.append(output) 

########################## for rf ##########################
feature_imp3 = feature_imp2

method = 'rf'
best_rf = joblib.load(filepath+'models/k'+k+'_'+kmer+'_'+predicting+'_'+method+'.joblib') # load model
best_rf.fit(X_train, y_train)
explainer = shap.Explainer(best_rf, readfile_by_kmer_df)
shap_values = explainer.shap_values(readfile_by_kmer_df)

def shapley_feature_ranking(j, X):
    feature_order = np.argsort(np.mean(np.abs(shap_values[j]), axis=0))
    return pd.DataFrame(
        {"class": best_rf.classes_[j],
         "model": "Random Forest",
            "features": [X.columns[i] for i in feature_order[::-1][:3]],
            "importance": [
                np.mean(np.abs(shap_values[j]), axis=0)[i] for i in feature_order[::-1][:3]
            ],
        }
    )
for i in range(len(best_rf.classes_)):
    output = shapley_feature_ranking(i, readfile_by_kmer_df)
    feature_imp3 = feature_imp3.append(output) 


########################## export ##########################
feature_imp3.rename(columns={"class": "ribotype"})
feature_imp3.to_csv('/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/collated_important_kmers_k'+str(kmer_size)+'.csv', index = True, header=True)

########################## add true library ribotypes ##########################
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
kmer_size = int(kmer) 

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
    if ribotype_list.count(',') == 0:
        library_id_2_ribotype[id_] = ribotype_list.replace("{'", '').replace("'}", '') 
        kmer_list = list(library_kmer_2_id.keys())
        id_list = list(library_kmer_2_id.values())
        kmer_position = id_list.index(id_)
        kmer_ = kmer_list[kmer_position]
        library_kmer_2_ribotype[kmer_] = ribotype_list.replace("{'", '').replace("'}", '') 

#### add to feature extraction list
library_ribotype = list()
for k in feature_imp3['features']:
    library_ribotype.append(library_kmer_2_ribotype[k])
Counter(library_ribotype) 
feature_imp3['feature_ribotype'] = library_ribotype

########################## export in prep for network analysis ##########################
feature_imp3.to_csv('/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/with_truth_collated_important_kmers_k'+str(kmer_size)+'.csv', index = True, header=True)
