import shap
import venn
import pandas as pd
import numpy as np
import joblib
import matplotlib.pyplot as plt
import os
from sklearn.linear_model import LogisticRegression

np.bool = np.bool_

k='31' # repeat w 21 and 17
kmer_size = int(k)
kmer = 'all' # repeat with 'kmer' aka unique kmer library
method = 'l1'

filepath = '/rds/general/user/dds122/projects/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/'
readfile_by_kmer_df = pd.read_csv('/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/10ribo_X_test_readfile_by_all_absence_df_k'+str(kmer_size)+'.csv', index_col = 0)
feature_imp = pd.DataFrame(columns=['class', 'data', 'features', 'importance'])
X_train = pd.read_csv('/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/10ribo_X_train_readfile_by_all_absence_df_k'+str(kmer_size)+'.csv', index_col = 0)
y_train = pd.read_csv('/rds/general/project/hda-22-23/live/Summer_projects/dds122/data/05-25-2023/10ribo_y_train_readfile_by_all_kmer_trueribotype_df_k'+str(kmer_size)+'.csv', index_col = 0)

# store shap values for top three kmers per ribotype
predicting_list = ['binary', 'sumcount', 'mediancount', 'meancount']
for predicting in predicting_list:
    model = joblib.load(filepath+'models/k'+k+'_'+kmer+'_'+predicting+'_'+method+'.joblib') # load ovo model
    labels=[1, 2, 5, 14, 15, 17, 20, 27, 78, 106]
    scratchmodel = LogisticRegression(**model.estimator.get_params()) # feed optimized parameters into logreg
    scratchmodel.fit(X_train, y_train)
    explainer = shap.LinearExplainer(scratchmodel, readfile_by_kmer_df)
    shap_values = explainer.shap_values(readfile_by_kmer_df)

    def shapley_feature_ranking(j, X):
        feature_order = np.argsort(np.mean(np.abs(shap_values[j]), axis=0))
        return pd.DataFrame(
            {"class": labels[j],
             "data": predicting,
                "features": [X.columns[i] for i in feature_order[::-1][:3]], ## restrict to top three
                "importance": [
                    np.mean(np.abs(shap_values[j]), axis=0)[i] for i in feature_order[::-1][:3] ## restrict to top three
                ],
            }
        )

    for lol in range(len(labels)):
        output = shapley_feature_ranking(lol, readfile_by_kmer_df)
        feature_imp = feature_imp.append(output)

# print and save the venn diagrams
labels=[1, 2, 5, 14, 15, 17, 20, 27, 78, 106]
for ribotype in labels:
    num_ribotype = ribotype
    subset_data = feature_imp[feature_imp['class']==ribotype]
    binary_list = subset[subset['data']=='binary']['features'].tolist()
    sum_list = subset[subset['data']=='sumcount']['features'].tolist()
    median_list = subset[subset['data']=='mediancount']['features'].tolist()
    mean_list = subset[subset['data']=='meancount']['features'].tolist()

    labels = venn.get_labels([set(binary_list), set(sum_list), set(median_list), set(mean_list)], fill=['number'])
    fig, ax = venn.venn4(labels, names=['Binary', 'Sum', 'Median', 'Mean'])
    plt.title('Comparing Selected K-mers from L1 Logistic Regression for K='+k+" for Ribotype "+str(ribotype))
    plt.savefig(filepath+'venn/datacomp/k'+k+'_'+kmer+'_'+method+'_ribo'+str(ribotype)+'.png')
    plt.show()
