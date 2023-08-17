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
predicting = 'kmer' # or 'all'

readfile_by_kmer_df = 
########################## for l1 ##########################
method = 'l1'

########################## for l2 ##########################
feature_imp2 = feature_imp1
method = 'l2'
logreg2 = joblib.load(filepath+'models/k'+k+'_'+kmer+'_'+predicting+'_'+method+'.joblib') # load model




########################## for rf ##########################
feature_imp3 = feature_imp2

method = 'rf'
best_rf = joblib.load(filepath+'models/k'+k+'_'+kmer+'_'+predicting+'_'+method+'.joblib') # load model
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
