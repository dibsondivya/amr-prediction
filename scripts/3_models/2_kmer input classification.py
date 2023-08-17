import joblib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import RandomizedSearchCV, GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay, classification_report
from pprint import pprint
import os
from sklearn.svm import SVC
from sklearn.multiclass import OneVsOneClassifier
from sklearn.linear_model import LogisticRegression
import csv
from multiprocessing import Pool as ThreadPool  
import collections

k_list = ['17', '21', '31']

## initialize ## to change
filepath = os.getcwd() + '/'
kmer = 'all'
predicting = 'meancount'

### define
def results_gen(k):

    X_train = pd.read_csv(filepath+'10ribo_X_train_readfile_by_'+kmer+'_'+predicting+'_df_k'+k+'.csv', index_col=0)
    X_test= pd.read_csv(filepath+'10ribo_X_test_readfile_by_'+kmer+'_'+predicting+'_df_k'+k+'.csv', index_col=0)
    y_train = pd.read_csv(filepath+'10ribo_y_train_readfile_by_all_kmer_trueribotype_df_k'+k+'.csv', index_col=0)
    y_test = pd.read_csv(filepath+'10ribo_y_test_readfile_by_all_kmer_trueribotype_df_k'+k+'.csv', index_col=0)


    ##dont scale!

    ############### random forest #################
    # Tuning
    method = 'rf'
    # Number of trees in random forest
    n_estimators = [int(x) for x in np.linspace(start = 200, stop = 2000, num = 10)]
    # Number of features to consider at every split
    max_features = ['auto', 'sqrt']
    # Maximum number of levels in tree
    max_depth = [int(x) for x in np.linspace(10, 110, num = 11)]
    max_depth.append(None)
    # Minimum number of samples required to split a node
    min_samples_split = [2, 5, 10]
    # Minimum number of samples required at each leaf node
    min_samples_leaf = [1, 2, 4]
    # Method of selecting samples for training each tree
    bootstrap = [True, False]
    # Create the random grid
    random_grid = {'n_estimators': n_estimators,
                'max_features': max_features,
                'max_depth': max_depth,
                'min_samples_split': min_samples_split,
                'min_samples_leaf': min_samples_leaf,
                'bootstrap': bootstrap}

    # Use the random grid to search for best hyperparameters
    # First create the base model to tune
    rf = RandomForestClassifier()
    # Random search of parameters, using 3 fold cross validation, 
    # search across 100 different combinations, and use all available cores
    rf_random = RandomizedSearchCV(estimator = rf, param_distributions = random_grid, n_iter = 100, cv = 3, verbose=2, random_state=42, n_jobs = -1)
    rf_random.fit(X_train, y_train)

    ## print best params
    print('RF BEST PARAMS FOR K = '+k)
    print(rf_random.best_params_)

    ## best param model
    best_rf = RandomForestClassifier(**rf_random.best_params_)
    best_rf.fit(X_train, y_train)
    # predict
    y_tuned_pred = rf_random.predict(X_test)
    # classification report
    report = classification_report(y_test, y_tuned_pred, labels=[1, 2, 5, 14, 15, 17, 20, 27, 78, 106], output_dict=True)
    report.update({"accuracy": {"precision": None, "recall": None, "f1-score": report["accuracy"], "support": report['macro avg']['support']}})
    df = pd.DataFrame(report).transpose()
    df.to_csv(filepath+'models/classification_report/k'+k+'_'+kmer+'_'+predicting+'_'+method+'.csv', header=True)
    cm = confusion_matrix(y_test, y_tuned_pred, labels=[1, 2, 5, 14, 15, 17, 20, 27, 78, 106])
    cm = ConfusionMatrixDisplay(cm, display_labels = [1, 2, 5, 14, 15, 17, 20, 27, 78, 106])
    cm.plot()
    plt.savefig(filepath+'models/confusion_matrix/k'+k+'_'+kmer+'_'+predicting+'_'+method+'.png')

    joblib.dump(best_rf, filepath+'models/k'+k+'_'+kmer+'_'+predicting+'_'+method+'.joblib')

    ############## l1 ############
    ## ovo and l1
    # define model
    method='l1'
    model = LogisticRegression(penalty='l1', solver='liblinear') # default is l2
    # define ovo strategy
    ovo = OneVsOneClassifier(model)
    # hyperparam tuning

    LRparam_grid = {
        'estimator__C': [0.001, 0.01, 0.1, 1, 10, 100, 1000]
        }
    ovo_tuned = GridSearchCV(estimator = ovo, param_grid = LRparam_grid, refit=True, cv=10, verbose=1, scoring='balanced_accuracy', n_jobs = -1)

    # fit model
    ovo_tuned.fit(X_train, y_train)
    print('L1 BEST PARAMS FOR K = '+k)
    print(ovo_tuned.best_params_)
    # export results
    df = pd.DataFrame(ovo_tuned.cv_results_)
    df.to_csv(filepath+'models/hyperparams/k'+k+'_'+kmer+'_'+predicting+'_'+method+'.csv', header=True)
    # make predictions
    yhat = ovo_tuned.predict(X_test)

    report = classification_report(y_test, yhat, labels=[1, 2, 5, 14, 15, 17, 20, 27, 78, 106], output_dict=True)
    report.update({"accuracy": {"precision": None, "recall": None, "f1-score": report["accuracy"], "support": report['macro avg']['support']}})
    df = pd.DataFrame(report).transpose()
    df.to_csv(filepath+'models/classification_report/k'+k+'_'+kmer+'_'+predicting+'_'+method+'.csv', header=True)
    cm = confusion_matrix(y_test, yhat, labels=[1, 2, 5, 14, 15, 17, 20, 27, 78, 106])
    cm = ConfusionMatrixDisplay(cm, display_labels = [1, 2, 5, 14, 15, 17, 20, 27, 78, 106])
    cm.plot()
    plt.savefig(filepath+'models/confusion_matrix/k'+k+'_'+kmer+'_'+predicting+'_'+method+'.png')

    joblib.dump(ovo_tuned.best_estimator_, filepath+'models/k'+k+'_'+kmer+'_'+predicting+'_'+method+'.joblib')

    ############## l2 ############
    ## ovo and l2
    method='l2'
    # define model
    model = LogisticRegression() # default is l2
    # define ovo strategy
    ovo = OneVsOneClassifier(model)
    # hyperparam tuning
    LRparam_grid = {
        'estimator__C': [0.001, 0.01, 0.1, 1, 10, 100, 1000]
        }
    ovo_tuned = GridSearchCV(estimator = ovo, param_grid = LRparam_grid, refit=True, cv=10, verbose=1, scoring='balanced_accuracy', n_jobs = -1)

    # fit model
    ovo_tuned.fit(X_train, y_train)
    print('L2 BEST PARAMS FOR K = '+k)
    print(ovo_tuned.best_params_)
    # export results
    df = pd.DataFrame(ovo_tuned.cv_results_)
    df.to_csv(filepath+'models/hyperparams/k'+k+'_'+kmer+'_'+predicting+'_'+method+'.csv', header=True)
    # make predictions
    yhat = ovo_tuned.predict(X_test)

    report = classification_report(y_test, yhat, labels=[1, 2, 5, 14, 15, 17, 20, 27, 78, 106], output_dict=True)
    report.update({"accuracy": {"precision": None, "recall": None, "f1-score": report["accuracy"], "support": report['macro avg']['support']}})
    df = pd.DataFrame(report).transpose()
    df.to_csv(filepath+'models/classification_report/k'+k+'_'+kmer+'_'+predicting+'_'+method+'.csv', header=True)
    cm = confusion_matrix(y_test, yhat, labels=[1, 2, 5, 14, 15, 17, 20, 27, 78, 106])
    cm = ConfusionMatrixDisplay(cm, display_labels = [1, 2, 5, 14, 15, 17, 20, 27, 78, 106])
    cm.plot()
    plt.savefig(filepath+'models/confusion_matrix/k'+k+'_'+kmer+'_'+predicting+'_'+method+'.png')

    joblib.dump(ovo_tuned.best_estimator_, filepath+'models/k'+k+'_'+kmer+'_'+predicting+'_'+method+'.joblib')
    ############# svm ############
    # define model
    method = 'svm'
    SVCmodel = SVC()
    # define ovo strategy
    SVCovo = OneVsOneClassifier(SVCmodel)
    # hyperparam tuning
    SVCparam_grid = {
        'estimator__kernel':('sigmoid', 'rbf', 'poly'), 
        'estimator__gamma': [1, 0.1, 0.01, 0.001],
        'estimator__C':[0.1, 1, 10, 100]
    }
    SVCovo_tuned = RandomizedSearchCV(estimator = SVCovo, param_distributions = SVCparam_grid, refit=True, n_iter = 25, cv=3, verbose=1, random_state=42, n_jobs = -1)
    # fit model
    SVCovo_tuned.fit(X_train, y_train)
    # export hyperparams
    df = pd.DataFrame(SVCovo_tuned.cv_results_)
    df.to_csv(filepath+'models/hyperparams/k'+k+'_'+kmer+'_'+predicting+'_'+method+'.csv', header=True)
    # print best parameter after tuning
    print('SVM Best Params for k = '+k)
    print(SVCovo_tuned.best_params_)
    # print how our model looks after hyper-parameter tuning
    print('SVM Best Estimator for K = '+k)
    print(SVCovo_tuned.best_estimator_)
    # make predictions
    SVCyhat = SVCovo_tuned.predict(X_test)
  
    report = classification_report(y_test, SVCyhat, labels=[1, 2, 5, 14, 15, 17, 20, 27, 78, 106], output_dict=True)
    report.update({"accuracy": {"precision": None, "recall": None, "f1-score": report["accuracy"], "support": report['macro avg']['support']}})
    df = pd.DataFrame(report).transpose()
    df.to_csv(filepath+'models/classification_report/k'+k+'_'+kmer+'_'+predicting+'_'+method+'.csv', header=True)
    cm = confusion_matrix(y_test, SVCyhat, labels=[1, 2, 5, 14, 15, 17, 20, 27, 78, 106])
    cm = ConfusionMatrixDisplay(cm, display_labels = [1, 2, 5, 14, 15, 17, 20, 27, 78, 106])
    cm.plot()
    plt.savefig(filepath+'models/confusion_matrix/k'+k+'_'+kmer+'_'+predicting+'_'+method+'.png')
    joblib.dump(SVCovo_tuned.best_estimator_, filepath+'models/k'+k+'_'+kmer+'_'+predicting+'_'+method+'.joblib')


## process samples    
n_threads = 3    
pool = ThreadPool(n_threads) 
tmp_res = pool.map_async(results_gen, k_list)
output = tmp_res.get()
pool.close() 
pool.join() 

