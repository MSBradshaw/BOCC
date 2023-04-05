"""
The purpose of this script is to create an xgboost regressor and train it to model the emperical p-value of a BOCC cluster
Along the way it should also optimize hyperparameters and tune it's own features

cluster_id	cluster_size	gene_ratio	HPO_ratio	num_sig_go_enrichment_terms	sig_go_enrichment_p_vals	sig_go_enrichment_fdr_corrected_p_vals	sig_go_enrichment_terms	go_sig_threshold	max_norm_cell_type_specificity	max_norm_cell_type_comma_sep_string	num_of_diseases	max_norm_disease_specificity	max_norm_disease_comma_sep_string	mg2_pairs_count	mg2_not_pairs_count	mg2_portion_families_recovered	avg_embeddedness	avg_internal_degree	conductance	cut_ratio	normalized_cut	expansion	triangle_participation_ratio	surprise	significance	newman_girvan_modularity	internal_edge_density	edges_inside	hub_dominance	max_plof	mean_plof	median_plof	std_plof	sum_plof	snowballing_pvalue	num_new_edges_on_any_node
"""
import os
import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import cross_val_score
import sklearn
import pickle

NUM_JOBS = 8
HYPERPARAMS = {
        'learning_rate': list(x/1000 for x in range(1,1000,5)), # 0.0001 to 0.999 by increments of 0.005 
        'gamma': list(x/100 for x in range(1,100,5)) + list(range(1,11,1)), # 0.01 to 1 by .05 and 1 to 10 by 1
        'n_estimators': list(range(1,1001)), # 1-1000
        'max_depth': list(range(1,15)), # depth of trees 1-15
        'max_leaves': list(range(1,5)),
        'subsample': [ x/100 for x in range(1,101,2)], # .01 - 1 by .02
        'booster' : ['dart'],
        'n_jobs': [8]
        }

def load_data(file_path):
    data = pd.read_csv(file_path,sep='\t')
    # remove the column cluster_id
    to_drop = ['cluster_id','sig_go_enrichment_terms','go_sig_threshold','max_norm_cell_type_comma_sep_string','num_new_edges_on_any_node','sig_go_enrichment_p_vals','mg2_portion_families_recovered','mg2_not_pairs_count','mg2_pairs_count','max_norm_disease_comma_sep_string','sig_go_enrichment_fdr_corrected_p_vals']
    for name in to_drop:
        data = data.drop(name,axis=1)    
    # pop the snowballing_pvalue column
    labels = list(data.pop('snowballing_pvalue'))
    # for col in data.columns:
    #     print(col)
    #     print(data[col])
    #     print('\n\n')
    return data, labels

def load_files(files):
    X = None
    y = None
    for f in files:
        tmp_X, tmp_y = load_data(f)
        if X is None:
            X = tmp_X
            y = tmp_y
        else:
            X = pd.concat([X,tmp_X])
            y = y + tmp_y
    return X, y

def train_model(X, y):
    # create an xgboost regressor
    model = xgb.XGBRegressor()
    # lets do grid search hyperparameter optimization
    gs = GridSearchCV(model, HYPERPARAMS, cv=10, n_jobs=NUM_JOBS, verbose=1, scoring='neg_mean_squared_error')
    # fit the model
    gs.fit(X,y)
    # print the best parameters
    print(gs.best_params_)
    # print the best score
    print(gs.best_score_)
    # export the drig search results to a file
    pd.DataFrame(gs.cv_results_).to_csv('grid_search_results.tsv',sep='\t')
    # export the best model as a pickle
    with open('best_model.pkl','wb') as f:
        pickle.dump(gs.best_estimator_,f)
    




def main():
    #'FinalBOCCFeatures/2019/paris.infomap.2019.bocc_res.tsv'
    files_2019 = ['FinalBOCCFeatures/2019/' + f for f in os.listdir('FinalBOCCFeatures/2019/')]
    X19, y19 = load_files(files_2019)
    assert(X19.shape[0] == len(y19))
    train_model(X19,y19)
    print('Complete')


main()

# python Scripts/train_model.py