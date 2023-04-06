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
import matplotlib.pyplot as plt
import seaborn as sns


NUM_JOBS = 32
HYPERPARAMS = {
        'learning_rate': list(x/1000 for x in range(1,100,5)), # 0.0001 to 0.999 by increments of 0.005 
        'gamma': list(x/100 for x in range(1,100,5)), # 0.01 to 1 by .05 and 1 to 10 by 1
        'n_estimators': list(range(1,500,10)), # 1-1000
        'max_depth': list(range(1,15)), # depth of trees 1-15
        'max_leaves': list(range(1,5)),
        'subsample': [ x/100 for x in range(1,101,5)], # .01 - 1 by .05
        'booster' : ['dart']
        }

# HYPERPARAMS = {
#         'learning_rate': [.04,.05], # 0.0001 to 0.999 by increments of 0.005 
#         'gamma': [0.001,0.002], # 0.01 to 1 by .05 and 1 to 10 by 1
#         'n_estimators': [100,150], # 1-1000
#         'max_depth': [1,5,10], # depth of trees 1-15
#         'max_leaves': [1,2,3],
#         'subsample': [.01,.05,.2], # .01 - 1 by .02
#         'booster' : ['dart']
#         }
print(HYPERPARAMS)

def load_data(file_path):
    data = pd.read_csv(file_path,sep='\t')
    # remove the column cluster_id
    to_drop = ['cluster_id','sig_go_enrichment_terms','go_sig_threshold','max_norm_cell_type_comma_sep_string','num_new_edges_on_any_node','sig_go_enrichment_p_vals','mg2_portion_families_recovered','mg2_not_pairs_count','mg2_pairs_count','max_norm_disease_comma_sep_string','sig_go_enrichment_fdr_corrected_p_vals']
    for name in to_drop:
        data = data.drop(name,axis=1)    
    # pop the snowballing_pvalue column
    labels = list(data.pop('snowballing_pvalue'))
    with open('features.tsv','w') as f:
        f.write('Feature\tMin-Max\tMean\tMedian\n')
        for col in data.columns:
            minv = np.min(data[col])
            maxv = np.max(data[col])
            medianv = np.median(data[col])
            meanv = np.mean(data[col])
            # if min is not Nan round it
            if not np.isnan(minv):
                minv = round(minv,3)
            # if max is not Nan round it
            if not np.isnan(maxv):
                maxv = round(maxv,3)
            # if median is not Nan round it
            if not np.isnan(medianv):
                medianv = round(medianv,3)
            # if mean is not Nan round it
            if not np.isnan(meanv):
                meanv = round(meanv,3)
            f.write('{}\t{}-{}\t{}\t{}\n'.format(col,minv,maxv,meanv,medianv))
            # print(data[col])
            # print('\n\n')
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
    model = xgb.XGBRegressor(n_jobs=-1)
    # lets do grid search hyperparameter optimization
    print('Starting grid search')
    gs = GridSearchCV(model, HYPERPARAMS, cv=10, verbose=1, scoring='neg_mean_squared_error')
    print('Fitting the model')
    # fit the model
    gs.fit(X,y)
    # print the best parameters
    print(gs.best_params_)
    # print the best score
    print(gs.best_score_)
    # export the drig search results to a file
    pd.DataFrame(gs.cv_results_).to_csv('grid_search_results.smaller.tsv',sep='\t')
    # export the best model as a pickle
    with open('best_model.smaller.pkl','wb') as f:
        pickle.dump(gs.best_estimator_,f)
    

def plot_feature_correlation(X,outfile):
    # plot the correlation matrix
    corr = X.corr()
    fig, ax = plt.subplots(figsize=(10,10))
    # sns.heatmap(corr, 
    #         xticklabels=corr.columns,
    #         yticklabels=corr.columns,ax=ax,cmap='vlag')
    # plot the correlation matrix as a clustermap
    sns.clustermap(corr,
            xticklabels=corr.columns,
            yticklabels=corr.columns,
            cmap='vlag')
    plt.tight_layout()
    plt.savefig(outfile)
    # print all pairs in the matrix with correlation > 0.9 or < -0.9
    for i in range(corr.shape[0]):
        for j in range(i+1,corr.shape[1]):
            if corr.iloc[i,j] > 0.9 or corr.iloc[i,j] < -0.9:
                print(corr.index[i],corr.columns[j],corr.iloc[i,j])
    print()

def drop_correlated_features(X):
    to_drop = ['HPO_ratio','surprise','conductance','significance','edges_inside','expansion']
    for name in to_drop:
        X = X.drop(name,axis=1)
    return X

def drop_all_but(X):
    keepers = ['hub_dominance','avg_embeddedness','cluster_size','gene_ratio','triangle_participation_ratio','avg_internal_degree','max_norm_disease_specificity','internal_edge_density','cut_ratio','sum_plof','newman_girvan_modularity','num_of_diseases']
    to_drop = [c for c in X.columns if c not in keepers]
    for name in to_drop:
        X = X.drop(name,axis=1)
    return X

def main():
    #'FinalBOCCFeatures/2019/paris.infomap.2019.bocc_res.tsv'
    files_2019 = ['FinalBOCCFeatures/2019/' + f for f in os.listdir('FinalBOCCFeatures/2019/')]
    X19, y19 = load_files(files_2019)
    plot_feature_correlation(X19,'Figures/feature_correlation.png')
    X19 = drop_correlated_features(X19)
    plot_feature_correlation(X19,'Figures/feature_correlation.post_correlation_dropping.png')
    X19 = drop_all_but(X19)
    plot_feature_correlation(X19,'Figures/feature_correlation.final_selected_features.png')
    assert(X19.shape[0] == len(y19))
    print('Training model')
    train_model(X19,y19)
    print('Complete')


if __name__ == '__main__':
    main()


# python Scripts/train_model.py