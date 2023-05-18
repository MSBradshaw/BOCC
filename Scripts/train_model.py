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
from sklearn.metrics import mean_squared_error, roc_curve
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import cross_val_score
from sklearn.decomposition import PCA
import sklearn
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import random
from sklearn.preprocessing import RobustScaler
from sklearn.feature_selection import SequentialFeatureSelector
from sklearn_genetic import GASearchCV
from sklearn_genetic.space import Categorical, Integer, Continuous
from sklearn.metrics import make_scorer
from scipy.stats import pearsonr
import shap

# import friedmanchisquared
from scipy import stats

# imoprt precision_score, recall_score, f1_score, roc_auc_score, accuracy_score
from sklearn.metrics import precision_score, recall_score, f1_score, roc_auc_score, accuracy_score
# import confusion matrix
from sklearn.metrics import confusion_matrix


NUM_JOBS = 32
HYPERPARAMS = {
        'learning_rate': list(x/1000 for x in range(1,100,5)), # 0.0001 to 0.999 by increments of 0.005 
        'gamma': list(x/100 for x in range(1,100,5)), # 0.01 to 1 by .05 and 1 to 10 by 1
        'n_estimators': list(range(1,500,10)), # 1-1000
        'max_depth': list(range(1,15)), # depth of trees 1-15
        'max_leaves': list(range(1,10)),
        'subsample': [ x/100 for x in range(1,101,5)], # .01 - 1 by .05
        'booster' : ['dart']
        }

# print the number of combinations of hyperparameters to search
print('Number of combinations of hyperparameters to search: {}'.format(np.prod([len(x) for x in HYPERPARAMS.values()])))

HYPERPARAMS_Genetic_Algo = {
        'learning_rate': Continuous(1e-3, 1e-1, distribution='uniform'),
        'gamma': Continuous(1e-2, 1, distribution='uniform'),
        'n_estimators': Integer(1, 500),
        'max_depth': Integer(1, 30),
        'max_leaves': Integer(1, 10),
        'subsample': Continuous(1e-2, 1, distribution='uniform'),
        'booster' : Categorical(['dart'])
        }

HYPERPARAMS_p35 = {'learning_rate': 0.04043818815600482, 'gamma': 0.5713083222302612, 'n_estimators': 196, 'max_depth': 25, 'max_leaves': 2, 'subsample': 0.20159801360534982, 'booster': 'dart'}
HYPERPARAMS_p1 = {'learning_rate': 0.05272395382072624, 'gamma': 0.36083958112876147, 'n_estimators': 35, 'max_depth': 6, 'max_leaves': 10, 'subsample': 0.7744273119528029, 'booster': 'dart'}
HYPERPARAMS_p05 = {'learning_rate': 0.055507261946723785, 'gamma': 0.5737961358822451, 'n_estimators': 214, 'max_depth': 6, 'max_leaves': 7, 'subsample': 0.7520797264062326, 'booster': 'dart'}
HYPERPARAMS_p35_regression = {'learning_rate': 0.06608502330035836, 'gamma': 0.012905252618403522, 'n_estimators': 416, 'max_depth': 8, 'max_leaves': 5, 'subsample': 0.9514669623440667, 'booster': 'dart'}
# classification params
classification_params = {'learning_rate': 0.07305940290241641, 
                         'gamma': 0.3795158439477443, 
                         'n_estimators': 98, 
                         'max_depth': 13, 
                         'max_leaves': 4, 
                         'subsample': 0.05772500174333884, 
                         'booster': 'dart'}

# classification params
regression_params = {'learning_rate': 0.0768765080451186, 
                     'gamma': 0.9691138788518666, 
                     'n_estimators': 356, 
                     'max_depth': 11, 
                     'max_leaves': 4, 
                     'subsample': 0.7529052269196642, 
                     'booster': 'dart'}

just_regression_params = {'learning_rate': 0.011440851730240585, 'gamma': 0.020537510220500793, 'n_estimators': 331, 'max_depth': 7, 'max_leaves': 10, 'subsample': 0.11060470434528673, 'booster': 'dart'}

def kendals_tau(ranks1,ranks2):
    assert len(ranks1) == len(ranks2)
    # source https://www.youtube.com/watch?v=zVufp7cJ8S4
    tau, p_value = stats.kendalltau(ranks1,ranks2)
    return tau, p_value

def load_data(file_path):
    print(file_path)
    data = pd.read_csv(file_path,sep='\t')
    data['algo'] = '.'.join(file_path.split('/')[-1].split('.')[0:2])
    # remove the column cluster_id
    to_drop = ['significance','sig_go_enrichment_terms','go_sig_threshold','max_norm_cell_type_comma_sep_string','num_new_edges_on_any_node','sig_go_enrichment_p_vals','mg2_portion_families_recovered','mg2_not_pairs_count','mg2_pairs_count','max_norm_disease_comma_sep_string','sig_go_enrichment_fdr_corrected_p_vals']
    for name in to_drop:
        if name in data.columns:
            data = data.drop(name,axis=1)
        else:
            print('WARNING: {} not in data'.format(name))
    # pop the snowballing_pvalue column
    labels = list(data.pop('snowballing_pvalue'))
    with open('features.tsv','w') as f:
        f.write('Feature\tMin-Max\tMean\tMedian\n')
        for col in data.columns:
            if col == 'algo':
                continue
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

def find_hyperparams(X, y):
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

def grid_search(X, y, downsample=True):
    # downsample the majority class
    if downsample:
        # downsample the majority class
        # get the indices of the majority class
        majority_indices = [i for i,x in enumerate(y) if x == 1]
        # get the indices of the minority class
        minority_indices = [i for i,x in enumerate(y) if x != 1]
        # get the number of majority class samples
        majority_count = len(majority_indices)
        # get the number of minority class samples
        minority_count = len(minority_indices)
        # get the number of samples to downsample the majority class to
        downsample_count = minority_count
        # get the indices of the majority class to downsample
        print('Downsampling {} samples to {}'.format(majority_count, downsample_count))
        downsample_indices = random.sample(majority_indices, downsample_count)
        # get the indices of the minority class
        upsample_indices = minority_indices
        # get the indices of the samples to keep
        keep_indices = downsample_indices + upsample_indices
        # downsample the majority class
        X = X.iloc[keep_indices]
        y = [y[i] for i in keep_indices]
    # test train split of data
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    # create an xgboost regressor
    model = xgb.XGBRegressor(n_jobs=-1)
    # lets do grid search hyperparameter optimization
    # print('Starting grid search')
    gs = GridSearchCV(model, HYPERPARAMS, cv=10, verbose=1, scoring='neg_mean_squared_error')
    # print('Fitting the model')
    # fit the model
    gs.fit(X_train,y_train)
    # print the best parameters
    print('Best parameters')
    print(gs.best_params_)
    # print the best score
    print('Best score')
    print(gs.best_score_)
    #pickle the model
    with open('best_model.pkl','wb') as f:
        pickle.dump(gs.best_estimator_,f)

def grid_search_classification(X, y, downsample=True):
    # downsample the majority class
    if downsample:
        majority_indices = [i for i,x in enumerate(y) if x != 1]
        minority_indices = [i for i,x in enumerate(y) if x == 1]
        majority_count = len(majority_indices)
        minority_count = len(minority_indices)
        downsample_count = minority_count
        print('Downsampling {} samples to {}'.format(majority_count, downsample_count))
        downsample_indices = random.sample(majority_indices, downsample_count)
        upsample_indices = minority_indices
        keep_indices = downsample_indices + upsample_indices
        X = X.iloc[keep_indices]
        y = [y[i] for i in keep_indices]
    # test train split of data
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    # create an xgboost regressor
    model = xgb.XGBClassifier(n_jobs=-1)
    # lets do grid search hyperparameter optimization
    # print('Starting grid search')
    gs = GridSearchCV(model, HYPERPARAMS, cv=10, verbose=1, scoring='precision')
    # print('Fitting the model')
    # fit the model
    gs.fit(X_train,y_train)
    # print the best parameters
    print('Best parameters')
    print(gs.best_params_)
    # print the best score
    print('Best score')
    print(gs.best_score_)
    #pickle the model
    with open('best_model.classification.pkl','wb') as f:
        pickle.dump(gs.best_estimator_,f)

def genetic_optimization_classification(X, y, downsample=True):
    # downsample the majority class
    if downsample:
        majority_indices = [i for i,x in enumerate(y) if x != 1]
        minority_indices = [i for i,x in enumerate(y) if x == 1]
        majority_count = len(majority_indices)
        minority_count = len(minority_indices)
        downsample_count = minority_count
        print('Downsampling {} samples to {}'.format(majority_count, downsample_count))
        downsample_indices = random.sample(majority_indices, downsample_count)
        upsample_indices = minority_indices
        keep_indices = downsample_indices + upsample_indices
        X = X.iloc[keep_indices]
        y = [y[i] for i in keep_indices]
    # test train split of data
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    # create an xgboost regressor
    model = xgb.XGBClassifier()
    # lets do grid search hyperparameter optimization
    print('Starting genetic algorithm search')
    def precision_plus_recall(ground_truth, predictions):
        precision = precision_score(ground_truth, predictions)
        recall = recall_score(ground_truth, predictions)
        return precision + recall

    # loss_func will negate the return value of my_custom_loss_func,
    #  which will be np.log(2), 0.693, given the values for ground_truth
    #  and predictions defined below.
    pr_sum  = make_scorer(precision_plus_recall, greater_is_better=True)

    go = GASearchCV(estimator=model, cv=10, param_grid=HYPERPARAMS_Genetic_Algo, verbose=True, scoring='f1', n_jobs=10, population_size=40, generations=100)
    # print('Fitting the model')
    # fit the model
    go.fit(X_train,y_train)
    # print the best parameters
    print('Best parameters')
    print(go.best_params_)
    # print the best score
    print('Best score')
    print(go.best_score_)
    #pickle the model
    with open('best_model.classification.f1.go.pkl','wb') as f:
        pickle.dump(go.best_estimator_,f)
    return go.best_params_, go.best_score_

def genetic_optimization_regression(X, y, downsample=False, threshold=1):
    # downsample the majority class
    if downsample:
        majority_indices = [i for i,x in enumerate(y) if x >= threshold]
        minority_indices = [i for i,x in enumerate(y) if x < threshold]
        majority_count = len(majority_indices)
        minority_count = len(minority_indices)
        downsample_count = minority_count
        print('Downsampling {} samples to {}'.format(majority_count, downsample_count))
        downsample_indices = random.sample(majority_indices, downsample_count)
        upsample_indices = minority_indices
        keep_indices = downsample_indices + upsample_indices
        X = X.iloc[keep_indices]
        y = [y[i] for i in keep_indices]
    # test train split of data
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    # create an xgboost regressor
    model = xgb.XGBRegressor()
    # lets do grid search hyperparameter optimization
    print('Starting genetic algorithm search')
    go = GASearchCV(estimator=model, cv=10, param_grid=HYPERPARAMS_Genetic_Algo, verbose=True, scoring='neg_mean_squared_error', n_jobs=10, population_size=40, generations=100)
    # print('Fitting the model')
    # fit the model
    go.fit(X_train,y_train)
    # print the best parameters
    print('Best parameters')
    print(go.best_params_)
    # print the best score
    print('Best score')
    print(go.best_score_)
    #pickle the model
    with open('best_model.regression.go.pkl','wb') as f:
        pickle.dump(go.best_estimator_,f)

    return go.best_params_, go.best_score_


def train_model_classifier(X, y, normaize=True, downsample=True):
    # turn it into a binary classification problem
    y = [ 1 if x < 0.05 else 0 for x in y]
    if downsample:
        # downsample the majority class
        # get the indices of the majority class
        majority_indices = [i for i,x in enumerate(y) if x == 0]
        # get the indices of the minority class
        minority_indices = [i for i,x in enumerate(y) if x != 0]
        # get the number of majority class samples
        majority_count = len(majority_indices)
        # get the number of minority class samples
        minority_count = len(minority_indices)
        # get the number of samples to downsample the majority class to
        downsample_count = minority_count
        # get the indices of the majority class to downsample
        downsample_indices = random.sample(majority_indices, downsample_count)
        # get the indices of the minority class
        upsample_indices = minority_indices
        # get the indices of the samples to keep
        keep_indices = downsample_indices + upsample_indices
        # downsample the majority class
        X = X.iloc[keep_indices]
        y = [y[i] for i in keep_indices]
    # print the number of samples in each class
    print('# y == 0: {}'.format(len([x for x in y if x == 0])))
    print('# y == 1: {}'.format(len([x for x in y if x != 0])))

    # make a test train split
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    if normaize:
        train_transformer = RobustScaler().fit(X_train)
        X_train = train_transformer.transform(X_train)

        test_transformer = RobustScaler().fit(X_test)
        X_test = test_transformer.transform(X_test)

    model = xgb.XGBClassifier(n_jobs=-1)
    print('Fitting the model')
    scores = cross_val_score(model, X, y, cv=10, scoring='precision')
    print('Precision', scores)
    scores = cross_val_score(model, X, y, cv=10, scoring='recall')
    print('Recall', scores)

    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    print('Precision', precision_score(y_test, y_pred))
    print('Recall', recall_score(y_test, y_pred))
    print('F1', f1_score(y_test, y_pred))
    print('Accuracy', accuracy_score(y_test, y_pred))
    print('Confusion matrix')  
    tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
    print('tn: {}, fp: {}\n fn: {}, tp: {}'.format(tn, fp, fn, tp))
    print(confusion_matrix(y_test, y_pred))
    return tn, fp, fn, tp, precision_score(y_test, y_pred), recall_score(y_test, y_pred)
      
def train_model_regressor(X, y, normaize=True, downsample=True):
    # turn it into a binary classification problem
    # y = [ 1 if x < 0.05 else 0 for x in y]
    if downsample:
        # downsample the majority class
        # get the indices of the majority class
        majority_indices = [i for i,x in enumerate(y) if x == 1]
        # get the indices of the minority class
        minority_indices = [i for i,x in enumerate(y) if x != 1]
        # get the number of majority class samples
        majority_count = len(majority_indices)
        # get the number of minority class samples
        minority_count = len(minority_indices)
        # get the number of samples to downsample the majority class to
        downsample_count = minority_count
        # get the indices of the majority class to downsample
        downsample_indices = random.sample(majority_indices, downsample_count)
        # get the indices of the minority class
        upsample_indices = minority_indices
        # get the indices of the samples to keep
        keep_indices = downsample_indices + upsample_indices
        # downsample the majority class
        X = X.iloc[keep_indices]
        y = [y[i] for i in keep_indices]
    # print the number of samples in each class
    print('# p == 1: {}'.format(len([x for x in y if x == 1])))
    print('# p != 1: {}'.format(len([x for x in y if x != 1])))

    # make a test train split
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    if normaize:
        train_transformer = RobustScaler().fit(X_train)
        X_train = train_transformer.transform(X_train)

        test_transformer = RobustScaler().fit(X_test)
        X_test = test_transformer.transform(X_test)

    model = xgb.XGBRegressor(n_jobs=-1)
    print('Fitting the model')
    scores = cross_val_score(model, X, y, cv=10, scoring='neg_mean_squared_error')
    print('negative MSE', scores)
    
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    mse = mean_squared_error(y_test, y_pred)
    print('Mean Squred Error', mse)
    return mse, y_pred, y_test

def plot_feature_correlation(X,outfile):
    # plot the correlation matrix
    corr = X.corr()
    fig, ax = plt.subplots(figsize=(10,10))
    sns.clustermap(corr,
            xticklabels=corr.columns,
            yticklabels=corr.columns,
            cmap='vlag',
            annot=True,
            fmt='.2f',
            linewidths=.75,
            center=0)

    plt.tight_layout()
    plt.savefig(outfile)
    plt.clf()
    # print all pairs in the matrix with correlation > 0.9 or < -0.9
    for i in range(corr.shape[0]):
        for j in range(i+1,corr.shape[1]):
            if corr.iloc[i,j] > 0.9 or corr.iloc[i,j] < -0.9:
                print(corr.index[i],corr.columns[j],corr.iloc[i,j])
    print()

def drop_correlated_features(X):
    to_drop = ['HPO_ratio','surprise','conductance','edges_inside','expansion']
    for name in to_drop:
        X = X.drop(name,axis=1)
    return X

def drop_all_but(X):
    # keepers = ['algo', 'hub_dominance','avg_embeddedness','cluster_size','gene_ratio','triangle_participation_ratio','avg_internal_degree','max_norm_disease_specificity','internal_edge_density','cut_ratio','sum_plof','newman_girvan_modularity','num_of_diseases']
    keepers = ['algo', 'hub_dominance','avg_embeddedness','cluster_size','gene_ratio','triangle_participation_ratio','avg_internal_degree','max_norm_disease_specificity','internal_edge_density','cut_ratio','newman_girvan_modularity','num_of_diseases']
    to_drop = [c for c in X.columns if c not in keepers]
    for name in to_drop:
        X = X.drop(name,axis=1)
    return X

def drop_algo(X):
    to_drop = ['algo']
    for name in to_drop:
        X = X.drop(name,axis=1)
    return X

def plot_pca(X,y):
    y = [ 0 if i != 1 else 1 for i in y]
    pca = PCA(n_components=4)
    
    fig, ax = plt.subplots(2,2,figsize=(10,10))
    for i,algo in enumerate(X['algo'].unique()):
        xi = i//2
        yi = i%2
        sub = X[X['algo'] == algo]
        # sub set y too
        suby = [j for j,v in zip(y,list(X['algo'] == algo)) if v]
        subsub = sub.drop('algo',axis=1)
        X_pca = pca.fit_transform(subsub)
        ax[xi,yi].scatter(X_pca[:,0],X_pca[:,1],c=suby)
        ax[xi,yi].set_xlabel('PC1')
        ax[xi,yi].set_ylabel('PC2')
        ax[xi,yi].set_title(algo)
        plt.tight_layout()
        plt.savefig('Figures/pca_1_2.png')
        plt.clf()
    
    # do the same for PC2 and 3
    fig, ax = plt.subplots(2,2,figsize=(10,10))
    for i,algo in enumerate(X['algo'].unique()):
        xi = i//2
        yi = i%2
        sub = X[X['algo'] == algo]
        suby = [j for j,v in zip(y,list(X['algo'] == algo)) if v]
        subsub = sub.drop('algo',axis=1)
        X_pca = pca.fit_transform(subsub)
        ax[xi,yi].scatter(X_pca[:,1],X_pca[:,2],c=suby)
        ax[xi,yi].set_xlabel('PC2')
        ax[xi,yi].set_ylabel('PC3')
        ax[xi,yi].set_title(algo)
        plt.tight_layout()
        plt.savefig('Figures/pca_2_3.png')
        plt.clf()
    
    # and again for PC3 and 4
    fig, ax = plt.subplots(2,2,figsize=(10,10))
    for i,algo in enumerate(X['algo'].unique()):
        xi = i//2
        yi = i%2
        sub = X[X['algo'] == algo]
        suby = [j for j,v in zip(y,list(X['algo'] == algo)) if v]
        subsub = sub.drop('algo',axis=1)
        X_pca = pca.fit_transform(subsub)
        ax[xi,yi].scatter(X_pca[:,2],X_pca[:,3],c=suby)
        ax[xi,yi].set_xlabel('PC3')
        ax[xi,yi].set_ylabel('PC4')
        ax[xi,yi].set_title(algo)
        plt.tight_layout()
        plt.savefig('Figures/pca_3_4.png')
        plt.clf()

def find_nan(X):
    # for each column
    for col in X.columns:
        # find the number of nan values
        if X[col].isna().any():
            print(col)
            print(X[col].isna().sum())
        
def plot_pval_dist(y,year,figname):
    # plot the distribution of p values
    plt.hist(y,bins=20)
    plt.xlabel('p value')
    plt.ylabel('count')
    plt.savefig(figname)
    plt.close()

def select_features(X,y,n,downsample=True,normaize=True):
    cols = X.columns
    y = [ 1 if x < 0.05 else 0 for x in y]
    if downsample:
        majority_indices = [i for i,x in enumerate(y) if x == 0]
        minority_indices = [i for i,x in enumerate(y) if x != 0]
        minority_count = len(minority_indices)
        downsample_count = minority_count
        downsample_indices = random.sample(majority_indices, downsample_count)
        upsample_indices = minority_indices
        keep_indices = downsample_indices + upsample_indices
        X = X.iloc[keep_indices]
        y = [y[i] for i in keep_indices]

    # make a test train split
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    if normaize:
        train_transformer = RobustScaler().fit(X_train)
        X_train = train_transformer.transform(X_train)
        X_train = pd.DataFrame(X_train,columns=X.columns)

        test_transformer = RobustScaler().fit(X_test)
        X_test = test_transformer.transform(X_test)
        X_test = pd.DataFrame(X_test,columns=X.columns)


    model = xgb.XGBClassifier(n_jobs=-1)
    sfs = SequentialFeatureSelector(estimator=model, n_features_to_select=n)
    sfs.fit(X_train, y_train)

    # train the model with the selected features
    print(sfs.get_support())
    print(type(X_train))
    print(X_train.shape)

    # subset X_train and X_test to only the selected features - assuming they are dataframes
    X_train = X_train.iloc[:,sfs.get_support()]
    X_test = X_test.iloc[:,sfs.get_support()]


    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    precision = precision_score(y_test, y_pred)
    recall = recall_score(y_test, y_pred)
    return precision, recall, [c for v,c in zip(sfs.get_support(),cols) if v]

def select_all_features(X,y):
    selection_results = {"i":[],"precision":[],"recall":[],"features":[],"normalized":[],'downsampled':[]}
    for norm in [True,False]:
        for downsample in [True,False]:
            for i in range(1,len(X.columns)-1):
                # print(i)
                print(i, norm, downsample)
                precision, recall, features = select_features(X,y,n=i,normaize=norm,downsample=downsample)
                print(i, precision, recall, features, norm, downsample)
                print()
                selection_results['i'].append(i)
                selection_results['precision'].append(precision)
                selection_results['recall'].append(recall)
                selection_results['features'].append(str(features))
                selection_results['normalized'].append(norm)
                selection_results['downsampled'].append(downsample)
    selection_results = pd.DataFrame(selection_results)
    selection_results.to_csv('ClassificationResults/selection_results.csv',index=False)

# find the best feature set from 'selection_results.csv'
def find_best_features():
    df = pd.read_csv('ClassificationResults/selection_results.csv')
    # sort by recall
    df = df.sort_values(by='recall',ascending=False)
    # assign a recall rank
    df['recall_rank'] = df['recall'].rank(method='first',ascending=False)
    # sort by precision
    df = df.sort_values(by='precision',ascending=False)
    # assign a precision rank
    df['precision_rank'] = df['precision'].rank(method='first',ascending=False)
    # assign a rank
    df['rank'] = df['recall_rank'] + df['precision_rank']
    # sort by rank
    df = df.sort_values(by='rank',ascending=True)
    # print the top 10
    print(df.head(5))
    # print all the information in the top 5
    for i in range(5):
        print(df.iloc[i])
    # save the ranked and sorted results
    df.to_csv('ClassificationResults/selection_results_ranked.csv',index=False)

def select_features_regression(X,y,n,downsample=True,normaize=True):
    cols = X.columns
    if downsample:
        majority_indices = [i for i,x in enumerate(y) if x == 1]
        minority_indices = [i for i,x in enumerate(y) if x != 1]
        minority_count = len(minority_indices)
        downsample_count = minority_count
        downsample_indices = random.sample(majority_indices, downsample_count)
        upsample_indices = minority_indices
        keep_indices = downsample_indices + upsample_indices
        X = X.iloc[keep_indices]
        y = [y[i] for i in keep_indices]

    # make a test train split
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    if normaize:
        train_transformer = RobustScaler().fit(X_train)
        X_train = train_transformer.transform(X_train)
        X_train = pd.DataFrame(X_train,columns=X.columns)

        test_transformer = RobustScaler().fit(X_test)
        X_test = test_transformer.transform(X_test)
        X_test = pd.DataFrame(X_test,columns=X.columns)


    model = xgb.XGBRegressor(n_jobs=-1)
    sfs = SequentialFeatureSelector(estimator=model, n_features_to_select=n, scoring='neg_mean_squared_error')
    sfs.fit(X_train, y_train)

    # train the model with the selected features
    print(sfs.get_support())
    print(type(X_train))
    print(X_train.shape)

    # subset X_train and X_test to only the selected features - assuming they are dataframes
    X_train = X_train.iloc[:,sfs.get_support()]
    X_test = X_test.iloc[:,sfs.get_support()]


    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    mse = mean_squared_error(y_test, y_pred)
    print('Trained MSE', mse)
    return mse, [c for v,c in zip(sfs.get_support(),cols) if v]

def select_all_features_regression(X,y,outdir='RegressionResults',norm_only=False,downsample_only=False):
    selection_results = {"i":[], "MSE":[], "features":[], "normalized":[], 'downsampled':[]}
    for norm in [True,False]:
        if norm_only and not norm:
            continue
        for downsample in [True]:
            if downsample_only and not downsample:
                continue
            for i in range(1,len(X.columns)-1):
                # print(i)
                print(i, norm, downsample)
                mse, features = select_features_regression(X,y,n=i,normaize=norm,downsample=downsample)
                print(i, mse, features, norm, downsample)
                print()
                selection_results['i'].append(i)
                selection_results['MSE'].append(mse)
                selection_results['features'].append(str(features))
                selection_results['normalized'].append(norm)
                selection_results['downsampled'].append(downsample)
    selection_results = pd.DataFrame(selection_results)
    # check if the output directory exists
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    selection_results.to_csv('{}/selection_results.csv'.format(outdir),index=False)

def find_best_features_regression(outdir='RegressionResults'):
    df = pd.read_csv('{}/selection_results.csv'.format(outdir))
    # sort my MSE
    df = df.sort_values(by='MSE',ascending=True)
    # print the top 10
    print(df.head(5))
    # print all the information in the top 5
    for i in range(5):
        print(df.iloc[i])
    # save the ranked and sorted results
    df.to_csv('{}/selection_results_ranked.csv'.format(outdir),index=False)

# take two list of features and plot the correlation between each pair
def plot_2sets_of_feature_correlation(X,set1,set2):
    # get the correlation matrix
    cors = np.zeros((len(set1),len(set2)))
    for i,col1 in enumerate(set1):
        for j,col2 in enumerate(set2):
            corr = np.corrcoef(X[col1],X[col2])[0,1]
            cors[i,j] = corr
    # plot the correlation matrix
    fig, ax = plt.subplots(figsize=(10,10))
    im = ax.imshow(cors, cmap='RdBu')
    # We want to show all ticks...
    ax.set_xticks(np.arange(len(set2)))
    ax.set_yticks(np.arange(len(set1)))
    # ... and label them with the respective list entries
    ax.set_xticklabels(set2)
    ax.set_yticklabels(set1)
    # rotate x ticks
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")
    # add the color scale
    fig.colorbar(im, ax=ax)

    plt.savefig('Figures/feature_correlation.png')
    plt.clf()

# function to remove all columns whose name contains the substring 'plof'
def drop_plof(X):
    cols = [c for c in X.columns if 'plof' not in c]
    return X[cols]

# function to subset the data to only the columns given in a list
def subset_data(X,cols):
    return X[cols]

# plot y pred vs y true and save with a given name
def plot_y_pred_vs_y_true(y_pred,y_true,title,name):
    fig, ax = plt.subplots(figsize=(10,10))
    ax.scatter(y_pred,y_true)
    ax.set_xlabel('y_pred')
    ax.set_ylabel('y_true')
    ax.set_title(title)
    plt.savefig(name)
    plt.clf()

# FinalBOCCFeaturesLOF
# do the same as main but with LOF
def main_lof():
    #'FinalBOCCFeatures/2019/paris.infomap.2019.bocc_res.tsv'
    files_2019 = ['FinalBOCCFeaturesLOF/2019/' + f for f in os.listdir('FinalBOCCFeaturesLOF/2019/')]
    X19, y19 = load_files(files_2019)

    print(X19.shape)
    print('num p < .05',len([i for i in y19 if i < .05]))
    plot_feature_correlation(X19,'Figures/feature_correlation.png')
    
    assert(X19.shape[0] == len(y19))
    find_nan(X19)
    # plot_pca(X19,y19)
    X19 = drop_algo(X19)
    X19 = drop_plof(X19)
    
   
    # # Select All Features for regression
    select_all_features_regression(X19,y19,'RegressionResultsLOF')
    find_best_features_regression('RegressionResultsLOF')

# FinalBOCCFeatures
# do the same as main but with
def pick_just_regression_features():
    #'FinalBOCCFeatures/2019/paris.infomap.2019.bocc_res.tsv'
    files_2019 = ['FinalBOCCFeatures/2019/' + f for f in os.listdir('FinalBOCCFeatures/2019/')]
    X19, y19 = load_files(files_2019)

    print(X19.shape)
    print('num p < .05',len([i for i in y19 if i < .05]))
    plot_feature_correlation(X19,'Figures/feature_correlation.png')
    
    assert(X19.shape[0] == len(y19))
    find_nan(X19)
    # plot_pca(X19,y19)
    X19 = drop_algo(X19)
    X19 = drop_plof(X19)
    
    # # Select All Features for regression
    select_all_features_regression(X19,y19,'JustRegressionResults')
    find_best_features_regression('JustRegressionResults')

def train_and_rank(X_train,y_train,X_test,y_test,outfile,normalize=False,downsample=False):
    # normalize the data
    if normalize:
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X_train)
        X_test = scaler.transform(X_test)
    # downsample the data
    if downsample:
        majority_indices = [i for i,x in enumerate(y_train) if x == 1]
        minority_indices = [i for i,x in enumerate(y_train) if x != 1]
        minority_count = len(minority_indices)
        downsample_count = minority_count
        downsample_indices = random.sample(majority_indices, downsample_count)
        upsample_indices = minority_indices
        keep_indices = downsample_indices + upsample_indices
        X_train = X_train.iloc[keep_indices]
        y_train = [y_train[i] for i in keep_indices]
    # train the model
    model = LinearRegression()
    model.fit(X_train,y_train)
    # get the predictions
    y_pred = model.predict(X_test)
    # get the MSE
    mse = mean_squared_error(y_test,y_pred)
    # get the correlation
    X_test['emprical_p'] = y_test
    X_test['predicted_p'] = y_pred
    # give a rank based on empirical p value
    X_test = X_test.sort_values(by='emprical_p',ascending=True)
    X_test['rank'] = range(1,len(X_test)+1)
    # give a rank based on predicted p value
    X_test = X_test.sort_values(by='predicted_p',ascending=True)
    X_test['rank_pred'] = range(1,len(X_test)+1)
    # plot the ranks
    fig, ax = plt.subplots(figsize=(10,10))
    ax.scatter(X_test['rank'],X_test['rank_pred'])
    ax.set_xlabel('rank')
    ax.set_ylabel('rank_pred')
    ax.set_title('Rank vs Rank Pred')
    plt.savefig(outfile)
    
def rank_2020_with_2019():
    files_2019 = ['FinalBOCCFeatures/2019/' + f for f in os.listdir('FinalBOCCFeatures/2019/')]
    X19, y19 = load_files(files_2019)
    files_2020 = ['FinalBOCCFeatures/2020/' + f for f in os.listdir('FinalBOCCFeatures/2020/')]
    X20, y20 = load_files(files_2020)
    find_nan(X19)
    find_nan(X20)
    X19 = drop_algo(X19)
    X20 = drop_algo(X20)
    
    list1 = ['cluster_size', 'num_sig_go_enrichment_terms', 'num_of_diseases', 'avg_internal_degree', 'normalized_cut', 'triangle_participation_ratio', 'newman_girvan_modularity', 'edges_inside', 'hub_dominance', 'sum_plof']
    list2 = ['max_norm_cell_type_specificity', 'num_of_diseases', 'max_norm_disease_specificity', 'avg_embeddedness', 'conductance', 'cut_ratio', 'normalized_cut', 'expansion', 'triangle_participation_ratio', 'internal_edge_density', 'hub_dominance', 'mean_plof', 'median_plof', 'sum_plof']
    
    X19_r1 = subset_data(X19.copy(),list1)
    X20_r1 = subset_data(X20.copy(),list1)
    X19_r2 = subset_data(X19.copy(),list2)
    X20_r2 = subset_data(X20.copy(),list2)

    train_and_rank(X19_r1,y19,X20_r1,y20,'Figures/rank_2020_with_2019.r1.png',normalize=False,downsample=False)
    train_and_rank(X19_r2,y19,X20_r2,y20,'Figures/rank_2020_with_2019.r2.png',normalize=True,downsample=True)

# function that trains a classifier with optimal features and returns the predictions for ALL the data
def classifier_train_and_predict_all(X,y,normalize=False,downsample=False,params=None):
    random.seed(3)
    Xog = X.copy()
    yog = y.copy()
    # convert y to discrete values
    y = [1 if x != 1 else 0 for x in y]

    # downsample the data
    if downsample:
        majority_indices = [i for i,x in enumerate(y) if x != 1]
        minority_indices = [i for i,x in enumerate(y) if x == 1]
        minority_count = len(minority_indices)
        downsample_count = minority_count
        downsample_indices = random.sample(majority_indices, downsample_count)
        upsample_indices = minority_indices
        keep_indices = downsample_indices + upsample_indices
        X = X.iloc[keep_indices]
        y = [y[i] for i in keep_indices]
    
    # test train split
    X_train, X_test, y_train, y_test = train_test_split(X,y,test_size=0.2,random_state=3)

    # normalize the data
    scaler = None
    if normalize:
        scaler = RobustScaler()
        scaler.fit(X_train)
        X_train = scaler.transform(X_train)
        X_test = scaler.transform(X_test)
    
    # train the model
    if params is not None:
        print('using params for classifier')
        model = xgb.XGBClassifier(**params)
    else:
        model = xgb.XGBClassifier()
    model.fit(X_train,y_train)

    # print precioson and recall
    y_pred = model.predict(X_test)
    print('precision',precision_score(y_test,y_pred))
    print('recall',recall_score(y_test,y_pred))

    # normalize the Xog
    if normalize:
        Xog = scaler.transform(Xog)
    
    # get the predictions for all the data
    y_pred_all = model.predict(Xog)
    y_prob_all = model.predict_proba(Xog)
    return y_pred_all, model, scaler, y_prob_all

# function that trains a regressor
def train_and_predict_regressor(X,y,normalize=False,downsample=False):
    Xog = X.copy()
    yog = y.copy()

    # test train split
    X_train, X_test, y_train, y_test = train_test_split(X,y,test_size=0.2,random_state=42)

    # normalize the data
    scaler = None
    if normalize:
        scaler = RobustScaler()
        scaler.fit(X_train)
        X_train = scaler.transform(X_train)
        X_test = scaler.transform(X_test)
    
    # train the model
    model = xgb.XGBRegressor()
    model.fit(X_train,y_train)

    # print precioson and recall
    y_pred = model.predict(X_test)
    print('mean_squared_error',mean_squared_error(y_test,y_pred))

    # normalize the Xog
    if normalize:
        Xog = scaler.transform(Xog)
    
    # get the predictions for all the data
    y_pred_all = model.predict(Xog)
    return y_pred_all, model, scaler

def train_regressor_with_classifier():
    # load the data
    files_2019 = ['FinalBOCCFeatures/2019/' + f for f in os.listdir('FinalBOCCFeatures/2019/')]
    X19, y19 = load_files(files_2019)
    # drop the algorithm column
    X19 = drop_algo(X19)
    # drop plof columns
    X19 = drop_plof(X19)
    # sebset for classification
    X19_c = subset_data(X19.copy(),['gene_ratio', 'HPO_ratio', 'num_sig_go_enrichment_terms', 'num_of_diseases', 'max_norm_disease_specificity', 'cut_ratio', 'expansion', 'newman_girvan_modularity', 'edges_inside'])

    y_preds = train_and_predict_all(X19_c,y19,normalize=False,downsample=True)

    # subset X and y based on the predictions
    X19 = X19.iloc[y_preds == 1]
    y19 = [y19[i] for i in range(len(y19)) if y_preds[i] == 1]

    # # find best feature set, the results of this are hard coded below since this takes a long time to run
    # select_all_features_regression(X19,y19,'RegressionTrainedOnClassification',norm_only=False,downsample_only=True)
    # find_best_features_regression('RegressionTrainedOnClassification')

    # # subset for regression (these feature are those found to be the best from the previous step)
    X19_r = subset_data(X19.copy(),['cluster_size', 'gene_ratio', 'HPO_ratio', 'num_of_diseases', 'conductance', 'cut_ratio', 'triangle_participation_ratio', 'newman_girvan_modularity', 'internal_edge_density', 'hub_dominance'])
    y_r_pred = train_and_predict_regressor(X19_r, y19, normalize=False, downsample=True)

    # scatter plot of the predictions vs truth
    plt.scatter(y19,y_r_pred,alpha=0.5,s=1)
    plt.xlabel('True')
    plt.ylabel('Predicted')
    plt.savefig('Figures/RegressionTrainedOnClassification.png')
    plt.clf()

# def train a regressor
def train_regressor(X,y,normalize=False,downsample=False,params=None):
    Xog = X.copy()
    yog = y.copy()

    # test train split
    # X_train, X_test, y_train, y_test = train_test_split(X,y,test_size=0.2,random_state=42)

    # normalize the data
    scaler = None
    if normalize:
        scaler = RobustScaler()
        scaler.fit(X)
        X = scaler.transform(X)
        X = pd.DataFrame(X)
    
    # downsample the data
    if downsample:
        majority_indices = [i for i,x in enumerate(y) if x == 1]
        minority_indices = [i for i,x in enumerate(y) if x != 1]
        minority_count = len(minority_indices)
        downsample_count = minority_count
        downsample_indices = random.sample(majority_indices, downsample_count)
        upsample_indices = minority_indices
        keep_indices = downsample_indices + upsample_indices
        X = X.iloc[keep_indices]
        y = [y[i] for i in keep_indices]

    
    # train the model
    if params is None:
        model = xgb.XGBRegressor()
    else:
        print('Using regression params')
        model = xgb.XGBRegressor(**params)
    model.fit(X,y)

    return model, scaler

def train_classifier(X,y,normalize=False,downsample=False,params=None):
    random.seed(0)
    # normalize the data
    scaler = None
    if normalize:
        scaler = RobustScaler()
        scaler.fit(X)
        X = scaler.transform(X)
        X = pd.DataFrame(X)
    
    # downsample the data
    if downsample:
        majority_indices = [i for i,x in enumerate(y) if x != 1]
        minority_indices = [i for i,x in enumerate(y) if x == 1]
        minority_count = len(minority_indices)
        downsample_count = minority_count
        downsample_indices = random.sample(majority_indices, downsample_count)
        upsample_indices = minority_indices
        keep_indices = downsample_indices + upsample_indices
        X = X.iloc[keep_indices]
        y = [y[i] for i in keep_indices]

    
    # train the model
    if params is not None:
        print('using params')
        model = xgb.XGBClassifier(**params)
    else:
        model = xgb.XGBClassifier()
    model.fit(X,y)

    return model, scaler

def rank_plot(x,y,output,title):
    fig, ax = plt.subplots(2,2)
    # set fig size
    fig.set_size_inches(10,10)

    ax[1,0].scatter(x,y,alpha=0.5,s=1)
    ax[1,0].set_xlabel('Rank')
    ax[1,0].set_ylabel('Predicted Rank')
    # remove top and right borders
    ax[1,0].spines['top'].set_visible(False)
    ax[1,0].spines['right'].set_visible(False)
    # plot density of rank at 0,0
    ax[0,0].hist(x,bins=100)
    # plot density of predicted rank at 1,1
    ax[1,1].hist(y,bins=100,orientation='horizontal')
    # remove top and right borders
    ax[1,1].spines['top'].set_visible(False)
    ax[1,1].spines['right'].set_visible(False)
    # remove left and bottom borders
    ax[0,0].spines['top'].set_visible(False)
    ax[0,0].spines['right'].set_visible(False)
    # hide plot 0,1
    ax[0,1].axis('off')
    ax[0,0].set_title(title)
    plt.tight_layout
    plt.savefig(output)
    plt.clf()

def practical_application(train_paths,test_paths,features_c,features_r,c_norm,c_downsample,r_norm,r_downsample,c_params,r_params,threshold=1):
    # load the data
    X_train, y_train = load_files(train_paths)
    X_test, y_test = load_files(test_paths)
    
    # subset the data
    X_train_c = subset_data(X_train.copy(),features_c)
    X_test_c = subset_data(X_train.copy(),features_c)
    
    # train the classifier
    c_y_pred, c_model, c_scaler, c_y_prob_all = classifier_train_and_predict_all(X_train_c,y_train,normalize=c_norm,downsample=c_downsample,params=c_params)

    # change c_y_prob_all to be a list of the probabilities of the positive class
    c_y_prob_all = [x[1] for x in c_y_prob_all]
    c_rank = {'cluster_ID':list(range(len(y_train))),'emperical_p': list(y_train), 'predicted_prob': list(c_y_prob_all)}
    c_rank = pd.DataFrame(c_rank)
    c_rank['emperical_p_rank'] = c_rank['emperical_p'].rank(ascending=True)
    c_rank['predicted_prob_rank'] = c_rank['predicted_prob'].rank(ascending=False)
    # sort by predicted prob rank smaller first
    c_rank = c_rank.sort_values(by=['predicted_prob_rank'])
    print(c_rank.head)

    # kendals tau
    c_tau, c_tau_p = kendals_tau(c_rank['emperical_p_rank'],c_rank['predicted_prob_rank'])
    print('Classification Kendals Tau: ',c_tau)
    print('Classification Kendals Tau p-value: ',c_tau_p)
    rank_plot(c_rank['emperical_p_rank'],c_rank['predicted_prob_rank'],'Figures/rank_plot.classification.19.png','')
    # kendals tau without p=1
    c_rank = c_rank[c_rank['emperical_p'] != 1]
    c_tau, c_tau_p = kendals_tau(c_rank['emperical_p_rank'],c_rank['predicted_prob_rank'])
    print('without p=1 Classification Kendals Tau: ',c_tau)
    print('without p=1 Classification Kendals Tau p-value: ',c_tau_p)
    rank_plot(c_rank['emperical_p_rank'],c_rank['predicted_prob_rank'],'Figures/rank_plot.classification.filtered.19.png','')


    # subset train based on predictions
    X_train_r = X_train.iloc[c_y_pred == 1]
    X_train_r = subset_data(X_train_r.copy(),features_r)
    y_train = [y_train[i] for i in range(len(y_train)) if c_y_pred[i] == 1]

    plot_feature_correlation(X_train_r,'Figures/feature_correlation.regression.19.png')

    # train the regressor
    r_model, r_scaler = train_regressor(X_train_r,y_train,normalize=r_norm,downsample=r_downsample,params=r_params)

    # Predict the test data
    if c_scaler is not None:
        X_test_c = c_scaler.transform(X_test)
    else:
        X_test_c = X_test.copy()
    X_test_c = subset_data(X_test_c.copy(),features_c)
    c_y_pred_test = c_model.predict(X_test_c)

    # print the class break down of c_y_pred_test
    print('Class Breakdown of c_y_pred_test')
    # print(np.bincount(c_y_pred_test.astype(int)))

    # subset c_y_pred_test based on the predictions
    print('prefiltering',X_test.shape)
    X_test_filtered = X_test.iloc[c_y_pred_test == 1]
    y_test_filtered = [y_test[i] for i in range(len(y_test)) if c_y_pred_test[i] == 1]
    print('postfiltering',X_test_filtered.shape)
    
    # predict the test data
    if r_scaler is not None:
        X_test_r = r_scaler.transform(X_test_filtered)
    else:
        X_test_r = X_test_filtered.copy()

    X_test_r = subset_data(X_test_r.copy(),features_r)
    r_y_pred_test = r_model.predict(X_test_r)

    # plot y_test_filtered vs r_y_pred_test
    plt.scatter(y_test_filtered,r_y_pred_test)
    plt.xlabel('emperical p-value')
    plt.ylabel('predicted p-value')
    plt.title('Predicted vs Emperical p-value')
    # remove top and right borders
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig('Figures/classification_regression_predicted_vs_emperical.png')
    plt.clf()
        
    # put r_y_pred_test and c_y_pred_test back together and return
    results = []
    print(len(c_y_pred_test))
    r_index = 0
    for i,x in enumerate(c_y_pred_test):
        if x == 0:
            results.append(np.inf)
        else:
            # print(r_y_pred_test[i])
            results.append(r_y_pred_test[r_index])
            r_index += 1
    assert len(results) == len(c_y_pred_test)
    # how mnay results are not inf
    print('num not inf',len([i for i in results if i != np.inf]))
    # how many total results   
    print('num total',len(results))
    # assign a rank to the results, if the result is 0, assign it the rank of the lowest value in the test set
    results_rank = pd.DataFrame({'result':results,'index':list(range(len(results)))})
    # sort the df in increasing order
    results_rank = results_rank.sort_values(by='result',ascending=True)
    print('increasing order')
    # print(results_rank)
    # assign a rank to each result allow for ties
    results_rank['rank'] = results_rank['result'].rank(method='average')
    print('ranked')
    # print(results_rank)
    # sort the df in the original order
    results_rank = results_rank.sort_values(by='index',ascending=True)
    print('original order')
    # print(results_rank)
    # return the rank
    return results_rank['rank'].tolist()

def main():
    #'FinalBOCCFeatures/2019/paris.infomap.2019.bocc_res.tsv'
    files_2019 = ['FinalBOCCFeatures/2019/' + f for f in os.listdir('FinalBOCCFeatures/2019/')]
    X19, y19 = load_files(files_2019)

    print(X19.shape)
    print('num p < .05',len([i for i in y19 if i < .05]))
    plot_feature_correlation(X19,'Figures/feature_correlation.png')
    
    assert(X19.shape[0] == len(y19))
    find_nan(X19)
    # plot_pca(X19,y19)
    X19 = drop_algo(X19)
    X19 = drop_plof(X19)
    
   
    # # Select All Features:
    select_all_features(X19,y19)
    find_best_features()

    # Regression
    # train_model_regressor(X19,y19,normaize=True,downsample=True)
    select_all_features_regression(X19,y19)
    print()
    find_best_features_regression()
    # best_regression1 = ['gene_ratio', 'HPO_ratio', 'num_sig_go_enrichment_terms', 'max_norm_cell_type_specificity', 'num_of_diseases', 'max_norm_disease_specificity', 'avg_embeddedness', 'conductance', 'cut_ratio', 'normalized_cut', 'triangle_participation_ratio', 'newman_girvan_modularity', 'edges_inside', 'hub_dominance', 'mean_plof']
    # best_regression2 = ['cluster_size', 'gene_ratio', 'HPO_ratio', 'num_sig_go_enrichment_terms', 'num_of_diseases', 'conductance', 'cut_ratio', 'expansion', 'triangle_participation_ratio', 'surprise', 'newman_girvan_modularity', 'internal_edge_density', 'hub_dominance', 'max_plof', 'mean_plof', 'std_plof']
    
    # # subset the data to only the best features
    # r1_X19 = subset_data(X19.copy(),best_regression1)
    # r2_X19 = subset_data(X19.copy(),best_regression2)

    # # train the model on the best features
    # mse1, ypred1, ytest1 = train_model_regressor(r1_X19,y19,normaize=False,downsample=False)
    # mse2, ypred2, ytest2 = train_model_regressor(r2_X19,y19,normaize=True,downsample=True)
    
    # plot_y_pred_vs_y_true(ypred1,ytest1,'Not Downsampled & Not Normalized','Figures/y_pred_vs_y_true_regression1.png')
    # plot_y_pred_vs_y_true(ypred2,ytest2,'Downsampled & Normalized','Figures/y_pred_vs_y_true_regression2.png')
    
    # list1 = ['cluster_size', 'num_sig_go_enrichment_terms', 'num_of_diseases', 'avg_internal_degree', 'normalized_cut', 'triangle_participation_ratio', 'newman_girvan_modularity', 'edges_inside', 'hub_dominance', 'sum_plof']
    # list2 = ['max_norm_cell_type_specificity', 'num_of_diseases', 'max_norm_disease_specificity', 'avg_embeddedness', 'conductance', 'cut_ratio', 'normalized_cut', 'expansion', 'triangle_participation_ratio', 'internal_edge_density', 'hub_dominance', 'mean_plof', 'median_plof', 'sum_plof']
    # plot_2sets_of_feature_correlation(X19,list1,list2)

def read_file(filename, fields_to_add):
    data = []
    with open(filename, 'r') as file:
        for line in file:
            fields = line.strip().split('\t')
            tuple_fields = []
            for i in fields_to_add:
                tuple_fields.append(fields[i])
            data.append(tuple(tuple_fields))
    return data

def get_top_n(scores, n):
    top_scores = []
    current_score = None
    count = 0
    for i in range(len(scores)):
        if count >= n:
            if scores[i][1] == current_score: top_scores.append(scores[i])
            else: break
        else:
            if current_score is None or scores[i][1] > current_score:
                current_score = scores[i][1]
            top_scores.append(scores[i])
            count += 1
    return top_scores

def intersection_count(list1, list2):
    set1 = set(list1)
    set2 = set(list2)
    return len(set1.intersection(set2))

def get_cdf(data):
    count, bins_count = np.histogram(data, bins=100)
    pdf = count / sum(count)
    cdf = np.cumsum(pdf)
    return cdf

def top_x_roc(e_file,p_file,outfile):
    E = read_file(e_file, [0,1])
    E.sort(key=lambda x: x[1])
    P = read_file(p_file, [0,1])
    P.sort(key=lambda x: x[1])

    top_1 = get_top_n(E, 1) 
    tp_fp = [] 
    for i in range(1,len(E)):
        top_E = [x[0] for x in get_top_n(E, i)]
        top_P = [x[0] for x in get_top_n(P, i)]
        TP = intersection_count(top_E, top_P)
        FP = len(top_P) - TP
        # print(TP,FP)
        tp_fp.append((TP,FP))
    # print(tp_fp)

    fp = [x[1] for x in tp_fp]
    tp = [x[0] for x in tp_fp]

    cdf_tp = get_cdf(tp)
    cdf_fp = get_cdf(fp)

    fig, ax  = plt.subplots()
    ax.plot(cdf_fp,cdf_tp)
    ax.set_ylabel('CDF TP')
    ax.set_xlabel('CDF FP')
    # remove top and right border
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.savefig(outfile)
    plt.clf()

def do_19v20_practical_application():
    files_2019 = ['FinalBOCCFeatures/2019/' + f for f in os.listdir('FinalBOCCFeatures/2019/')]
    files_2020 = ['FinalBOCCFeatures/2020/' + f for f in os.listdir('FinalBOCCFeatures/2020/')]
    c_features = ['gene_ratio', 'HPO_ratio', 'num_sig_go_enrichment_terms', 'num_of_diseases', 'max_norm_disease_specificity', 'cut_ratio', 'expansion', 'newman_girvan_modularity', 'edges_inside']
    r_features = ['cluster_size', 'gene_ratio', 'HPO_ratio', 'num_of_diseases', 'conductance', 'cut_ratio', 'triangle_participation_ratio', 'newman_girvan_modularity', 'internal_edge_density', 'hub_dominance']
    ranks_20 = practical_application(files_2019,
                          files_2020,
                          c_features,
                          r_features,
                          c_norm=False,
                          c_downsample=True,
                          r_norm=False,
                          r_downsample=True,
                          c_params=classification_params,
                          r_params=regression_params)
    
    # load the 2020 files and get the p values
    X20, y20 = load_files(files_2020)
    print('X2-shape',X20.shape)
    print('# predicted',len(ranks_20))
    # plot the ranks vs their actual p values
    # plot_ranks_vs_pvalues(ranks_20,y20,'Figures/ranks_vs_pvalues.19v20.png')
    # create a df of ranks_20 with clusterIDs
    p_df = pd.DataFrame({'clusterID':list(range(len(y20))),'emprical_p':y20})
    ranks_df = pd.DataFrame({'clusterID':list(range(len(y20))),'rank':ranks_20})
    df = pd.DataFrame({'clusterID':list(range(len(y20))),'predicted_p_rank':ranks_20,'emprical_p':y20})
    # create a rank columns based on emprical_p
    df['rank'] = df['emprical_p'].rank(method='min')
    print(df)
    # sort by clusterID
    df.sort_values(by=['clusterID'],inplace=True)
    # write p_df to file
    p_df.to_csv('.pvalues.19v20.tsv',index=False,sep='\t')
    # write ranks_df to file
    ranks_df.to_csv('.ranks.19v20.tsv',index=False,sep='\t')

def do_19v20_practical_application_p35():
    # set seed
    random.seed(42)
    files_2019 = ['FinalBOCCFeatures/2019/' + f for f in os.listdir('FinalBOCCFeatures/2019/')]
    files_2020 = ['FinalBOCCFeatures/2020/' + f for f in os.listdir('FinalBOCCFeatures/2020/')]
    c_features = ['num_sig_go_enrichment_terms', 'num_of_diseases', 'avg_embeddedness', 'conductance', 'normalized_cut', 'triangle_participation_ratio', 'newman_girvan_modularity', 'edges_inside']
    r_features = ['num_sig_go_enrichment_terms', 'num_of_diseases', 'avg_embeddedness', 'conductance', 'normalized_cut', 'triangle_participation_ratio', 'newman_girvan_modularity', 'edges_inside']
    # load the files
    X19, y19 = load_files(files_2019)
    X20, y20 = load_files(files_2020)
    # subset the features
    X19_c = X19[c_features]
    X20_c = X20[c_features]
    X19_r = X19[r_features]
    X20_r = X20[r_features]
    # subset X19 to be only those with y < .35
    threshold = .35
    X19_r = X19_r.iloc[[i for i,x in enumerate(y19) if x < threshold]]
    y19_r = [y19[i] for i in range(len(y19)) if y19[i] < threshold]
    # train the classifier
    # threshold y19 for training classifier
    y19_c = [1 if p < threshold else 0 for p in y19]
    # downsample
    majority_indices = [i for i,x in enumerate(y19_c) if x == 0]
    minority_indices = [i for i,x in enumerate(y19_c) if x == 1]
    minority_count = len(minority_indices)
    downsample_count = minority_count
    downsample_indices = random.sample(majority_indices, downsample_count)
    upsample_indices = minority_indices
    keep_indices = downsample_indices + upsample_indices
    X19_c = X19_c.iloc[keep_indices]
    y19_c = [y19_c[i] for i in keep_indices]
    # create classifier
    c_model, c_scaler = train_classifier(X19_c,y19_c,normalize=False,downsample=False,params=HYPERPARAMS_p35)
    # train the regressor
    r_model, r_scaler = train_regressor(X19_r,y19_r,normalize=False,downsample=False,params=HYPERPARAMS_p35_regression)
    # predict X20_c with classifier
    if c_scaler is not None:
        X20_c = c_scaler.transform(X20_c)
    y20_c_pred = c_model.predict(X20_c)
    # subset X20_r based on y20_c_pred
    X20_c_r = X20_r.iloc[[i for i,x in enumerate(y20_c_pred) if x == 1]]
    # predict X20_c_r with regressor
    if r_scaler is not None:
        X20_c_r = r_scaler.transform(X20_c_r)
    y20_c_r_pred = r_model.predict(X20_c_r)
    # create a df of clusterID and predicted p values
    predictions = []
    r_i = 0
    for x in y20_c_pred:
        if x == 0:
            # indivitive that the p-value is bad
            predictions.append(1)
        else:
            predictions.append(y20_c_r_pred[r_i])
            r_i += 1
    df = pd.DataFrame({'clusterID':list(range(len(y20))),'predicted_p':predictions,'emprical_p':y20})
    # rank y20 is ascending order
    df['real_rank'] = df['emprical_p'].rank(method='min')
    # rank predicted_p is ascending order too
    df['predicted_rank'] = df['predicted_p'].rank(method='min')
    # subset df to include only those with predicted_p != 1
    df = df[df['emprical_p'] != 1]
    t,p = kendals_tau(list(df['real_rank']), list(df['predicted_rank']))
    print('Kendals Tau',t)
    print('Kendals Tau p-value',p)
    # plot rank vs rank
    rank_plot(df['real_rank'],df['predicted_rank'],'Figures/rank_plot.19v20.p35.png','p<.35 classification -> regression\ntau={}'.format(t))
    


def do_19v20_regression_only_roc():
    files_2019 = ['FinalBOCCFeatures/2019/' + f for f in os.listdir('FinalBOCCFeatures/2019/')]
    files_2020 = ['FinalBOCCFeatures/2020/' + f for f in os.listdir('FinalBOCCFeatures/2020/')]
    
    # load the files
    X19, y19 = load_files(files_2019)
    X20, y20 = load_files(files_2020)
    
    r_features = ['cluster_size', 'gene_ratio', 'num_of_diseases', 'max_norm_disease_specificity', 'conductance', 'edges_inside', 'hub_dominance']
    # sub set the features
    X19 = X19[r_features]
    X20 = X20[r_features]
    model, scaler = train_regressor(X19,y19,normalize=True,downsample=True)
    # predict X20
    X20 = scaler.transform(X20)
    y20_pred = model.predict(X20)
    df = pd.DataFrame({'clusterID':list(range(len(y20))),'predicted_p':y20_pred,'emprical_p':y20})
    # save df of just clusterID and predicted p values
    df1 = df[['clusterID','predicted_p']]
    df1.to_csv('.predicted_p.regression_only.19v20.tsv',index=False,sep='\t')
    # save df of just clusterID and emprical p values
    df2 = df[['clusterID','emprical_p']]
    df2.to_csv('.emprical_p.regression_only.19v20.tsv',index=False,sep='\t')
    # plot the roc curve
    top_x_roc('.predicted_p.regression_only.19v20.tsv','.emprical_p.regression_only.19v20.tsv','Figures/roc_curve.regression_only.19v20.png')

    # create columns of ranks based on predicted p values
    df['rank'] = df['predicted_p'].rank(ascending=False)
    # create a rank for actual p values
    df['actual_rank'] = df['emprical_p'].rank(ascending=False)
    # sort by clusterID
    df = df.sort_values(by=['clusterID'])
    # w = kendals_w(list(df['rank']),list(df['actual_rank']))
    # print('Regression Only', w)

def train_classifier_with_threshold_and_params(X_train,X_test,y_train,y_test,threshold,params):
    y_train = [1 if p < threshold else 0 for p in y_train]
    y_test = [1 if p < threshold else 0 for p in y_test]
    # print the number p == 1 in y19
    print('y19',sum(y_train),len(y_train))
    print('y20',sum(y_test),len(y_test))

    # train the classifier
    model, scaler = train_classifier(X_train,y_train,normalize=False,downsample=True,params=params)
    # predict X20
    if scaler is not None:
        X_test = scaler.transform(X_test)
        X_test = pd.DataFrame(X_test,columns=features)
    
    y_test_pred = model.predict(X_test)
    y_test_probs = model.predict_proba(X_test)
    # print(y20_probs)
    y_test_probs = y_test_probs[:,1]
    return y_test, y_test_pred, y_test_probs

def classifier_roc():
    features = ['gene_ratio', 'HPO_ratio', 'num_sig_go_enrichment_terms', 'num_of_diseases', 'max_norm_disease_specificity', 'cut_ratio', 'expansion', 'newman_girvan_modularity', 'edges_inside']
    # load the 2019 files
    files_2019 = ['FinalBOCCFeatures/2019/' + f for f in os.listdir('FinalBOCCFeatures/2019/')]
    X19, y19 = load_files(files_2019)
    # sub set the features
    X19 = X19[features]
    plot_feature_correlation(X19,'Figures/feature_correlation.classification_9.2019.png')
    # load the 2020 files
    files_2020 = ['FinalBOCCFeatures/2020/' + f for f in os.listdir('FinalBOCCFeatures/2020/')]
    X20, y20 = load_files(files_2020)
    og_y20 = y20.copy()
    # sub set the features
    X20 = X20[features]

    y20_real_p5, y20_pred_p5, y20_probs_p5 = train_classifier_with_threshold_and_params(X19,X20,y19,y20,0.05,None)
    y20_real_1, y20_pred_1, y20_probs_1 = train_classifier_with_threshold_and_params(X19,X20,y19,y20,1,None)
    y20_real_1_params, y20_pred_1_params, y20_probs_1_params = train_classifier_with_threshold_and_params(X19,X20,y19,y20,1,classification_params)

    fpr_1_params, tpr_1_params, thresh_1_params = roc_curve(y20_real_1_params,y20_probs_1_params)
    fpr_1, tpr_1, thresh_1 = roc_curve(y20_real_1,y20_probs_1)
    fpr_p5, tpr_p5, thresh_p5 = roc_curve(y20_real_p5,y20_probs_p5)
    roc_auc_1_params = roc_auc_score(y20_real_1_params,y20_probs_1_params)
    roc_auc_1 = roc_auc_score(y20_real_1,y20_probs_1)
    roc_auc_p5 = roc_auc_score(y20_real_p5,y20_probs_p5)

    print('ROC AUC 1 params',roc_auc_1_params)
    print('ROC AUC 1',roc_auc_1)
    print('ROC AUC p5',roc_auc_p5)
    
    # print confusion matrix for each of the thresholds
    print('Confusion Matrix t=1 with params')
    print(confusion_matrix(y20_real_1_params,y20_pred_1_params))
    print('Confusion Matrix t=1')
    print(confusion_matrix(y20_real_1,y20_pred_1))
    print('Confusion Matrix t=0.05')
    print(confusion_matrix(y20_real_p5,y20_pred_p5))

    # plot the roc curve
    plt.clf()
    fig, ax = plt.subplots()
    ax.plot(fpr_1_params, tpr_1_params, label='p<1.00 with hyperparams')
    ax.plot(fpr_1, tpr_1, label='p<1.00')
    ax.plot(fpr_p5, tpr_p5, label='p<0.05')
    ax.plot([0, 1], [0, 1], 'k--')
    ax.set_xlabel('False Positive Rate')    
    ax.set_ylabel('True Positive Rate')
    ax.set_title('ROC Curve')
    ax.legend(loc="lower right")
    # remove top and right border
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.savefig('Figures/roc_curve.classifier.19v20.png')
    plt.clf()

def start_grid_search():
    # load the 2019 files
    files_2019 = ['FinalBOCCFeatures/2019/' + f for f in os.listdir('FinalBOCCFeatures/2019/')]
    X19, y19 = load_files(files_2019)
    r_features = ['cluster_size', 'gene_ratio', 'HPO_ratio', 'num_of_diseases', 'conductance', 'cut_ratio', 'triangle_participation_ratio', 'newman_girvan_modularity', 'internal_edge_density', 'hub_dominance']
    # sub set the features
    X19 = X19[r_features]
    grid_search(X19,y19)

def start_grid_search():
    # load the 2019 files
    files_2019 = ['FinalBOCCFeatures/2019/' + f for f in os.listdir('FinalBOCCFeatures/2019/')]
    X19, y19 = load_files(files_2019)
    features = ['gene_ratio', 'HPO_ratio', 'num_sig_go_enrichment_terms', 'num_of_diseases', 'max_norm_disease_specificity', 'cut_ratio', 'expansion', 'newman_girvan_modularity', 'edges_inside']
    # sub set the features
    X19 = X19[features]
    # threshold y
    threshold = 1
    y19 = [1 if p < threshold else 0 for p in y19]
    grid_search_classification(X19,y19)

def start_genetic_algo_search():
    # load the 2019 files
    files_2019 = ['FinalBOCCFeatures/2019/' + f for f in os.listdir('FinalBOCCFeatures/2019/')]
    X19, y19 = load_files(files_2019)
    features = ['gene_ratio', 'HPO_ratio', 'num_sig_go_enrichment_terms', 'num_of_diseases', 'max_norm_disease_specificity', 'cut_ratio', 'expansion', 'newman_girvan_modularity', 'edges_inside']
    # sub set the features
    X19 = X19[features]
    # threshold y
    threshold = 1
    y19 = [1 if p < threshold else 0 for p in y19]
    genetic_optimization_classification(X19, y19, downsample=True)

def start_genetic_algo_search_regression():
    # load the 2019 files
    files_2019 = ['FinalBOCCFeatures/2019/' + f for f in os.listdir('FinalBOCCFeatures/2019/')]
    X19, y19 = load_files(files_2019)
    og_y19 = y19.copy()
    # regression features
    features = ['cluster_size', 'gene_ratio', 'num_of_diseases', 'max_norm_disease_specificity', 'conductance', 'edges_inside', 'hub_dominance']    
    # sub set the features
    X19 = X19[features]
    # threshold y
    threshold = 1
    y19 = [1 if p < threshold else 0 for p in y19]
    # train the classifier
    model, scaler = train_classifier(X19,y19,normalize=False,downsample=True,params=classification_params)
    # predict X19
    if scaler is not None:
        X19 = scaler.transform(X19)
        X19 = pd.DataFrame(X19,columns=features)
    y19_pred = model.predict(X19)
    # subset X19 for that that are predicted to be 1
    X19_for_regression = X19[y19_pred == 1]
    # y19_for_regression = og_y19[y19_pred == 1]
    y19_for_regression = [og_y19[i] for i,y in enumerate(y19_pred) if y == 1]
    print('X19_for_regression',X19_for_regression.shape)
    genetic_optimization_regression(X19_for_regression, y19_for_regression)

def start_genetic_algo_for_just_regression():
    # load the 2019 files
    files_2019 = ['FinalBOCCFeatures/2019/' + f for f in os.listdir('FinalBOCCFeatures/2019/')]
    # load 2020
    files_2020 = ['FinalBOCCFeatures/2020/' + f for f in os.listdir('FinalBOCCFeatures/2020/')]
    X19, y19 = load_files(files_2019)
    X20, y20 = load_files(files_2020)
    # combine X19 and X20
    X = pd.concat([X19,X20])
    y = y19 + y20
    og_y19 = y19.copy()
    # regression features
    features = ['num_sig_go_enrichment_terms', 'num_of_diseases', 'avg_embeddedness', 'conductance', 'normalized_cut', 'triangle_participation_ratio', 'newman_girvan_modularity', 'edges_inside']
    # sub set the features
    X = X[features]
    genetic_optimization_regression(X, y,downsample=True)

def start_genetic_algo_for_p35_regression():
    # load the 2019 files
    files_2019 = ['FinalBOCCFeatures/2019/' + f for f in os.listdir('FinalBOCCFeatures/2019/')]
    X19, y19 = load_files(files_2019)
    og_y19 = y19.copy()
    # regression features
    features = ['num_sig_go_enrichment_terms', 'num_of_diseases', 'avg_embeddedness', 'conductance', 'normalized_cut', 'triangle_participation_ratio', 'newman_girvan_modularity', 'edges_inside']
    # sub set the features
    X19 = X19[features]
    # threshold y
    threshold = 0.35
    # remove samples from X first
    X19 = X19[[x < threshold for x in y19]]
    y19 = [p for p in y19 if p < threshold]
    genetic_optimization_regression(X19, y19,downsample=False)

def normalize_rank(xs):
    # normalize list x to range from 0-1
    return [x / max(xs) for x in xs]

def filter_and_flip(train_files,test_files,features,plot_prefix,params=None,downsample=True,threshold=1,plot=True):
    # load the train files
    X_train, y_train = load_files(train_files)
    # load the test files
    X_test, y_test = load_files(test_files)
    # sub set the features
    X_train = X_train[features]
    X_test = X_test[features]
    # threshold y
    y_train_og = y_train.copy()
    y_test_og = y_test.copy()
    y_train = [1 if p < threshold else 0 for p in y_train]
    y_test = [1 if p < threshold else 0 for p in y_test]
   
    # train the classifier
    model, scaler = train_classifier(X_train,y_train,normalize=False,downsample=False,params=params)
    # predict X_test
    if scaler is not None:
        X_test = scaler.transform(X_test)
        X_test = pd.DataFrame(X_test,columns=features)
    y_test_pred = model.predict(X_test)
    # get test probs
    y_test_probs = model.predict_proba(X_test)

    # print the precision and recall
    print('precision',precision_score(y_test, y_test_pred))
    print('recall',recall_score(y_test, y_test_pred))

    # ROC curve for the classifier
    fpr, tpr, thresholds = roc_curve(y_test, y_test_probs[:,1])
    auc = roc_auc_score(y_test, y_test_probs[:,1])

    if not plot:
        return fpr, tpr, auc

    plt.plot(fpr, tpr, label='ROC curve (area = %0.2f)' % auc)
    plt.plot([0, 1], [0, 1], 'k--')  # random predictions curve
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    # remove top and right border
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # set figure size
    plt.gcf().set_size_inches(6, 4)
    plt.legend(loc="lower right")
    plt.savefig(plot_prefix + '_classifier_roc.png',dpi=300)
    plt.clf()

    # make a df with the ranks
    print(model.classes_)
    y_test_probs = [x[1] for x in y_test_probs]
    c_rank = {'cluster_ID':list(range(len(y_test))),'emperical_p': list(y_test_og), 'predicted_prob': list(y_test_probs)}
    c_rank = pd.DataFrame(c_rank)
    c_rank['emperical_p_rank'] = c_rank['emperical_p'].rank(ascending=True)
    c_rank['predicted_prob_rank'] = c_rank['predicted_prob'].rank(ascending=True)
    c_rank['emperical_p_rank_normalized'] = normalize_rank(c_rank['emperical_p_rank'])
    c_rank['predicted_prob_rank_normalized'] = normalize_rank(c_rank['predicted_prob_rank'])

    c_rank = c_rank[c_rank['emperical_p'] != 1]

    # kendall tau
    tau, p_value = kendals_tau(c_rank['emperical_p_rank'], c_rank['predicted_prob_rank'])
    print('Kendall Tau:',tau)
    print('p-value:',p_value)

    # pearson correlation
    corr, p_value = pearsonr(c_rank['emperical_p_rank'], c_rank['predicted_prob_rank'])
    print('Pearson Correlation:',corr)
    print('p-value:',p_value)

    # remove p = 1
    rank_plot(c_rank['emperical_p_rank'], c_rank['predicted_prob_rank'], plot_prefix + '_classifier_rank.png',title='Kendall\'s Tau = ' + str(round(tau,2)))
    return fpr, tpr, auc


    

def do_all_filter_and_flip():
    print('Filter and Flip 2019 - 2022')
    # list 2019 files
    files_2019 = ['FinalBOCCFeatures/2019/' + f for f in os.listdir('FinalBOCCFeatures/2019/')]
    # list 2020 files
    files_2020 = ['FinalBOCCFeatures/2020/' + f for f in os.listdir('FinalBOCCFeatures/2020/')]
    # list 2021 files
    files_2021 = ['FinalBOCCFeatures/2021/' + f for f in os.listdir('FinalBOCCFeatures/2021/')]

    features = ['gene_ratio', 'HPO_ratio', 'num_sig_go_enrichment_terms', 'num_of_diseases', 'max_norm_disease_specificity', 'cut_ratio', 'expansion', 'newman_girvan_modularity', 'edges_inside']

    print('With Threshold = 1.00')
    print('2019 2020')
    filter_and_flip(files_2019,
                    files_2020, 
                    features,
                    plot_prefix='Figures/2019v2020_',
                    params=classification_params,
                    downsample=True)
    print('2019 2021')
    filter_and_flip(files_2019,
                    files_2021,
                    features,
                    plot_prefix='Figures/2019v2021_',
                    params=classification_params,
                    downsample=True)
    
    print('With Threshold = .1')
    print('2019 2020')
    filter_and_flip(files_2019,
                    files_2020, 
                    features,
                    plot_prefix='Figures/2019v2020_tp1_',
                    params=HYPERPARAMS_p1,
                    downsample=True,
                    threshold=.1)
    print('2019 2021')
    filter_and_flip(files_2019,
                    files_2021,
                    features,
                    plot_prefix='Figures/2019v2021_tp1_',
                    params=HYPERPARAMS_p1,
                    downsample=True,
                    threshold=.1)
    
    print('With Threshold = .35')
    print('2019 2020')
    filter_and_flip(files_2019,
                    files_2020, 
                    features,
                    plot_prefix='Figures/2019v2020_tp35_',
                    params=HYPERPARAMS_p35,
                    downsample=True,
                    threshold=.35)
    print('2019 2021')
    filter_and_flip(files_2019,
                    files_2021,
                    features,
                    plot_prefix='Figures/2019v2021_tp35_',
                    params=HYPERPARAMS_p35,
                    downsample=True,
                    threshold=.35)
    
    print('With Threshold = .05')
    print('2019 2020')
    filter_and_flip(files_2019,
                    files_2020, 
                    features,
                    plot_prefix='Figures/2019v2020_tp05_',
                    params=HYPERPARAMS_p05,
                    downsample=True,
                    threshold=.05)
    print('2019 2021')
    filter_and_flip(files_2019,
                    files_2021,
                    features,
                    plot_prefix='Figures/2019v2021_tp05_',
                    params=HYPERPARAMS_p05,
                    downsample=True,
                    threshold=.05)
    
def threshold_rocs():
    print('Filter and Flip 2019 - 2022')
    # list 2019 files
    files_2019 = ['FinalBOCCFeatures/2019/' + f for f in os.listdir('FinalBOCCFeatures/2019/')]
    # list 2020 files
    files_2020 = ['FinalBOCCFeatures/2020/' + f for f in os.listdir('FinalBOCCFeatures/2020/')]
    # list 2021 files
    files_2021 = ['FinalBOCCFeatures/2021/' + f for f in os.listdir('FinalBOCCFeatures/2021/')]

    features = ['gene_ratio', 'HPO_ratio', 'num_sig_go_enrichment_terms', 'num_of_diseases', 'max_norm_disease_specificity', 'cut_ratio', 'expansion', 'newman_girvan_modularity', 'edges_inside']
 
    res = {'threshold':[],'auc':[]}
    res_tp_fp = {'threshold':[],'tpr':[],'fpr':[]}
    for t in [0.01,.05,.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.75,.8,.85,.9,.95,1.00]:
        print(t)
        fpr, tpr, auc = filter_and_flip(files_2019,
                        files_2020, 
                        features,
                        plot_prefix='Figures/2019v2020_',
                        params=classification_params,
                        downsample=True,
                        threshold=t,
                        plot=False)
        print(t, auc)
        res['threshold'].append(t)
        res['auc'].append(auc)
        res_tp_fp['threshold'].append(t)
        res_tp_fp['tpr'].append(tpr)
        res_tp_fp['fpr'].append(fpr)
    # plot t vs auc
    fig, ax = plt.subplots()
    ax.plot(res['threshold'],res['auc'])
    ax.set_xlabel('Threshold')
    ax.set_ylabel('AUC')
    # remove top and right spines from the ax
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # set x ticks
    ax.set_xticks([0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1])
    plt.tight_layout()
    plt.savefig('Figures/threshold_auc.png')
    plt.clf()

    
def just_regression_train_and_test(train_files,test_files,features,params,downsample=True):
    # load the train files
    X_train, y_train = load_files(train_files)
    # load the test files
    X_test, y_test = load_files(test_files)
    # sub set the features
    X_train = X_train[features]
    X_test = X_test[features]
    # downsample
    if downsample:
        majority_indices = [i for i,x in enumerate(y_train) if x == 1]
        minority_indices = [i for i,x in enumerate(y_train) if x != 1]
        majority_count = len(majority_indices)
        minority_count = len(minority_indices)
        downsample_count = minority_count
        print('Downsampling {} samples to {}'.format(majority_count, downsample_count))
        downsample_indices = random.sample(majority_indices, downsample_count)
        upsample_indices = minority_indices
        keep_indices = downsample_indices + upsample_indices
        X_train = X_train.iloc[keep_indices]
        y_train = [y_train[i] for i in keep_indices]
    
    # train the model
    model = xgb.XGBRegressor(**params)
    model.fit(X_train, y_train)
    # predict
    y_test_probs = model.predict(X_test)
    # create a df with id, real rank, predicted rank
    c_rank = {'cluster_ID':list(range(len(y_test))),'emperical_p': list(y_test), 'predicted_prob': list(y_test_probs)}
    c_rank = pd.DataFrame(c_rank)
    c_rank['emperical_p_rank'] = c_rank['emperical_p'].rank(ascending=True)
    c_rank['predicted_prob_rank'] = c_rank['predicted_prob'].rank(ascending=True)
    return list(c_rank['emperical_p_rank']), list(c_rank['predicted_prob_rank']), y_test

def just_regression():
    features = ['num_sig_go_enrichment_terms', 'num_of_diseases', 'avg_embeddedness', 'conductance', 'normalized_cut', 'triangle_participation_ratio', 'newman_girvan_modularity', 'edges_inside']
    # list 2019 files
    files_2019 = ['FinalBOCCFeatures/2019/' + f for f in os.listdir('FinalBOCCFeatures/2019/')]
    # list 2020 files
    files_2020 = ['FinalBOCCFeatures/2020/' + f for f in os.listdir('FinalBOCCFeatures/2020/')]
    # list 2021 files
    files_2021 = ['FinalBOCCFeatures/2021/' + f for f in os.listdir('FinalBOCCFeatures/2021/')]
    real_ranks20, pred_ranks_20, y_test_20 = just_regression_train_and_test(files_2019,files_2020,features,regression_params)
    tau20, p20 = kendals_tau(real_ranks20, pred_ranks_20)
    rank_plot(real_ranks20, pred_ranks_20, 'Figures/2019v2020_just_regression_rank.png',title='Kendall\'s Tau = ' + str(round(tau20,2)))
    # filter real_ranks20 and pred_ranks_20 to only include points with real rank < 100
    pred_ranks_20 = [pred_ranks_20[i] for i,x in enumerate(y_test_20) if x < 1]
    real_ranks20 = [real_ranks20[i] for i,x in enumerate(y_test_20) if x < 1]
    tau20, p20 = kendals_tau(real_ranks20, pred_ranks_20)
    rank_plot(real_ranks20, pred_ranks_20, 'Figures/2019v2020_just_regression_rank_non_p1.png',title='Kendall\'s Tau = ' + str(round(tau20,2)))


def do_genetic_optimization_for_p1_and_p35():
    # list 2019 files
    files_2019 = ['FinalBOCCFeatures/2019/' + f for f in os.listdir('FinalBOCCFeatures/2019/')]
    # these are the features determined from using just regression, JustRegressionResults/
    features = ['num_sig_go_enrichment_terms', 'num_of_diseases', 'avg_embeddedness', 'conductance', 'normalized_cut', 'triangle_participation_ratio', 'newman_girvan_modularity', 'edges_inside']
    # load files
    X, y = load_files(files_2019)
    # sub set the features
    X = X[features]

    # print('-----------------Threshold = 0.1-----------------')
    # # threshold the ys as .1
    yp1 = [1 if x < .1 else 0 for x in y]
    genetic_optimization_classification(X,yp1,downsample=True)
    # print('\n\n\n\n\n')
    # print('-----------------Threshold = 0.35-----------------')
    yp35 = [1 if x < .35 else 0 for x in y]
    genetic_optimization_classification(X,yp35,downsample=True)
    print('-----------------Threshold = 0.05-----------------')
    yp05 = [1 if x < .05 else 0 for x in y]
    genetic_optimization_classification(X,yp05,downsample=True)


def do_genetic_optimization_for_p1_and_p35_with_2019_and_2020():
    # list 2019 files
    files_2019 = ['FinalBOCCFeatures/2019/' + f for f in os.listdir('FinalBOCCFeatures/2019/')]
    files_2020 = ['FinalBOCCFeatures/2020/' + f for f in os.listdir('FinalBOCCFeatures/2020/')]
    files_2021 = ['FinalBOCCFeatures/2021/' + f for f in os.listdir('FinalBOCCFeatures/2021/')]
    # these are the features determined from using just regression, JustRegressionResults/
    features = ['num_sig_go_enrichment_terms', 'num_of_diseases', 'avg_embeddedness', 'conductance', 'normalized_cut', 'triangle_participation_ratio', 'newman_girvan_modularity', 'edges_inside']
    # load files
    X, y = load_files(files_2019)
    X20, y20 = load_files(files_2020)
    X21, y21 = load_files(files_2021)

    # create list of names for the samples in X19 based on algo and the row index
    X19_cluster_ids = ['{}.{}:{}'.format(a,'2019',str(i)) for a,i in zip(X['algo'],X['cluster_id'])]
    X20_cluster_ids = ['{}.{}:{}'.format(a,'2020',str(i)) for a,i in zip(X20['algo'],X20['cluster_id'])]
    X21_cluster_ids = ['{}.{}:{}'.format(a,'2021',str(i)) for a,i in zip(X21['algo'],X21['cluster_id'])]

    X = pd.concat([X,X20])
    y = y + y20
    # sub set the features
    X_og = X.copy()
    X21_og = X21.copy()

    X = X[features]
    X21 = X21[features]

    # print('-----------------Threshold = 0.1-----------------')
    # # # threshold the ys as .1
    yp1 = [1 if x < .1 else 0 for x in y]
    # p1_params, p1_score = genetic_optimization_classification(X,yp1,downsample=True)
    # # print('\n\n\n\n\n')
    # print('-----------------Threshold = 0.35-----------------')
    yp35 = [1 if x < .35 else 0 for x in y]
    # p35_params, p35_score = genetic_optimization_classification(X,yp35,downsample=True)
    # print('-----------------Threshold = 0.05-----------------')
    yp05 = [1 if x < .05 else 0 for x in y]
    # p05_params, p05_score = genetic_optimization_classification(X,yp05,downsample=True)
    # print('-----------------Threshold = 1.00-----------------')
    y_1 = [1 if x < 1 else 0 for x in y]
    # _1_params, _1_score = genetic_optimization_classification(X,y_1,downsample=True)

    _1_params = {'learning_rate': 0.05686513701078857, 'gamma': 0.3690725162929928, 'n_estimators': 211, 'max_depth': 3, 'max_leaves': 2, 'subsample': 0.20829085379990156, 'booster': 'dart'}
    p35_params = {'learning_rate': 0.021628224103092883, 'gamma': 0.3221929530606452, 'n_estimators': 87, 'max_depth': 3, 'max_leaves': 7, 'subsample': 0.32206709706671505, 'booster': 'dart'}
    p1_params = {'learning_rate': 0.005681916432964979, 'gamma': 0.3422791197286503, 'n_estimators': 209, 'max_depth': 6, 'max_leaves': 6, 'subsample': 0.25470023288701704, 'booster': 'dart'}
    p05_params = {'learning_rate': 0.03456314296918745, 'gamma': 0.9629066305213869, 'n_estimators': 48, 'max_depth': 7, 'max_leaves': 7, 'subsample': 0.37969973950855684, 'booster': 'dart'}

    # train the classifiers
    model_p1, scaler_p1 = train_classifier(X,yp1,normalize=False,downsample=True,params=p1_params)
    model_p35, scaler_p35 = train_classifier(X,yp35,normalize=False,downsample=True,params=p35_params)
    model_p05, scaler_p05 = train_classifier(X,yp05,normalize=False,downsample=True,params=p05_params)
    model_1, scaler_1 = train_classifier(X,y_1,normalize=False,downsample=True,params=_1_params)

    # predict X21
    y21_pred_p1 = model_p1.predict_proba(X21)[:,1]
    y21_pred_p35 = model_p35.predict_proba(X21)[:,1]
    y21_pred_p05 = model_p05.predict_proba(X21)[:,1]
    y21_pred_1 = model_1.predict_proba(X21)[:,1]

    # predict X21
    y21_pred_binary_p1 = model_p1.predict(X21)
    y21_pred_binary_p35 = model_p35.predict(X21)
    y21_pred_binary_p05 = model_p05.predict(X21)
    y21_pred_binary_1 = model_1.predict(X21)

    # get roc curve
    y21_p1 = [1 if x < .1 else 0 for x in y21]
    fpr_p1, tpr_p1, thresholds_p1 = roc_curve(y21_p1, y21_pred_p1)
    auc_p1 = roc_auc_score(y21_p1, y21_pred_p1)

    y21_p35 = [1 if x < .35 else 0 for x in y21]
    fpr_p35, tpr_p35, thresholds_p35 = roc_curve(y21_p35, y21_pred_p35)
    auc_p35 = roc_auc_score(y21_p35, y21_pred_p35)

    y21_p05 = [1 if x < .05 else 0 for x in y21]
    fpr_p05, tpr_p05, thresholds_p05 = roc_curve(y21_p05, y21_pred_p05)
    auc_p05 = roc_auc_score(y21_p05, y21_pred_p05)

    y21_1 = [1 if x < 1 else 0 for x in y21]
    fpr_1, tpr_1, thresholds_1 = roc_curve(y21_1, y21_pred_1)
    auc_1 = roc_auc_score(y21_1, y21_pred_1)

    # plot the roc curves
    plt.plot(fpr_1, tpr_1, label='p < 1.00 (area = %0.2f)' % auc_1)
    plt.plot(fpr_p35, tpr_p35, label='p < 0.35 (area = %0.2f)' % auc_p35)
    plt.plot(fpr_p1, tpr_p1, label='p < 0.10 (area = %0.2f)' % auc_p1)
    plt.plot(fpr_p05, tpr_p05, label='p < 0.05 (area = %0.2f)' % auc_p05)
    plt.plot([0, 1], [0, 1], 'k--')  # random predictions curve
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    # remove top and right border
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # set figure size
    plt.gcf().set_size_inches(5, 5)
    plt.legend(loc="lower right",frameon=False)
    plt.savefig('Figures/trained_19_20_test_2021_roc.png',dpi=300)
    plt.clf()

    # write results to TSV with the following columns:
    # cluster_id, gene, 1.00, 0.35, 0.10, 0.05
    all_coms = {}
    g21 = load_clusters('SubComs/2021/paris.greedy.2021.coms.txt','paris.greedy.2021:')
    # add g21 to all_coms
    all_coms.update(g21)
    # do it for cesna, walktrap and infomap now
    c21 = load_clusters('SubComs/2021/paris.cesna.2021.coms.txt','paris.cesna.2021:')
    all_coms.update(c21)
    w21 = load_clusters('SubComs/2021/paris.walktrap.2021.coms.txt','paris.walktrap.2021:')
    all_coms.update(w21)
    i21 = load_clusters('SubComs/2021/paris.infomap.2021.coms.txt','paris.infomap.2021:')
    all_coms.update(i21)

    print('Writing results to Results/2021_predictions.tsv')
    outfile = open('Results/2021_predictions.tsv','w')
    outfile.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format('cluster_id','gene','p < 1.00','p < 0.35','p < 0.10','p < 0.05'))
    for cid, _1, p35, p1, p05 in zip(X21_cluster_ids, y21_pred_binary_1, y21_pred_binary_p35, y21_pred_binary_p1, y21_pred_binary_p05):
        if cid not in all_coms:
            continue
        for node in all_coms[cid]:
            if 'HP:' in node:
                continue
            outfile.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(cid,node,_1,p35,p1,p05))
        

def regress_rank_tau(X_train,y_train,X_test,y_test,y_test_pred,params,prefix=''):
    model_p1 = xgb.XGBRegressor(**params)
    model_p1.fit(X_train, y_train)
    # get the X21_r samples that where predicted 1 by the classifier
    X21_r_predicted = X_test[[x == 1 for x in y_test_pred]]
    y21_r_predicted = [ y_test[i] for i,x in enumerate(y_test_pred) if x == 1]
    y21_pred = model_p1.predict(X21_r_predicted)
    # create dataframe with cluster_id, real rank, predicted rank
    c_rank = {'cluster_ID':list(range(len(y21_r_predicted))),'emperical_p': list(y21_r_predicted), 'predicted_prob': list(y21_pred)}
    c_rank = pd.DataFrame(c_rank)
    c_rank['emperical_p_rank'] = c_rank['emperical_p'].rank(ascending=True)
    c_rank['predicted_prob_rank'] = c_rank['predicted_prob'].rank(ascending=True)
    # get kendall's tau
    tau, p_value = kendals_tau(c_rank['emperical_p_rank'], c_rank['predicted_prob_rank'])
    print(prefix,'Kendall Tau:',tau)



    # p35_params, p35_score = genetic_optimization_regression(X_p35,y_p35,downsample=True)
    # p05_params, p05_score = genetic_optimization_regression(X_p05,y_p05,downsample=True)
    # _1_params, _1_score = genetic_optimization_regression(X_1,y_1,downsample=True)

# function that trains on 2019 and 2020 and tests on 2021, using the 0.35 threshold and the optimzed parameters, then calculate SHAP values for each feature
def do_shap_analysis():
    # list 2019 files
    files_2019 = ['FinalBOCCFeatures/2019/' + f for f in os.listdir('FinalBOCCFeatures/2019/')]
    files_2020 = ['FinalBOCCFeatures/2020/' + f for f in os.listdir('FinalBOCCFeatures/2020/')]
    files_2021 = ['FinalBOCCFeatures/2021/' + f for f in os.listdir('FinalBOCCFeatures/2021/')]

    # these are the features determined from using just regression, JustRegressionResults/
    features = ['num_sig_go_enrichment_terms', 'num_of_diseases', 'avg_embeddedness', 'conductance', 'normalized_cut', 'triangle_participation_ratio', 'newman_girvan_modularity', 'edges_inside']
    # load files
    X, y = load_files(files_2019)
    X20, y20 = load_files(files_2020)
    X21, y21 = load_files(files_2021)
    X = pd.concat([X,X20])
    y = y + y20
    # sub set the features
    X_og = X.copy()
    X21_og = X21.copy()

    # for each columns in X report the min and max
    for c in X.columns:
        print(c,':',min(X[c]),max(X[c]))

    X = X[features]
    X21 = X21[features]

    yp35 = [1 if x < .35 else 0 for x in y]

    p35_params = {'learning_rate': 0.021628224103092883, 'gamma': 0.3221929530606452, 'n_estimators': 87, 'max_depth': 3, 'max_leaves': 7, 'subsample': 0.32206709706671505, 'booster': 'dart'}
    model_p35, scaler_p35 = train_classifier(X,yp35,normalize=False,downsample=True,params=p35_params)

    # check if shap_values.p exists
    if os.path.exists('shap_values.p'):
        print('Loading shap_values.p and explainer.p')
        shap_values = pickle.load(open('shap_values.p','rb'))
        explainer = pickle.load(open('explainer.p','rb'))
    else:
        print('Calculating SHAP values')
        explainer = shap.Explainer(model_p35.predict, X21)
        shap_values = explainer(X)
        print('Saving shap_values.p and explainer.p')
        # pickle the shap_values
        pickle.dump(shap_values, open('shap_values.p','wb'))
        # picle explainer
        pickle.dump(explainer, open('explainer.p','wb'))
    
    shap.plots.beeswarm(shap_values)
    shap.plots.bar(shap_values)



    
    
    


    

def load_clusters(filename,prefix):
    com_dict = {}
    for line in open(filename,'r'):
        row = line.strip().split('\t')
        c_name = prefix + row[0]
        com_dict[c_name] = row[1:]
    return com_dict
    
    


def make_classification_results():
    # load files 2019-2022
    files_2019 = ['FinalBOCCFeatures/2019/' + f for f in os.listdir('FinalBOCCFeatures/2019/')]
    files_2020 = ['FinalBOCCFeatures/2020/' + f for f in os.listdir('FinalBOCCFeatures/2020/')]
    files_2021 = ['FinalBOCCFeatures/2021/' + f for f in os.listdir('FinalBOCCFeatures/2021/')]
    # files_2022 = ['FinalBOCCFeatures/2022/' + f for f in os.listdir('FinalBOCCFeatures/2022/')]
    # load files 2019-2022
    X19, y19 = load_files(files_2019)
    X20, y20 = load_files(files_2020)
    X21, y21 = load_files(files_2021)
    # create list of names for the samples in X19 based on algo and the row index
    X19_cluster_ids = ['{}.{}:{}'.format(a,'2019',str(i)) for a,i in zip(X19['algo'],X19['cluster_id'])]
    X20_cluster_ids = ['{}.{}:{}'.format(a,'2020',str(i)) for a,i in zip(X20['algo'],X20['cluster_id'])]
    X21_cluster_ids = ['{}.{}:{}'.format(a,'2021',str(i)) for a,i in zip(X21['algo'],X21['cluster_id'])]
    
    # X22, y22 = load_files(files_2022)
    features = ['num_sig_go_enrichment_terms', 'num_of_diseases', 'avg_embeddedness', 'conductance', 'normalized_cut', 'triangle_participation_ratio', 'newman_girvan_modularity', 'edges_inside']
    # sub set the features
    X19 = X19[features]
    X20 = X20[features]
    X21 = X21[features]
    thresholds = [1.00,.35,.1]
    results_2021 = {'cluster_id':[],}
    classifiers = {}
    for t in thresholds:
        # threshold y for 2019
        y19_c = [1 if p < t else 0 for p in y19]
        # train the classifier
        model, scaler = train_classifier(X19,y19_c,normalize=False,downsample=True,params=classification_params)
        classifiers[t] = model
    m1 = classifiers[1.00]
    results_2021['p < 1.00'] = m1.predict(X21)
    m35 = classifiers[.35]
    results_2021['p < 0.35'] = m35.predict(X21)
    mp1 = classifiers[.1]
    results_2021['p < 0.10'] = mp1.predict(X21)
    results_2021['cluster_id'] = X21_cluster_ids

    # roc curve of each threshold
    y21_1 = [1 if y < 1 else 0 for y in y21]
    y21_35 = [1 if y < .35 else 0 for y in y21]
    y21_01 = [1 if y < .1 else 0 for y in y21]
    fpr_1, tpr_1, thresholds_1 = roc_curve(y21_1, m1.predict_proba(X21)[:,1])
    auc_1 = roc_auc_score(y21_1, results_2021['p < 1.00'])
    fpr_35, tpr_35, thresholds_35 = roc_curve(y21_35, m35.predict_proba(X21)[:,1])
    auc_35 = roc_auc_score(y21_35, results_2021['p < 0.35'])
    fpr_10, tpr_10, thresholds_10 = roc_curve(y21_01, mp1.predict_proba(X21)[:,1])
    auc_10 = roc_auc_score(y21_01, results_2021['p < 0.10'])
    # plot the roc curves
    plt.plot(fpr_1, tpr_1, label='p < 1.00 (area = %0.2f)' % auc_1)
    plt.plot(fpr_35, tpr_35, label='p < 0.35 (area = %0.2f)' % auc_35)
    plt.plot(fpr_10, tpr_10, label='p < 0.10 (area = %0.2f)' % auc_10)
    plt.plot([0, 1], [0, 1], 'k--')  # random predictions curve
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    # remove top and right border
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # set figure size
    plt.gcf().set_size_inches(5, 5)
    plt.legend(loc="lower right")
    plt.savefig('Figures/2021_classifier_roc.png',dpi=300)
    plt.clf()

    # to tsv
    results_2021 = pd.DataFrame(results_2021)
    results_2021.to_csv('Results/2021.classification.tsv',sep='\t',index=False)

    all_coms = {}
    g21 = load_clusters('SubComs/2021/paris.greedy.2021.coms.txt','paris.greedy.2021:')
    # add g21 to all_coms
    all_coms.update(g21)
    # do it for cesna, walktrap and infomap now
    c21 = load_clusters('SubComs/2021/paris.cesna.2021.coms.txt','paris.cesna.2021:')
    all_coms.update(c21)
    w21 = load_clusters('SubComs/2021/paris.walktrap.2021.coms.txt','paris.walktrap.2021:')
    all_coms.update(w21)
    i21 = load_clusters('SubComs/2021/paris.infomap.2021.coms.txt','paris.infomap.2021:')
    all_coms.update(i21)
    # create a new dictionary with the colulmns, cluster_id, gene p < 1.00, p < 0.35, p < 0.10
    results_2021_gene_wise = {'cluster_id':[],'gene':[],'p < 1.00':[],'p < 0.35':[],'p < 0.10':[],'p < 0.05':[]}
    missing_cluster = set()
    for i,row in results_2021.iterrows():
        cluster_id = row['cluster_id']
        if cluster_id not in all_coms:
            missing_cluster.add(cluster_id)
            continue
        for gene in all_coms[cluster_id]:
            if 'HP:' in gene:
                continue
            results_2021_gene_wise['cluster_id'].append(cluster_id)
            results_2021_gene_wise['gene'].append(gene)
            results_2021_gene_wise['p < 1.00'].append(row['p < 1.00'])
            results_2021_gene_wise['p < 0.35'].append(row['p < 0.35'])
            results_2021_gene_wise['p < 0.10'].append(row['p < 0.10'])
            results_2021_gene_wise['p < 0.05'].append(row['p < 0.05'])
    results_2021_gene_wise = pd.DataFrame(results_2021_gene_wise)
    results_2021_gene_wise.to_csv('Results/2021.classification.gene_wise.tsv',sep='\t',index=False)
    print('Missing Clusters:',missing_cluster)
    print('Missing Clusters Count:',len(missing_cluster))
    # print number of p < 1.00, p < 0.35, p < 0.10
    print('gene-wise p < 1.00',len([x for x in results_2021_gene_wise['p < 1.00'] if x == 1]))
    print('gene-wise p < 0.35',len([x for x in results_2021_gene_wise['p < 0.35'] if x == 1]))
    print('gene-wise p < 0.10',len([x for x in results_2021_gene_wise['p < 0.10'] if x == 1]))
    # print the same but from results_2021
    print('p < 1.00',len([x for x in results_2021['p < 1.00'] if x == 1]))
    print('p < 0.35',len([x for x in results_2021['p < 0.35'] if x == 1]))
    print('p < 0.10',len([x for x in results_2021['p < 0.10'] if x == 1]))
    # print the same for results_2021 but as a %
    print('p < 1.00',len([x for x in results_2021['p < 1.00'] if x == 1]) / len(results_2021['p < 1.00']))
    print('p < 0.35',len([x for x in results_2021['p < 0.35'] if x == 1]) / len(results_2021['p < 0.35']))
    print('p < 0.10',len([x for x in results_2021['p < 0.10'] if x == 1]) / len(results_2021['p < 0.10']))




if __name__ == '__main__':
    # classifier_roc()
    # start_grid_search()
    # do_19v20_regression_only_roc()
    # start_grid_search()
    # do_19v20_practical_application()
    # train_regressor_with_classifier()
    # main_lof()
    # Do feature selection for if we use ONLY regression from the beginning
    # pick_just_regression_features()
    # main()
    # print('-----------------Classification-----------------')
    # start_genetic_algo_search()
    # print('-----------------Regression-----------------')
    # start_genetic_algo_search_regression()
    # do_all_filter_and_flip()
    # start_genetic_algo_for_just_regression()
    # just_regression()
    # threshold_rocs()
    # do_genetic_optimization_for_p1_and_p35()
    # start_genetic_algo_for_p35_regression()
    # do_19v20_practical_application_p35()
    # make_classification_results()

    do_genetic_optimization_for_p1_and_p35_with_2019_and_2020()
    # do_shap_analysis()
    


# python Scripts/train_model.py