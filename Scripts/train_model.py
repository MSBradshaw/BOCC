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
        'max_leaves': list(range(1,5)),
        'subsample': [ x/100 for x in range(1,101,5)], # .01 - 1 by .05
        'booster' : ['dart']
        }

print(HYPERPARAMS)

def kendals_w(ranks1,ranks2):
    assert len(ranks1) == len(ranks2)
    # source https://www.youtube.com/watch?v=zVufp7cJ8S4
    fstat = stats.friedmanchisquare(ranks1,ranks2)
    q = fstat[0]
    n = len(ranks1)
    w = q / (n * (k-1))
    return w

def load_data(file_path):
    data = pd.read_csv(file_path,sep='\t')
    data['algo'] = '.'.join(file_path.split('/')[-1].split('.')[0:2])
    # remove the column cluster_id
    to_drop = ['significance','cluster_id','sig_go_enrichment_terms','go_sig_threshold','max_norm_cell_type_comma_sep_string','num_new_edges_on_any_node','sig_go_enrichment_p_vals','mg2_portion_families_recovered','mg2_not_pairs_count','mg2_pairs_count','max_norm_disease_comma_sep_string','sig_go_enrichment_fdr_corrected_p_vals']
    for name in to_drop:
        data = data.drop(name,axis=1)    
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
        majority_indices = [i for i,x in enumerate(y) if x == 1]
        minority_indices = [i for i,x in enumerate(y) if x != 1]
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
    with open('best_model.pkl','wb') as f:
        pickle.dump(gs.best_estimator_,f)

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
        for downsample in [True,False]:
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
def classifier_train_and_predict_all(X,y,normalize=False,downsample=False):
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
    X_train, X_test, y_train, y_test = train_test_split(X,y,test_size=0.2,random_state=42)

    # normalize the data
    scaler = None
    if normalize:
        scaler = RobustScaler()
        scaler.fit(X_train)
        X_train = scaler.transform(X_train)
        X_test = scaler.transform(X_test)
    
    # train the model
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
    return y_pred_all, model, scaler

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
def train_regressor(X,y,normalize=False,downsample=False):
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
    model = xgb.XGBRegressor()
    model.fit(X,y)

    return model, scaler

def train_classifier(X,y,normalize=False,downsample=False):
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
    model = xgb.XGBClassifier()
    model.fit(X,y)

    return model, scaler

def practical_application(train_paths,test_paths,features_c,features_r,c_norm,c_downsample,r_norm,r_downsample):
    # load the data
    X_train, y_train = load_files(train_paths)
    X_test, y_test = load_files(test_paths)
    
    # subset the data
    X_train_c = subset_data(X_train.copy(),features_c)
    X_test_c = subset_data(X_train.copy(),features_c)
    
    # train the classifier
    c_y_pred, c_model, c_scaler = classifier_train_and_predict_all(X_train_c,y_train,normalize=c_norm,downsample=c_downsample)

    # subset train based on predictions
    X_train_r = X_train.iloc[c_y_pred == 1]
    X_train_r = subset_data(X_train_r.copy(),features_r)
    y_train = [y_train[i] for i in range(len(y_train)) if c_y_pred[i] == 1]
    # train the regressor
    r_model, r_scaler = train_regressor(X_train_r,y_train,normalize=r_norm,downsample=r_downsample)

    # Predict the test data
    if c_scaler is not None:
        X_test_c = c_scaler.transform(X_test)
    else:
        X_test_c = X_test.copy()
    X_test_c = subset_data(X_test_c.copy(),features_c)
    c_y_pred_test = c_model.predict(X_test_c)

    # print the class break down of c_y_pred_test
    print('Class Breakdown of c_y_pred_test')
    print(np.bincount(c_y_pred_test.astype(int)))

    # subset c_y_pred_test based on the predictions
    # X_test_filtered = X_test.iloc[c_y_pred_test == 1]
    # y_test_filtered = [y_test[i] for i in range(len(y_test)) if c_y_pred_test[i] == 1]
    
    # predict the test data
    if r_scaler is not None:
        X_test_r = r_scaler.transform(X_test)
    else:
        X_test_r = X_test.copy()

    X_test_r = subset_data(X_test_r.copy(),features_r)
    r_y_pred_test = r_model.predict(X_test_r)
    
    # put r_y_pred_test and c_y_pred_test back together and return
    results = []
    for i,x in enumerate(c_y_pred_test):
        if x == 0:
            results.append(np.inf)
        else:
            print(r_y_pred_test[i])
            results.append(r_y_pred_test[i])
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
    print(results_rank)
    # assign a rank to each result allow for ties
    results_rank['rank'] = results_rank['result'].rank(method='average')
    print('ranked')
    print(results_rank)
    # sort the df in the original order
    results_rank = results_rank.sort_values(by='index',ascending=True)
    print('original order')
    print(results_rank)
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
    print(tp_fp)

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
                          r_downsample=True)
    
    # load the 2020 files and get the p values
    X20, y20 = load_files(files_2020)
    # plot the ranks vs their actual p values
    # plot_ranks_vs_pvalues(ranks_20,y20,'Figures/ranks_vs_pvalues.19v20.png')
    # create a df of ranks_20 with clusterIDs
    p_df = pd.DataFrame({'clusterID':list(range(len(y20))),'emprical_p':y20})
    ranks_df = pd.DataFrame({'clusterID':list(range(len(y20))),'rank':ranks_20})
    df = pd.DataFrame({'clusterID':list(range(len(y20))),'predicted_p_rank':ranks_20,'emprical_p':y20})
    # create a rank columns based on emprical_p
    df['rank'] = df['emprical_p'].rank(method='min')
    # sort by clusterID
    df.sort_values(by=['clusterID'],inplace=True)
    # write p_df to file
    p_df.to_csv('.pvalues.19v20.tsv',index=False,sep='\t')
    # write ranks_df to file
    ranks_df.to_csv('.ranks.19v20.tsv',index=False,sep='\t')
    top_x_roc('.ranks.19v20.tsv','.pvalues.19v20.tsv','Figures/top_x_roc.19v20.png')
    # calc kendall's W
    # w = kendals_w(list(df['rank']),list(df['predicted_p_rank']))
    # print('Classification & Regression', w)

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


def classifier_roc():
    features = ['gene_ratio', 'HPO_ratio', 'num_sig_go_enrichment_terms', 'num_of_diseases', 'max_norm_disease_specificity', 'cut_ratio', 'expansion', 'newman_girvan_modularity', 'edges_inside']
    # load the 2019 files
    files_2019 = ['FinalBOCCFeatures/2019/' + f for f in os.listdir('FinalBOCCFeatures/2019/')]
    X19, y19 = load_files(files_2019)
    # sub set the features
    X19 = X19[features]
    # load the 2020 files
    files_2020 = ['FinalBOCCFeatures/2020/' + f for f in os.listdir('FinalBOCCFeatures/2020/')]
    X20, y20 = load_files(files_2020)
    og_y20 = y20.copy()
    # sub set the features
    X20 = X20[features]

    threshold = 1
    y19 = [1 if p < threshold else 0 for p in y19]
    y20 = [1 if p < threshold else 0 for p in y20]
    # print the number p == 1 in y19
    print('y19',sum(y19),len(y19))
    print('y20',sum(y20),len(y20))

    # train the classifier
    model, scaler = train_classifier(X19,y19,normalize=False,downsample=True)
    # predict X20
    if scaler is not None:
        X20 = scaler.transform(X20)
        X20 = pd.DataFrame(X20,columns=features)
    
    y20_pred = model.predict(X20)
    y20_probs = model.predict_proba(X20)
    print(y20_probs)
    y20_probs = y20_probs[:,1]
    print(y20_probs)

    print(len(y20_pred))
    fpr, tpr, thresh = roc_curve(y20,y20_probs)

    print('fpr',fpr)
    print('tpr',tpr)

    # plot the roc curve
    plt.plot(fpr,tpr)
    # plot a diagonal line
    plt.plot([0,1],[0,1],'k--',c='red')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    # remove top and right spines
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.title('19,20 ROC Curve threshold = {}'.format(threshold))
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



if __name__ == '__main__':
    # classifier_roc()
    start_grid_search()
    # do_19v20_regression_only_roc()
    # start_grid_search()
    # do_19v20_practical_application()
    # train_regressor_with_classifier()
    # main_lof()
    # main()
    


# python Scripts/train_model.py