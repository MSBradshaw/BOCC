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
from sklearn.decomposition import PCA
import sklearn
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import random
from sklearn.preprocessing import RobustScaler
from sklearn.feature_selection import SequentialFeatureSelector

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


def train_model(X, y):
    # create an xgboost regressor
    model = xgb.XGBRegressor(n_jobs=-1)
    # lets do grid search hyperparameter optimization
    # print('Starting grid search')
    # gs = GridSearchCV(model, HYPERPARAMS, cv=10, verbose=1, scoring='neg_mean_squared_error')
    print('Fitting the model')
    # fit the model
    # model.fit(X,y, scoring='neg_mean_squared_error')
    # cross validate the model
    scores = cross_val_score(model, X, y, cv=10, scoring='neg_mean_squared_error')
    print(scores)

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
    return mse

    

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

def select_all_features_regression(X,y):
    selection_results = {"i":[], "MSE":[], "features":[], "normalized":[], 'downsampled':[]}
    for norm in [True,False]:
        for downsample in [True,False]:
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
    selection_results.to_csv('RegressionResults/selection_results.csv',index=False)

def find_best_features_regression():
    df = pd.read_csv('RegressionResults/selection_results.csv')
    # sort my MSE
    df = df.sort_values(by='MSE',ascending=True)
    # print the top 10
    print(df.head(5))
    # print all the information in the top 5
    for i in range(5):
        print(df.iloc[i])
    # save the ranked and sorted results
    df.to_csv('RegressionResults/selection_results_ranked.csv',index=False)

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

def main():
    #'FinalBOCCFeatures/2019/paris.infomap.2019.bocc_res.tsv'
    files_2019 = ['FinalBOCCFeatures/2019/' + f for f in os.listdir('FinalBOCCFeatures/2019/')]
    X19, y19 = load_files(files_2019)

    print(X19.shape)
    print('num p < .05',len([i for i in y19 if i < .05]))
    plot_feature_correlation(X19,'Figures/feature_correlation.png')
    
    assert(X19.shape[0] == len(y19))
    find_nan(X19)
    plot_pca(X19,y19)
    X19 = drop_algo(X19)
    
   
    # # Select All Features:
    # select_all_features(X19,y19)
    # find_best_features()

    # Regression
    # train_model_regressor(X19,y19,normaize=True,downsample=True)
    select_all_features_regression(X19,y19)
    print()
    find_best_features_regression()
    list1 = ['cluster_size', 'num_sig_go_enrichment_terms', 'num_of_diseases', 'avg_internal_degree', 'normalized_cut', 'triangle_participation_ratio', 'newman_girvan_modularity', 'edges_inside', 'hub_dominance', 'sum_plof']
    list2 = ['max_norm_cell_type_specificity', 'num_of_diseases', 'max_norm_disease_specificity', 'avg_embeddedness', 'conductance', 'cut_ratio', 'normalized_cut', 'expansion', 'triangle_participation_ratio', 'internal_edge_density', 'hub_dominance', 'mean_plof', 'median_plof', 'sum_plof']
    # plot_2sets_of_feature_correlation(X19,list1,list2)


    


if __name__ == '__main__':
    main()


# python Scripts/train_model.py