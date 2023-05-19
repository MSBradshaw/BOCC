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
from sklearn.metrics import roc_curve
from sklearn.model_selection import GridSearchCV
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
import shap
from sklearn.metrics import precision_score, recall_score, f1_score, roc_auc_score, accuracy_score
from sklearn.metrics import confusion_matrix

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
# classification params
classification_params = {'learning_rate': 0.07305940290241641, 
                         'gamma': 0.3795158439477443, 
                         'n_estimators': 98, 
                         'max_depth': 13, 
                         'max_leaves': 4, 
                         'subsample': 0.05772500174333884, 
                         'booster': 'dart'}

def load_data(file_path):
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

def genetic_optimization_classification(X, y, downsample=True):
    # set seed
    random.seed(42)
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

def select_features(X,y,n,downsample=True,normaize=True):
    # set seed
    random.seed(42)
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

# function to remove all columns whose name contains the substring 'plof'
def drop_plof(X):
    cols = [c for c in X.columns if 'plof' not in c]
    return X[cols]

def feature_selection():
    files_2019 = ['FinalBOCCFeatures/2019/' + f for f in os.listdir('FinalBOCCFeatures/2019/')]
    X19, y19 = load_files(files_2019)

    print(X19.shape)
    print('num p < .05',len([i for i in y19 if i < .05]))
    plot_feature_correlation(X19,'Figures/feature_correlation.png')
    
    assert(X19.shape[0] == len(y19))
    find_nan(X19)
    # remove the columns cluster_id
    X19 = X19.drop('cluster_id', axis=1)
    # plot_pca(X19,y19)
    X19 = drop_algo(X19)
    X19 = drop_plof(X19)
    
    # # Select All Features:
    select_all_features(X19,y19)
    find_best_features()

def get_tpr_fpr_auc(train_files,test_files,features,plot_prefix,params=None,downsample=True,threshold=1,plot=True):
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

    return fpr, tpr, auc

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
        fpr, tpr, auc = get_tpr_fpr_auc(files_2019,
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

def optimize_train_test_report():
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

    # subset the features
    X = X[features]
    X21 = X21[features]

    # print('-----------------Threshold = 0.1-----------------')
    # # # threshold the ys as .1
    yp1 = [1 if x < .1 else 0 for x in y]
    p1_params, p1_score = genetic_optimization_classification(X,yp1,downsample=True)
    # # print('\n\n\n\n\n')
    # print('-----------------Threshold = 0.35-----------------')
    yp35 = [1 if x < .35 else 0 for x in y]
    p35_params, p35_score = genetic_optimization_classification(X,yp35,downsample=True)
    # print('-----------------Threshold = 0.05-----------------')
    yp05 = [1 if x < .05 else 0 for x in y]
    p05_params, p05_score = genetic_optimization_classification(X,yp05,downsample=True)
    # print('-----------------Threshold = 1.00-----------------')
    y_1 = [1 if x < 1 else 0 for x in y]
    _1_params, _1_score = genetic_optimization_classification(X,y_1,downsample=True)

    # there are the hard coded best params so the genetic algorithm doesn't have to be run every time
    # _1_params = {'learning_rate': 0.05686513701078857, 'gamma': 0.3690725162929928, 'n_estimators': 211, 'max_depth': 3, 'max_leaves': 2, 'subsample': 0.20829085379990156, 'booster': 'dart'}
    # p35_params = {'learning_rate': 0.021628224103092883, 'gamma': 0.3221929530606452, 'n_estimators': 87, 'max_depth': 3, 'max_leaves': 7, 'subsample': 0.32206709706671505, 'booster': 'dart'}
    # p1_params = {'learning_rate': 0.005681916432964979, 'gamma': 0.3422791197286503, 'n_estimators': 209, 'max_depth': 6, 'max_leaves': 6, 'subsample': 0.25470023288701704, 'booster': 'dart'}
    # p05_params = {'learning_rate': 0.03456314296918745, 'gamma': 0.9629066305213869, 'n_estimators': 48, 'max_depth': 7, 'max_leaves': 7, 'subsample': 0.37969973950855684, 'booster': 'dart'}

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
    outfile.close()

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

if __name__ == '__main__':
    feature_selection()
    threshold_rocs()
    optimize_train_test_report()
    do_shap_analysis()

# python Scripts/train_model.py
