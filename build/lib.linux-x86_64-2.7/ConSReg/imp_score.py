'''
Functions for ranking features by stability selection
'''

# Author: (Alex) Qi Song <alexsong@vt.edu>

from joblib import Parallel, delayed
from sklearn.model_selection import KFold
from train_test import train_and_val

from rpy2.robjects.packages import importr # Function for loading R library
from rpy2.robjects import numpy2ri
from utils import get_max_abs

import numpy as np
import pandas as pd
import warnings

# Load R libraries
BASE = importr("base")
GGL = importr("gglasso")

# Activate python<->R objects conversion
numpy2ri.activate()

# Suppress scikit learn future warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

'''
Tune the best lambda using 10 fold cross validation
'''
def get_best_lambda(features, labels, fold = 10, verbose = 0, n_jobs = 1):
    '''
    Parameters
    ----------
    features : numpy 2d array. feature matrix
    labels: numpy 1d array. labes
    fold : int. for cross-validation
    verbose : int. verbose level of Parallel()
    n_jobs : int. Number of jobs
    
    Returns
    -------
    best_lambda : float. the best lambda
    '''
    
    k_fold = KFold(n_splits = fold, shuffle = True)
    res = Parallel(n_jobs = n_jobs, verbose = verbose)(delayed(train_and_val)(
                                                                features_train = features[train_index,],
                                                                labels_train = labels[train_index,],
                                                                features_val = features[test_index,],
                                                                labels_val = labels[test_index,],
                                                                ml_engine = "lrlasso"
                                                            )
                                       for train_index, test_index in k_fold.split(features))
    best_lambdas = map(lambda x: x[1], res)
    best_lambda = np.mean(best_lambdas)
    return(best_lambda)

'''
def get_coef(features, labels, reg_lambda, return_model = False):
    
    # Train model using best lambda and get selected features
    args = {'x':features, 'y':labels, 'loss':'logit', 'lambda':reg_lambda}
    model = GGL.gglasso(**args)

    coef = np.array(model[model.names.index("beta")])
    
    if return_model:
        return(coef, model)
    else:
        return(coef)
'''

'''
Train a randomized LRLASSO model and get model coefficients
'''
def run_rand_LRLASSO(features, labels, reg_lambda, sample_fraction):
    '''
    Parameters
    ----------
    features : numpy 2d array. feature matrix
    labels: numpy 1d array. labes
    reg_lambda : float. lambda value for regularization
    sample_fraction : float. what fraction of sample 
    should be used for stability selection?
    
    Returns
    -------
    coef: numpy 1d array. binary vector of non-zero-weight features
    '''
    
    # Re-initialize random seeds
    np.random.seed(None)

    # Randomly sample the training examples
    index = np.random.choice(features.shape[0],int(features.shape[0]*sample_fraction),replace = False)
    features_now = features[index,]
    labels_now = labels[index]
    
    # Min-max transformation
    features_now = features_now/get_max_abs(np.array(features_now))
    
    # Randomly scale the features
    scaler = np.random.uniform(low=1e-4,high=1,size = features_now.shape[1])
    features_now = np.multiply(features_now,scaler)

    # Convert to R objects
    features_now = BASE.as_matrix(features[index,])
    labels_now = BASE.as_matrix(labels[index,])

    # Train model and get coefficients
    args = {'x':features_now, 'y':labels_now, 'loss':'logit', 'lambda':reg_lambda}
    model = GGL.gglasso(**args)

    coef = np.array(model[model.names.index("beta")])
    non_zero_coef = (coef != 0)
    
    return(non_zero_coef)

'''
LRLASSO importance score by stability selection
'''
def LRLASSO_imp_score(features, labels, reg_lambda, n_resampling = 200, verbose = 0, sample_fraction = 0.5, n_jobs = 1):
    '''
    Parameters
    ----------
    features : numpy 2d array. feature matrix
    labels: numpy 1d array. labes
    reg_lambda : float. lambda value for regularization
    n_resampling : int. number of resampling times.
    verbose : int. verbose level of Parallel()
    sample_fraction : float. what fraction of sample 
    should be used for stability selection?
    
    n_jobs : int. Number of jobs
    
    Returns
    -------
    scores: numpy 1d array. importance scores
    '''
    
    result = Parallel(n_jobs = n_jobs, verbose = verbose)(delayed(run_rand_LRLASSO)(
                                                                            features = features,
                                                                            labels = labels, 
                                                                            reg_lambda = reg_lambda,
                                                                            sample_fraction = sample_fraction)
                                                for i in range(n_resampling))
    
    # Compute the ranking scores from stability selection result
    scores = np.sum(np.concatenate(result,axis=1),axis = 1, keepdims = True)
    scores = scores / float(n_resampling)
    
    return(scores)

'''
Find best lambda and use this best lamda to perform stability selection
'''
def imp_score(feature_mat, n_resampling=200, verbose = 0, fold = 10, n_jobs = 1):
    '''
    Parameters
    ----------
    feature_mat : pandas dataframe. feature matrix containing the labels
    n_resampling : int. number of resampling times.
    verbose : int. verbose level of Parallel()
    fold : int. for cross-validation
    n_jobs : int. Number of jobs
    
    Returns
    -------
    scores: numpy 1d array. importance scores
    '''
    
    # Convert to numpy array
    features = np.array(feature_mat.iloc[:,1:])
    labels = np.array(feature_mat.iloc[:,0])
    
    # feature names
    feature_names = feature_mat.iloc[:, 1:].columns
    
    # gglasso requires that negative label should be -1
    labels[labels == 0] = -1
    
    # Get best lambda by cv
    best_lambda = get_best_lambda(features = features, labels = labels, verbose = verbose, fold = fold, n_jobs = n_jobs)

    '''
    # Train model with best lambda and get coefficients
    if return_model:
        coef, model = get_coef(features = features, return_model = return_model, labels = labels, reg_lambda = best_lambda)
    else:
        coef = get_coef(features = features, return_model = return_model, labels = labels, reg_lambda = best_lambda)
    '''
    
    # compute importance scores from LRLASSO
    args = {
                'features':features,
                'labels':labels,
                'reg_lambda':best_lambda, 
                'n_resampling':n_resampling,
                'verbose':verbose,
                'sample_fraction':0.5,
                'n_jobs':n_jobs
            }
    feature_scores = LRLASSO_imp_score(**args)
    
    # format result
    # ranks = np.concatenate([coef, feature_scores],axis=1)
    ranks = pd.DataFrame(feature_scores, columns = ["scores"], index = feature_names)
    ranks = ranks.sort_values(by = "scores",ascending = False)

    return(ranks)
    

