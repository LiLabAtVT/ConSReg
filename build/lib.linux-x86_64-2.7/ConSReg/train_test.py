'''
Functions for cross-validating the feature matrix using
different machine learning engine
'''

# Author: (Alex) Qi Song <alexsong@vt.edu>

from sklearn.svm import LinearSVC,l1_min_c
from sklearn.linear_model import SGDClassifier,LogisticRegression
from sklearn.calibration import CalibratedClassifierCV
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score,precision_recall_curve,auc,confusion_matrix,recall_score
from scipy.stats import pearsonr

from rpy2.robjects.packages import importr # This is for loading R library
from rpy2.robjects import numpy2ri
from rpy2.robjects import FloatVector
from utils import get_max_abs

import warnings
import numpy as np

# Suppress scikit learn future warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# Load R libraries
BASE = importr("base")
GGL = importr("gglasso")
RRF = importr("RRF")

# Activate python<->R objects conversion
numpy2ri.activate()

'''
Train using training feature matrix and get best parameter from
validation feature matrix 
'''
def train_and_val(features_train,labels_train,features_val,labels_val, ml_engine, feature_group = None):
    '''
    Parameters
    ----------
    features_train : numpy 2d array. training feature matrix
    labels_train : numpy 1d array. training labels. 0 for negative and 1 for positive
    features_val : numpy 2d array. validation feature matrix
    labels_val : numpy 1d array. validation labels. 0 for negative and 1 for positive
    ml_engine: string. machine learning engine to use. Possible choice : 1) 'lrlasso': logistic lasso
                                                                          2) 'lglasso' logistic group lasso
                                                                          3) 'lgen': logistic elastic net
                                                                          4) 'grrf' : guided regularized random forest
                                                                          5) 'lsvm': linear support vector machine
                                                                          6) 'lgpcc': logistic regression + pearson correlation coefficient
    
    Returns
    -------
    best_model : the best model selected from validation
    best_param : the best hyper-parameters from validation
    feature_num : number of selected features from best model
    '''
    
    # Train model and get predicted score
    if ml_engine == "lrlasso" or ml_engine == "lglasso":

        # Convert to R objects
        features_train = BASE.as_matrix(features_train)
        features_val = BASE.as_matrix(features_val)

        # gglasso requires negative class to have label '-1'. Replace 0 with -1
        # For convenience of comparison, labels_val will not be converted to R object
        labels_train[labels_train == 0] = -1
        labels_val[labels_val == 0] = -1
        labels_train = BASE.as_vector(labels_train)
        labels_val = labels_val.reshape(labels_val.shape[0],1)

        # To avoid name conflict with python keyword 'lambda', use dictionary to pass function argument
        if ml_engine == "lrlasso":
            args = {'x':features_train, 'y':labels_train, 'loss':'logit', 'lambda.factor':0.01}
        elif ml_engine == "lglasso":
            args = {'x':features_train, 'y':labels_train, 'group':BASE.as_vector(feature_group), 'loss':'logit', 'lambda.factor':0.01}

        # Train model on training data set
        best_model = GGL.gglasso(**args)

        # Predict on validation data set
        pred = GGL.predict_gglasso(best_model,type = 'class',newx = features_val)
        pred = np.array(pred)
        
        # Get sequence of lambdas
        lambda_seq = np.array(best_model[best_model.names.index('lambda')])
        
        # Get the lambda which gives highest accuracy on validation dataset
        best_idx = np.argmax(np.sum(pred == labels_val,axis = 0))
        best_param = lambda_seq[best_idx]
        
        # Get number of selected features from the best model
        coef = np.array(best_model[best_model.names.index("beta")])[:,best_idx]
        feature_num = coef[coef != 0].shape[0]
        
    elif ml_engine == "lgen":
        
        # Generate sequence of two parameters: alpha and l1_ratio
        alpha_list = np.logspace(-3,3,5)
        l1_ratio_list = np.linspace(0,1,5)
        
        best_acc = 0
        best_param = None
        best_model = None
        
        # Test the different alpha and l1_ratio on validation dataset
        for alpha in alpha_list:
            for l1_ratio in l1_ratio_list:
                lgen = SGDClassifier(loss = 'log',penalty = 'elasticnet', alpha = alpha, l1_ratio = l1_ratio).fit(features_train, labels_train)
                pred = lgen.predict(features_val)
                acc = np.sum(labels_val == pred)
                if acc > best_acc:
                        best_acc = acc
                        best_param = (alpha,l1_ratio)
                        best_model = lgen
                        
        # Get number of selected features from the best model
        feature_num = best_model.coef_[best_model.coef_ != 0].shape[0]
        
    elif ml_engine == "grrf":
        
        # Convert to R objects
        features_train = BASE.as_matrix(features_train)
        features_val = BASE.as_matrix(features_val)
        labels_train = BASE.as_vector(labels_train)
        
        rf = RRF.RRF(features_train,BASE.as_factor(labels_train), flagReg = 0, ntree = 100) # build an ordinary RF 

        # Get importance score
        imp = rf[rf.names.index("importance")] 
        imp = np.array(imp)
        imp = imp/(np.max(imp))
        best_acc = 0
        best_param = None
        best_model = None
        
        # Test different gamma on validation dataset
        for gamma in (0,0.5,1):
            coefReg = (1-gamma) + gamma*imp
            coefReg = FloatVector(coefReg)
            grrf = RRF.RRF(features_train,BASE.as_factor(labels_train), flagReg=1, coefReg=coefReg, ntree = 100)
            pred = np.array(RRF.predict_RRF(grrf, features_val)) - 1
            acc = np.sum(labels_val == pred)
            if acc > best_acc:
                best_acc = acc
                best_param = gamma
                best_model = grrf
        
        # Get number of selected features from the best model
        feature_num = np.array(best_model[best_model.names.index("feaSet")]).shape[0]
    
    elif ml_engine == "lsvm":
        '''
        Calculate the lower bound of C for a null model
        If C goes smaller than this value, the model would 
        end up selecting no features
        '''
        min_c = l1_min_c(features_train, labels_train)

        # log spaced list of C parameters
        c_list = np.logspace(np.log10(min_c),3,10)
        best_acc = 0
        best_param = None
        best_model = None

        # Train on training dataset, validate each C on validation dataset
        for C in c_list:
            svm = LinearSVC(C = C,penalty = 'l1',dual=False).fit(features_train,labels_train)
            pred = svm.predict(features_val)
            acc = np.sum(labels_val == pred)
            if acc > best_acc:
                best_acc = acc
                best_param = C
                best_model = svm
        
        # Get number of selected features from the best model
        feature_num = best_model.coef_[best_model.coef_ != 0].shape[0]
        
    elif ml_engine == "lgpcc":
        
        cor_coef = []
        for i in range(features_train.shape[1]):
            x = features_train[:,i].astype(np.float32)
            y = labels_train.astype(np.float32)
            cor_coef.append(pearsonr(x,y)[0])
        
        cor_coef = np.array(cor_coef)
        order = np.argsort(np.abs(cor_coef))[-1:0:-1]
        to_include = np.arange(50,order.shape[0],50)
        best_acc = 0
        best_param = None
        best_model = None
        
        for i in to_include:
            lg = LogisticRegression().fit(features_train[:,order[:i]],labels_train)
            pred = lg.predict(features_val[:,order[:i]])
            acc = np.sum(pred == labels_val)
            if acc > best_acc:
                best_acc = acc
                best_param = order[:i]
                best_model = lg
                
        # Get number of selected features from the best model
        feature_num = best_model.coef_[best_model.coef_ != 0].shape[0]
        
    return(best_model,best_param,feature_num)

'''
This function first splits feature matrix into training matrix, validation matrix
and test matrix. Then it calls train_and_val() to get the best model from training matrix
and validation matrix. This best model will be tested by test matrix.
'''
def test_by_cv(feature_mat, diff_tab, ml_engine, feature_group):
    '''
    Parameters
    ----------
    feature_mat : pandas dataframe. feature matrix
    diff_tab : pandas dataframe. A table of fold change, basemean padj and gene names,
    following the format of deseq2 output
    
    ml_engine : string. machine learning engine to use. Possible choice : 1) 'lrlasso': logistic lasso
                                                                          2) 'lglasso' logistic group lasso
                                                                          3) 'lgen': logistic elastic net
                                                                          4) 'grrf' : guided regularized random forest
                                                                          5) 'lsvm': linear support vector machine
                                                                          6) 'lgpcc': logistic regression + pearson correlation coefficient
    
    feature_group : pandas dataframe. feature group from CoReg, or other analysis. Only valid when ml_engine = 'lrlasso'
    rep : int. number of replicates to have for each datset
    
    Returns
    -------
    auroc : float. auc of roc
    feature_num: int. number of selected features
    sn: float. sensitivity
    sp: float. specificity
    auprc : float. auc of prc curve
    best_param : best hyper parameter tunned from validation dataset
    '''
    features = np.array(feature_mat.iloc[:,1:],dtype = np.float128)
    labels = np.array(feature_mat.iloc[:,0])

    features_train_val, features_test, labels_train_val, labels_test = train_test_split(features,labels,train_size = 0.8)
    features_train, features_val, labels_train, labels_val = train_test_split(features_train_val, labels_train_val, train_size = 0.75)
    
    # Min-max transformation(normalization). Apply the same scaler to all feature matrices
    max_abs_fc = get_max_abs(features_train)
    features_train = features_train / max_abs_fc
    features_val = features_val / max_abs_fc
    features_test = features_test / max_abs_fc

    # Get the best model trained on training set, with best hyperparameters decided from validation set
    (best_model, best_param, feature_num) = train_and_val(features_train, labels_train, features_val, labels_val, ml_engine, feature_group)
    
    # Predict on test dataset.
    if ml_engine == "lrlasso":
        pred_raw = GGL.predict_gglasso(best_model,type = 'link',s = best_param, newx = BASE.as_matrix(features_test))
        pred_proba = 1/(1+np.exp(-np.array(pred_raw))) # Convert to numpy array. Calculate sigmoid output.
    
    elif ml_engine == "lglasso":
        pred_raw = GGL.predict_gglasso(best_model,type = 'link',s = best_param, newx = BASE.as_matrix(features_test))
        pred_proba = 1/(1+np.exp(-np.array(pred_raw))) # Convert to numpy array. Calculate sigmoid output
    
    elif ml_engine == "lgen":
        pred_proba = best_model.predict_proba(features_test)[:,1]
    
    elif ml_engine == "grrf":
        pred_raw = RRF.predict_RRF(best_model, features_test,type = "prob")
        pred_proba = np.array(pred_raw)[:,1]
    
    elif ml_engine == "lsvm":
        cal_svm = CalibratedClassifierCV(best_model, cv = 5).fit(features_train,labels_train)
        pred_proba = cal_svm.predict_proba(features_test)[:,1]
    
    elif ml_engine == "lgpcc":
        pred_proba = best_model.predict_proba(features_test[:,best_param])[:,1]
    
    # Get evaluation metrics, auc-roc, auc-prc
    labels_pred = np.round(pred_proba).astype(np.int32)
    precision,recall = precision_recall_curve(labels_test,pred_proba)[:2]
    auroc = roc_auc_score(labels_test,pred_proba)
    auprc = auc(recall,precision)
    best_param = np.mean(best_param)
    
    return(auroc, auprc, feature_num, best_param)