import pandas as pd
import numpy as np
import re

from joblib import Parallel, delayed
from collections import defaultdict
from comb_peak_files import comb_peak_files
from peak_weight import get_weight
from edge_list import edge_list_from_peak, edge_list_to_nx_graph
from feature_mat import get_all_feature_mat
from train_test import test_by_cv
from imp_score import imp_score
from coreg import get_coreg_module

'''
Main object of ConSReg
'''
class ConSReg:
    
    def __init__(self):
        self.weight_adj = None
        self.TF_to_target_graph = None
        self.target_to_TF_graph = None
        self.dap_seq_el = None
        self.feature_mat_list = []
        self.feature_group = []
        self.imp_scores_UR = None
        self.imp_scores_DR = None
        self.diff_tab_list = []
        self.diff_name_list = []
        self.networks_UR = []
        self.networks_DR = []
        self.auroc = []
        self.auprc = []
    
    '''
    process the binding site data (DAP-seq)
    '''
    def preprocess(self, dap_file, diff_file, atac_file,  gff_file, upTSS = 3000, downTSS = 500, verbose = True):
        '''
        Parameters
        ----------
        dap_folder : string. Folder of DAP-seq peak files (bed format)
        diff_file : list. File names of differential contrasts, following
        the format of DESeq2 output file
        
        atac_file : string. File name of atac peak files (bed format)
        gff_file : string. File name of genome annotation gff file

        upTSS : positions relative to upstream region of TSS. This is used
        for finding nearest gene for each binding site

        downTSS: positions relative to downstream region of TSS. This is used
        for finding nearest gene for each binding site

        Returns
        -------
        self
        '''
        
        if verbose:
            print("Merging DAP-seq peaks...")
            
        # Merge DAP-seq peaks
        peak_tab = comb_peak_files(dap_file)

        # TF->target edge list identified using peaks within the specified regions
        el_anno = edge_list_from_peak(peak_tab,gff_file,upTSS = 3000, downTSS = 500)

        # This reduces the number of peaks to be calculated for weights. Use the peaks that have genes annotated within the specified region
        # Note that repeated edges can exist because we want to keep all peaks for one TF->target interaction so as to overlap with all possible
        # open chromatin regions
        peak_tab = peak_tab.iloc[el_anno.index,]

        # Generate edge list without repeated edges
        el_nr = el_anno.groupby(el_anno.columns.tolist(), as_index = False).size()
        el_nr = pd.DataFrame(el_nr.values, index = el_nr.index)
        el_nr.reset_index(inplace = True)
        el_nr.columns = ["TFID","targetID","count"]
        self.dap_seq_el = el_nr.loc[:,["TFID","targetID"]]

        if verbose:
            print("Done")
            print("Overlapping DAP-seq with ATAC-seq...")
            
        # Compute the overlaps between DAP-seq and ATAC-seq. Use the results as weights
        weight = get_weight(peak_tab, atac_file)

        # Sum the weights for the same TF-target interaction
        weight_dict = defaultdict(int)
        for i,row in enumerate(el_anno.itertuples()):
            weight_dict[(row.TFID,row.targetID)] += weight[i]

        weight_df = pd.DataFrame(weight_dict.values())
        weight_df.index = pd.MultiIndex.from_tuples(weight_dict.keys())

        # Generate the adj of weights
        self.weight_adj = weight_df.unstack(fill_value = 0).transpose()
        self.weight_adj.index = self.weight_adj.index.droplevel(level=0)

        # Convert edge list to networkx graph objects
        self.TF_to_target_graph, self.target_to_TF_graph = edge_list_to_nx_graph(el_nr)
        
        if verbose:
            print("Done")
            print("Reading diff tables...")
        
        for file in diff_file:
            self.diff_name_list.append(re.match(".*/([^/]*)$",file).group(1))
            self.diff_tab_list.append(pd.read_csv(file,index_col = 0))
        
        if verbose:
            print("Done")
            
        return(self)   
    
    '''
    Generate feature matrices for all conditions
    '''
    def gen_feature_mat(self, verbose = True):
        if self.weight_adj is None or self.TF_to_target_graph is None or self.target_to_TF_graph is None or len(self.diff_tab_list) == 0:
            print("Analysis aborted. Data is not preprocssed yet")
            return(self)
            
        if verbose:
            print("Generating feature matrices...")

        for diff_tab in self.diff_tab_list:
            feature_mat_up,feature_mat_down = get_all_feature_mat(diff_tab, self.TF_to_target_graph, self.target_to_TF_graph, self.weight_adj)
            self.feature_mat_list.append((feature_mat_up,feature_mat_down))

        if verbose:
            print("Done")
        
        return(self)
    
    '''
    Get importance scores for each feature matrix
    '''
    def compute_imp_score(self, n_resampling = 200, n_jobs = 1, verbose = True):
        '''
        Parameters
        ----------
        feature_mat_list : a list of numpy 2d array. feature matrices

        Returns
        -------
        self
        '''

        if verbose:
            print("Performing stability selection and compute importance score for each TF...")
        
        if len(self.feature_mat_list) == 0:
            print("Analysis aborted. There are no feature matrices generated yet.")
            return(self)
        
        # Generate importance scores for each feature matrix
        imp_score_list_up = []
        imp_score_list_down = []
        
        for (feature_mat_up, feature_mat_down) in self.feature_mat_list:
            if feature_mat_up is not None:
                imp_score_up = imp_score(feature_mat_up, n_resampling = n_resampling, verbose = 0, fold = 10, n_jobs = n_jobs)
            else:
                imp_score_up = None
            
            if feature_mat_down is not None:
                imp_score_down = imp_score(feature_mat_down, n_resampling = n_resampling, verbose = 0, fold = 10, n_jobs = n_jobs)
            else:
                imp_score_down = None
                
            imp_score_list_up.append(imp_score_up)
            imp_score_list_down.append(imp_score_down)

        self.imp_scores_UR = pd.concat(imp_score_list_up,axis = 1).fillna(0)
        self.imp_scores_DR = pd.concat(imp_score_list_down,axis = 1).fillna(0)
        
        self.imp_scores_UR.columns = self.diff_name_list
        self.imp_scores_DR.columns = self.diff_name_list
        
        if verbose:
            print("Done")
            
        return(self)
    
    '''
    Get feature group, by running coreg
    '''
    def get_feature_group(self, deep_split, verbose = True):
        '''
        Parameters
        ----------
        deep_split : int. Control the split of hierarchical tree. Larger value will result
        in smaller coreg modules. possible values : {0,1,2,3,4}
        
        Returns
        -------
        self
        '''
        
        if verbose:
            print("Generating CoReg modules...")
        
        if len(self.diff_tab_list) == 0:
            print("Analysis aborted. There are no diff tables generated yet.")
            return(self)
        
        if self.dap_seq_el is None:
            print("Analysis aborted. DAP-seq edge list is not generated yet")
            return(self)
            
        for diff_tab in self.diff_tab_list:
            feature_group = get_coreg_module(self.dap_seq_el, diff_tab, deepSplit = deep_split)
            self.feature_group.append(feature_group)
        
        if verbose:
            print("Done")
        
        return(self)
    
    '''
    Infer networks
    '''
    def gen_networks(self, imp_cutoff = 0.5, verbose = True):
        '''
        Parameters
        ----------
        imp_score_list : a list of pandas dataframe. a list of dataframe, each corresponds to 
        importance scores for one condition

        feature_mat_list : a list of numpy 2d array. feature matrices
        imp_cutoff : float. TFs will be selected if its importance socre > imp_cutoff

        Returns
        -------
        self
        '''
        
        if self.imp_scores_UR is None and self.imp_scores_DR is None:
            print("Analysis aborted. There are no importance scores generated yet.")
            return(self)
            
        if len(self.feature_mat_list) == 0:
            print("Analysis aborted. There are no feature matrices generated yet.")
            return(self)
        
        imp_scores_UR_list = map(lambda x: self.imp_scores_UR.loc[:,x], self.imp_scores_UR.columns)
        imp_scores_DR_list = map(lambda x: self.imp_scores_DR.loc[:,x], self.imp_scores_DR.columns)
        
        to_iter = zip(imp_scores_UR_list, imp_scores_DR_list, self.feature_mat_list)
        
        for imp_score_up, imp_score_down, (feature_mat_up, feature_mat_down) in to_iter:

            if feature_mat_up is not None:
                # Select TFs with imp score > imp_cutoff
                select_TF_up = imp_score_up[imp_score_up > imp_cutoff].index
                
                # Generate adjacency matrix then edge list
                adj_up = feature_mat_up.loc[(feature_mat_up.loc[:,"label"] == 1), select_TF_up]
                
                # Convert adj to edge list
                el_up = adj_up.stack().reset_index()
                el_up.columns = ["targetID","TFID","edge_weight"]
                el_up = el_up.loc[(el_up.loc[:,"edge_weight"] != 0),["TFID","targetID"]]
                
                # Sort edge list by TF
                el_up.sort_values(by = "TFID",axis = 0,inplace = True)
            else:
                el_up = None
            
            if feature_mat_down is not None:
                select_TF_down = imp_score_down[imp_score_down > imp_cutoff].index

                adj_down = feature_mat_down.loc[(feature_mat_down.loc[:,"label"] == 1), select_TF_down]

                el_down = adj_down.stack().reset_index()           
                el_down.columns = ["targetID","TFID","edge_weight"]
                el_down = el_down.loc[(el_down.loc[:,"edge_weight"] != 0),["TFID","targetID"]]
                el_down.sort_values(by = "TFID",axis = 0,inplace = True)
            else:
                el_down = None
            
            self.networks_UR.append(el_up)
            self.networks_DR.append(el_down)
        return(self)

    '''
    Run cross validation. Tune best parameters from validation data; Train model on training data
    with best parameters. Test model on test data
    '''
    def eval_by_cv(self, ml_engine, n_jobs = 1, verbose = True, rep = 5):
        '''
        Parameters
        ----------
        ml_engine : string. machine learning engine to use. Possible choice : 1) 'lrlasso': logistic lasso
                                                                          2) 'lglasso' logistic group lasso
                                                                          3) 'lgen': logistic elastic net
                                                                          4) 'grrf' : guided regularized random forest
                                                                          5) 'lsvm': linear support vector machine
                                                                          6) 'lgpcc': logistic regression + pearson correlation coefficient
        
        feature_group : pandas dataframe. feature group from CoReg, or other analysis. Only valid when ml_engine = 'lrlasso'
        
        Returns
        -------
        self
        '''
    
        if verbose:
            print("Performing cross-validation for each feature matrix using {} engine...".format(ml_engine))

        if len(self.feature_mat_list) == 0:
            print("Analysis aborted. There are no feature matrices generated yet.")
            return(self)
        
        if len(self.diff_tab_list) == 0:
            print("Analysis aborted. There are no diff tables generated yet.")
            return(self)
        
        # Generate feature group if necessary
        if ml_engine != 'lglasso':
            self.feature_group = [None] * len(self.diff_tab_list)
        else:
            if len(self.feature_group) == 0:
                print("Analysis aborted. 'lglasso' is selected but there are no feature groups generated yet")
                return(self)
        
        iter_list = zip(self.diff_name_list, self.diff_tab_list, self.feature_mat_list, self.feature_group)
        
        out = Parallel(n_jobs = n_jobs, verbose=15)(delayed(parallel_eval)(
                                            diff_name,
                                            diff_tab,
                                            feature_mat_up,
                                            feature_mat_down,
                                            feature_group,
                                            ml_engine = ml_engine,
                                            rep = rep
                                            ) 
                          for diff_name, diff_tab,(feature_mat_up, feature_mat_down),feature_group in iter_list)
            
        # Format output result
        self.auroc = pd.DataFrame(map(lambda x: x[0], out),columns = ["diff_name","auroc_mean_UR", "auroc_std_UR", "auroc_mean_DR", "auroc_std_DR"])
        self.auprc = pd.DataFrame(map(lambda x: x[1], out),columns = ["diff_name","auprc_mean_UR", "auprc_std_UR", "auprc_mean_DR", "auprc_std_DR"])
        
        if verbose:
            print("Done")
        return(self)

'''
Helper function
'''
def parallel_eval(diff_name, diff_tab, feature_mat_up, feature_mat_down, feature_group, ml_engine, rep):
            
    if feature_mat_up is None:
        auroc_mean_up,auroc_std_up,auprc_mean_up,auprc_std_up = None,None,None,None
    else:
        auroc_up = []
        auprc_up = []
        for i in range(rep):
            auroc,auprc,_,_ = test_by_cv(feature_mat_up, diff_tab, ml_engine = ml_engine, feature_group = feature_group)
            auroc_up.append(auroc)
            auprc_up.append(auprc)
        auroc_mean_up = np.mean(auroc_up)
        auroc_std_up = np.std(auroc_up)
        auprc_mean_up = np.mean(auprc_up)
        auprc_std_up = np.std(auprc_up)

    if feature_mat_down is None:
        auroc_mean_down,auroc_std_down,auprc_mean_down,auprc_std_down = None,None,None,None
    else:
        auroc_down = []
        auprc_down = []
        for i in range(rep):
            auroc,auprc,_,_ = test_by_cv(feature_mat_down, diff_tab, ml_engine = ml_engine, feature_group = feature_group)
            auroc_down.append(auroc)
            auprc_down.append(auprc)
        auroc_mean_down = np.mean(auroc_down)
        auroc_std_down = np.std(auroc_down)
        auprc_mean_down = np.mean(auprc_down)
        auprc_std_down = np.std(auprc_down)    

    auroc_res = (diff_name, auroc_mean_up, auroc_std_up, auroc_mean_down, auroc_std_down)
    auprc_res = (diff_name, auprc_mean_up, auprc_std_up, auprc_mean_down, auprc_std_down)
    return(auroc_res, auprc_res)