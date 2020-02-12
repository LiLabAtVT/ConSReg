import pandas as pd
import numpy as np
import re
import pickle
from joblib import Parallel, delayed
from collections import defaultdict,namedtuple
from os import path
from itertools import chain

from ConSReg.comb_peak_files import comb_peak_files
from ConSReg.peak_weight import get_weight
from ConSReg.edge_list import edge_list_from_peak, edge_list_to_nx_graph
from ConSReg.feature_mat import get_all_feature_mat
from ConSReg.train_test import test_by_cv
from ConSReg.imp_score import imp_score
from ConSReg.coreg import get_coreg_module

feature_mat_list = namedtuple('feature_mat_list','comp_names UR_feature_mat_list DR_feature_mat_list')
'''
Function to read a previously saved ConSReg object
'''
def load_obj(file_name):
    
    # Check if file exists
    if not path.exists(file_name):
        raise IOError(file_name + " does not exists")
        
    with open(file_name, 'rb') as input:
        obj = pickle.load(input)

    if obj.__class__.__name__ is not 'ConSReg':
        print("Sorry...loaded object is not a ConSReg object!")
        return(None)
    else:  
        return(obj)
    
'''
Helper function for parallelization
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

'''
Helper function
'''
def check_file_list(file_list):
    to_use = []
    for i,file in enumerate(file_list):
        if not path.exists(file):
            print(file + " does not exists. Skipped...")
        else:
            to_use.append(file)
    return(to_use)

'''
Main object of ConSReg
'''
class ConSReg:
    
    def __init__(self):

        # These attributes are used and modified internally
        self.__weight_adj = None
        self.__TF_to_target_graph = None
        self.__target_to_TF_graph = None
        self.__dap_seq_el = None
        self.__diff_tab_list = []
        
        # Modifying these attributes from outside of class is not recommended
        self._feature_mat_list_final = []
        self._feature_mat_list_dap = []
        self._feature_mat_list_reweighted = []
        self._feature_group = []
        self._diff_name_list = []
        
        # These are output results
        self.imp_scores_UR = None
        self.imp_scores_DR = None
        self.networks_UR = []
        self.networks_DR = []
        self.auroc = []
        self.auprc = []
        
        # Set up some flags for conveniently checking the status
        self.__preprocess = False
        self.__feature_final_generated = False
        self.__feature_dap_generated = False
        self.__feature_reweighted_generated = False
        self.__imp_generated = False
        self.__feature_group_generated = False
        self.__networks_generated = False
        self.__evaluated = False
        self.__use_peak_signal = False
    
    '''
    Function to save current ConSReg object
    '''
    def save_obj(self, file_name):
        
        with open(file_name, 'wb') as output:
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)
        
        return(self)
    '''
    process the binding site data (DAP-seq)
    '''
    def preprocess(self, dap_files, diff_files, atac_file,  gff_file, dap_chr_col = 0, dap_chr_start_col = 1, dap_chr_end_col = 2, dap_strand_col = None, dap_signal_col = None, atac_chr_col = 0, atac_chr_start_col = 1, atac_chr_end_col = 2, atac_signal_col = None, up_tss = 3000, down_tss = 500, up_type = 'all', down_type = 'all', use_peak_signal = False, n_jobs = 1, verbose = True):
        '''
        Parameters
        ----------
        dap_file : list. File names of DAP-seq peak files (bed format)
        diff_file : list. File names of differential contrasts, following
        the format of DESeq2 output file
        
        atac_file : string. File name of atac peak files (bed format). None if no atac-seq file is available
        gff_file : string. File name of genome annotation gff file
        
        dap_chr_col: int: column number for dap-seq chromosome information, 0 indexed.
        dap_chr_start_col: int, column number for dap-seq peak start position, 0 indexed.
        dap_chr_end_col: int, column number for dap-seq peak end position, 0 indexed.
        dap_strand_col: int/None, column number for dap-seq peak strand information, 0 indexed.
        dap_signal_col: int/None, column number for dap-seq peak signal value, 0 indexed.
        
        atac_chr_col: column number for atac-seq chromosome information, 0 indexed.
        atac_chr_start_col: column number for atac-seq peak start position, 0 indexed.
        atac_chr_end_col: column number for atac-seq peak end position, 0 indexed.
        atac_signal_col: column number for atac-seq peak signal value, 0 indexed.
        
        up_tss : positions relative to upstream region of TSS. This is used
        for finding nearest gene for each binding site

        down_tss: positions relative to downstream region of TSS. This is used
        for finding nearest gene for each binding site
        
        up_type: type of binding sites. 'all' or 'intergenic'
        down_type: type of binding sites. 'all' or 'intron' or 'non_intron'
        
        use_peak_signal: True/False. Whether to use peak signal for ATAC-seq and DAP-seq?
        use_atac_peak_signal: True/False. 
        
        n_jobs: int, number of jobs (for parallelization)
        verbose: bool, whether to print out details?
        Returns
        -------
        self
        '''
        ################# Check input arguments ################
        # Check value of file names
        if type(dap_files) != list:
            raise ValueError("Argument 'dap_file' should be a list")
        if type(diff_files) != list:
            raise ValueError("Argument 'diff_file' should be a list")
        if type(atac_file) != str and atac_file is not None:
            raise ValueError("Argument 'atac_file' should be a string")
        if type(gff_file) != str:
            raise ValueError("Argument 'gff_file' should be a string")
        
        # Check column number of dap-seq and atac-seq
        for name,param in zip(['dap_chr_col', 'dap_chr_start_col', 'dap_chr_end_col', 'atac_chr_col', 'atac_chr_start_col', 'atac_chr_end_col'],[dap_chr_col, dap_chr_start_col, dap_chr_end_col, atac_chr_col, atac_chr_start_col, atac_chr_end_col]):
            if type(param) != int:
                raise ValueError("Argument '"+ name +"' should be an integer")
        for name,param in zip(['dap_strand_col','dap_signal_col','atac_signal_col'],[dap_strand_col,dap_signal_col,atac_signal_col]):
            if param is not None and type(param) != int:
                raise ValueError("Argument '"+ name +"' should be an integer or None")
        
        # Check value of TSS positions and binding site type
        if type(up_tss) != int:
            raise ValueError("Argument 'up_tss' should be an integer")
        if type(down_tss) != int:
            raise ValueError("Argument 'down_tss' should be an integer")
        if up_type is not 'all' and up_type is not "intergenic":
            raise ValueError("Argument 'up_type' should be 'all' or 'intergenic'")
        if down_type not in ['all', 'intron', 'non_intron']:
            raise ValueError("Argument 'down_type' should be 'all' or 'intron' or 'non_intron'")
        
        # Check 'verbose'
        if type(verbose) != bool:
            raise ValueError("Argument 'verbose' should be a boolean variable (True/False)")
        if type(use_peak_signal) != bool:
            raise ValueError("Argument 'use_peak_signal' should be a boolean variable (True/False)")
        
        dap_files = check_file_list(dap_files)
        diff_files = check_file_list(diff_files)
        
        if len(dap_files) == 0:
            raise IOError("None of files specified in argument 'dap_files' exists")
            
        if len(diff_files) == 0:
            raise IOError("None of files specified in argument 'diff_files' exists")
            
        if atac_file is not None:
            if not path.exists(atac_file):
                raise IOError(atac_file + " does not exists")
        if not path.exists(gff_file):
            raise IOError(gff_file + " does not exists")
        ################# End checking ################
                             
        if verbose:
            print("Merging DAP-seq peaks...")
        
        if self.__preprocess:
            print("Analysis aborted. There exists a preprocessed dataset.")
            return(self)
        
        self.__use_peak_signal = use_peak_signal
        
        # Merge DAP-seq peaks
        peak_tab = comb_peak_files(dap_files, dap_chr_col = dap_chr_col, dap_chr_start_col = dap_chr_start_col, dap_chr_end_col = dap_chr_end_col, dap_strand_col = dap_strand_col, dap_signal_col = dap_signal_col)

        # TF->target edge list identified using peaks within the specified regions
        if verbose:
            print("Done")
            print("Assigning CREs to nearest genes...")
        
        peak_tab_list = np.array_split(peak_tab, n_jobs)
        results = Parallel(n_jobs = n_jobs)(delayed(edge_list_from_peak)(
                                                tab,
                                                gff_file,
                                                up_tss = up_tss,
                                                down_tss = down_tss,
                                                up_type = up_type,
                                                down_type = down_type
                                                ) 
                              for tab in peak_tab_list)
        el_anno = pd.concat(map(lambda x: x[0],results,),axis = 0, ignore_index = True)
        peak_tab = pd.concat(map(lambda x: x[1],results),axis = 0, ignore_index = True)
        
        # Add dap-seq signal values to edge list
        el_anno = pd.concat([el_anno,peak_tab.loc[:,"signalValue"]],axis = 1)
        
        '''
        Note that repeated edges can exist because we want to
        keep all peaks for one TF->target interaction so as to 
        overlap with all possible open chromatin regions
        '''
        # Generate edge list without repeated edges. 
        peak_count = el_anno.groupby(['TFID','targetID'], as_index = False).size()
        peak_sum = el_anno.groupby(['TFID','targetID'], as_index = False).sum()
        peak_sum['count'] = peak_count.values
        el_nr = peak_sum
        self.__dap_seq_el = el_nr.loc[:,["TFID","targetID"]]

        if verbose:
            print("Done")
            print("Overlapping CREs with ATAC-seq...")
            
        # Compute the overlaps between DAP-seq and ATAC-seq. Use the results as weights
        if atac_file is not None:
            peak_tab_list = np.array_split(peak_tab, n_jobs)
            results = Parallel(n_jobs=n_jobs)(delayed(get_weight)(
                                                tab,
                                                atac_file,
                                                atac_chr_col = atac_chr_col,
                                                atac_chr_start_col = atac_chr_start_col,
                                                atac_chr_end_col = atac_chr_end_col,
                                                atac_signal_col = atac_signal_col,
                                                use_peak_signal = use_peak_signal
                                              )
                                              for tab in peak_tab_list)
            weight = list(chain(*results))
            
            # Sum the weights for the same TF-target interaction
            weight_dict = defaultdict(int)
            for i,row in enumerate(el_anno.itertuples()):
                weight_dict[(row.TFID,row.targetID)] += weight[i]
                
            weight_df = pd.DataFrame(list(weight_dict.values()))
            weight_df.index = pd.MultiIndex.from_tuples(list(weight_dict.keys()))

            # Generate the adj of weights
            self.__weight_adj = weight_df.unstack(fill_value = 0).transpose()
            self.__weight_adj.index = self.__weight_adj.index.droplevel(level=0)
        else:
            # If atac-seq data is not provided:
            self.__weight_adj = None 
            
        # Convert edge list to networkx graph objects
        self.__TF_to_target_graph, self.__target_to_TF_graph = edge_list_to_nx_graph(el_nr)
        
        if verbose:
            print("Done")
            print("Reading diff tables...")
        
        for file in diff_files:
            self._diff_name_list.append(re.match(".*/([^/]*)$",file).group(1))
            self.__diff_tab_list.append(pd.read_csv(file,index_col = 0))
        
        if verbose:
            print("Done")
        
        # Mark the flag. Preprocessing is done
        self.__preprocess = True
        return(self)   
    
    '''
    Generate feature matrices for all conditions
    '''
    def gen_feature_mat(self, neg_type, verbose = True):
                             
        ################# Check input arguments ################
        # Check neg_type
        if neg_type not in ["udg","ndeg","leg","high_mean"]:
            raise ValueError("Argument 'neg_type' can be only specified as one of 'udg','ndeg','leg' or 'high_mean'")
        
        # Check 'verbose'
        if type(verbose) != bool:
            raise ValueError("Argument 'verbose' should be a boolean variable")
        ################# End checking ################
                             
        if not self.__preprocess:
            print("Analysis aborted. Data is not preprocssed yet")
            return(self)
        
        if verbose:
            print("Existing feature matrices will be overwritten.")
            print("Generating feature matrices...")
        
        self._feature_mat_list_dap = []
        self._feature_mat_list_reweight = []
        self._feature_mat_list_final = []
        
        for diff_tab in self.__diff_tab_list:
            feature_mat_list_up,feature_mat_list_down = get_all_feature_mat(diff_tab, neg_type, self.__TF_to_target_graph, self.__target_to_TF_graph, self.__weight_adj, self.__use_peak_signal)
            self._feature_mat_list_dap.append((feature_mat_list_up[0],feature_mat_list_down[0]))
            self._feature_mat_list_reweight.append((feature_mat_list_up[1],feature_mat_list_down[1]))
            self._feature_mat_list_final.append((feature_mat_list_up[2],feature_mat_list_down[2]))

        if verbose:
            print("Done")
        
        # Mark the flag. Feature matrices have been generated
        self.__feature_final_generated = True
        self.__feature_dap_generated = True
        self.__feature_reweighted_generated = True
        
        return(self)
    
    '''
    Functions to extract different types of feature matrices.
    '''
    def get_feature_mat_dap(self):
        '''
        Parameters
        ----------
        self
        
        Returns
        -------
        self
        '''
        if not self.__feature_dap_generated:
            print("Analysis aborted. Feature matrices have not been generated yet.")
            
        comp_names = []
        ur_mat_list = []
        dr_mat_list = []
        for comp_name,(ur_mat,dr_mat) in zip(self._diff_name_list, self._feature_mat_list_dap):
            comp_names.append(comp_name)
            ur_mat_list.append(ur_mat)
            dr_mat_list.append(dr_mat)
        
        to_return = feature_mat_list(comp_names, ur_mat_list, dr_mat_list)
        return(to_return)
            
    def get_feature_mat_reweight(self):
        '''
        Parameters
        ----------
        self
        
        Returns
        -------
        self
        '''
        if not self.__feature_reweighted_generated:
            print("Analysis aborted. Feature matrices have not been generated yet.")
            
        comp_names = []
        ur_mat_list = []
        dr_mat_list = []
        for comp_name,(ur_mat,dr_mat) in zip(self._diff_name_list, self._feature_mat_list_reweight):
            comp_names.append(comp_name)
            ur_mat_list.append(ur_mat)
            dr_mat_list.append(dr_mat)
        
        to_return = feature_mat_list(comp_names, ur_mat_list, dr_mat_list)
        return(to_return)
        
    def get_feature_mat_final(self):
        '''
        Parameters
        ----------
        self
        
        Returns
        -------
        self
        '''
        if not self.__feature_final_generated:
            print("Analysis aborted. Feature matrices have not been generated yet.")
            
        comp_names = []
        ur_mat_list = []
        dr_mat_list = []
        for comp_name,(ur_mat,dr_mat) in zip(self._diff_name_list, self._feature_mat_list_final):
            comp_names.append(comp_name)
            ur_mat_list.append(ur_mat)
            dr_mat_list.append(dr_mat)
        
        to_return = feature_mat_list(comp_names, ur_mat_list, dr_mat_list)
        return(to_return)
    
    '''
    Get importance scores for each feature matrix
    '''
    def compute_imp_score(self, n_resampling = 200, n_jobs = 1, verbose = True):
        '''
        Parameters
        ----------
        n_resampling : number of resampling times
        n_jobs: number of jobs (for parallelization)
        verbose: output detailed results or not?

        Returns
        -------
        self
        '''
        ################# Check input arguments ################
        if type(n_resampling) != int:
            raise ValueError("Argument 'n_resampling' should be an integer")
        elif n_resampling <= 1:
            raise ValueError("Argument 'n_resampling' should be a positive integer")
        
        if type(n_jobs) != int:
            raise ValueError("Argument 'n_jobs' should be an integer")
        elif n_jobs <= 1:
            raise ValueError("Argument 'n_jobs' should be a positive integer")
                             
        # Check 'verbose'
        if type(verbose) != bool:
            raise ValueError("Argument 'verbose' should be a boolean variable")
        ################# End checking ################
                             
        if verbose:
            print("Existing importance scores will be overwritten.")
            print("Performing stability selection and compute importance score for each TF...")
        
        if self.__imp_generated:
            print("Analysis aborted. Importance scores have been already generated")
            return(self)
        
        # Generate importance scores for each feature matrix
        imp_score_list_up = []
        imp_score_list_down = []
        
        for (feature_mat_up, feature_mat_down) in self._feature_mat_list_final:
            if feature_mat_up is not None:
                imp_score_up = imp_score(feature_mat_up, n_resampling = n_resampling, verbose = verbose, fold = 10, n_jobs = n_jobs)
            else:
                imp_score_up = None
            
            if feature_mat_down is not None:
                imp_score_down = imp_score(feature_mat_down, n_resampling = n_resampling, verbose = verbose, fold = 10, n_jobs = n_jobs)
            else:
                imp_score_down = None
                
            imp_score_list_up.append(imp_score_up)
            imp_score_list_down.append(imp_score_down)

        self.imp_scores_UR = pd.concat(imp_score_list_up,axis = 1).fillna(0)
        self.imp_scores_DR = pd.concat(imp_score_list_down,axis = 1).fillna(0)
        
        # Only use names of contrasts that did not generate None feature matrix
        self.imp_scores_UR.columns = [self._diff_name_list[i] for i,_ in enumerate(imp_score_list_up) if _ is not None]
        self.imp_scores_DR.columns = [self._diff_name_list[i] for i,_ in enumerate(imp_score_list_down) if _ is not None] 
        
        if verbose:
            print("Done")
        
        # Mark the flag. Importance scores have been generated
        self.__imp_generated = True
        
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
        ################# Check input arguments ################
        if deep_split not in [0,1,2,3,4]:
            raise ValueError("Argument 'deep_split' can only be chosen from integer values {0,1,2,3,4}")
                             
        # Check 'verbose'
        if type(verbose) != bool:
            raise ValueError("Argument 'verbose' should be a boolean variable")
        ################# End checking ################
                             
        if verbose:
            print("Existing feature groups will be overwritten.")
            print("Generating CoReg modules...")
        
        if not self.__preprocess:
            print("Analysis aborted. Data is not preprocssed yet")
            return(self)
            
        self._feature_group = []
        for diff_tab in self.__diff_tab_list:
            feature_group = get_coreg_module(self.__dap_seq_el, diff_tab, deepSplit = deep_split)
            self._feature_group.append(feature_group)
        
        if verbose:
            print("Done")
        
        # Mark the flag
        self.__feature_group_generated = True
        
        return(self)
    
    '''
    Infer networks
    '''
    def gen_networks(self, imp_cutoff = 0.5, verbose = True):
        '''
        Parameters
        ----------
        imp_cutoff : float. TFs will be selected if its importance socre > imp_cutoff
        verbose: output detailed results or not？

        Returns
        -------
        self
        '''
        ################# Check input arguments ################
        if type(imp_cutoff) != float:
            raise ValueError("Argument 'imp_cutoff' should be a float number")
        elif imp_cutoff < 0 or imp_cutoff > 1:
            raise ValueError("Argument 'imp_cutoff' should be within [0,1]")
                             
        # Check 'verbose'
        if type(verbose) != bool:
            raise ValueError("Argument 'verbose' should be a boolean variable")
        ################# End checking ################                     
        
        if not self.__imp_generated or not self.__feature_final_generated:
            print("Analysis aborted. There are no importance scores/feature matrices generated yet.")
            return(self)
        
        imp_scores_UR_list = list(map(lambda x: self.imp_scores_UR.loc[:,x], self.imp_scores_UR.columns))
        imp_scores_DR_list = list(map(lambda x: self.imp_scores_DR.loc[:,x], self.imp_scores_DR.columns))
        
        to_iter = zip(imp_scores_UR_list, imp_scores_DR_list, self._feature_mat_list_final)
        
        if verbose:
            print("Existing networks will be overwritten.")
            print("Generating networks...")
            
        self.networks_UR = []
        self.networks_DR = []
        
        # None feature matrix will be skipped
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
            
            if feature_mat_down is not None:
                select_TF_down = imp_score_down[imp_score_down > imp_cutoff].index

                adj_down = feature_mat_down.loc[(feature_mat_down.loc[:,"label"] == 1), select_TF_down]

                el_down = adj_down.stack().reset_index()           
                el_down.columns = ["targetID","TFID","edge_weight"]
                el_down = el_down.loc[(el_down.loc[:,"edge_weight"] != 0),["TFID","targetID"]]
                el_down.sort_values(by = "TFID",axis = 0,inplace = True)
            
            self.networks_UR.append(el_up)
            self.networks_DR.append(el_down)
        
        if verbose:
            print("Done")
            
        # Mark the flag. Networks were already generated
        self.__networks_generated = True
        
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
        ################# Check input arguments ################
        if ml_engine not in ['lrlasso','lglasso','lgen','grrf','lsvm','lgpcc']:
            raise ValueError("Argument 'ml_engine' can only be chosen from {'lrlasso','lglasso','lgen','grrf','lsvm','lgpcc'}")
        
        if type(n_jobs) != int:
            raise ValueError("Argument 'n_jobs' should be an integer")
        elif n_jobs <=0:
            raise ValueError("Argument 'n_jobs' should be a positive integer")
                             
        if type(rep) != int:
            raise ValueError("Argument 'rep' should be an integer")
        elif rep <=0:
            raise ValueError("Argument 'rep' should be a positive integer")
                             
        # Check 'verbose'
        if type(verbose) != bool:
            raise ValueError("Argument 'verbose' should be a boolean variable")
        ################# End checking ################       
                             
        if verbose:
            print("Performing cross-validation for each feature matrix using {} engine...".format(ml_engine))
            print("Old evaluation results will be ovewritten")

        if not self.__feature_final_generated:
            print("Analysis aborted. There are no feature matrices generated yet.")
            return(self)
        
        # Generate feature group if necessary
        if ml_engine != 'lglasso':
            self._feature_group = [None] * len(self.__diff_tab_list)
        else:
            if not self._feature_group_generated:
                print("Analysis aborted. 'lglasso' is selected but there are no feature groups generated yet")
                return(self)
        
        iter_list = zip(self._diff_name_list, self.__diff_tab_list, self._feature_mat_list_final, self._feature_group)
        
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
        self.auroc = pd.DataFrame(list(map(lambda x: x[0], out)),columns = ["diff_name","auroc_mean_UR", "auroc_std_UR", "auroc_mean_DR", "auroc_std_DR"])
        self.auprc = pd.DataFrame(list(map(lambda x: x[1], out)),columns = ["diff_name","auprc_mean_UR", "auprc_std_UR", "auprc_mean_DR", "auprc_std_DR"])
        
        if verbose:
            print("Done")
        
        # Mark the flag
        self.__evaluated = True
        return(self)
