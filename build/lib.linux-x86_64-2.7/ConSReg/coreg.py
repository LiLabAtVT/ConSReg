'''
Function for getting coreg module
'''

# Author: (Alex) Qi Song <alexsong@vt.edu>

from rpy2.robjects import pandas2ri,StrVector
from pandas import read_csv,DataFrame,concat
from rpy2.robjects.packages import importr # This is for loading R library

def get_coreg_module(edge_list, diff_tab, deepSplit):
    
    if diff_tab is not None:
        
        # Get expressed genes
        expr_genes = diff_tab.loc[diff_tab.loc[:,"baseMean"] > 0,].index
    
    else:
        expr_genes = concat([edge_list.loc[:,"TFID"], edge_list.loc[:,"targetID"]],axis = 0).drop_duplicates()
    
    # Filter interactions by expressed TFs and targets
    expr_TFs_idx = edge_list.loc[:,"TFID"].isin(expr_genes)
    expr_TFs_ID = edge_list.loc[expr_TFs_idx,"TFID"].unique()
    expr_targets_idx = edge_list.loc[:,"targetID"].isin(expr_genes)
    expr_edge_list = edge_list.loc[expr_TFs_idx & expr_targets_idx,]
    gene_names = StrVector(expr_TFs_ID)

    # Prepare R environment
    base = importr('base')
    CR = importr('CoReg')
    pandas2ri.activate()
    
    '''
    Run CoReg in R environment, get the module assignment
    To speed up computation, CoReg calculation will be only 
    performed on the subset specified by gene_names
    '''
    graph = CR.networkFromEdgeList(expr_edge_list)
    coreg_module = pandas2ri.ri2py(CR.CoReg(graph,gene_names,deepSplit=deepSplit)[0])
    coreg_module.index = coreg_module.loc[:,"ID"]
    
    # Get module info for TFs, return as a dictionary
    expr_TFs_ID = edge_list.loc[expr_TFs_idx,"TFID"].unique()
    TF_module = coreg_module.loc[expr_TFs_ID,"module"].astype(int)
    
    return(TF_module)