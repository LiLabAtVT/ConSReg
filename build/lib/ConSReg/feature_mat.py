'''
Functions for construcing feature matrix
'''

# Author: (Alex) Qi Song <alexsong@vt.edu>

import pandas as pd
import numpy as np
from networkx.convert_matrix import to_pandas_adjacency
from itertools import chain

'''
Get the TFs for the target genes from the graph
'''
def get_TF(genes, target_to_TF_graph):
    '''
    Parameters
    ----------
    genes : list like object of strings. gene name

    target_to_TF_graph : networkx DiGraph object. Digraph object of target->TF graph

    Returns
    -------
    TFs : a set of strings. TF names for the target genes
    '''
    TFs = set()
    _ = list(map(lambda x: TFs.update(target_to_TF_graph.neighbors(x)), genes))
    return(TFs)

'''
Filter TF->target and target->TF graph. Modify diff_tab. Assign log2FoldChange = 0 
to 'baseMean' of negative genes
'''
def filter(diff_tab, direction, neg_type, TF_to_target_graph, target_to_TF_graph):
    '''
    Parameters
    ----------
    diff_tab : pandas dataframe. A table of fold change, basemean padj and gene names,
    following the format of deseq2 output

    direction: string. 'up' for up-regulated 'down' for down-regulated
    neg_type: string. See the paper for more details
        "udg" -- undetected genes (UDGs), which have a mean expression value equal to zero.
    "leg" -- low-expressed genes (LEGs), which have mean expression between 0 and 0.5
    "ndeg" -- non-significantly differentially expressed genes (NDEGs), which have p-value > 0.05
    "high_mean" -- high mean genes, ndegs which have mean expression <= 8.
    TF_to_target_graph : networkx DiGraph object. Digraph object of TF->target graph
    target_to_TF_graph : networkx DiGraph object. Digraph object of target->TF graph

    Returns
    -------
    DEGs : pandas index. gene names of DEGs
    neg_genes : pandas index. gene names of negative genes
    DEG_TFs : a set of strings. gene names of the TFs of DEGs
    neg_TFs : a set of strings. gene names of the TFs of negative genes
    diff_tab : pandas dataframe. Modified table of fold change, basemean padj and gene names,
    following the format of deseq2 output

    TF_to_target_graph : networkx DiGraph object. Filtered graph. Digraph object of TF->target graph
    target_to_TF_graph : networkx DiGraph object. Filtered graph. Digraph object of target->TF graph
    '''

    '''
    For undetected genes (baseMean = 0) in both conditions, assign log2FoldChange = 0 to them
    '''
    all_genes = set(TF_to_target_graph.nodes).intersection(set(diff_tab.index))
    TF_to_target_graph = TF_to_target_graph.subgraph(all_genes)
    target_to_TF_graph = target_to_TF_graph.subgraph(all_genes)

    # Filter the genes by FC and p-values, get up DEGs, down DEGs and negative genes.
    if direction == "up":
        DEGs = set(diff_tab.loc[(diff_tab.loc[:,"padj"] < 0.05) & (diff_tab.loc[:,"log2FoldChange"] > 0),].index)

        if neg_type == "udg":
            neg_genes = diff_tab.loc[(diff_tab.loc[:,"baseMean"] == 0),].index
            diff_tab.loc[neg_genes,"log2FoldChange"] = 1
        elif neg_type == "ndeg":
            neg_genes = set(diff_tab.loc[(diff_tab.loc[:,"padj"] >= 0.05) & (diff_tab.loc[:,"log2FoldChange"] >= 0),].index)	
        elif neg_type == "leg":
            neg_genes = diff_tab.loc[(diff_tab.loc[:,"baseMean"] <= 0.5) & (diff_tab.loc[:,"baseMean"] > 0) & (diff_tab.loc[:,"log2FoldChange"] >= 0) & (diff_tab.loc[:,"log2FoldChange"] <= 0.5),].index
        elif neg_type == "high_mean":
            neg_genes = diff_tab.loc[(diff_tab.loc[:,"baseMean"] <= 8) & (diff_tab.loc[:,"padj"] >= 0.05) & (diff_tab.loc[:,"log2FoldChange"] > 0),].index
    
    # If direction is "down"
    else:
        DEGs = set(diff_tab.loc[(diff_tab.loc[:,"padj"] < 0.05) & (diff_tab.loc[:,"log2FoldChange"] < 0),].index)
    
        if neg_type == "udg":
            neg_genes = diff_tab.loc[(diff_tab.loc[:,"baseMean"] == 0),].index
            diff_tab.loc[neg_genes,"log2FoldChange"] = 1
        elif neg_type == "ndeg":
            neg_genes = set(diff_tab.loc[(diff_tab.loc[:,"padj"] >= 0.05) & (diff_tab.loc[:,"log2FoldChange"] < 0),].index)	
        elif neg_type == "leg":
            neg_genes = diff_tab.loc[(diff_tab.loc[:,"baseMean"] <= 0.5) & (diff_tab.loc[:,"baseMean"] > 0) & (diff_tab.loc[:,"log2FoldChange"] < 0) & (diff_tab.loc[:,"log2FoldChange"] > -0.5),].index
        elif neg_type == "high_mean":
            neg_genes = diff_tab.loc[(diff_tab.loc[:,"baseMean"] <= 8) & (diff_tab.loc[:,"padj"] >= 0.05) & (diff_tab.loc[:,"log2FoldChange"] < 0),].index        
        
    # Some genes might not have TF connections now. We only keep the genes with
    # at least one TF connection
    DEGs = [gene for gene,deg in TF_to_target_graph.in_degree(DEGs) if deg > 0]
    neg_genes = [gene for gene,deg in TF_to_target_graph.in_degree(neg_genes) if deg > 0]
    
    # If there are no DEG genes or background genes
    if len(DEGs) == 0 or len(neg_genes) == 0:
        return(None)
    
    # Sort genes based on FC
    DEGs = diff_tab.loc[DEGs,].sort_values('log2FoldChange', axis=0, ascending = False, inplace = False).index
    neg_genes = diff_tab.loc[neg_genes,['log2FoldChange']].abs().sort_values('log2FoldChange', axis=0, inplace = False).index
    
    # Balance the number of DEGs and neg_genes
    # Randomly select the same number of neg_genes
    if len(DEGs) > len(neg_genes):
        DEGs = DEGs[:len(neg_genes)]
    else:
        neg_genes = np.random.choice(neg_genes,size=len(DEGs),replace = False)
        neg_genes = pd.Index(neg_genes)
    
    # Get the TFs of top ranked genes. TFs should be expressed i.e. baseMean > 0
    DEG_TFs = get_TF(DEGs,target_to_TF_graph)         # Get a set
    neg_TFs = get_TF(neg_genes,target_to_TF_graph)     # Get a set
    
    # Remove the non expressed TFs and TFs with NA fold change
    non_exp_TFs = set(diff_tab.loc[diff_tab.loc[:,"baseMean"] == 0,].index)
    NA_TFs = set(diff_tab.loc[diff_tab.loc[:,"log2FoldChange"].isna(),].index)
    TFs_to_remove = non_exp_TFs.union(NA_TFs)
    DEG_TFs = DEG_TFs.difference(TFs_to_remove)
    neg_TFs = neg_TFs.difference(TFs_to_remove)
    
    return(DEGs, neg_genes, DEG_TFs, neg_TFs, diff_tab, TF_to_target_graph, target_to_TF_graph)

'''
Generate a single feature matrix
'''
def gen_feature_mat(diff_tab, target_to_TF_graph, weight_adj, DEGs, neg_genes, DEG_TFs, neg_TFs, use_peak_signal):
    '''
    Parameters
    ----------
    diff_tab : pandas dataframe. A table of fold change, basemean padj and gene names,
    following the format of deseq2 output
    
    target_to_TF_graph : networkx DiGraph object. Digraph object of target->TF graph
    weight_adj : pandas dataframe. Weight matrix computed from overlapping regions
    between DAP-seq peaks and ATAC-seq peaks
    
    DEGs:  pandas index. gene names of DEGs
    neg_genes : pandas index. gene names of negative genes
    DEG_TFs : a set of strings. gene names of the TFs of DEGs
    neg_TFs : a set of strings. gene names of the TFs of negative genes
    
    Returns
    -------
    feature_mat_dap : a pandas datarame, dap-seq feature matrix (with 1,0 as values)
    feature_mat_reweight : a pandas dataframe, atac-seq reweighted feature matrix
    feature_mat_final : a pandas datarame, final feature matrix
    '''
    
    row_genes = DEGs | neg_genes
    col_genes = DEG_TFs | neg_TFs # Make sure there are no duplicate TFs
    
    '''
    Generate the feature matrix by multiply subset of adj matrix and fold change column
    When weight_adj is not None, weights are used (overlap from atac-seq and dap-seq peak signals are included).
    Otherwise, only use the count of dap-seq peaks.
    '''
    if weight_adj is not None:
        adj = weight_adj.loc[row_genes,col_genes]
    else:
        if use_peak_signal:
            adj = to_pandas_adjacency(target_to_TF_graph.subgraph(row_genes | col_genes),weight = 'signalValue').loc[row_genes,col_genes]
        else:
            adj = to_pandas_adjacency(target_to_TF_graph.subgraph(row_genes | col_genes),weight = 'count').loc[row_genes,col_genes]
    
    fc = diff_tab.loc[col_genes,'log2FoldChange']
    feature_mat_dap = adj.where(adj == 0, 1) 
    feature_mat_reweight = adj
    feature_mat_final = pd.DataFrame(adj.values * fc.values.reshape(1,fc.shape[0]),index = row_genes, columns = col_genes)
    
    # Re-order the featureMat
    feature_mat_dap = feature_mat_dap.loc[list(DEGs)+list(neg_genes)]
    feature_mat_reweight = feature_mat_reweight.loc[list(DEGs)+list(neg_genes)]
    feature_mat_final = feature_mat_final.loc[list(DEGs)+list(neg_genes),]
    
    # Add label column
    label = pd.DataFrame(int(feature_mat_final.shape[0]/2)*[1] + int(feature_mat_final.shape[0]/2)*[0], index = feature_mat_final.index, columns = ['label'])
    feature_mat_dap = pd.concat([label, feature_mat_dap], axis = 1)
    feature_mat_reweight = pd.concat([label, feature_mat_reweight], axis = 1)
    feature_mat_final = pd.concat([label, feature_mat_final], axis = 1)
    
    return([feature_mat_dap, feature_mat_reweight, feature_mat_final])

'''
Get up-regulated (UR) feature matrix and down-regulated (DR) feature matrix
'''
def get_all_feature_mat(diff_tab, neg_type, TF_to_target_graph, target_to_TF_graph, weight_adj, use_peak_signal):
    '''
    Parameters
    ----------
    diff_tab : pandas dataframe. A table of fold change, basemean padj and gene names,
    following the format of deseq2 output
    neg_type: string. See the paper for more details
        "udg" -- undetected genes (UDGs), which have a mean expression value equal to zero.
    "leg" -- low-expressed genes (LEGs), which have mean expression between 0 and 0.5
    "ndeg" -- non-significantly differentially expressed genes (NDEGs), which have p-value > 0.05
    "high_mean" -- high mean genes, ndegs which have mean expression <= 8. 
    
    TF_to_target_graph : networkx DiGraph object. Digraph object of TF->target graph
    target_to_TF_graph : networkx DiGraph object. Digraph object of target->TF graph
    weight_adj -- weight adjacency matrix computed from open chromatin region data (e.g. ATAC-seq)
    
    Returns
    -------
    feature_mat_list_up : A list of pandas dataframes, UR feature matrices (dap-seq feature matrix, reweighted feature matrix, and final feature matrix).In each matrix, the first column is the class label and other columns are values for each TF(feature)
    feature_mat_list_down : A list of pandas dataframes, DR feature matrices (dap-seq feature matrix, reweighted feature matrix, and final feature matrix).In each matrix, the first column is the class label and other columns are values for each TF(feature)
    '''
    filter_up = filter(diff_tab,"up", neg_type, TF_to_target_graph, target_to_TF_graph)
    
    if filter_up is not None:
        DEGs_up, neg_genes_up, DEGs_TFs_up, neg_TFs_up, diff_tab_up, TF_to_target_graph_up, target_to_TF_graph_up = filter_up
        feature_mat_list_up = gen_feature_mat(diff_tab_up, target_to_TF_graph_up, weight_adj, DEGs_up, neg_genes_up, DEGs_TFs_up, neg_TFs_up, use_peak_signal)
    else:
        feature_mat_list_up = [None,None,None]

    filter_down = filter(diff_tab, "down", neg_type, TF_to_target_graph, target_to_TF_graph)
    
    if filter_down is not None:
        DEGs_down, neg_genes_down, DEGs_TFs_down, neg_TFs_down, diff_tab_down, TF_to_target_graph_down, target_to_TF_graph_down = filter_down
        feature_mat_list_down = gen_feature_mat(diff_tab_down, target_to_TF_graph_down, weight_adj, DEGs_down, neg_genes_down, DEGs_TFs_down, neg_TFs_down, use_peak_signal)
    else:
        feature_mat_list_down = [None,None,None]

    return(feature_mat_list_up,feature_mat_list_down)
