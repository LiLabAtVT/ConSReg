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
    _ = map(lambda x: TFs.update(target_to_TF_graph.neighbors(x)), genes)
    return(TFs)

'''
Filter TF->target and target->TF graph. Modify diff_tab. Assign log2FoldChange = 0 
to 'baseMean' of negative genes
'''
def filter(diff_tab, direction, TF_to_target_graph, target_to_TF_graph):
    '''
    Parameters
    ----------
    diff_tab : pandas dataframe. A table of fold change, basemean padj and gene names,
    following the format of deseq2 output
    
    direction: string. 'up' for up-regulated 'down' for down-regulated
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
        neg_genes = diff_tab.loc[(diff_tab.loc[:,"baseMean"] == 0),].index
        diff_tab.loc[neg_genes,"log2FoldChange"] = 1
    else:
        DEGs = set(diff_tab.loc[(diff_tab.loc[:,"padj"] < 0.05) & (diff_tab.loc[:,"log2FoldChange"] < 0),].index)
        neg_genes = diff_tab.loc[(diff_tab.loc[:,"baseMean"] == 0),].index
        diff_tab.loc[neg_genes,"log2FoldChange"] = 1
        
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
def gen_feature_mat(diff_tab, target_to_TF_graph, weight_adj, DEGs, neg_genes, DEG_TFs, neg_TFs):
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
    feature_mat : a pandas datarame. feature matrix
    '''
    
    row_genes = DEGs | neg_genes
    col_genes = DEG_TFs | neg_TFs # Make sure there are no duplicate TFs
    
    # Generate the feature matrix by multiply subset of adj matrix and fold change column
    if weight_adj is not None:
        adj = weight_adj.loc[row_genes,col_genes]
    else:
        adj = to_pandas_adjacency(Target2TFGraph.subgraph(row_genes | col_genes),weight = 'count').loc[row_genes,col_genes]
    fc = diff_tab.loc[col_genes,'log2FoldChange']
    feature_mat = pd.DataFrame(adj.values * fc.values.reshape(1,fc.shape[0]),index = row_genes, columns = col_genes)
    
    # Re-order the featureMat
    feature_mat = feature_mat.loc[list(DEGs)+list(neg_genes),]
    
    # Add label column
    label = pd.DataFrame(feature_mat.shape[0]/2*[1] + feature_mat.shape[0]/2*[0], index = feature_mat.index, columns = ['label'])
    feature_mat = pd.concat([label, feature_mat], axis = 1)
    return(feature_mat)

'''
Get up-regulated (UR) feature matrix and down-regulated (DR) feature matrix
'''
def get_all_feature_mat(diff_tab, TF_to_target_graph, target_to_TF_graph, weight_adj):
    '''
    Parameters
    ----------
    diff_tab : pandas dataframe. A table of fold change, basemean padj and gene names,
    following the format of deseq2 output
    
    TF_to_target_graph : networkx DiGraph object. Digraph object of TF->target graph
    target_to_TF_graph : networkx DiGraph object. Digraph object of target->TF graph
    weight_adj -- weight adjacency matrix computed from open chromatin region data (e.g. ATAC-seq)
    
    Returns
    -------
    feature_mat_up : A pandas dataframe. UR feature matrix. The first column is the class label and other columns are values for each TF(feature)
    feature_mat_down : A pandas dataframe. DR feature matrix. The first column is the class label and other columns are values for each TF(feature)
    '''
    filter_up = filter(diff_tab,"up", TF_to_target_graph, target_to_TF_graph)
    
    if filter_up is not None:
        DEGs_up, neg_genes_up, DEGs_TFs_up, neg_TFs_up, diff_tab_up, TF_to_target_graph_up, target_to_TF_graph_up = filter_up
        feature_mat_up = gen_feature_mat(diff_tab_up, target_to_TF_graph_up, weight_adj, DEGs_up, neg_genes_up, DEGs_TFs_up, neg_TFs_up)
    else:
        feature_mat_up = None

    filter_down = filter(diff_tab, "down", TF_to_target_graph, target_to_TF_graph)
    
    if filter_down is not None:
        DEGs_down, neg_genes_down, DEGs_TFs_down, neg_TFs_down, diff_tab_down, TF_to_target_graph_down, target_to_TF_graph_down = filter_down
        feature_mat_down = gen_feature_mat(diff_tab_down, target_to_TF_graph_down, weight_adj, DEGs_down, neg_genes_down, DEGs_TFs_down, neg_TFs_down)
    else:
        feature_mat_down = None

    return(feature_mat_up,feature_mat_down)
