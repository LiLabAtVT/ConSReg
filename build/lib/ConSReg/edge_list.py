'''
Functions for generate DAP-seq edge list and convert edge list
'''

# Author: (Alex) Qi Song <alexsong@vt.edu>

from networkx.convert_matrix import from_pandas_edgelist # Using networkx 2.1
from networkx import DiGraph
from rpy2.robjects import pandas2ri
from pandas import read_csv,DataFrame
from rpy2.robjects.packages import importr # This is for loading R library

# import R libraries
base = importr('base')
GF = importr('GenomicFeatures')
GR = importr('GenomicRanges')
BG = importr("BiocGenerics")
IR = importr('IRanges')
S4V = importr('S4Vectors')
CSK = importr('ChIPseeker')
meth = importr('methods')
    
"""
Process the peak table and use ChIPseeker R package to find
nearest genes. Nearest genes were considered as target gene
affected by the TF. TF ID is defined in the 'locusID' column
"""
def edge_list_from_peak(peak_tab, gff_file, upTSS = 3000, downTSS = 500):
    """
    Parameters
    ----------
	peak_tab : pandas dataFrame. Table of peaks
	
    gff_file: string. File name of gff file
	
    upTSS : positions relative to upstream region of TSS
    
    downTSS: positions relative to downstream region of TSS
    
    Returns
    -------
    el: pandas dataFrame. edge list of TF->target
	"""
    pandas2ri.activate()
    txdb = GF.makeTxDbFromGFF(gff_file)

    # Find the nearest genes for all peaks
    peaks = GR.GRanges(
          seqnames = S4V.Rle(base.unlist(peak_tab['chr'])),
          ranges = IR.IRanges(base.unlist(peak_tab['chrStart']),base.unlist(peak_tab['chrEnd'])),
          strand = S4V.Rle(BG.strand(base.gsub('\\.','*',base.unlist(peak_tab['strand'])))),
          TFID = base.unlist(peak_tab['TFID'])
        )
    
    # annotate peaks
    anno_peak_tab = pandas2ri.ri2py(base.as_data_frame(meth.slot(CSK.annotatePeak(peaks,TxDb = txdb),"anno")))
    index = (anno_peak_tab['distanceToTSS'] > -upTSS) & (anno_peak_tab['distanceToTSS'] <= downTSS)
    
    # Make edge list
    el = anno_peak_tab.loc[index,('TFID','geneId')]
    el.columns = ('TFID','targetID')
    el.index -= 1 # R starts the index with 1. To make it consistent with python, minus the index
    return(el)

# Convert edgelist dataframe to networkx directed graph object
def edge_list_to_nx_graph(edge_list):
    """
    Parameters
    ----------
	edge_list : pandas dataFrame. edge list of TF->target
    
    Returns
    -------
    TF_to_target_graph : networkx DiGraph object. Digraph object of TF->target graph
    
    target_to_TF_grap : networkx DiGraph object. Digraph object of target-> graph
	"""
    TF_to_target_graph = from_pandas_edgelist(df = edge_list,source = 'TFID',target = 'targetID', edge_attr = "count", create_using = DiGraph())
    target_to_TF_graph = from_pandas_edgelist(df = edge_list,source = 'targetID',target = 'TFID', edge_attr = "count", create_using = DiGraph())
    return(TF_to_target_graph, target_to_TF_graph)