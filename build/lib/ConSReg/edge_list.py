'''
Functions for generate DAP-seq edge list and convert edge list
'''

# Author: (Alex) Qi Song <alexsong@vt.edu>

from networkx.convert_matrix import from_pandas_edgelist # Using networkx 2.1
from networkx import DiGraph
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr # This is for loading R library
import pandas as pd 

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
def edge_list_from_peak(peak_tab, gff_file, up_tss, down_tss, up_type, down_type):
    """
    Parameters
    ----------
    peak_tab : pandas dataFrame. Table of peaks
    
    gff_file: string. File name of gff file
    
    up_tss : positions relative to upstream region of TSS. This is used
    for finding nearest gene for each binding site

    down_tss: positions relative to downstream region of TSS. This is used
    for finding nearest gene for each binding site

    up_type: type of binding sites. 'all' or 'intergenic'
    down_type: type of binding sites. 'all' or 'intron' or 'no_intron'. None includes all CREs and 'intron' includes only CREs in introns and 'no_intron' includes all CREs other than intron CREs.
    
    Returns
    -------
    el: pandas dataFrame. edge list of TF->target
    """
    pandas2ri.activate()
    txdb = GF.makeTxDbFromGFF(gff_file) 
    peak_tab = peak_tab.reset_index(drop=True)
    
    # Find the nearest genes for all peaks
    peaks = GR.GRanges(
          seqnames = S4V.Rle(base.unlist(peak_tab['chr'])),
          ranges = IR.IRanges(base.unlist(peak_tab['chrStart']),base.unlist(peak_tab['chrEnd'])),
          strand = S4V.Rle(BG.strand(base.gsub('\\.','*',base.unlist(peak_tab['strand'])))),
          TFID = base.unlist(peak_tab['TFID'])
        )
    
    # annotate peaks
    anno_priority = base.c("Intergenic","Intron","Exon","5UTR","3UTR","Downstream","Promoter")
    anno_peak_tab = pandas2ri.ri2py(base.as_data_frame(meth.slot(CSK.annotatePeak(peaks,TxDb = txdb,genomicAnnotationPriority = anno_priority, level = "gene"),"anno")))
    anno_peak_tab.reset_index(inplace=True, drop=True) # R object is 1-based. Make it 0-based here.
    '''
    Extract binding sites located in the upstream of the TSS of the nearest gene or 
    binding sites located within up_tss bp to the upstream of the TSS of the nearest gene
    '''
    if up_type is 'all':
        up_index = (anno_peak_tab['distanceToTSS'] < 0) & (anno_peak_tab['distanceToTSS'] > -up_tss)
    else:
        up_index = (anno_peak_tab['distanceToTSS'] < 0) & (anno_peak_tab['distanceToTSS'] > -up_tss) & pd.Series(pandas2ri.ri2py(base.grepl("Intergenic",anno_peak_tab['annotation'])),dtype = bool)
    
    '''
    Extract binding sites located in the first intron of nearest gene or
    binding sites located within down_tss bp to the downstream of the TSS of the nearest genes
    '''
    if down_type is 'all':
        down_index = (anno_peak_tab['distanceToTSS'] >= 0) & (anno_peak_tab['distanceToTSS'] <= down_tss)
    elif down_type is 'intron':
        intron_index = pd.Series(pandas2ri.ri2py(base.grepl("Intron",anno_peak_tab['annotation'])),dtype = bool)
        down_index = (anno_peak_tab['distanceToTSS'] >= 0) & (anno_peak_tab['distanceToTSS'] <= down_tss) & intron_index
    else:
        no_intron_index = pd.Series(pandas2ri.ri2py(base.grep("Intron",anno_peak_tab['annotation'],invert = True)),dtype = bool)
        down_index = (anno_peak_tab['distanceToTSS'] >= 0) & (anno_peak_tab['distanceToTSS'] <= down_tss) & no_intron_index

    index = up_index | down_index
    
    # Make edge list
    el = anno_peak_tab.loc[index,('TFID','geneId')]
    el.columns = ('TFID','targetID')
    peak_tab = peak_tab.loc[index,:]
    
    return(el,peak_tab)

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
    TF_to_target_graph = from_pandas_edgelist(df = edge_list,source = 'TFID',target = 'targetID', edge_attr = ["count","signalValue"], create_using = DiGraph())
    target_to_TF_graph = from_pandas_edgelist(df = edge_list,source = 'targetID',target = 'TFID', edge_attr = ["count","signalValue"], create_using = DiGraph())
    return(TF_to_target_graph, target_to_TF_graph)
