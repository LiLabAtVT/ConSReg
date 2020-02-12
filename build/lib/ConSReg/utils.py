'''
Some helper functions
'''
from intervaltree import IntervalTree
from pandas import read_table,unique,concat,DataFrame
from numpy import abs
'''
Read bed file. Convert to a dictionary of interval trees:
{chrName:[[start1,end1],[start2,end2]...],chrName:[[start1,end1],[start2,end2]...]}
and a dictonary of signal values
{chrName:{(start1,end1):signalValue...},...chrName:chrName:{(start1,end1):signalValue...}}
'''
def bed2Int(atac_file, atac_chr_col, atac_chr_start_col, atac_chr_end_col, atac_signal_col):
    """
    Parameters
    ----------
    atac_file: string. File name of ATAC-seq bed file
    atac_chr_col: column number for atac-seq chromosome information, 0 indexed.
    atac_chr_start_col: column number for atac-seq peak start position, 0 indexed.
    atac_chr_end_col: column number for atac-seq peak end position, 0 indexed.
    atac_signal_col: column number for atac-seq peak signal, 0 indexed.
    
    Returns
    -------
    chrom_int: dictionary of IntervalTree.
    signal_val: dictionary
    """
    chrom_int = dict()
    signal_val = dict()
    
    if atac_signal_col is not None:
        intervals = read_table(atac_file).iloc[:,[atac_chr_col, atac_chr_start_col, atac_chr_end_col, atac_signal_col]]
        intervals.columns = ["chr","chrStart","chrEnd","signalValue"]
    else:
        intervals = read_table(atac_file).iloc[:,[atac_chr_col, atac_chr_start_col, atac_chr_end_col]]
        
        # If not provided, ATAC-seq signal is considered to be 1 for all peaks.
        intervals = concat([intervals,DataFrame([1] * intervals.shape[0], index = intervals.index)], axis = 1)
        intervals.columns = ["chr","chrStart","chrEnd","signalValue"]
    
    intervals.loc[:,'chr'] = intervals.loc[:,'chr'].astype(str)
    
    # Construct an interval tree for each chromosome
    # And a dictionary to store signal values of ATAC-seq
    for chr_name in unique(intervals.iloc[:,0]):
        idx = intervals.iloc[:,0] == chr_name
        int_tuples = list(intervals.loc[idx,["chrStart","chrEnd"]].itertuples(index = False,name = None))
        signalValues = intervals.loc[idx,"signalValue"].tolist()
        chrom_int[chr_name] = IntervalTree().from_tuples(int_tuples)
        signal_val[chr_name] = {int_tuple:val for int_tuple,val in zip(int_tuples,signalValues)}
    return(chrom_int,signal_val)
'''
Get the maximum absolute value from a feature matrix
'''
def get_max_abs(features):
    '''
    Parameters:
    ----------
    features : a numpy 2d array of features
    
    Returns:
    ----------
    max_abs : a float number. Maximum absolute fold change value
    '''
    if features[features > 0].size > 0:
        max_abs_pos = abs(features[features > 0].max())
    else:
        max_abs_pos = 0

    if features[features < 0].size > 0:
        max_abs_neg = abs(features[features < 0].min())
    else:
        max_abs_neg = 0
    
    max_abs = max(max_abs_pos, max_abs_neg)
    return(max_abs)