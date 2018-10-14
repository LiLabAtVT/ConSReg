'''
Some helper functions
'''
from intervaltree import IntervalTree
from pandas import read_table,unique
from numpy import abs

'''
Read bed file. Convert to a dictionary of interval trees:
{chrName:[[start1,end1],[start2,end2]...],chrName:[[start1,end1],[start2,end2]...]}
'''
def bed2Int(bed_file):
    """
    Parameters
    ----------
    bed_file : string. bed file file name
    
    Returns
    -------
    chrom_int: dictionary of IntervalTree. 
    """
    chrom_int = dict()
    intervals = read_table(bed_file,header=None).iloc[:,:3]
    
    # Construct an interval tree for each chromosome
    for chr_name in unique(intervals.iloc[:,0]):
        int_tuples = list(intervals.loc[intervals.iloc[:,0] == chr_name,1:].itertuples(index = False,name = None))
        chrom_int[chr_name] = IntervalTree().from_tuples(int_tuples)
    return(chrom_int)

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