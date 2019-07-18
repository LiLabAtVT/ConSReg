"""
Finding weight for each peak.
"""

# Author: (Alex) Qi Song <alexsong@vt.edu>

from utils import bed2Int
from numpy import array

'''
Get weight for each peak
'''
def get_weight(peak_tab, chrom_bed_file):
    """
    Parameters
    ----------
    peak_tab : pandas DataFrame. a merged table of peaks
    
    chrom_bed_file: string. File name of chromatin feature (ATAC-seq bed file)
    
    Returns
    -------
    weight_tab: pandas DataFrame. table of weights 
    """
    
    chrom_int = bed2Int(chrom_bed_file)
    weight = []
    
    # Iterate over each peak and compute weight by computing overlapping regions
    for row in peak_tab.itertuples():
        found_intervals = chrom_int[row.chr].search(row.chrStart,row.chrEnd)
        found_intervals = sorted(found_intervals,key = lambda x: x[0])
        
        if len(found_intervals) > 0:
            for i,each_interval in enumerate(found_intervals):
                found_intervals[i] = [max(each_interval[0],row.chrStart),min(each_interval[1],row.chrEnd)]

            end_pos = found_intervals[0][1]
            overlap_len = found_intervals[0][1] - found_intervals[0][0] + 1

            for each_interval in found_intervals:
                if each_interval[1] > end_pos:
                    if each_interval[0] <= end_pos:
                        overlap_len += each_interval[1] - end_pos        
                    else:
                        overlap_len += each_interval[1] - each_interval[0] + 1   
                    end_pos = each_interval[1]
            weight.append(float(overlap_len) / float(row.chrEnd - row.chrStart + 1))

        else:
            weight.append(0)
            
    return(weight)