"""
Finding weight for each peak.
"""

# Author: (Alex) Qi Song <alexsong@vt.edu>

from ConSReg.utils import bed2Int
from numpy import array,multiply
'''
Get atac-seq weight for each dap-seq peak
'''
def get_weight(peak_tab, atac_file, atac_chr_col, atac_chr_start_col, atac_chr_end_col, atac_signal_col, use_peak_signal):
    """
    Parameters
    ----------
    peak_tab : pandas DataFrame. a merged table of peaks
    atac_file: string. File name of ATAC-seq bed file
    atac_chr_col: column number for atac-seq chromosome information, 0 indexed.
    atac_chr_start_col: column number for atac-seq peak start position, 0 indexed.
    atac_chr_end_col: column number for atac-seq peak end position, 0 indexed.
    use_peak_signal: whether to use dap-seq peak signal.
    
    Returns
    -------
    weight_tab: pandas DataFrame. table of weights 
    """
    
    chrom_int, atac_signal_val = bed2Int(atac_file, atac_chr_col = atac_chr_col, atac_chr_start_col = atac_chr_start_col, atac_chr_end_col = atac_chr_end_col, atac_signal_col = atac_signal_col)
    weight = []
    
    # Iterate over each peak and compute weight by computing overlapping regions
    for row in peak_tab.itertuples():
        found_intervals = chrom_int[row.chr].search(row.chrStart,row.chrEnd)
        found_intervals = sorted(found_intervals,key = lambda x: x[0])
        
        if len(found_intervals) > 0:
            overlap_intervals = []
            for each_interval in found_intervals:
                overlap_intervals.append((max(each_interval[0],row.chrStart),min(each_interval[1],row.chrEnd)))

            end_pos = -1
            overlap_len = []
            atac_signal = []
            for each_f_int,each_o_int in zip(found_intervals, overlap_intervals):
                if each_o_int[1] > end_pos:
                    if each_o_int[0] <= end_pos:
                        overlap_len.append(each_o_int[1] - end_pos)
                    else:
                        overlap_len.append(each_o_int[1] - each_o_int[0] + 1)
                    end_pos = each_o_int[1]
                    atac_signal.append(atac_signal_val[row.chr][(each_f_int[0],each_f_int[1])])
            
            overlap_len = array(overlap_len)
            atac_signal_weight = multiply(array(atac_signal),overlap_len).sum()
            if use_peak_signal:
                weight.append(atac_signal_weight / float(row.chrEnd - row.chrStart + 1) * row.signalValue)
            else:
                weight.append(overlap_len.sum() / float(row.chrEnd - row.chrStart + 1))
        else:
            weight.append(0)
            
    return(weight)
