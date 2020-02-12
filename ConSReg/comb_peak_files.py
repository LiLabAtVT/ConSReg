'''
Combine peaks files in a folder
'''

# Author: (Alex) Qi Song <alexsong@vt.edu>

import re
import pandas as pd
import numpy as np

def comb_peak_files(peak_file_list, dap_chr_col, dap_chr_start_col, dap_chr_end_col, dap_strand_col, dap_signal_col):
    """
    Parameters
    ----------
    peak_file_list : list of strings. File names for all peak files. File should be named as TFID.narrowPeak
    dap_chr_col: int: column number for dap-seq chromosome information, 0 indexed.
    dap_chr_start_col: int, column number for dap-seq peak start position, 0 indexed.
    dap_chr_end_col: int, column number for dap-seq peak end position, 0 indexed.
    dap_strand_col: int/None, column number for dap-seq peak strand information, 0 indexed.
    dap_signal_col: int/NOne, column number for dap-seq peak score, 0 indexed.
    
    Returns
    -------
    peak_tab_final: pandas dataFrame. Merged peak table, with TFID added to the last column
    """
    peak_tab_final = []
    
    # Iterate through peak files
    for peak_file in peak_file_list:
        
        # Read bed files as pandas dataframe
        peak_tab = pd.read_table(peak_file)
        
        # Get TFID
        TFID = re.match(".*?([^/]*)\.narrowPeak$",peak_file).group(1)
        TFID_col = pd.DataFrame(np.array([TFID]*peak_tab.shape[0]).reshape(peak_tab.shape[0],1))
        peak_tab = pd.concat([peak_tab,TFID_col],axis = 1)
        
        if dap_strand_col is not None:
            dap_strand = peak_tab.iloc[:,[dap_strand_col]]
        else:
            dap_strand = pd.DataFrame(['.']*peak_tab.shape[0])
            
        if dap_signal_col is not None:
            dap_signal = peak_tab.iloc[:,[dap_signal_col]]
            
            # Normalize the peak signal values within each library
            dap_signal /= dap_signal.max()
        else:
            # If there are no signal values provided, signal values will be assigned as ones
            dap_signal = pd.DataFrame([1]*peak_tab.shape[0])
            
        peak_tab = pd.concat([peak_tab.iloc[:,[dap_chr_col, dap_chr_start_col, dap_chr_end_col]],dap_strand,dap_signal,peak_tab.iloc[:,-1]],axis = 1)
        
        # Put column names
        peak_tab.columns = ("chr","chrStart","chrEnd","strand","signalValue","TFID")
        peak_tab_final.append(peak_tab)
    
    # Merge all dataframes
    peak_tab_final = pd.concat(peak_tab_final,ignore_index=True)
    peak_tab_final.loc[:,'chr'] = peak_tab_final.loc[:,'chr'].astype(str)
    return(peak_tab_final)