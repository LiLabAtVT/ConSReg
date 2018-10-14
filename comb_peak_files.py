'''
Combine peaks files in a folder
'''

# Author: (Alex) Qi Song <alexsong@vt.edu>

import os
import re
import pandas as pd
from numpy import array

def comb_peak_files(peak_file_folder):
    """
    Parameters
    ----------
    peak_file_folder : string. Folder containing all the peak files. File should be named as TFID.bed
    
    Returns
    -------
    peak_tab_final: pandas dataFrame. Merged peak table, with TFID added to the last column
    """
    os.chdir(peak_file_folder)
    peak_tab_final = []
    
    # Iterate through peak files
    for file in os.listdir("."):
        if os.path.isfile(file):
            
            # Read bed files as pandas dataframe
            peak_tab = pd.read_table(file)
            TFID = re.match("(.*)\..*?$",file).group(1)
            TFID_col = pd.DataFrame(array([TFID]*peak_tab.shape[0]).reshape(peak_tab.shape[0],1))
            peak_tab = pd.concat([peak_tab,TFID_col],axis = 1)
            peak_tab.columns = ("chr","chrStart","chrEnd","name","score","strand","signalValue","pValue","qValue","peak","TFID")
            peak_tab_final.append(peak_tab)
    
    # Merge all dataframes
    peak_tab_final = pd.concat(peak_tab_final,ignore_index=True)
    return(peak_tab_final)