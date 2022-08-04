# ConSReg 1.1.7
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
Condition-specific regulations
- [Getting Started](#getting-started)
  * [1. Installation](#1-installation)
    + [1.1 Required packages](#11-required-packages)
      - [1.1.1 Python](#111-python)
      - [1.1.2 R](#112-r)
    + [1.2 Installation choices](#12-installation-choices)
    + [1.3 Installation instructions](#13-installation-instructions)
      - [1.3.1 installing ConSReg by Anaconda (with environment.yml)](#131-installing-consreg-by-anaconda--with-environmentyml-)
      - [1.3.2 installing ConSReg by Aanconda (manual installation)](#132-installing-consreg-by-aanconda--manual-installation-)
      - [1.3.3 installing ConSReg by singularity image](#133-installing-ConSReg-by-singularity-image)
      - [1.3.4 Manual installation for all dependencies](#134-manual-installation-for-all-dependencies)
        * [install R](#install-r)
        * [install R packages](#install-r-packages)
        * [Python installation](#python-installation)
  * [2. Sample datasets](#2-sample-datasets)
  * [3. Input files and format](#3-input-files-and-format)
  * [4. Analysis steps](#4-analysis-steps)
  * [4. Publication](#4-publication)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>
<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>
# Getting Started
## 1. Installation
### 1.1 Required packages
#### 1.1.1 Python
- python = 3.6
- numpy == 1.16.2
- scipy == 1.1.0
- pandas == 0.21.1
- joblib >= 0.11
- rpy2==2.8.6
- networkx >= 2
- sklearn >= 0.19.1
- intervaltree == 2.1.0
#### 1.1.2 R
- ChIPSeeker == 1.16.1
- CoReg == 1.0.1
- gglasso == 1.4
- RRF == 1.9
- R >= 3.5.1
### 1.2 Installation choices
We provide four ways to install ConSReg and its dependencies:
1. installing ConSReg by Anaconda (with environment.yml)
2. install ConSReg by Anaconda (manual installation)
3. singularity image
4. manual installation for all dependecies
Users only need to choose one of them that works. We would suggest to starting with `1` or `3`, as these two options are very easy and straightforward.
### 1.3 Installation instructions
#### 1.3.1 installing ConSReg by Anaconda (with environment.yml)
Since ConSReg is dependent on both Python and R packages, we recommend installing ConSReg by Anaconda to easily set up the running environment. You may retrive Anaconda from [here](https://www.anaconda.com/) and install the version corresponding to your OS.   
Once Anaconda is installed in your os, create a new conda environment using the [`environment.yml`](https://raw.githubusercontent.com/LiLabAtVT/ConSReg/master/environment.yml) file in this repository:
```shell
conda env create -f environment.yml
```
Alternatively, you can type in
```shell
conda env create -f https://raw.githubusercontent.com/LiLabAtVT/ConSReg/master/environment.yml
```
This will create a new conda environment named `consreg` which contains `ConSReg` package and all its dependecies. You can then type in `conda activate consreg` to activate this environment or `conda deactivate` to deactivate this environment. For more information, please refer to official documentation of Ananconda.
#### 1.3.2 installing ConSReg by Aanconda (manual installation)
Alternatively, after Anaconda is installed in your OS (see **1.3.1**), run the following commands to create an new environment and install ConSReg and all its dependencies into the new environment:
```bash
conda create -y -n consreg python=3.6 # The new environment name is 'consreg'. You may use other name instead.
conda activate consreg
conda install -y -c bioconda --no-channel-priority bioconductor-chipseeker
conda install -y --no-channel-priority r-base r-essentials
conda install -y --no-channel-priority -c conda-forge r-gglasso r-rrf r-devtools
pip install ConSReg
```
Then ConSReg environment can be activated by `conda activate consreg` and disabled by `conda deactivate`
#### 1.3.3 installing ConSReg by singularity image
Singularity is a container system which creasts lightweight container that hosts all system dependencies and environment for a given software package. Users may simply pull container image from the cloud and then run the program inside the container without having to installing it in their own machines. You may install Singularity following the instructions [here](https://docs.sylabs.io/guides/3.5/user-guide/quick_start.html).  
To install ConSReg using Singularity, you may simply pull our prebuilt singularity image for ConSReg：
```bash
singularity pull -U library://alexsong0374/consreg/consreg_singularity_ubuntu20.04
```
This will create an image file called "consreg_singularity_ubuntu20.04_latest.sif" locally. To run python environment with ConSReg installed, you may run the local container by:
```bash
singularity run consreg_singularity_ubuntu20.04_latest.sif python3
```
Alternatively, you may run jupyter notebook inside the container:
```bash
singularity run consreg_singularity_ubuntu20.04_latest.sif jupyter notebook
```
Additonally, we provide the singularity definination file (**consreg_singularity.def**) in this repository for rebuilding the container for ConSReg. You may rebuild the container by yourself and add other packages you want.
#### 1.3.4 Manual installation for all dependencies
If all of the above steps fail, you may manually install all dependencies by following the steps below.
##### install R
If R is not already installed， you may follow these steps to build R from source code. Otherwise, you may skip this section and start from 1.2.2

First, disable any conda environment, if there is an active one.
```shell
conda deactivate
```
Download R source code from CRAN (https://cran.r-project.org/). You may use any version you like. It is recommended to use R version > 3.0.0. This ensures that rpy2 works correctly with R.
```shell
# Download R 3.6.1
wget https://cran.r-project.org/src/base/R-3/R-3.6.1.tar.gz
```
Decompress the downloaded file
```shell
tar -zvxf R-3.6.1
```
In the decompressed folder, configure R by:
```shell
./configure prefix=path_to_install_R --enable-R-shlib
```
`--prefix=` specifies a writeable directory to install R into. `--enable-R-shlib` flag was added to build R shared libraries.

In the decompressed folder, compile R
```shell
make
```
Install R into the specified directory:
```shell
make install
```
Add a line to ~/.bashrc to tell the OS where to look for R 
```shell
export PATH=path_to_R_bin_directory:$PATH
```
Add the following line to ~/.bashrc. This is for telling rpy2 where to look for dynamic libraries.
```shell
export LD_LIBRARY_PATH=/home/alexsong/R/3.6.1/lib64/R/lib:$LD_LIBRARY_PATH
```
Apply the changes to environment variables `PATH` and `LD_LIBRARY_PATH`:
```shell
source ~/.bashrc
```
##### install R packages
ConSReg requires several R packages: `ChIPseeker`, `CoReg`, `gglasso` and `RRF`.

It is recommended to deactivate any conda environment when installing R packages, as it may add the environment-specific path which may fail the installation. If any conda environment is active, you may deactivate it by:
```shell
conda deactivate
```
To install `ChIPSeeker` from bioconductor, type the following commands in R (for R 3.6 or higher version):
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ChIPseeker")
```
For older version of R, type the following commands in R:
```R
source("https://bioconductor.org/biocLite.R")
biocLite("ChIPseeker")
```
Please refer to the instructions described [here](https://www.bioconductor.org/packages/release/bioc/html/ChIPseeker.html) for more details.

To install `CoReg` pakcage from GitHub, type the following commands in R environment:
```R
install.packages("devtools")
library(devtools)
install_github("LiLabAtVT/CoReg")
```
Please refer to the GitHub page of `CoReg` project for more details: 
[link](https://github.com/LiLabAtVT/CoReg)

To install `gglasso` package from CRAN, type the following commands in R environment:
```R
install.pacakges("gglasso")
```
Please refer to the link [here](https://cran.r-project.org/web/packages/gglasso/index.html) for more details.

To install `RRF` package from CRAN, type the following commands in R environment:
```R
install.pacakges("RRF")
```
Please refer to the link [here](https://cran.r-project.org/web/packages/RRF/index.html) for more details.
##### Python installation
ConSReg can be installed by pip:
```shell
pip install ConSReg
```
Sometime rpy2 may throw out error message when imported in Python. This problem may arise because rpy2 was built with the R version that is different from the one it is linked to when imported in Python. To fix this, you may remove rpy2 package then reinstall it with 'no-cache-dir' flag:

```shell
pip install ConSReg --no-cache-dir
```

Alternatively, you may want to install ConSReg in development mode to be able to edit the package by yourself. To do so, simply `git clone` this repository and then under the directory that contains `setup.py`, type in:
```shell
pip install -e .
```
## 2. Sample datasets
Sample datasets can be found in `data` folder.

## 3. Input files and format
ConSReg can take three types of genomic data as inputs:
1. Open chromatin position information from ATAC-seq data, represented as a bed file, which should include three columns: Column1 represents chromosome name (e.g, chr1, chr2); Column2 represents start position of an open chromatin region; Column3 represents end position of an open chromatin region (See also the sample data files under `data/atac_seq_all_peaks` in this repository)

2. TF binding location information from DAP-seq/ChIP-seq data, represented as narrowPeak files, which should contain the same types of columns stated for ATAC-seq data (See also the sample data files under `data/dap_seq_all_peaks` in this repository). Note that the binding information for different TFs should be in separate files and each file is named by the name of the corresponding TF. For example, file “AT1G01060.narrowPeak” contains only the TF binding locations for the TF “AT1G01060”.

3. Differentially expressed gene (DEG) information from RNA-seq/microarray data, represented as **comma-separated values (CSV) table files**, which should contain four columns: The **first column** represents gene name; a column named **“baseMean”** that represents mean expression for the gene; a column named **“log2FoldChange”** that represents log2-scaled fold change; a column named **“padj”** that represents the adjusted p-values from statistical test for differential expressions. This type of DEG information usually could be generated from DESeq2 package (See DESeq2 page for more information: https://bioconductor.org/packages/release/bioc/html/DESeq2.html). ConSReg can take multiple DEG tables as inputs and each table corresponds to a DEG analysis from one experiment.

## 4. Analysis steps
We provide code for analyzing the sample datasets in two jupyter notebooks located in the root folder of this project: **bulk_analysis.ipynb** (for bulk RNA-seq data) and **single_cell_analysis.ipynb** (for single cell RNA-seq data).

## 4. Publication
Please cite the followint paper if you use ConSReg in your research:  
Qi Song, Jiyoung Lee, Shamima Akter, Ruth Grene, Song Li. "Prediction of condition-specific regulatory genes using machine learning." Nucleic acids research 48.11 (2020): e62-e62.
