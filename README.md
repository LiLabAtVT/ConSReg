# ConSReg 1.1.4
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
Condition-specific regulations
- [ConSReg 1.1.4](#consreg-114)
- [Getting Started](#getting-started)
  * [1. Installation](#1-installation)
    + [1.1 Required packages](#11-required-packages)
      - [1.1.1 Python](#111-python)
      - [1.1.2 R](#112-r)
    + [1.2 Easy installation by Anaconda (recommended)](#12-easy-installation-by-anaconda--recommended-)
    + [1.3 Manual installation (Skip this section if 1.2 is successful)](#13-manual-installation--skip-this-section-if-12-is-successful-)
    + [1.3.1 R installation](#131-r-installation)
      - [install R](#install-r)
      - [install R packages](#install-r-packages)
    + [1.3.2 Python installation](#132-python-installation)
  * [2. Sample datasets](#2-sample-datasets)
  * [3. Analysis](#3-analysis)
  * [4. Publication](#4-publication)

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
### 1.2 Easy installation by Anaconda (recommended)
Since ConSReg is dependent on both Python and R packages, we recommend installing ConSReg by Anaconda to easily set up the running environment. You may retrive Anaconda from [here](https://www.anaconda.com/) and install the version corresponding to your OS.  
Once Anaconda is installed in your OS, run the following commands to create an new environment and install ConSReg and all its dependencies into the new environment:
```bash
conda create -y -n consreg python=3.6 # The new environment name is 'consreg'. You may use other name instead.
conda activate consreg
conda install -y -c bioconda --no-channel-priority bioconductor-chipseeker
conda install -y --no-channel-priority r-base r-essentials
conda install -y --no-channel-priority -c conda-forge r-gglasso r-rrf r-devtools
pip install ConSReg
```
Then ConSReg environment can be activated by `conda activate consreg` and disabled by `conda deactivate`
### 1.3 Manual installation (Skip this section if 1.2 is successful)
### 1.3.1 R installation
#### install R
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
#### install R packages
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
### 1.3.2 Python installation
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

## 3. Analysis
We provide code for analyzing the sample datasets in two jupyter notebooks located in the root folder of this project: **bulk_analysis.ipynb** (for bulk RNA-seq data) and **single_cell_analysis.ipynb** (for single cell RNA-seq data).

## 4. Publication
Please cite the followint paper if you use ConSReg in your research:  
Qi Song, Jiyoung Lee, Shamima Akter, Ruth Grene, Song Li. "Prediction of condition-specific regulatory genes using machine learning." Nucleic acids research 48.11 (2020): e62-e62.
