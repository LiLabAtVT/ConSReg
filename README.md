# ConSReg
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
Condition-specific regulations

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
### 1.2 R installation
#### 1.2.1 install R
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
#### 1.2.2 install R packages
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
### 1.3 Python installation
We recommend the users to create a new Python environment for ConSReg using Anaconda and install ConSReg in this environment. This can guarantee ConSReg work with correct dependencies. However, installing ConSReg without conda environment is also welcome.

To create a new environment using conda:
```shell
conda create --name consreg python=3.6
```
Activate the new environment
```shell
conda activate consreg
```
Then ConSReg can be then installed using pip:
```shell
pip install ConSReg
```
Sometime rpy2 may throw out error message when imported in Python. This problem may arise because rpy2 was built with the R version that is different from the one it is linked to when imported in Python. To fix this, you may remove rpy2 package then reinstall it with 'no-cache-dir' flag:

```shell
pip install ConSReg --no-cache-dir
```
You may refer to https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html for more information about installation and usage of Anaconda.

## 2. Sample datasets
Sample datasets can be found in `data` folder.

## 3. Analysis
We provide code for analyzing the sample datasets in two jupyter notebooks located in the root folder of this project: **bulk_analysis.ipynb** (for bulk RNA-seq data) and **single_cell_analysis.ipynb** (for single cell RNA-seq data).

## 4. Publication
ConSReg is currently in review at **Genome Research**. We will soon provide a pre-print version of our manuscript. 

Qi Song, Jiyoung Lee, Shamima Akter, Ruth Grene, Song Li.  Accurate prediction of condition-specific regulatory maps in Arabidopsis using integrated genomic data (in review)
