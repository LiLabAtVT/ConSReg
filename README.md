# ConSReg
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
Condition-specific regulations

# Getting Started
## 1. Installation
### 1.1 Required packages
#### 1.1.1 Python
- python = 2.7
- numpy >= 1.9.0
- scipy >= 1.1.0
- pandas == 0.21.1
- joblib >= 0.12.5
- rpy2==2.8.6
- networkx >= 2
- sklearn >= 0.18.1
- intervaltree == 2.1.0
#### 1.1.2 R
- ChIPSeeker == 1.16.1
- CoReg == 1.0.1
- gglasso == 1.4
- RRF == 1.9
- R >= 3.5.1
### 1.2 Python installation
ConSReg can be installed using pip:
```shell
pip install --user ConSReg
```
### 1.3 R installation
ConSReg requires several R packages: `ChIPseeker`, `CoReg`, `gglasso` and `RRF`.

To install `ChIPSeeker` from bioconductor, type the following commands in R environment:
```R
source("https://bioconductor.org/biocLite.R")
biocLite("ChIPseeker")
```
Please refer to the instructions described [here](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html) for more details.

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

## 2. Sample datasets
Sample datasets can be found in `data` folder.

## 3. Analysis
We provide code for analyzing the sample datasets in two jupyter notebooks located in the root folder of this project: **bulk_analysis.ipynb** (for bulk RNA-seq data) and **single_cell_analysis.ipynb** (for single cell RNA-seq data).

## 4. Publication
ConSReg is currently in review at **Genome Research**. We will soon provide a pre-print version of our manuscript. 

Qi Song, Jiyoung Lee, Shamima Akter, Ruth Grene, Song Li.  Accurate prediction of condition-specific regulatory maps in Arabidopsis using integrated genomic data (in review)
