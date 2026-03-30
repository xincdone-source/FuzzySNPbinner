# FuzzySNPbinner

Fuzzy-Optimized Crosspoints and BIN Mapping for Multi-crop SNP Data.
This R package modifies the original SNPbinner's crosspoint detection and binning by incorporating fuzzy mathematics for adaptive parameters across different crops.

## 1. Prerequisites (环境要求)
- **R** (>= 4.1)
- **Python** (>= 3.7) 
本包的底层算法由 Python 实现，通过 R 包 `reticulate` 进行跨语言调用。请确保您的系统中已安装 Python，并且 `reticulate` 能够正常发现您的 Python 环境。

## 2. Installation (安装)
您可以使用 `devtools` 直接从 GitHub 安装本包：

```R
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("xincdone-source/FuzzySNPbinner")