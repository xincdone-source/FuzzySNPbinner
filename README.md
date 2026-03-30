FuzzySNPbinner
Fuzzy-Optimized Crosspoints and BIN Mapping for Multi-crop SNP Data.
This R package modifies the original SNPbinner's crosspoint detection and binning by incorporating fuzzy mathematics for adaptive parameters across different crops.

1. Prerequisites (环境要求)
R (>= 4.1)

Python (>= 3.7)
本包的底层算法由 Python 实现，通过 R 包 reticulate 进行跨语言调用。请确保您的系统中已安装 Python，并且 reticulate 能够正常发现您的 Python 环境。

2. Installation (安装)
您可以使用 devtools 直接从 GitHub 安装本包：

R
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("xincdone-source/FuzzySNPbinner")

3. Quick Start (快速上手)
我们提供了一个极小的内置测试数据集，您可以通过以下代码快速体验完整的分析流程：

R
library(FuzzySNPbinner)

# 1. 获取内置测试数据路径 (极小测试片段)
test_tsv <- system.file("extdata", "test_data.tsv", package = "FuzzySNPbinner")

# 2. 手动设定物理阈值 (请根据您的群体类型和预期重组率调整)
# 针对内置的极小测试片段，我们使用较小的阈值进行演示：
my_min_length <- 10000      # 过滤掉小于 10 kb 的碎片状态
my_min_bin_size <- 100000   # 划分 BIN 的最小单位为 100 kb

# 3. 运行模糊优化交叉点检测 (HMM概率自动适应，物理阈值手动控制)
out_crosp <- tempfile(fileext = ".crosp.csv")
fuzzy_run_crosspoints(
  input_path = test_tsv, 
  output_path = out_crosp, 
  min_length = my_min_length, 
  min_bin_size = my_min_bin_size,
  seq_depth = 10 
)

# 4. 运行 BIN 划分
out_bin <- tempfile(fileext = ".bin.csv")
fuzzy_run_bins(
  input_path = out_crosp, 
  output_path = out_bin, 
  min_bin_size = my_min_bin_size
)
