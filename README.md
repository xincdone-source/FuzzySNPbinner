FuzzySNPbinner

Fuzzy-Optimized Crosspoints and BIN Mapping for Multi-crop SNP Data.

This R package modifies the original SNPbinner's crosspoint detection and binning by incorporating fuzzy mathematics for adaptive parameters across different crops.

1. Prerequisites
R (>= 4.1)
Python (>= 3.7)

This package relies on Python for its underlying algorithm and uses the reticulate R package to interface between R and Python. Please ensure that Python is installed on your system and that reticulate can access your Python environment.

2. Installation

To install this package directly from GitHub, use the following commands:

# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}

# Install FuzzySNPbinner from GitHub
devtools::install_github("xincdone-source/FuzzySNPbinner")
3. Quick Start

We provide a small built-in test dataset for you to quickly experience the full analysis workflow:

# Load the package
library(FuzzySNPbinner)

# 1. Get the path to the built-in test data (a small test fragment)
test_tsv <- system.file("extdata", "test_data.tsv", package = "FuzzySNPbinner")

# 2. Manually set physical thresholds (adjust according to your population type and expected recombination rate)
# For the built-in small test fragment, we use smaller thresholds for demonstration:
my_min_length <- 10000      # Filter out fragments smaller than 10 kb
my_min_bin_size <- 100000   # Set the minimum bin size to 100 kb

# 3. Run fuzzy-optimized crosspoint detection (HMM probabilities are automatically adapted, physical thresholds are manually controlled)
out_crosp <- tempfile(fileext = ".crosp.csv")
fuzzy_run_crosspoints(
  input_path = test_tsv, 
  output_path = out_crosp, 
  min_length = my_min_length, 
  min_bin_size = my_min_bin_size,
  seq_depth = 10 
)

# 4. Run BIN mapping
out_bin <- tempfile(fileext = ".bin.csv")
fuzzy_run_bins(
  input_path = out_crosp, 
  output_path = out_bin, 
  min_bin_size = my_min_bin_size
)
