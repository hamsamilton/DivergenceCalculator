# DivergenceCalculator

## Overview

`DivergenceCalculator` is an R package for calculating transcriptional divergence for RNA-seq experiments. It is a small, focused package, consisting of a handful of functions for preprocessing RNA-seq data and calculating transcriptomic divergence in a way that has been experimentally validated.

## Installation

You can install the development version of DivergenceCalculator from GitHub:

```r
# install.packages("devtools")
devtools::install_github("hamsamilton/DivergenceCalculator")
```

## Features

### Read Depth Standardization

Mapping depth is the number of exons detected by RNA-seq, or the column sum in a typical count table. We have shown that mapping depth is a key confounding technical factor when calculating transcriptional divergence, but downsampling samples to a shared mapping depth addresses this issue. Given a matrix of raw counts, this function downsamples the reads from each sample (column) to level numreads, and returns the table. Samples with a mapping depth below numreads are ignored. We have also shown that there is a key "resolution limit" of mapping depth below which there are too few reads to estimate transcriptional divergence accurately, so be cautious when setting numreads below (1e6 - 3e6).

```r
standardized_counts <- StandardizeReadDepth(
  expression_matrix, 
  numreads = 1e6,
  ntimes = 1
)
```

### Transcriptomic Divergence Calculation

Once you have normalized the read table, the function CalcDivergence is used to estimate transcriptional divergence. This returns a dataframe with 2 columns, the divergence values and the name of each sample (column)

```r
TranscriptionalDivergence <- CalcDivergence(standardized_counts)
```

### Gene Expression Distribution Visualization

Visualizing the gene expression distributions used to calculate transcriptional divergence helps to check to ensure the calculation is valid and to communicate the results with others. This function generates plots of gene expression distributions of all samples with the same axes, so you can see the differences in divergence between samples easily.

```r
# Create distribution plots
plots <- mkGeneDistPlots(standardized_counts)

# Display the first plot
plots[[1]]
```

## Workflow Example

```r
library(DivergenceCalculator)

# Load your expression data
# expression_data <- read.csv("path/to/expression_data.csv", row.names=1)

# Standardize read depth
standardized_counts <- StandardizeReadDepth(expression_data, numreads = 5e5)

# Visualize expression distributions
dist_plots <- mkGeneDistPlots(standardized_counts)

# Calculate divergence metrics
DivergenceDf <- CalcDivergence(standardized_counts)


## Citation

If you use this package in your research, please cite this github page for now:


## License

This package is released under the MIT License.
