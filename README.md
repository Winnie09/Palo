Single-cell Unbiased Visualization with SCUBI
====

## Overview
Single-cell sequencing measures thousands of features for each cell, which is hard to perceive directly. Dimension reduction techniques such as PCA, t-SNE, and UMAP thus become essential to map the cells with high-dimensional information to a low dimensional space. The reduced dimension representations are then visualized with scatterplots to understand the latent structure of the data, to identify cell subpopulations or developmental trajectories, or to compare across different samples. However, we observe significant biases in such visualization procedure which could lead to problematic interpretations of the data in many real applications. The first bias arises when visualizing the gene expression levels or cell identities. The scatterplot only shows a subset of cells plotted last and the cells plotted earlier are masked and unseen. The second bias arises when comparing the cell type compositions from different types of samples. The scatterplot is biased by the unbalanced total numbers of cells across samples.

To this end, we develop a new visualization method, SCUBI, that tackles these biases. To address the first bias, SCUBI partitions the coordinate space into small non-overlapping squares and visualizes the aggregated information of cells within each square. To address the second bias, SCUBI calculates the proportions of cells falling in non-overlapping squares and visualizes the differences of cell proportions across samples. SCUBI more faithfully reveals the real biological information in the single-cell genomic datasets. 

SCUBI comes with six modes designed for six different scenarios in single-cell data visualization. The first three modes are designed for cell-level visualizations and the last three modes are designed for comparisons across different types of samples.

## SCUBI Installation

SCUBI software can be installed via Github.
Users should have R installed on their computer before installing SCUBI. R can be downloaded here: http://www.r-project.org/.
To install the latest version of SCUBI package via Github, run following commands in R:
```{r }
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("Winnie09/SCUBI")
```

## User Manual
Please visit this webpage for the user manual: https://winnie09.github.io/Wenpin_Hou/pages/SCUBI.html

## Contact the Author
Author: Wenpin Hou, Zhicheng Ji

Report bugs and provide suggestions by sending email to:

Maintainer: Wenpin Hou (whou10@jhu.edu)

Or open a new issue on this Github page
