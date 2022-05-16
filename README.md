Palo: Spatially-aware color palette optimization for single-cell and spatial data
====

## Overview
In the exploratory data analysis of single-cell or spatial genomic data, single cells or spatial spots are often visualized using a two-dimensional plot where each cluster is marked with a different color. With tens of clusters, current visualization methods will often result in visually similar colors assigned to spatially neighbouring clusters, making it hard to distinguish and identify the boundary between clusters. To address this issue, we developed `Palo` that optimizes the color palette assignment for single-cell and spatial data in a spatially aware manner. `Palo` identifies pairs of clusters that are spatially neighbouring to each other, and assigns visually different colors to those neighbouring clusters. We demonstrate that `Palo` results in better visualization in real single-cell and spatial genomic datasets. 


## Palo Installation

`Palo` software can be installed via Github.
Users should have R installed on their computer before installing `Palo. R` can be downloaded here: http://www.r-project.org/.
To install the latest version of `Palo` package via Github, run following commands in R:
```{r }
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("Winnie09/Palo")
```

## User Manual
Please visit this webpage for the user manual: https://winnie09.github.io/Wenpin_Hou/pages/Palo.html

## Contact the Author
Author: Wenpin Hou, Zhicheng Ji

Report bugs and provide suggestions by sending email to:

Maintainer: Wenpin Hou (whou10@jhu.edu)

Or open a new issue on this Github page

