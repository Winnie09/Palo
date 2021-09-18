#' Suggest resolution for multi-sample regression
#'
#' Suggest the resolution to be used by SCUBI multi-sample regression mode
#'
#' This function suggests the resolution to be used in SCUBI function of multi-sample regression. For each resolution, the averaged number of cells within each square is calculated. To determine the suggested resolution, a plot will be generated to study the change the averaged number of cells with resolution, and the suggested resolution is picked at the elbow location.
#' @param dim1 Numeric vector for the first dimension, such as the first dimension of UMAP or PCA.
#' @param dim2 Numeric vector for the second dimension.
#' @param sample Character vector of cell-sample assignment (each cell belong to which sample).
#' @param design A numeric data.frame or matrix indicating the covariates of the samples. It should not contain the intercept term. The row names of the matrix should contain the unique values in 'sample' variable.
#' @param resolution A numeric vector of resolution values to be tested.
#' @param plot Either TRUE or FALSE. If TRUE, the function will directly generate a plot. If FALSE, the function will return the data frame used for plotting.
#' @import ggplot2
#' @export
#' @return If plot=TRUE, a ggplot object. If plot=FALSE, a vector .
#' @author Wenpin Hou<whou10@@jhu.edu>, Zhicheng Ji
#' @examples 
#' resolution_sample_expression(dim1=rnorm(1000),dim2=rnorm(1000),sample=rep(c('a','b'),500),design=data.frame(c(0,1),row.names=c('a','b')))

resolution_sample_expression <- function(dim1, dim2, sample, design, resolution=c(1:50)*20/max(diff(range(dim1)),diff(range(dim2))),plot=TRUE) {
  design <- as.matrix(design)
  m <- 1 - sapply(resolution,function(factor) {
    g <- paste0(round(dim1 * factor), "_", round(dim2 * factor))
    mean(tapply(design[sample,1],list(g),function(i) length(unique(i))==1))
  })
  if (plot==TRUE) {
    ggplot(data.frame(x=resolution,y=m),aes(x=x,y=y)) + geom_point() + geom_line() + theme_classic() + xlab('Resolution') + ylab('Proportion of valid squares') + coord_cartesian(ylim = c(0,1))
  } else {
    data.frame(resolution=resolution,percentage=m)
  }
}
