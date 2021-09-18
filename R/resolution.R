#' Suggest resolution
#'
#' Suggest the resolution to be used by SCUBI
#'
#' This function suggests the resolution to be used in various SCUBI functions. For each resolution, the averaged number of cells within each square is calculated. To determine the suggested resolution, a plot will be generated to study the change the averaged number of cells with resolution, and the suggested resolution is picked at the elbow location.
#' @param dim1 Numeric vector for the first dimension, such as the first dimension of UMAP or PCA.
#' @param dim2 Numeric vector for the second dimension.
#' @param resolution A numeric vector of resolution values to be tested.
#' @param plot Either TRUE or FALSE. If TRUE, the function will directly generate a plot. If FALSE, the function will return the data frame used for plotting.
#' @import ggplot2
#' @export
#' @return If plot=TRUE, a ggplot object. If plot=FALSE, a vector .
#' @author Wenpin Hou<whou10@@jhu.edu>, Zhicheng Ji
#' @examples 
#' resolution(dim1=rnorm(1000),dim2=rnorm(1000))

resolution <- function(dim1, dim2, resolution=c(1:50)*20/max(diff(range(dim1)),diff(range(dim2))),plot=TRUE) {
  m <- sapply(resolution,function(factor) {
    g <- paste0(round(dim1 * factor), "_", round(dim2 * factor))
    #mean(as.vector(table(g)) > 1)
    mean(as.vector(table(g)))
  })
  if (plot==TRUE) {
    ggplot(data.frame(x=resolution,y=m),aes(x=x,y=y)) + geom_point() + geom_line() + theme_classic() + xlab('Resolution') + ylab('Average number of data points') + coord_cartesian(ylim = c(0,200))
  } else {
    m
  }
}
