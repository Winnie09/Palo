#' SCUBI for continuous variable
#'
#' SCUBI for plotting a continuous variable
#'
#' This function partitions the whole coordinate space into multiple non-overlapping squares, and calculates and visualizes the averaged value of a user-defined continuous variable for cells falling in each square It takes as input the coordinates of the cells (e.g. UMAP or PCA) and the values of the continuous variable.
#' @param dim1 Numeric vector for the first dimension, such as the first dimension of UMAP or PCA.
#' @param dim2 Numeric vector for the second dimension.
#' @param feature Numeric vector of a continuous variable (e.g. mitochondrial read proportion).
#' @param resolution Numeric value of resolution. A higher value will result in smaller squares, thus higher resolution.
#' @param plot Either TRUE or FALSE. If TRUE, the function will directly generate a plot. If FALSE, the function will return the data frame used for plotting. This is useful for more complicated operations.
#' @param palette A vector of color codes for plotting palette.
#' @import ggplot2
#' @export
#' @return If plot=TRUE, a ggplot object. If plot=FALSE, a data.frame with the following columns: xmin,xmax,ymin,ymax of the square, and the averaged expresison value of that square.
#' @author Wenpin Hou<whou10@@jhu.edu>, Zhicheng Ji
#' @examples 
#' scubi_continuous(dim1=rnorm(10000),dim2=rnorm(10000),feature=rnorm(10000))

scubi_continuous <- function(dim1,dim2,feature,resolution=10*20/max(diff(range(dim1)),diff(range(dim2))),plot=TRUE,palette=rainbow(15)[c(11:1,15)]) {
  resolution <- resolution * as.numeric(formatC(20/max(diff(range(dim1)),diff(range(dim2))),digits = 1))
  g <- paste0(round(dim1*resolution),'_',round(dim2*resolution))
  tabd <- table(g)
  tab <- as.vector(tabd)
  names(tab) <- names(tabd)
  me <- rowsum(feature,g)[,1]
  me <- me/tab[names(me)]
  center <- cbind(as.numeric(sub('_.*','',names(me))),as.numeric(sub('.*_','',names(me))))
  pd <- data.frame(xmin=(center[,1]-0.5)/resolution,xmax=(center[,1]+0.5)/resolution,ymin=(center[,2]-0.5)/resolution,ymax=(center[,2]+0.5)/resolution,feature=me)
  if (plot) {
    ggplot() + geom_rect(data=pd,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=feature)) + theme_classic() + scale_fill_gradientn(colors=palette) + theme(legend.position = 'bottom',legend.key.height = unit(0.5, "cm"),legend.margin=margin(t = 0, unit='cm'))
  } else {
    pd
  }
}
