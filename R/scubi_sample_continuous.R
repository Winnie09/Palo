#' SCUBI for multi-sample continuous variable comparison
#'
#' SCUBI for comparing a continuous variable across multiple samples
#'
#' This function partitions the whole coordinate space into multiple non-overlapping squares, calculates the proportions of discrete categories for cells falling in each square, and calculates and visualizes the contrast between sample categories. It takes as input the coordinates of the cells (e.g. UMAP or PCA), the cell-sample assignment (each cell belong to which sample), and a numeric design matrix for the samples. An intercept will be automatically added to the design matrix so do not include any intercept term. The first column of the design matrix will always be the contrast of interest and will be visualized.
#' @param dim1 Numeric vector for the first dimension, such as the first dimension of UMAP or PCA.
#' @param dim2 Numeric vector for the second dimension.
#' @param feature Numeric vector of a continuous variable (e.g. mitochondrial read proportion).
#' @param sample Character vector of cell-sample assignment (each cell belong to which sample).
#' @param design A numeric data.frame or matrix indicating the covariates of the samples. It should not contain the intercept term. The row names of the matrix should contain the unique values in 'sample' variable.
#' @param smooth Either TRUE or FALSE. If TRUE, the coefficients will be smoothed across the coordinate space with loess.
#' @param clip Either TRUE or FALSE. If TRUE, will enforce values larger than 0.975 quantile to equal 0.975 quantile, and enforce values smaller than 0.025 quantile to equal 0.025 quantile.
#' @param scale Either TRUE or FALSE. If TRUE, the values will be scaled to have mean 0 and standard deviation of 1.
#' @param resolution Numeric value of resolution. A higher value will result in smaller squares, thus higher resolution.
#' @param plot Either TRUE or FALSE. If TRUE, the function will directly generate a plot. If FALSE, the function will return the data frame used for plotting. This is useful for more complicated operations.
#' @param palette A vector of color codes for plotting palette.
#' @import ggplot2
#' @export
#' @return If plot=TRUE, a ggplot object. If plot=FALSE, a data.frame with the following columns: xmin,xmax,ymin,ymax of the square, and the averaged expresison value of that square.
#' @author Wenpin Hou<whou10@@jhu.edu>, Zhicheng Ji 
#' @examples 
#' scubi_sample_continuous(dim1=rnorm(10000),dim2=rnorm(10000),feature=rnorm(10000),sample=sample(c('S1','S2'),10000,replace=TRUE),design=data.frame(group=c(1,0),row.names=c('S1','S2')))

scubi_sample_continuous <- function(dim1,dim2,feature,sample,design,smooth=TRUE,clip=TRUE,scale=TRUE,resolution=10*20/max(diff(range(dim1)),diff(range(dim2))),plot=TRUE,palette=rainbow(15)[c(11:1,15)]) {
  resolution <- resolution * as.numeric(formatC(20/max(diff(range(dim1)),diff(range(dim2))),digits = 1))
  design <- as.matrix(design)
  g <- paste0(round(dim1*resolution),'_',round(dim2*resolution))
  comb <- paste0(g,' ',sample)
  expr <- tapply(feature,comb,mean)
  
  square <- sub(' .*','',names(expr))
  x <- cbind(1,design[sub('.* ','',names(expr)),,drop=F])
  me <- sapply(unique(square),function(s) {
    id <- which(square==s)
    if (length(unique(x[id,2]))==1) {
      NA
    } else {
      tx <- x[id,]
      ty <- expr[id]
      (chol2inv(chol(crossprod(tx))) %*% t(tx) %*% ty)[2,]
    }
  })
  center <- cbind(as.numeric(sub('_.*','',names(me))),as.numeric(sub('.*_','',names(me))))
  if (clip) {
    q1 <- quantile(me,0.025,na.rm=T) 
    q2 <- quantile(me,0.975,na.rm=T) 
    me[me <= q1] <- q1
    me[me >= q2] <- q2
  }
  if (smooth) {
    aid <- which(!is.na(me))
    me[aid] <- fitted(loess(me[aid]~center[aid,1]+center[aid,2],span=0.1,degree = 1))
  }
  if (scale) me <- scale(me)
  pd <- data.frame(xmin=(center[,1]-0.5)/resolution,xmax=(center[,1]+0.5)/resolution,ymin=(center[,2]-0.5)/resolution,ymax=(center[,2]+0.5)/resolution,feature=me)
  if (plot) {
    ggplot() + geom_rect(data=pd,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=feature)) + theme_classic() + scale_fill_gradientn(colors=palette) + theme(legend.position = 'bottom',legend.key.height = unit(0.5, "cm"),legend.margin=margin(t = 0, unit='cm'))
  } else {
    pd
  }
}
