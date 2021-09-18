#' SCUBI for discrete variable
#'
#' SCUBI for plotting a discrete variable
#'
#' This function partitions the whole coordinate space into multiple non-overlapping squares, and calculates and visualizes the proportions of discrete categories for cells falling in each square. It takes as input the coordinates of the cells (e.g. UMAP or PCA) and the values of the discrete variable.
#' @param dim1 Numeric vector for the first dimension, such as the first dimension of UMAP or PCA.
#' @param dim2 Numeric vector for the second dimension.
#' @param feature Character vector of a discrete variable (e.g. gender of the sample the cell was collected).
#' @param normalize Either TRUE or FALSE. If TRUE, for each discrete category, number of cells in each region will be normalized by the total number of cells. This is to address the bias induced by the unbalanced total number of cells for discrete categories.
#' @param resolution Numeric value of resolution. A higher value will result in smaller squares, thus higher resolution.
#' @param plot Either TRUE or FALSE. If TRUE, the function will directly generate a plot. If FALSE, the function will return the data frame used for plotting. This is useful for more complicated operations.
#' @param palette A vector of color codes for plotting palette.
#' @import ggplot2
#' @export
#' @return If plot=TRUE, a ggplot object. If plot=FALSE, a data.frame with the following columns: xmin,xmax,ymin,ymax of the square, and the averaged expresison value of that square.
#' @author Wenpin Hou<whou10@@jhu.edu>, Zhicheng Ji
#' @examples 
#' scubi_discrete(dim1=rnorm(10000),dim2=rnorm(10000),feature=sample(c('A','B'),10000,replace=TRUE))

scubi_discrete <- function(dim1,dim2,feature,normalize=FALSE,resolution=10*20/max(diff(range(dim1)),diff(range(dim2))),plot=TRUE,palette=rainbow(length(unique(feature)))) {
  resolution <- resolution * as.numeric(formatC(20/max(diff(range(dim1)),diff(range(dim2))),digits = 1))
  feature <- as.character(feature)
  sumfeature <- table(feature)
  g <- paste0(round(dim1*resolution),'_',round(dim2*resolution))
  tabd <- table(g)
  tab <- as.vector(tabd)
  names(tab) <- names(tabd)
  me <- table(paste0(g,':-:',feature))
  meg <- sub(':-:.*','',names(me))
  pd <- do.call(rbind,lapply(unique(meg),function(i) {
    tmp <- me[which(meg==i)]
    if (normalize) {
      tmp <- round(tmp/sumfeature[sub('.*:-:','',names(tmp))] * max(sumfeature))
    }
    center <- cbind(as.numeric(sub('_.*','',i)),as.numeric(sub('.*_','',i)))
    if (length(tmp) == 1) {
      mn <- names(tmp)
      data.frame(xmin=(center[,1]-0.5)/resolution,xmax=(center[,1]+0.5)/resolution,ymin=(center[,2]-0.5)/resolution,ymax=(center[,2]+0.5)/resolution,feature=sub('.*:-:','',mn),stringsAsresolutions = F)
    } else {
      mn <- names(which.max(tmp))
      totaln <- sum(tmp)
      facid <- 2:(totaln-1)
      cand <- facid[sapply(facid,function(i) totaln %% i)==0]
      if (length(cand)==0) {
        tmp <- tmp * 3
        totaln <- sum(tmp)
        facid <- 2:(totaln-1)
        cand <- facid[sapply(facid,function(i) totaln %% i)==0]
      }
      f1 <- cand[which.min(abs(cand-sqrt(totaln)))]
      f2 <- totaln/f1
      
      tmp <- tmp[names(tmp)!=mn]
      df1 <- data.frame(xmin=(center[,1]-0.5)/resolution,xmax=(center[,1]+0.5)/resolution,ymin=(center[,2]-0.5)/resolution,ymax=(center[,2]+0.5)/resolution,feature=sub('.*:-:','',mn),stringsAsresolutions = F)
      df2 <- data.frame(xmin=rep((center[,1]-0.5+seq(0,1,length.out=f1+1)[-f1-1])/resolution,f2),xmax=rep((center[,1]-0.5+seq(0,1,length.out=f1+1)[-1])/resolution,f2),ymin=rep((center[,2]-0.5+seq(0,1,length.out=f2+1)[-f2-1])/resolution,each=f1),ymax=rep((center[,2]-0.5+seq(0,1,length.out=f2+1)[-1])/resolution,each=f1),stringsAsresolutions = F)
      df2 <- df2[order(colSums((t(df2)-unlist(df2[1,]))^2)),]
      df2 <- df2[1:sum(tmp),]
      df2$feature <- sub('.*:-:','',rep(names(tmp),tmp))
      rbind(df1,df2)
    }
  }))
  if (plot) {
    ggplot() + geom_rect(data=pd,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=feature)) + theme_classic() + scale_fill_manual(values=palette) + theme(legend.position = 'bottom',legend.key.height = unit(0.5, "cm"),legend.margin=margin(t = 0, unit='cm'))
  } else {
    pd
  }
}
