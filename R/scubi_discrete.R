#' SCUBI for discrete variable
#'
#' SCUBI for plotting a discrete variable
#'
#' This function partitions the whole coordinate space into multiple non-overlapping squares, and calculates and visualizes the proportions of discrete categories for cells falling in each square. It takes as input the coordinates of the cells (e.g. UMAP or PCA) and the values of the discrete variable.
#' @param dim1 Numeric vector for the first dimension, such as the first dimension of UMAP or PCA.
#' @param dim2 Numeric vector for the second dimension.
#' @param feature Character vector of a discrete variable (e.g. gender of the sample the cell was collected).
#' @param resolution Numeric value of resolution. A higher value will result in smaller squares, thus higher resolution.
#' @param plot Either TRUE or FALSE. If TRUE, the function will directly generate a plot. If FALSE, the function will return the data frame used for plotting. This is useful for more complicated operations.
#' @param palette A vector of color codes for plotting palette.
#' @import ggplot2 reshape2 pdist clue
#' @export
#' @return If plot=TRUE, a ggplot object. If plot=FALSE, a data.frame with the following columns: xmin,xmax,ymin,ymax of the square, and the averaged expresison value of that square.
#' @author Wenpin Hou<whou10@@jhu.edu>, Zhicheng Ji
#' @examples 
#' scubi_discrete(dim1=rnorm(10000),dim2=rnorm(10000),feature=sample(c('A','B'),10000,replace=TRUE))

scubi_discrete <- function(dim1,dim2,feature,resolution=10*20/max(diff(range(dim1)),diff(range(dim2))),plot=TRUE,palette=rainbow(length(unique(feature)))) {
  feature <- as.character(feature)
  sumfeature <- table(feature)
  g <- paste0(round(dim1*resolution),'_',round(dim2*resolution))
  tabg <- table(g)
  singleg <- names(tabg)[tabg==1]
  center <- cbind(as.numeric(sub('_.*','',singleg)),as.numeric(sub('.*_','',singleg)))
  singlepd <- data.frame(xmin=(center[,1]-0.5)/resolution,xmax=(center[,1]+0.5)/resolution,ymin=(center[,2]-0.5)/resolution,ymax=(center[,2]+0.5)/resolution,feature=feature[match(singleg,g)],stringsAsFactors = F)
  
  pd <- do.call(rbind,lapply(names(tabg)[tabg > 1],function(i) {
    tmpid <- which(g==i)
    tmpfeat <- feature[tmpid]
    center <- cbind(as.numeric(sub('_.*','',i)),as.numeric(sub('.*_','',i)))
    if (length(unique(tmpfeat)) == 1) {
      data.frame(xmin=(center[,1]-0.5)/resolution,xmax=(center[,1]+0.5)/resolution,ymin=(center[,2]-0.5)/resolution,ymax=(center[,2]+0.5)/resolution,feature=unique(tmpfeat),stringsAsFactors = F)
    } else {
      n <- length(tmpid)
      seq <- n:(n+2)
      l <- sapply(seq,function(i) {
        max(which(sapply(1:(sqrt(i)),function(j) i%%j==0)))
      })
      n2 <- seq[which.max(l)]
      f1 <- max(l)
      f2 <- n2/max(l)
      df <- data.frame(xmin=rep((center[,1]-0.5+seq(0,1,length.out=f1+1)[-f1-1])/resolution,f2),xmax=rep((center[,1]-0.5+seq(0,1,length.out=f1+1)[-1])/resolution,f2),ymin=rep((center[,2]-0.5+seq(0,1,length.out=f2+1)[-f2-1])/resolution,each=f1),ymax=rep((center[,2]-0.5+seq(0,1,length.out=f2+1)[-1])/resolution,each=f1),stringsAsFactors = F)
      dfc <- cbind(rowMeans(df[,c('xmin','xmax')]),rowMeans(df[,c('ymin','ymax')]))
      pointpos <- cbind(dim1[tmpid],dim2[tmpid])
      dm <- as.matrix(pdist(pointpos,dfc))
      if (nrow(dm) > 500) {
        # use suboptimal but faster routine
        dmm <- melt(dm)
        dmm <- dmm[order(dmm[,3]),]
        assign <- rep(NA,nrow(dm))
        while(sum(is.na(assign)) > 0) {
          ndid <- which(!duplicated(dmm[,1]) & !duplicated(dmm[,2]))
          assign[dmm[ndid,1]] <- dmm[ndid,2]
          dmm <- dmm[!dmm[,1]%in%dmm[ndid,1] & !dmm[,2]%in%dmm[ndid,2],]
        }
        #sum(sapply(1:length(assign),function(i) dm[i,assign[i]]))
      } else {
        assign <- as.vector(solve_LSAP(dm))
      }
      df <- df[assign,]
      df$feature <- tmpfeat
      df
    }
  }))
  pd <- rbind(pd,singlepd)
  if (plot) {
    ggplot() + geom_rect(data=pd,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=feature)) + theme_classic() + scale_fill_manual(values=palette) + theme(legend.position = 'bottom',legend.key.height = unit(0.5, "cm"),legend.margin=margin(t = 0, unit='cm'))
  } else {
    pd
  }
}
