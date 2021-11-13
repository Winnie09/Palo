#' SCUBI for gene expression
#'
#' SCUBI for plotting the expression of one gene
#'
#' This function partitions the whole coordinate space into multiple non-overlapping squares, and calculates and visualizes the averaged expression of cells falling in each square. It takes as input the coordinates of the cells (e.g. UMAP or PCA), the gene expression count matrix, and the name of the gene to be plotted.
#' @param dim1 Numeric vector for the first dimension, such as the first dimension of UMAP or PCA.
#' @param dim2 Numeric vector for the second dimension.
#' @param count Count matrix for the gene expression data. Each row is a gene and each column is a cell. A sparse matrix is accepted.
#' @param gene A character value, a character vector, or a list of gene names. The gene names must exist in the count matrix. For a character value, the expression of that gene will be plotted. For a character vector, the read counts will be summed across these genes. For a list of gene names, each element is a gene or a vector of genes. A subplot will be generated for each element in the list.
#' @param resolution Numeric value of resolution. A higher value will result in smaller squares, thus higher resolution.
#' @param plot Either TRUE or FALSE. If TRUE, the function will directly generate a plot. If FALSE, the function will return the data frame used for plotting. This is useful for more complicated operations.
#' @param palette A vector of color codes for plotting palette.
#' @import ggplot2 Matrix
#' @export
#' @return If plot=TRUE, a ggplot object. If plot=FALSE, a data.frame with the following columns: xmin,xmax,ymin,ymax of the square, and the averaged expresison value of that square.
#' @author Wenpin Hou<whou10@@jhu.edu>, Zhicheng Ji
#' @examples 
#' scubi_expression(dim1=rnorm(10000),dim2=rnorm(10000),count=matrix(rbinom(100000,10,0.1),nrow=10,dimnames=list(LETTERS[1:10],NULL)),gene='A')
#' scubi_expression(dim1=rnorm(10000),dim2=rnorm(10000),count=matrix(rbinom(100000,10,0.1),nrow=10,dimnames=list(LETTERS[1:10],NULL)),gene=list('A','B'))

scubi_expression <- function(dim1,dim2,count,gene,resolution=10*20/max(diff(range(dim1)),diff(range(dim2))),plot=TRUE,palette=rainbow(15)[c(11:1,15)]) {
  g <- paste0(round(dim1*resolution),'_',round(dim2*resolution))
  rc <- rowsum(colSums(count),g)[,1]
  if (!is.list(gene)) {
    if (length(gene)==1) {
      me <- rowsum(count[gene,],g)[,1]  
    } else {
      me <- rowsum(colSums(count[gene,]),g)[,1]
    }
    me <- log2(me/rc*1e5 + 1)
    center <- cbind(as.numeric(sub('_.*','',names(me))),as.numeric(sub('.*_','',names(me))))
    pd <- data.frame(xmin=(center[,1]-0.5)/resolution,xmax=(center[,1]+0.5)/resolution,ymin=(center[,2]-0.5)/resolution,ymax=(center[,2]+0.5)/resolution,expression=me)
    if (plot) {
      ggplot() + geom_rect(data=pd,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=expression)) + theme_classic() + scale_fill_gradientn(colors=palette) + theme(legend.position = 'bottom',legend.key.height = unit(0.5, "cm"),legend.margin=margin(t = 0, unit='cm'))
    } else {
      pd
    }  
  } else {
    gl <- gene
    pd <- do.call(rbind,sapply(gl,function(gene) {
      if (length(gene)==1) {
        me <- rowsum(count[gene,],g)[,1]  
      } else {
        me <- rowsum(colSums(count[gene,]),g)[,1]
      }
      me <- log2(me/rc*1e5 + 1)
      center <- cbind(as.numeric(sub('_.*','',names(me))),as.numeric(sub('.*_','',names(me))))
      data.frame(Gene=paste0(gene,collapse = '_'),xmin=(center[,1]-0.5)/resolution,xmax=(center[,1]+0.5)/resolution,ymin=(center[,2]-0.5)/resolution,ymax=(center[,2]+0.5)/resolution,expression=me)
    },simplify = F))
    if (plot) {
      ggplot() + geom_rect(data=pd,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=expression)) + theme_classic() + scale_fill_gradientn(colors=palette) + theme(legend.position = 'bottom',legend.key.height = unit(0.5, "cm"),legend.margin=margin(t = 0, unit='cm')) + facet_wrap(~Gene)
    } else {
      pd
    }
  }
}
