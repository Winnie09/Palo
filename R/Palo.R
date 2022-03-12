#' Palo: Color palette optimization
#'
#' Palo: Color palette optimization for single-cell and spatial data
#'
#' @param position A matrix of low-dimensional coordinates for single-cell data or 2-D spatial positions for spatial transcriptomics data. The number of rows of the matrix is the number of cells and the matrix has two columns.
#' @param cluster A vector of cell or spots clusters. Should have the same length as the number of cells or spots. Should have the same order as the columns of \code{position}.
#' @param palette A user-defined character vector of palette.
#' @param iter A numeric value of optimization iterations.
#' @import MASS
#' @export
#' @return A named vector of colors.
#' @author Wenpin Hou<whou10@@jhu.edu>, Zhicheng Ji
#' @examples 
#' Palo(matrix(rnorm(200),ncol=2),rep(c('c1','c2'),50),c('red','blue'))

Palo <- function(position,cluster,palette,iter=1000) {
  cluster <- as.character(cluster)
  uc <- unique(cluster)
  sampcluster <- sample(cluster)
  k <- sapply(uc,function(i) {
    kde2d(position[cluster==i,1],position[cluster==i,2],n=100,lims=c(range(position[,1]),range(position[,2])))[[3]]
  })
  
  sampk <- sapply(uc,function(i) {
    kde2d(position[sampcluster==i,1],position[sampcluster==i,2],n=100,lims=c(range(position[,1]),range(position[,2])))[[3]]
  })
  
  uf <- expand.grid(1:ncol(k),1:ncol(k))
  uf <- uf[uf[,1]<uf[,2],]
  
  cut <- apply(sampk,2,function(i) quantile(i,0.95))

  ss <- apply(uf,1,function(i) {
    sum(k[,i[1]] > cut[i[1]] & k[,i[2]] > cut[i[2]])/sum(k[,i[1]] > cut[i[1]] | k[,i[2]] > cut[i[2]])
  })
  
  ss <- data.frame(t1=uc[uf[,1]],t2=uc[uf[,2]],score=ss,stringsAsFactors = F)
  
  scorefunc <- function(cv) {
    cdist <- as.matrix(dist(t(col2rgb(cv))))
    cscore <- apply(ss,1,function(i) {
      cdist[which(colnames(k)==i[1]),which(colnames(k)==i[2])]
    })
    sum(ss[,3] * cscore)  
  }
  sampscore <- sapply(1:iter,function(seed) {
    set.seed(seed)
    scorefunc(sample(palette))  
  })
  set.seed(which.max(sampscore))
  palette <- sample(palette)
  score <- scorefunc(palette)
  
  for (siter in 1:iter) {
    tc <- palette
    id <- sample(1:length(palette),2)
    tmp <- tc[id[2]]
    tc[id[2]] <- tc[id[1]]
    tc[id[1]] <- tmp
    tscore <- scorefunc(tc)
    if (tscore > score) {
      score <- tscore
      palette <- tc
    }
  }
  names(palette) <- uc
  palette  
}
