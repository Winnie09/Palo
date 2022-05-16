#' Palo: Color palette optimization
#'
#' Palo: Color palette optimization for single-cell and spatial data
#'
#' @param position A matrix of low-dimensional coordinates for single-cell data or 2-D spatial positions for spatial transcriptomics data. The number of rows of the matrix is the number of cells and the matrix has two columns.
#' @param cluster A vector of cell or spots clusters. Should have the same length as the number of cells or spots. Should have the same order as the columns of \code{position}.
#' @param palette A user-defined character vector of palette.
#' @param rgb_weight A numeric vector of weights (lengths of three) for RGB values when calculating color differences. Common choices of weights include (1,1,1), (3,4,2), and (2,4,3).
#' @param color_blind_fun A character value indicating if the color palette will be transformed before calculating the color distances for colorblind-friendly visualizations. The value has to be one of 'deutan','protan','tritan','desaturate', or NULL. If NULL, no transformation will be done.
#' @param init_iter A numeric value of number of optimization iterations for  initial optimization.
#' @param refine_iter A numeric value of number of optimization iterations for  refined optimization.
#' @param early_stop A numeric value to early stop the optimization process if the color score remains unchanged for early_stop consecutive exchanges.
#' @import MASS colorspace
#' @export
#' @return A named vector of colors.
#' @author Wenpin Hou<whou10@@jhu.edu>, Zhicheng Ji
#' @examples 
#' x <- matrix(rnorm(2000),ncol=2)
#' x[1:200,1] <- x[1:200,1] + 5
#' x[200 + 1:200,1] <- x[200 + 1:200,1] + 10
#' x[400 + 1:200,2] <- x[400 + 1:200,2] + 5
#' x[600 + 1:200,2] <- x[600 + 1:200,2] + 10
#' cluster <- paste0('cluster',rep(1:5,each=200))
#' palette <- c("orange","red","green","blue","yellow")
#' pal <- Palo(x,cluster,palette)
#' qplot(x[,1],x[,2],col=cluster) + scale_color_manual(values=pal)

Palo <- function(position,cluster,palette,rgb_weight=c(1,1,1),color_blind_fun=NULL,init_iter=1000,refine_iter=1000,early_stop=100) {
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
    if (!is.null(color_blind_fun)) {
      eval(parse(text=paste0('cv <- ',color_blind_fun,'(cv)')))
    }
    cdist <- as.matrix(dist(t(sqrt(rgb_weight)*col2rgb(cv))))
    cscore <- apply(ss,1,function(i) {
      cdist[which(colnames(k)==i[1]),which(colnames(k)==i[2])]
    })
    sum(ss[,3] * cscore)  
  }
  sampscore <- sapply(1:init_iter,function(seed) {
    set.seed(seed)
    scorefunc(sample(palette))  
  })
  set.seed(which.max(sampscore))
  palette <- sample(palette)
  score <- scorefunc(palette)
  
  idlecount <- 0
  for (siter in 1:refine_iter) {
    idlecount <- idlecount + 1
    tc <- palette
    id <- sample(1:length(palette),2)
    tmp <- tc[id[2]]
    tc[id[2]] <- tc[id[1]]
    tc[id[1]] <- tmp
    tscore <- scorefunc(tc)
    if (tscore > score) {
      score <- tscore
      palette <- tc
      idlecount <- 0
    }
    if (idlecount > early_stop) {
      break
    }
  }
  names(palette) <- uc
  palette  
}
