#' @author Nate Pope
#' @title Scale rasters 
#' @description Scale raster to range from 0-1
#' 
#' @param x Raster to be scaled
#' @export
scale_to_0_1 <- function(x) { #scales a given vector so that min(x)==0 and max(x)==1
  (x - min(x,na.rm=TRUE))/(max(x,na.rm=TRUE) - min(x,na.rm=TRUE))
}

#' @author Nate Pope
#' @title Lower triangle 
#' @description Get lower triangle of square matrix, return vector
#' 
#' @param x Square matrix
#' @export
lower <- function(x){ #returns lower triangular part of a matrix
  x[lower.tri(x)]
}

# PCA dist -------------------------------------------------------

#' @author Bill Peterman
#' @title PCA distance
#' @description Calculate Euclidean distance between principle component axes
#' @param gi `genind` object from adegenet
#' @param n_axes Number of principle component axes to retain and calculate distance between
#' @param scale (Default = TRUE) Calculated distance will be rescaled 0-1
pca_dist <- function(gi, 
                     n_axes = 64,
                     scale = T){
  a_tab <- adegenet::tab(gi)
  pc <- prcomp(a_tab)
  
  if(isTRUE(scale)){
    pc_dist <- as.matrix(scale_to_0_1(dist(pc$x[,1:n_axes])))
    
  } else {
    pc_dist <- as.matrix(dist(pc$x[,1:n_axes]))
    
  }
  return(pc_dist)
}