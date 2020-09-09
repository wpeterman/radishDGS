#' @author Bill Peterman
#' @title  Make a binary raster continuous
#' @description Function to create a continuous surface from a binary surface
#' 
#' @param bin_rast Binary (1/0) raster
#' @param window The size of the moving window to calculate mean habitat within (Default = 3)
#' 
#' @export

bin2cont <- function(bin_rast,
                     window = 3){

  cont_rast <- raster::focal(bin_rast,
                     w = matrix(1 / window^2, nc = window, nr = window))
  return(cont_rast)
}