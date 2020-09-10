#' @author Bill Peterman
#' @title Specify habitat threshold
#' @description Function to create binary raster from continuous habitat surface
#' 
#' @param hab_rast The simulated habitat surface
#' @param p The proportion of the landscape that should be habitat (0-1)
#' 
#' @export

habitat_threshold <- function(hab_rast,
                              p){
  
  hab_bin <- hab_rast
  hab_bin[] <- as.numeric(hab_rast[] >= quantile(hab_rast[], 1-p))
  hab_ <- hab_rast * hab_bin
  
  hab_[] <- (hab_[] - min(hab_[])) / (max(hab_[]) - min(hab_[]))
  
  hab_bin_rat <- ratify(hab_bin)
  
  hab_stack <- stack(list(hab_cont = hab_,
                          hab_bin = hab_bin_rat))
  names(hab_stack) <- c('hab_cont', 'hab_bin')
  return(hab_stack)
}