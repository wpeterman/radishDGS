#' @author Bill Peterman
#' @description Function to create binary raster from Fractal Brownian Motion surface
#' 
#' @param fbm The simulated fbm surface
#' @param p The proportion of the landscape that should be habitat (0-1)
#' 
#' @export

fbm_habitat <- function(fbm,
                        p){
  library(raster)
  
  fbm_bin <- fbm
  fbm_bin[] <- as.numeric(fbm[] >= quantile(fbm[], 1-p))
  fbm_hab <- fbm * fbm_bin
  
  fbm_hab[] <- (fbm_hab[] - min(fbm_hab[])) / (max(fbm_hab[]) - min(fbm_hab[]))
  
  return(stack(list(fbm_cont = fbm_hab,
                    fbm_bin = fbm_bin)))
}