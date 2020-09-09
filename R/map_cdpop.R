#' @title Map simulations
#' @author Bill Peterman
#' @description Function to map simulated demes onto conductance surface
#' @param params Object generated from `get_params` OR results from `cdpop_sim`
#' @param ... Not used
#' 
#' @return Plot of conductance surface with sampled demes
#' @export
#' 
#' @usage XXX
#' 
#' @examples 
#' ## Not Run:
#' ## map_cdpop(params) ##  
#' 
#' ## End (Not run)
#' 
map_cdpop <- function(params){

  conduct <- params$conductance_surface
  pts <- params$pts

  raster::plot(conduct, main = 'Sampled demes on conductance')
  raster::plot(pts, add = TRUE, pch = 19)
  
  
}