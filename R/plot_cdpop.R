#' @title Plot simulations
#' @author Bill Peterman
#' @description Function to plot genetic distances generated from CDPOP simulations
#' @param params Object generated from `get_params` OR results from `cdpop_sim`
#' @param ... Not used
#' 
#' @return 4-panel plot of two genetic distances X resistance distance and Euclidean distance
#' @export
#' 
#' @usage XXX
#' 
#' @examples 
#' ## Not Run:
#' ## plot_cdpop(params) ##  
#' 
#' ## End (Not run)
#' 
plot_cdpop <- function(params){
  par(mfrow = c(2,2))
  
  resist <- lower(params$trueRes)
  geoD <- lower(params$geoD)
  
  if(is.null(params$dc)){
    pca <- lower(params$pca)
    dps <- lower(params$dps)
    
    plot(dps ~ resist,
         xlab = "Resistance distance",
         ylab = "Proportion shared alleles")
    plot(pca ~ resist,
         xlab = "Resistance distance",
         ylab = "PCA distance")
    
    plot(dps ~ geoD,
         xlab = "Geographic distance",
         ylab = "Proportion shared alleles")
    plot(pca ~ geoD,
         xlab = "Geographic distance",
         ylab = "PCA distance")

    
  } else {
    dc <- lower(params$dc)
    fst <- lower(params$fst)
    
    plot(dc ~ resist,
         xlab = "Resistance distance",
         ylab = "Chord distance")
    plot(fst ~ resist,
         xlab = "Resistance distance",
         ylab = "Fst")

    plot(dc ~ geoD,
         xlab = "Geographic distance",
         ylab = "Chord distance")
    plot(fst ~ geoD,
         xlab = "Geographic distance",
         ylab = "Fst")

  }
  par(mfrow = c(1,1))
}