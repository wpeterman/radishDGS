#' @author Bill Peterman
#' @description Function to create correlated Fractal Brownian Motion surfaces
#' 
#' @param corr Desired correlation level between fbm surfaces
#' @param dim Dimension of simulated raster (Default = 400)
#' @param H Fractal dimension (Hurst parameter, H) controlling clustering in Fractal Brownian Motion landscape simulation model (Default = 0.5)
#' @param user_seed Seed can be set for reproducibility (Default = NA)
#' 
#' @export

fbm_corr <- function(corr,
                     dim,
                     H = 0.5,
                     user_seed = NA){
  library(NLMR)
  library(raster)
  
  sim_rast <- nlm_fbm(ncol = dim,
                      nrow = dim, 
                      fract_dim = H, 
                      user_seed = user_seed)
  
  rast_corr <- 0
  
  while(rast_corr > corr + 0.015 | rast_corr < corr - 0.015){
    rep_sim <- nlm_fbm(ncol = dim,
                       nrow = dim, 
                       fract_dim = H, 
                       user_seed = NA)
    
    corr_ <-  sqrt((1 / (corr^2)) - 1) 
    mat <- as.matrix(sim_rast)
    
    rep_sim_ <- rep_sim * corr_
    corr_rast <- sim_rast + rep_sim_
    
    rast_corr <- layerStats(stack(sim_rast, corr_rast), 'pearson')$`pearson correlation coefficient`[1,2]
    cat(paste0(round(rast_corr, digits = 4), '\n \n'))
  }
  
  cat(paste0("Surface correlation is ", round(rast_corr, 4), '\n \n'))
  
  
  corr_rast[] <- (corr_rast[] - min(corr_rast[])) / (max(corr_rast[]) - min(corr_rast[]))
  
  
  return(stack(list(sim_true = sim_rast,
                    corr_sim = corr_rast)))
}