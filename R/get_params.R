# Get parameters ------------------------------------------------------
#' @author Bill Peterman
#' @title Retrieve data
#' @description Function to retrieve CDPOP simulation files for `radish` analysis
#' @param dir Full path to directory where simulation results saved
#' @param ... Not used
#' 
#' @return Data frame or named list with loglikelihood table and AIC summary table for provided models
#' @export
#' 
#' @usage XXX
#' 
#' @examples 
#' ## Not Run:
#' ## ** TO BE COMPLETED ** ##  
#' 
#' ## End (Not run)
#' 
get_params <- function(dir,
                       ...){
  stopifnot(dir.exists(dir))
  dir_files <- list.files(dir, 
                          full.names = T)
  
  AllRes_list <- readRDS(dir_files[grep('AllResults_list.rds', 
                                        dir_files)])
  
  ## Import genind
  gi <- AllRes_list$sim_genind
  
  ## Conductance surface
  conductance_surface <- AllRes_list$conductance_surface
  names(conductance_surface) <- 'Conductance'
 
  
  ## Import locations
  pts_ <- AllRes_list$pts

  ## Import genetic distances
  if(any(grep("Dps_dist",dir_files))){
    Dps <- as.matrix(read.csv(dir_files[grep('Dps_dist', 
                                             dir_files)])[,-1])
    pca <- as.matrix(read.csv(dir_files[grep('pca_dist', 
                                             dir_files)])[,-1])
    geo <- as.matrix(read.csv(dir_files[grep('geographic', 
                                             dir_files)])[,-1])
    true_res <- as.matrix(read.csv(dir_files[grep('true_Resist', 
                                                  dir_files)])[,-1])
    
    out_list <- list(pts = pts_,
                     covariates = AllRes_list$covariates,
                     dps = Dps,
                     pca = pca,
                     geoD = geo,
                     trueRes = true_res,
                     conductance_surface = conductance_surface,
                     sim_genind = gi,
                     effect_size = AllRes_list$effect_size,
                     out_dir = dir,
                     all_dirs = dir_files)
  }
  
  if(any(grep("Fst_dist",dir_files))){
    Fst <- as.matrix(read.csv(dir_files[grep('Fst_dist', 
                                             dir_files)])[,-1])
    Dc <- as.matrix(read.csv(dir_files[grep('Dc_dist', 
                                            dir_files)])[,-1])
    geo <- as.matrix(read.csv(dir_files[grep('geographic', 
                                             dir_files)])[,-1])
    true_res <- as.matrix(read.csv(dir_files[grep('true_Resist', 
                                                  dir_files)])[,-1])
    
    out_list <- list(pts = pts_,
                     covariates = AllRes_list$covariates,
                     dc = Dc,
                     fst = Fst,
                     geoD = geo,
                     trueRes = true_res,
                     conductance_surface = conductance_surface,
                     sim_genind = gi,
                     effect_size = AllRes_list$effect_size,
                     out_dir = dir,
                     all_dirs = dir_files)
  }
  return(out_list)
}