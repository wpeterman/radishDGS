# Get parameters from fitted radish models ------------------------------------------------------
#' @author Bill Peterman
#' @title Get parameter estimates from fitted `radish` models
#' @description Function to fit `radish` models to CDPOP simulation data
#' @param Results_dir Full path to the top 'Results' directory where simulation results are saved
#' @param radish_model (Default = 'wishart) Specify which `radish` model you want to get parameter estimates from. Should be one of c('wishart', 'mlpe', 'ls'). Assumes there is only a single fitted model object for each measurement model in the directory.
#' @param save_table (Default = TRUE) The parameter estimate table will be saved as a CSV file
#' @param conv Extra parameter to be used with other functions when checking radish convergence
#' @param ... Not used
#' 
#' @return List of fitted `radish` models
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
radish_parameters <- function(Results_dir,
                              radish_model = 'wishart',
                              save_table = TRUE,
                              conv = NULL,
                              ...){
  stopifnot(dir.exists(Results_dir))
  
  params <- get_params(Results_dir)
  
  if(length(params$all_dirs[grep(paste0('--',radish_model), 
                                 params$all_dirs)]) != 1){
    stop("First fit a radish model OR Specify correct radish model to evaluate!")
  }
  
  mod <- readRDS(params$all_dirs[grep(paste0('--',radish_model), 
                                      params$all_dirs)])
  
  est_eff <- coef(mod)
  
  all_res <- readRDS(params$all_dirs[grep('AllResults', 
                                          params$all_dirs)])
  eff_size <- all_res$effect_size
  SE <- sqrt(diag(summary(mod)$vcov))
  
  if(!is.null(conv)){
    df <- data.frame(parameter = names(est_eff),
                     true_eff = eff_size,
                     est_eff = est_eff,
                     est_SE = SE,
                     converged = conv
    )
  } else {
    df <- data.frame(parameter = names(est_eff),
                     true_eff = eff_size,
                     est_eff = est_eff,
                     est_SE = SE)
  }
  
  if(isTRUE(save_table)){
    write.table(df,
                paste0(Results_dir, '--parameter_est.csv'),
                sep = ",",
                row.names = F,
                col.names = T) 
  }
  
  
  return(df)
}