# Fit radish models ------------------------------------------------------
#' @author Bill Peterman
#' @title Fit radish models to simulated data
#' @description Function to fit `radish` models to CDPOP simulation data
#' @param Results_dir Full path to the top 'Results' directory where simulation results are saved. If only a single iteration of a simulation is to be analyzed, point to that directory within the 'Results' directory.
#' @param conductance_model Default = `loglinear_conductance`. c('loglinear_conductance', 'linear_conductance')
#' @param measurement_model Default = `generalized_wishart`. c('generalized_wishart', 'mlpe', 'leastsquares')
#' @param optimizer (Default = 'newton')
#' @param gd Genetic distance to use. Specify 'dps' or 'pca' if using individual-based simulations; 'dc' or 'fst' for population-based simulations.
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
fit_radish2sims <- function(Results_dir,
                            conductance_model = 'loglinear_conductance',
                            measurement_model = 'generalized_wishart',
                            optimizer = 'newton',
                            gd,
                            ...){
  stopifnot(dir.exists(Results_dir))
  stopifnot(gd %in% c('dps', 'pca', 'fst', 'dc'))
  
  sim_dirs <- list.dirs(Results_dir, recursive = F)
  
  corr_list <- 
    param_list <- 
    mod_list <- vector('list', length(sim_dirs))
  
  if(length(sim_dirs) == 1){ ## Single run
    if(basename(sim_dirs) != 'data'){
      stop("Specify a correct path to simulation results!")
    }
    sim_dirs <- Results_dir
    
    cat(paste0("\n \n        ****  Processing ", basename(sim_dirs), "   ****", '\n \n'))
    
    params <- get_params(sim_dirs)
    gen_dist <- as.matrix(params[[gd]])
    
    ## Measurement model
    if(measurement_model == 'mlpe'){
      measurement_model <- radish::mlpe
      mm <- 'mlpe'
    } else if(measurement_model == 'generalized_wishart' |
              measurement_model == 'wishart'){
      measurement_model <- radish::generalized_wishart
      mm <- 'wishart'
    } else {
      measurement_model <- radish::leastsquares
      mm <- 'ls'
    }
    
    ## Conductance model
    cm <- conductance_model
    if(conductance_model == 'loglinear_conductance'){
      conductance_model <- radish::loglinear_conductance
      
    } else {
      conductance_model <- radish::linear_conductance
    }
    
    ## radish data
    radish_cond <- conductance_surface(params$covariates, 
                                       params$pts,
                                       directions = 8)
    
    rhs <- names(params$covariates)
    fmla_radish <- formula(paste0('gen_dist ~ ', paste0(rhs, collapse = " + ")))
    
    nloci <- adegenet::nLoc(params$sim_genind)
    
    ## Use testthat::capture_messages to capture and report when convergence fails????
    mod <- radish(formula = fmla_radish, 
                  data = radish_cond, 
                  conductance_model = conductance_model,
                  measurement_model = measurement_model,
                  nu = nloci, # Number of loci
                  optimizer = 'newton')
    
    mod_list[[1]] <- mod
    
    
    saveRDS(mod, paste0(sim_dirs,'/radish--',cm,'--', mm, '.rds'))
    
    radish_params <- radish_parameters(Results_dir = sim_dirs,
                                       radish_model = mm)
    param_list[[1]] <- radish_params
    
    # Resistance Corr ---------------------------------------------------------
    
    resistance_distance <- radish::radish_distance(theta = matrix(radish_params[[3]], 
                                                                  1, 
                                                                  length(radish_params[[3]])),
                                                   formula = fmla_radish, 
                                                   data = radish_cond, 
                                                   conductance_model = radish::loglinear_conductance)$distance[,,1] 
    
    resCorr <- cor(lower(params$trueRes),
                   lower(resistance_distance))
    
    
    # Raster Corr -------------------------------------------------------------
    
    trueCond <- params$conductance_surface 
    
    surf_ind <- rep(1, length(params$effect_size))
    optCond <-   raster::stackApply(params$covariates,
                                    indices = surf_ind, function(x, ...) exp(x %*% radish_params[[3]]))
    
    rastCorr <- layerStats(stack(trueCond, optCond),
                           stat = 'pearson')$`pearson correlation coefficient`[1,2]
    
    corr_df <- data.frame(resist_corr = resCorr,
                          conductRast_corr = rastCorr)
    corr_list[[1]] <- corr_df
    
    
  } ## End sim dir loop
  else { ## Multiple sims
    for(i in 1:length(sim_dirs)){
      cat(paste0("\n \n        ****  Processing ", basename(sim_dirs[i]), "   ****", '\n \n'))
      
      params <- get_params(sim_dirs[i])
      gen_dist <- as.matrix(params[[gd]])
      
      ## Measurement model
      mm_ <- measurement_model
      if(is.function(mm_)){
        measurement_model <- mm
      }
      if(measurement_model == 'mlpe'){
        measurement_model <- radish::mlpe
        mm <- 'mlpe'
      } else if(measurement_model == 'generalized_wishart' |
                measurement_model == 'wishart'){
        measurement_model <- radish::generalized_wishart
        mm <- 'wishart'
      } else {
        measurement_model <- radish::leastsquares
        mm <- 'ls'
      }
      
      ## Conductance model
      cm_ <- conductance_model
      if(is.function(cm_)){
        conductance_model <- cm
      }
      if(conductance_model == 'loglinear_conductance'){
        conductance_model <- radish::loglinear_conductance
        cm <- 'loglinear_conductance'
        
      } else {
        conductance_model <- radish::linear_conductance
        cm <- 'linear_conductance'
      }
      
      ## radish data
      radish_cond <- conductance_surface(params$covariates, 
                                         params$pts,
                                         directions = 8)
      
      rhs <- names(params$covariates)
      fmla_radish <- formula(paste0('gen_dist ~ ', paste0(rhs, collapse = " + ")))
      
      nloci <- adegenet::nLoc(params$sim_genind)
      
      ## Use testthat::capture_messages to capture and report when convergence fails????
      mod <- radish(formula = fmla_radish, 
                    data = radish_cond, 
                    conductance_model = conductance_model,
                    measurement_model = measurement_model,
                    nu = nloci, # Number of loci
                    optimizer = 'newton')
      
      mod_list[[i]] <- mod
      
      saveRDS(mod, paste0(sim_dirs[i],'/radish--',cm,'--', mm, '.rds'))
      
      sim_names <- basename(sim_dirs)
      
      radish_params <- radish_parameters(Results_dir = sim_dirs[i],
                                         radish_model = mm,
                                         save_table = FALSE)
      param_list[[i]] <- radish_params
      
      
      # Resistance Corr ---------------------------------------------------------
      
      resistance_distance <- radish::radish_distance(theta = matrix(radish_params[[3]], 
                                                                    1, 
                                                                    length(radish_params[[3]])),
                                                     formula = fmla_radish, 
                                                     data = radish_cond, 
                                                     conductance_model = radish::loglinear_conductance)$distance[,,1] 
      
      resCorr <- cor(lower(params$trueRes),
                     lower(resistance_distance))
      
      
      # Raster Corr -------------------------------------------------------------
      
      trueCond <- params$conductance_surface 
      
      surf_ind <- rep(1, length(params$effect_size))
      optCond <-   raster::stackApply(params$covariates,
                                      indices = surf_ind, function(x, ...) exp(x %*% radish_params[[3]]))
      
      rastCorr <- layerStats(stack(trueCond, optCond),
                             stat = 'pearson')$`pearson correlation coefficient`[1,2]
      
      corr_df <- data.frame(resist_corr = resCorr,
                            conductRast_corr = rastCorr)
      corr_list[[i]] <- corr_df
      
      names(corr_list) <- names(mod_list) <- names(param_list) <- sim_names
      
    } ## End sim dir loop
  }
  
  effects <- do.call(rbind, param_list)
  effects <- data.frame(iter = row.names(effects),
                        effects,
                        row.names = NULL)
  
  corrTab <- do.call(rbind, corr_list)
  corrTab <- data.frame(iter = row.names(corrTab),
                        corrTab,
                        row.names = NULL)
  
  write.table(effects, paste0(Results_dir,'/all_sims_parameters_est.csv'),
              sep = ",", row.names = F, col.names = T)
  write.table(corrTab, paste0(Results_dir,'/correlationTable.csv'),
              sep = ",", row.names = F, col.names = T)
  out_list <- list(fitted_radish = mod_list,
                   effects = effects,
                   corrTab = corrTab)  
  return(out_list)
}
