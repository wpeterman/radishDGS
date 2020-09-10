# AIC table ---------------------------------------------------------------
#' @author Bill Peterman
#' @title radish model selection 
#' @description Function to create AIC table of fitted models
#' @param mod_list List containing fitted `radish` objects
#' @param AICc Use second order AIC in ranking models (Default = FALSE). See Details
#' @param BIC Use BIC in ranking models (Default = FALSE). See Details
#' @param mod_names Optional. Specify names for fitted model objects. By default, the right hand side of the fitted `radish` model will be used as the model name.
#' @param verbose (Default = FALSE) Should the table be pronted to the console
#' @param ... Not used
#' 
#' @return Data frame with AIC summary table for provided models
#' @export
#' 
#' @usage XXX
#' 
#' @examples 
#' ## Not Run:
#' ## ** TO BE COMPLETED ** ##  
#' 
#' ## End (Not run)

aic_tab <- function(mod_list,
                    AICc = FALSE,
                    BIC = FALSE,
                    mod_names = NULL,
                    verbose = FALSE,
                    ...){
  ## All models comparable
  mod_dims <- as.vector(lapply(mod_list, function(x) x$dim))
  if(length(unique.default(mod_dims)) != 1L) {
    stop(cat("\n You are attempting to compare models with different number of sample locations or fit on different landscape surfaces. These are not valid comparisons. \n"))
  }
  
  if(is.null(mod_names)){
    mod_names <- as.vector(sapply(mod_list, function(x) as.character(x$formula)[2]))
  }
  mod_loglik <- as.vector(sapply(mod_list, function(x) x$loglik))
  mod_AIC <- as.vector(sapply(mod_list, function(x) x$aic))
  mod_df <- as.vector(sapply(mod_list, function(x) x$df))
  
  
  # AIC ---------------------------------------------------------------------
  if(AICc == FALSE & BIC == FALSE) {
    Delta_AIC <- mod_AIC - min(mod_AIC)
    AIC_wt <- exp(-0.5 * Delta_AIC)
    
    tab <- data.frame(model = mod_names,
                      K = mod_df,
                      AIC = mod_AIC,
                      Delta_AIC = Delta_AIC,
                      AIC_wt = AIC_wt / sum(AIC_wt),
                      row.names = NULL
    )
    
    tab_ <- tab[order(tab$AIC),]
    tab_$Cum.wt <- cumsum(tab_$AIC_wt)
    tab_$loglik <- mod_loglik[order(tab$AIC)]
    tab_[,3:7] <- round(tab_[,3:7], digits = 4)
    
    if(isTRUE(verbose)){
      return(print(tab_, row.names = FALSE))
    } else {
      return(tab_)
    }
    
    # AICc --------------------------------------------------------------------
    
  } else if(isTRUE(AICc)) {
    if(isTRUE(BIC)){
      stop(cat("\n Only a single information criterion can be calculated for a table. Set either `AICc` or `BIC` = FALSE \n"))
    }
    
    mod_n <- as.vector(sapply(mod_list, function(x) x$dim[2]))
    mod_AICc <- -2 * (mod_loglik) + 2 * mod_df * (mod_n / (mod_n - mod_df - 1))
    
    Delta_AICc <- mod_AICc - min(mod_AICc)
    AICc_wt <- exp(-0.5 * Delta_AICc)
    
    tab <- data.frame(model = mod_names,
                      K = mod_df,
                      AIC = mod_AIC,
                      AICc = mod_AICc,
                      Delta_AICc = Delta_AICc,
                      AICc_wt = AICc_wt / sum(AICc_wt),
                      row.names = NULL
    )
    
    tab_ <- tab[order(tab$AICc),]
    tab_$Cum.wt <- cumsum(tab_$AICc_wt)
    tab_$loglik <- mod_loglik
    tab_[,3:8] <- round(tab_[,3:8], digits = 4)
    
    if(isTRUE(verbose)){
      return(print(tab_, row.names = FALSE))
    } else {
      return(tab_)
    }
    
    # BIC ---------------------------------------------------------------------
    
  } else {
    mod_n <- as.vector(sapply(mod_list, function(x) x$dim[2]))
    mod_pairs <- mod_n * (mod_n - 1) / 2
    mod_BIC <- -2 * mod_loglik + mod_df * log(mod_pairs)
    
    Delta_BIC <- mod_BIC - min(mod_BIC)
    BIC_wt <- exp(-0.5 * Delta_BIC)
    
    tab <- data.frame(model = mod_names,
                      K = mod_df,
                      AIC = mod_AIC,
                      BIC = mod_BIC,
                      Delta_BIC = Delta_BIC,
                      BIC_wt = BIC_wt / sum(BIC_wt),
                      row.names = NULL
    )
    
    tab_ <- tab[order(tab$BIC),]
    tab_$Cum.wt <- cumsum(tab_$BIC_wt)
    tab_$loglik <- mod_loglik
    tab_[,3:8] <- round(tab_[,3:8], digits = 4)
    
    if(isTRUE(verbose)){
      return(print(tab_, row.names = FALSE))
    } else {
      return(tab_)
    }
  } ## End BIC
} ## End function