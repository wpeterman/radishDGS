# CV Model Selection ------------------------------------------------------
#' @author Bill Peterman
#' @title Conduct model selection on cross validation models
#' @description Function to create model selection tables from fitted cross validation models
#' @param cv_list List containing fitted cross validation objects
#' @param cv_names Optional. Specify names for fitted cv objects. By default, the right hand side of the fitted `radish` model will be used as the name.
#' @param aic If the full radish model was fit during cross validation, an AIC table can also be created (Default = FALSE)
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

cv_model_selection <- function(cv_list,
                               cv_names = NULL,
                               aic = FALSE,
                               ...){
  
  if(is.null(cv_names)){
    # cv_names <- names(cv_list)
    cv_names <- as.vector(sapply(cv_list, 
                                 function(x) as.character(x$train_mod$formula)[2]))
  }
  
  mod_df <- as.vector(sapply(cv_list, function(x) x$train_mod$df))
  mod_loglik <- as.vector(sapply(cv_list, function(x) x$cv_loglik))
  Delta_LL <- mod_loglik - max(mod_loglik)
  
  
  cv_tab <- data.frame(model = cv_names,
                       K = mod_df,
                       loglik = mod_loglik,
                       Delta_LL = Delta_LL,
                       row.names = NULL)
  
  cv_tab_ <- cv_tab[order(-cv_tab$loglik),]
  
  
  # AIC table ---------------------------------------------------------------
  
  if(isTRUE(aic)){
    mod_list <- lapply(cv_list, function(x) x$full_mod)
    mod_names <- as.vector(sapply(mod_list, function(x) as.character(x$formula)[2]))
    
    tab <- aic_tab(mod_list,
                   mod_names = mod_names)
    
    out <- list(loglik_tab = cv_tab_,
                AIC_tab = tab)
  } else {
    
    # No AIC table ------------------------------------------------------------
    
    out <- cv_tab_   
  }
  return(out)
}