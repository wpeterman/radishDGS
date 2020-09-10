# Model Cross Validation --------------------------------------------------
#' @author Bill Peterman
#' @title radish cross validation
#' @description Function to conduct cross validation of `radish` models
#' @param pts SpatialPoints object of sampled demes
#' @param covariates RasterStack of raster layers for `radish` model
#' @param fmla Formula for model to be assessed. LHS should be the name of the distance matrix, RHS contains names of layers in `covariates`
#' @param model c('mlpe', 'wishart')
#' @param prop_train Proportion of sample points used for fitting a training model. The remainder of points will be used in assessment.
#' @param fit_full Should the full `radish` model with all sample demes be fit (Default = TRUE)
#' @param ... Not used
#' 
#' @return Named list containing the `radish` model fit to the training data, the loglikehood of the test data, and (optionally) the model fit to the complete data set.
#' 
#' @usage XXX
#' @export
#' 
#' @examples 
#' ## Not Run:
#' ## ** TO BE COMPLETED ** ##  
#' 
#' ## End (Not run)
radish_cv <- function(pts,
                      covariates,
                      fmla,
                      model = 'mlpe',
                      nu = NULL,
                      prop_train = 0.8,
                      seed = NULL,
                      fit_full = TRUE,
                      ...){
  if(is.null(seed)){
    seed <- floor(runif(1, 1, 1e9))
  }
  
  set.seed(seed)
  
  total_ <- length(pts)
  train_ <- sort(sample(total_, floor(total_ * prop_train)))
  test_ <- which(!(c(1:total_) %in% train_))
  
  fmla <- as.formula(fmla)
  terms <- terms(fmla)
  vars <- as.character(attr(terms, "variables"))[-1]
  resp_name <- vars[1]
  response <- attr(terms, "response")
  gd_mat <- if(response == 1){
    get(vars[attr(terms, "response")], parent.frame())
  } else {
    stop("'formula' must have genetic distance matrix on lhs")
  }
  
  gd_mat_ <- gd_mat
  
  train_pts <- pts[train_]
  test_pts <- pts[test_]
  
  # fmla_ <- reformulate(attr(terms, "term.labels"))
  # fmla_radish <- formula(paste('gd_mat', fmla_, collapse = " + "))
  
  rhs <- attr(terms, "term.labels")
  fmla_radish <- formula(paste0('gd_mat ~ ', paste0(rhs, collapse = " + ")))
  
  if(model == 'mlpe'){
    measurement_model <- radish::mlpe
  } else if(model == 'generalized_wishart' |
            model == 'wishart'){
    measurement_model <- radish::generalized_wishart
  } else {
    measurement_model <- radish::leastsquares
  }
  
  # Fit training data -------------------------------------------------------
  
  cat("\n\n\n ----- Fitting training model ----- \n\n")
  cat(          paste0(rhs, collapse = " + "))
  cat("\n\n    --------------------   \n\n")
  
  gd_mat <- as.matrix(gd_mat_[train_,train_])
  
  training_surface <- conductance_surface(covariates, train_pts, directions = 8)
  fit <- try(radish(fmla_radish, 
                    training_surface,
                    conductance_model = radish::loglinear_conductance,
                    measurement_model = measurement_model,
                    nu = nu,
                    optimizer = 'newton'))
  if(class(fit) == 'try-error'){
    fit <- try(radish(fmla_radish, 
                      training_surface,
                      conductance_model = radish::loglinear_conductance,
                      measurement_model = measurement_model,
                      nu = nu,
                      optimizer = 'bfgs'))
  }
  
  if(class(fit) == 'try-error'){
    stop("Could not successfully optimize training model!")
  }
  
  
  # loglik of test data -----------------------------------------------------
  gd_mat <- as.matrix(gd_mat_[test_,test_])
  test_surface <- conductance_surface(covariates, test_pts, directions = 8)
  
  ll <- radish_grid(theta = matrix(coef(fit), 1, length(coef(fit))), 
                    formula = fmla_radish, 
                    data = test_surface, 
                    conductance_model = radish::loglinear_conductance, 
                    measurement_model = measurement_model,
                    nu = nu)
  
  
  # Full Model --------------------------------------------------------------
  if(isTRUE(fit_full)){
    cat("\n\n\n ----- Fitting full model ----- \n\n")
    cat(          paste0(rhs, collapse = " + "))
    cat("\n\n    --------------------   \n\n")
    
    gd_mat <- as.matrix(gd_mat_)
    full_surface <- conductance_surface(covariates, pts, directions = 8)
    full_fit <- try(radish(fmla_radish, 
                           full_surface,
                           conductance_model = radish::loglinear_conductance,
                           measurement_model = measurement_model,
                           theta = fit$mle$theta,
                           nu = nu,
                           optimizer = 'newton'))
    if(class(full_fit) == 'try-error'){
      full_fit <- try(radish(fmla_radish, 
                             full_surface,
                             conductance_model = radish::loglinear_conductance,
                             measurement_model = measurement_model,
                             nu = nu,
                             optimizer = 'bfgs'))
    }
    
    if(class(full_fit) == 'try-error'){
      stop("Could not successfully optimize full model!")
    }
    out <- list(train_mod = fit,
                cv_loglik = ll$loglik,
                full_mod = full_fit,
                seed = seed)
  } else {
    out <- list(train_mod = fit,
                cv_loglik = ll$loglik,
                seed = seed)
  }
  return(out)
}