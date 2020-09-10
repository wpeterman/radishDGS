#' @author Bill Peterman
#' @title Internal PopGenReport function
#' @description Function to run PopGenReport from R
#' 
#' @param sim_name Name for simulation results. Defaults to 'output'
#' @param pts Spatial points object
#' @param sim_dir Directory where simulation results will be written
#' @param resist_mat Resistance distance matrix
#' @param n_ind Default = 25; Number of individuals in each deme
#' @param n.offspring Default = 2; Poisson draw around ‘mean fecundity’
#' @param mig.rate Default = 0.2; Proportion of dispersing individuals that achieve the maximum dispersal distance
#' @param disp.rate Default = 0.25; Proportion of offspring that disperse
#' @param disp.max Default = 0.125; What percentage of the maximum pairwise distance can a dispersing individual traverse?
#' @param sex.ratio Default = 0.5
#' @param n.loci Default = 20
#' @param n.allels Default = 20
#' @param steps Number of iterations to run simulation
#' @param mut.rate Default = 0.0005

#' 
#' @keywords internal
#' 
PopGen_sim <- function(sim_name = 'output_',
                       pts,
                       sim_dir,
                       resist_mat,
                       n_ind = 25,
                       n.offspring = 2,
                       mig.rate = 0.2,
                       disp.rate = 0.25,
                       disp.max = 0.125,
                       sex.ratio = 0.5,
                       n.loci = 20,
                       n.allels = 20,
                       steps = 200,
                       mut.rate = 0.0005){
  
  
  
  
  
  # Create directories ------------------------------------------------------
  
  if(!dir.exists(sim_dir)) dir.create(sim_dir, recursive = TRUE)
  suppressWarnings(
    dir.create(paste0(sim_dir,"/data/"), recursive = TRUE)
  )
  data_dir <- paste0(sim_dir,"/data/")
  
  
  # Fill NULL ---------------------------------------------------------------
  # Run PopGenReport  ------------------------------------------------------
  
  pop.sim_init <- init.popgensim(n.pops = length(pts),
                                 n.ind = n_ind, 
                                 sex.ratio = sex.ratio,
                                 n.loci = n.loci,
                                 n.allels = n.allels, 
                                 n.cov = 3)
  
  sim.out <-  run.popgensim(simpops = pop.sim_init,
                            steps = steps,
                            cost.mat = resist_mat,
                            n.offspring = n.offspring,
                            n.ind = n_ind,
                            mig.rate = mig.rate,
                            disp.max = disp.max * max(resist_mat),
                            disp.rate = disp.rate,
                            n.allels = n.allels,
                            mut.rate = mut.rate)
  
  pops.gi <- pops2genind(sim.out)
  

  
  
  # Wrap-up -----------------------------------------------------------------
  
  saveRDS(pops.gi,
          paste0(data_dir, "sim_genind.rds"))
  
  return(pops.gi)
  
}

# # PCA dist -------------------------------------------------------
# 
# pca_dist <- function(gi, 
#                      n_axes = 64){
#   a_tab <- adegenet::tab(gi)
#   pc <- prcomp(a_tab)
#   pc_dist <- as.matrix(dist(pc$x[,1:n_axes]))
#   return(pc_dist)
# }

# Random Samples ----------------------------------------------------------

## Randomly select populations and individuals from within populations
# 
# gi_samp <- function(gi,
#                     n_ind = 100) {
#   ind_samp <- sort(sample(1:nInd(gi), n_ind))
#   gi_s <- gi[ind_samp]
#   
#   out <- list(genind = gi_s,
#               pop_samp = ind_samp)
# }
