msprime_wrapper <- 
  function(master_seed, 
           reps, 
           covariates, 
           effect_size,
           buffer_size = 0.2,
           number_of_demes = 300,
           sampled_proportion_of_demes = 0.1,
           conductance_quantile_for_demes = 0.5,
           dispersal_kernel = function(dist) exp(-1./0.3 * dist^2),
           number_of_loci = 500,
           maf = 0.05)
{
  library(radish)

  reticulate::source_python(system.file("py/island_model.py",package="radishDGS"))

  stopifnot(reps > 0)
  stopifnot(length(effect_size) == dim(covariates)[3])
  stopifnot(buffer_size >= 0. & buffer_size < 0.5)
  stopifnot(number_of_demes > 0)
  stopifnot(sampled_proportion_of_demes > 0. & sampled_proportion_of_demes <= 1.)
  stopifnot(conductance_quantile_for_demes > 0.)
  stopifnot(maf >= 0. & maf < 0.5)

  set.seed(master_seed)
  random_seeds <- sample.int(1e8, reps)

  # create buffer surrounding study area
  buffer           <- covariates[[1]]
  buffer[]         <- 0
  buffer_size_row  <- floor(nrow(buffer) * buffer_size)
  buffer_size_col  <- floor(ncol(buffer) * buffer_size)
  if(buffer_size_row > 0.)
  {
    buffer[1:buffer_size_row,] <- 1
    buffer[nrow(buffer):(nrow(buffer)-buffer_size_row+1),] <- 1
  }
  if(buffer_size_col > 0.)
  {
    buffer[,1:buffer_size_col] <- 1
    buffer[,ncol(buffer):(ncol(buffer)-buffer_size_col+1)] <- 1
  }

  # log-linear conductance surface
  conductance <- covariates[[1]]
  if (length(covariates) > 1)
    conductance[] <- exp(effect_size * covariates[[1]][])
  else
    conductance[] <- exp(raster::getValues(covariates) %*% effect_size)

  output <- list()
  for (seed in random_seeds)
  {
    timer <- Sys.time()
    cat("\n---\nRunning simulation with random seed", seed, "on", as.character(timer), "\n")

    try({ #this ensures that failure of a given simulation will not stop the loop

      set.seed(seed)

      # simulate deme locations by sampling grid cells uniformly without replacement, where some condition is met
      threshold               <- quantile(conductance, conductance_quantile_for_demes, na.rm = TRUE)
      possible_cells          <- which(!is.na(conductance[]) & conductance[] >= threshold)
      all_deme_cells          <- sample(possible_cells, number_of_demes)
      all_deme_coords         <- xyFromCell(conductance, all_deme_cells, spatial=TRUE)

      # choose some number of demes inside study area (e.g. not in buffer) to sample
      number_of_sampled_demes <- floor(sampled_proportion_of_demes * number_of_demes)
      demes_not_in_buffer     <- which(!extract(buffer, all_deme_coords))
      sampled_demes           <- sort(sample(demes_not_in_buffer, number_of_sampled_demes)) #will fail if insufficient number

      # get "true" resistance distances among demes
      # NB: the effect sizes here are for CONDUCTANCE (the inverse of resistance)
      form <- as.formula(paste0("~", paste(names(covariates), sep="+")))
      resistance_model        <- radish::conductance_surface(covariates, all_deme_coords, directions=8)
      resistance_distance     <- radish::radish_distance(matrix(effect_size, 1, length(effect_size)), form, resistance_model, radish::loglinear_conductance)$distance[,,1]
      geographic_distance     <- radish::radish_distance(matrix(0, 1, length(effect_size)), form, resistance_model, radish::loglinear_conductance)$distance[,,1]

      # standardize resistance distance to [0,1] and calculate migration matrix
      resistance_distance    <- scale_to_0_1(resistance_distance)
      migration_matrix       <- dispersal_kernel(resistance_distance)
      diag(migration_matrix) <- 0 #diagonal must be 0 for msprime

      # diploid population sizes per deme: sample single individual per deme
      pop_size <- rep(1, number_of_demes)

      # track a diploid from each "sampled" deme
      # NB: should be less than or equal to 2*pop_size from previous step
      smp_size <- ifelse(1:number_of_demes %in% sampled_demes, 2, 0)

      # use msprime to simulate haploid genotypes
      # NB: a key point is that we only simulate haploids (e.g gametes)
      #     so diploid genotypes have to be formed from haplotypes. We also generally want to filter
      #     by minor allele frequency, as low frequency alleles are often less informative about population
      #     structure (and may violate model assumptions)
      SNPsim <- island_model(num_loci   = number_of_loci, #number of loci (e.g. SNPs)
                             split_time = 200, #generations in past where panmixia ended
                             anc_size   = 500, #ancestral population size, as number of DIPLOIDS
                             migr_mat   = migration_matrix, #migration matrix. element m_ij is proportion of population i that is replaced by individuals from population j, PER GENERATION. These should be small.
                             pop_size   = pop_size, #contemporary (DIPLOID) population sizes, as number of diploids
                             smp_size   = smp_size, #number of sampled HAPLOIDS per population (can be 0, in which case populations are modelled but not sampled)
                             ran_seed   = seed,
                             maf_filter = maf,
                             keep_trees = TRUE, #if true, keep coalescent tree for each locus as newick (could be used to simulate msats)
                             use_dtwf   = TRUE) #if true, use Wright-Fisher backward time simulations (discrete generations, better when population sizes are small/migr rates high)

      # create diploid genotypes by pooling two haploids from same spatial location
      genotypes          <- SNPsim$genotypes[,seq(1,ncol(SNPsim$genotypes),2)] + SNPsim$genotypes[,seq(2,ncol(SNPsim$genotypes),2)]
      demes_of_genotypes <- rep(1:number_of_demes, smp_size/2) #deme id for each sampled diploid
      was_sampled        <- demes_of_genotypes %in% sampled_demes
      demes_of_samples   <- demes_of_genotypes[was_sampled]

      # filter by minor allele frequency (only applies to SNPs)
      frequency  <- rowMeans(genotypes) / 2
      maf_filter <- frequency >= maf & frequency <= 1.-maf
      if (any(!maf_filter))
        warning("MAF rejection sampling failed, check python source")
      genotypes  <- genotypes[maf_filter,]
      frequency  <- frequency[maf_filter]
      cat("Simulation produced", nrow(genotypes), "SNPs passing MAF filter\n")

      # genetic covariance (for SNPs)
      normalized_genotypes <- (genotypes - frequency)/sqrt(2 * frequency * (1-frequency)) #"normalized" genotypes 
      genetic_covariance   <- t(normalized_genotypes) %*% normalized_genotypes / (nrow(normalized_genotypes)-1)

      # store results needed to fit models/reproduce simulation (won't store anything if simulation fails)
      timestamp     <- Sys.time()
      output[[paste0("seed", seed)]] <- 
                       list(covariates    = covariates, 
                            genotypes     = genotypes[,was_sampled], 
                            coords        = all_deme_coords[demes_of_samples,], 
                            rdistance     = resistance_distance[demes_of_samples,demes_of_samples],
                            gdistance     = geographic_distance[demes_of_samples,demes_of_samples],
                            migration     = migration_matrix[demes_of_samples,demes_of_samples],
                            covariance    = genetic_covariance[was_sampled,was_sampled],
                            frequency     = frequency,
                            random_seed   = seed,
                            timestamp     = timestamp)
    })
    timer <- Sys.time() - timer
    cat("Finished in", timer, attr(timer, "units"), "\n")
  }

  # track which ones failed (will help us troubleshoot)
  cat("\nFailed for", length(random_seeds)-length(output), "of", length(random_seeds), "simulations\n")
  failures <- random_seeds[!(random_seeds %in% names(output))]
  attr(output, "failed") <- failures

  output
}
