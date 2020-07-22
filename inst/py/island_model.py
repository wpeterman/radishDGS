import msprime
import numpy

def island_model (num_loci, split_time, anc_size, migr_mat, pop_size, smp_size, ran_seed, maf_filter = 0., keep_trees = True, use_dtwf = False, keep_variants = True, verbose = True):
    # time is in generations, sizes are absolute number DIPLOIDS
    # migration matrix entry m_ij is the proportion of i that comes FROM j
    num_loci = int(num_loci)
    split_time = float(split_time)
    anc_size = float(anc_size)
    migr_mat = numpy.array(migr_mat, dtype=float)
    pop_size = numpy.array(pop_size, dtype=float)
    smp_size = numpy.array(smp_size, dtype=int)
    ran_seed = int(ran_seed)
    assert(migr_mat.shape[0] == migr_mat.shape[1])
    assert(migr_mat.shape[0] == pop_size.shape[0])
    assert(migr_mat.shape[0] == smp_size.shape[0])
    assert(split_time > 0.)
    assert(anc_size > 0.)
    assert(num_loci > 0)
    assert(maf_filter > 0. and maf_filter < 0.5)
    num_popul = smp_size.shape[0]
    # define populations
    pop_confg = []
    dem_event = []
    for i in range(num_popul):
      pop_confg += [msprime.PopulationConfiguration(
        sample_size=smp_size[i], initial_size=pop_size[i])]
    # demographic events: split from ancestral lineage at same time
    for i in range(num_popul-1):
      dem_event += [msprime.MassMigration(
        time = split_time, source = i, 
        destination = num_popul-1, proportion = 1.)]
    dem_event += [msprime.PopulationParametersChange(
        time = split_time, initial_size = anc_size)]
    # simulate trees
    if use_dtwf:
        model = "dtwf"
    else:
        model = "hudson"
    mu = float(num_popul)/(numpy.mean(pop_size)*split_time) #heuristic; should generate mutations for most loci without wasting too much time
    if verbose:
        dd = msprime.DemographyDebugger(
                population_configurations=pop_confg,
                demographic_events=dem_event,
                migration_matrix=migr_mat)
        dd.print_history()
    sim = msprime.simulate(
            population_configurations = pop_confg,
            migration_matrix = migr_mat,
            demographic_events = dem_event,
            Ne = anc_size,
            mutation_rate = mu,
            num_replicates = num_loci,
            random_seed = ran_seed,
            model = model)
    # iterate over trees, extract genotypes and newick strings
    if verbose:
        print("[msprime_wrapper] mutating trees")
    trees = []
    genotypes = []
    for locus in sim:
        if keep_trees:
            trees += [locus.at(0.).newick()]
        if keep_variants:
            # if locus has no variants, or if it doesn't pass MAF filter, 
            # mutate it until it does
            if verbose:
                print("[msprime_wrapper] simulating locus ...")
            if locus.get_num_mutations() == 0:
                fail = True
            else:
                for site in locus.variants():
                    freq = float(numpy.sum(site.genotypes))/float(numpy.sum(smp_size)) 
                    if freq < maf_filter or freq > 1.-maf_filter:
                        fail = True
                        break
                    else:
                        fail = False
                        break
            while fail:
                ran_seed += 100
                locus = msprime.mutate(locus, rate = mu, keep = False, random_seed = ran_seed)
                # reject if first site doesn't pass MAF filter
                if locus.get_num_mutations() == 0:
                    fail = True
                else:
                    for site in locus.variants():
                        freq = float(numpy.sum(site.genotypes))/float(numpy.sum(smp_size)) 
                        if freq < maf_filter or freq > 1.-maf_filter:
                            fail = True
                            break
                        else:
                            fail = False
                            break
            # get first variant within locus
            for site in locus.variants():
                genotypes += [site.genotypes]
                break
    genotypes = numpy.array(genotypes)
    return {"genotypes" : genotypes, "trees" : trees}
