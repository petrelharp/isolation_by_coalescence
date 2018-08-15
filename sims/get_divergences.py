import msprime
import pyslim
import numpy as np

samples_per_cell = 40
num_subsamples = 2
num_chroms = 1

for outdir in ["run_015334", "run_003033", "run_010760"]:

    treefile = outdir + "/pop_100000.trees"
    decap = pyslim.load(treefile, slim_format=True)

    # recapitate
    recomb = msprime.RecombinationMap(positions = [0.0, decap.sequence_length], 
                                      rates = [1e-8, 0.0],
                                      num_loci = int(decap.sequence_length))

    pop_configs = [msprime.PopulationConfiguration() for _ in range(decap.num_populations)]

    recap = msprime.simulate(
                from_ts = decap, 
                population_configurations = pop_configs,
                recombination_map = recomb,
                start_time = decap.slim_generation,
                Ne=2e4)

    ts = msprime.mutate(recap, rate=1e-8)

    def grid_samples(ts, n, m, prob=1.0, shrink=0.75):
        '''
        Returns a list of lists of node IDs, the ones that fall in the middle
        shrink of each rectangle given by discretizing into an (nxm) grid.
        '''
        samples = [[[] for _ in range(n)] for _ in range(m)]
        locs = np.array(ts.tables.individuals.location)
        locs.resize((int(len(ts.tables.individuals.location,)/3), 3))
        max_x = np.ceil(100*max(locs[:,0]))/100
        max_y = np.ceil(100*max(locs[:,1]))/100
        dx = max_x / n
        dy = max_y / m
        for ind in ts.individuals():
            if np.random.uniform() < prob:
                i = int(np.floor(ind.location[0] / dx))
                j = int(np.floor(ind.location[1] / dy))
                # distance to center of cell
                x_diff = ind.location[0]/dx - (i + 0.5)
                y_diff = ind.location[1]/dy - (j + 0.5)
                if (abs(x_diff) < shrink / 2) and (abs(y_diff) < shrink / 2):
                    samples[i][j].extend(ind.nodes)
        # put these in random order
        for a in samples:
            for b in a:
                np.random.shuffle(b)
        return samples

    sample_grid = grid_samples(ts, 4, 4, prob=1.0, shrink=0.75)
    samples = [a for b in sample_grid for a in b]
    # extract nonoverlapping subsamples
    num_subs_per_cell = [min(samples_per_cell, int(len(a)/num_subsamples)) for a in samples]
    for a in num_subs_per_cell:
        assert(a > 0)

    for subs in range(num_subsamples):
        subsamples = [[a[j] for j in range(subs * k, (subs + 1) * k)] 
                      for a, k in zip(samples, num_subs_per_cell)]
        # breaks at the chromosomes
        windows = np.linspace(0.0, ts.sequence_length, num_chroms + 1)
        # bs = msprime.BranchLengthStatCalculator(ts)
        bs = msprime.SiteStatCalculator(ts)
        divs = np.array(bs.divergence(samples, windows=windows))
        samplefile = open(treefile + "." + str(subs) + ".samples.tsv", 'w')
        for a in samples:
            samplefile.write("\t".join(map(str,a)) + "\n")
        samplefile.close()
        outfile = treefile + "." + str(subs) + ".divs.tsv"
        header = '\t'.join(['chrom' + str(j) for j in range(1, num_chroms + 1)])
        np.savetxt(outfile, divs, delimiter='\t', header=header)
