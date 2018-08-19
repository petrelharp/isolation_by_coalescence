#!/usr/bin/env python3
description = '''Compute mean pairwise divergences across a discritization of the landscape.'''

import msprime
import pyslim
import numpy as np
import os, glob

### command-line parsing: execute this with '-h' for help
import argparse

parser = argparse.ArgumentParser(description=description)
parser.add_argument("--basedir", "-o", type=str, dest="basedir", required=True,
                    help="name of directory to save output files to.")
parser.add_argument("--logfile", "-g", type=str, dest="logfile", 
                    help="Name of log file [default: outdir/divergences.log]")
parser.add_argument("--samples_per_cell", "-s", type=int, dest="samples_per_cell", 
                    default=40, help="Number of samples per cell.")
parser.add_argument("--num_subsamples", "-k", type=int, dest="num_subsamples", 
                    default=2, help="Number of indpendent subsamples.")
parser.add_argument("--mutation_rate", "-u", type=float, dest="mutation_rate", 
                    help="Mutation rate.")
parser.add_argument("--grid_width", "-n", type=int, dest="grid_width", 
                    help="Number of cells across the landscape (x direction) in the discretization.")
parser.add_argument("--grid_height", "-m", type=int, dest="grid_height", 
                    help="Number of cells up the landscape (y direction) in the discretization.")

args = parser.parse_args()

if not os.path.isdir(args.basedir):
    os.mkdir(args.basedir)

if args.logfile is None:
    args.logfile = os.path.join(args.basedir, "divergences.log")

logfile = open(args.logfile, "w")

print(args, file=logfile)

outdir = args.basedir
samples_per_cell = args.samples_per_cell
num_subsamples = args.num_subsamples
mutation_rate = args.mutation_rate
num_chroms = 1
grid_width = args.grid_width
grid_height = args.grid_height

def grid_samples(ts, n, m, prob=1.0, shrink=0.75):
    '''
    Returns a list of lists of node IDs, the ones that fall in the middle
    shrink of each rectangle given by discretizing into an (nxm) grid.
    '''
    samples = [[[] for _ in range(m)] for _ in range(n)]
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


for treefile in glob.glob(os.path.join(outdir, "*.trees")):
    logfile.write("Reading {}".format(treefile))
    logfile.flush()
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

    sample_grid = grid_samples(recap, n=grid_width, m=grid_height, prob=1.0, shrink=0.75)
    samples = [a for b in sample_grid for a in b]
    # extract nonoverlapping subsamples
    num_subs_per_cell = [min(samples_per_cell, int(len(a)/num_subsamples)) for a in samples]
    for a in num_subs_per_cell:
        assert(a > 0)

    for subs in range(num_subsamples):
        subsamples = [[a[j] for j in range(subs * k, (subs + 1) * k)] 
                      for a, k in zip(samples, num_subs_per_cell)]
        samplefile = open(treefile + "." + str(subs) + ".samples.tsv", 'w')
        samplefile.write("\t".join(['id', 'population', 'x', 'y', 'individual']) + "\n")
        for k in range(len(subsamples)):
            for a in subsamples[k]:
                ind = recap.individual(recap.node(a).individual)
                samplefile.write("{}\t{}\t{}\t{}\t{}\n".format(a, k, ind.location[0], ind.location[1], recap.node(a).individual))
        samplefile.close()
        # we want individual divergences
        the_subsamples = [x for y in subsamples for x in y]
        # simplify first for speed
        sub_ts = recap.simplify(the_subsamples)
        mut_ts = msprime.mutate(sub_ts, rate=mutation_rate)
        # breaks at the chromosomes
        windows = np.linspace(0.0, mut_ts.sequence_length, num_chroms + 1)
        # bs = msprime.BranchLengthStatCalculator(sub_ts)
        bs = msprime.SiteStatCalculator(mut_ts)
        divs = np.array(bs.divergence([[x] for x in range(len(the_subsamples))], 
                                      windows=windows))
        outfile = treefile + "." + str(subs) + ".divs.tsv"
        header = '\t'.join(['chrom' + str(j) for j in range(1, num_chroms + 1)])
        np.savetxt(outfile, divs, delimiter='\t', header=header)

logfile.flush()
