import msprime
import pyslim
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

outdir = "test06/"
treefiles = [outdir + "pop_{:03d}.trees".format(x) for x in range(1, 21)]

for treefile in treefiles:
    # treefile = outdir + "pop_100.trees"
    outfile = treefile + ".locs.png"
    ts = pyslim.load(treefile, slim_format=True)
    locations = np.array([i.location for i in ts.individuals()])
    fig = plt.figure(figsize=(5,5))
    plt.scatter(locations[:,0], locations[:,1], marker='.')
    plt.savefig(outfile, dpi=288)

