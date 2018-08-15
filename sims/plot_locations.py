import msprime
import pyslim
import numpy as np
import matplotlib
import glob
matplotlib.use('Agg')
from matplotlib import pyplot as plt

outdir = "run_031087/"
# treefiles = [outdir + "pop_{:03d}.trees".format(x) for x in range(1, 21)]
treefiles = glob.glob(outdir + "pop_*.trees")

for treefile in treefiles:
    # treefile = outdir + "pop_100.trees"
    outfile = treefile + ".locs.png"
    ts = pyslim.load(treefile, slim_format=True)
    locations = np.array([i.location for i in ts.individuals()])
    fig = plt.figure(figsize=(max(locations[:,0]),max(locations[:,1])))
    plt.scatter(locations[:,0], locations[:,1], marker='.')
    plt.savefig(outfile, dpi=288)

