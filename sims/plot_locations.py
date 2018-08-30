import msprime
import pyslim
import numpy as np
import matplotlib
import glob
matplotlib.use('Agg')
from matplotlib import pyplot as plt

treefiles = glob.glob("run_0*/*.trees")

for treefile in treefiles:
    print(treefile)
    outfile = treefile + ".locs.png"
    locfile = treefile + ".locs.tsv"
    ts = msprime.load(treefile)
    locations = np.array([i.location for i in ts.individuals()])
    np.savetxt(locfile, locations, header="x\ty\tz")
    fig = plt.figure(figsize=(max(locations[:,0]),max(locations[:,1])))
    plt.scatter(locations[:,0], locations[:,1], marker='.')
    plt.savefig(outfile, dpi=288)
    print("... done!")

