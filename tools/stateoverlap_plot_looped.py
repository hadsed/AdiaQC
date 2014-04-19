'''
File: stateoverlap_plot_looped.py
Author: Hadayat Seddiqi
Date: 4.18.14
Description: Plots state overlaps with psi in time, but
             now the plots loop to show all possible
             memory-pair phase plots (and no 3D).
'''

import os, optparse, sys
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as pl
import itertools

# Command line options
if __name__=="__main__":
    parser = optparse.OptionParser("usage: %prog [options] arg1 arg2")
    parser.add_option("-b", "--binary", dest="binary", default=1,
                      type="int", 
                      help="Whether data file is in binary (1) or not (0).")
    (options, args) = parser.parse_args()
    binary = options.binary

# Get state overlap labels
axis = []
with open('state_overlap_labels.txt') as f:
    axis = f.readlines()

# Get probabilities
data = []
ftimes = []
if binary:
    try:
        for fname in os.listdir('.'):
            if fname[:15] == 'state_overlap_T' and fname[-4:] == '.npy':
                prob = np.load(fname)
                t = fname[15:].split('.txt')[0].split('T')[0]
                ftimes.append(float(t))
                data.append(prob)
    except IOError, e:
        print "IOError: No file eigenspectrum.dat.npy. \nYou probably forgot "+\
            "to specify the appropriate command-line argument: -b 0."
        sys.exit(0)
else:
    for fname in os.listdir('.'):
        if fname[:15] == 'state_overlap_T' and fname[-4:] == '.txt':
            prob = np.loadtxt(fname)
            t = fname[15:].split('.txt')[0].split('T')[0]
            ftimes.append(float(t))
            data.append(prob)

# Fix up the data
data = zip(*data)
sortidx = np.array(ftimes).argsort()
for icol, col in enumerate(data):
    data[icol] = np.array(col)[sortidx]
time = np.linspace(0,1,len(data[0]))

# Plot timeseries
pl.title('Six-var QUBO State Probabilities vs. Time')
pl.xlabel('Time')
pl.ylabel('Probability')
pl.ylim([0,1])
pl.xlim([0,1])
for icol, col in enumerate(data):
    pl.plot(time, col, label=axis[icol])
pl.legend(loc=2)
pl.show()

# Plot phase diagram
for pair in itertools.combinations(enumerate(data), 2):
    pl.xlabel(axis[pair[0][0]])
    pl.ylabel(axis[pair[1][0]])
    pl.title('Six-var QUBO Phase Diagram')
    pl.plot(pair[0][1], pair[1][1])
    pl.show()

# # Plot timeseries in 3D
# fig = pl.figure()
# from matplotlib import pyplot
# from mpl_toolkits.mplot3d import Axes3D
# ax = Axes3D(fig)
# ax.set_title('Phase Diagram in Time')
# ax.set_xlabel(axis[0])
# ax.set_ylabel(axis[1])
# ax.set_zlabel('Time')
# ax.scatter(xs=data[0], ys=data[1], zs=time)
# pl.show()
