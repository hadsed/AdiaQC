'''
File: eigspec_plot.py
Author: Hadayat Seddiqi
Date: 4.16.14
Description: Plots eigenspectrum from outputted data file.
'''

import os, optparse, sys
import itertools
import numpy as np
import pylab as pl

# Command line options
if __name__=="__main__":
    parser = optparse.OptionParser("usage: %prog [options] arg1 arg2")
    parser.add_option("-b", "--binary", dest="binary", default=1,
                      type="int", 
                      help="Whether data file is in binary (1) or not (0).")
    parser.add_option("-f", "--fname", dest="filename", default=None,
                      type="string", 
                      help="Output filename (without .png suffix).")
    parser.add_option("-e", "--eignum", dest="eignum", default=None,
                      type="int", 
                      help="Number of eigenvalues to plot (from smallest).")
    parser.add_option("-t", "--title", dest="title", default="Eigenspectrum",
                      type="string", 
                      help="Title for the graph.")
    (options, args) = parser.parse_args()
    binary = options.binary
    filename = options.filename
    eignum = options.eignum
    title = options.title

# Get data
if binary:
    try:
        eigs = np.load('eigenspectrum.dat.npy')
    except IOError, e:
        print "IOError: No file eigenspectrum.dat.npy. \nYou probably forgot "+\
            "to specify the appropriate command-line argument: -b 0."
        sys.exit(0)
else:
    eigs = np.loadtxt('eigenspectrum.dat')

time = eigs[:,0]/eigs[-1,0]
eigs = eigs[:,1:]
if not eignum:
    eignum = eigs.shape[1]

for col in xrange(eignum):
    pl.plot(time, eigs[:,col])

pl.title(title)
pl.xlabel('t/T')
pl.ylabel('Energy')
if filename is None:
    pl.savefig('eigenspectrum.png')
else:
    pl.savefig(filename+'.png')
