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
    (options, args) = parser.parse_args()
    binary = options.binary

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

# probs = np.load(fpath+str(num)+'/probsT'+str(t)+'.dat.npy')
# bstrs = [ line.rstrip('\n') for line in open(fpath+str(num)+
#                                              '/statelabels.txt') ]

time = eigs[:,0]/eigs[-1,0]
eigs = eigs[:,1:]

for col in xrange(eigs.shape[1]):
    pl.plot(time, eigs[:,col])

pl.title('Eigenspectrum')
pl.xlabel('Time')
pl.ylabel('Energy')
pl.savefig('eigenspectrum.png')
pl.show()
