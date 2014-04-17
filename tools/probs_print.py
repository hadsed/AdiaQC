'''
File: probs_print.py
Author: Hadayat Seddiqi
Date: 4.16.14
Description: Prints probabilities to screen.
'''

import os, optparse, sys
import numpy as np

# Command line options
if __name__=="__main__":
    parser = optparse.OptionParser("usage: %prog [options] arg1 arg2")
    parser.add_option("-f", "--filename", dest="filename", default=None,
                      type="string", 
                      help="Filename for output probabilities.")
    parser.add_option("-s", "--sort", dest="sort", default=0,
                      type="int", 
                      help="Sort by probabilities.")
    (options, args) = parser.parse_args()
    filename = options.filename
    sort = options.sort

# Get data
probs = np.load(filename)
bstrs = [ line.rstrip('\n') for line in open('statelabels.txt') ]
zipped = zip(bstrs, probs)

# Sort
if sort:
    zipped = sorted(zipped, key=lambda x: x[1])[::-1]

# Print
for pair in zipped:
    print pair
