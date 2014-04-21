'''
File: exp1_analysis.py
Author: Hadayat Seddiqi
Date: 4.20.14
Description: Plots and facts on hopfield experiment #1.
'''

import os, optparse, sys
import numpy as np
import json
import pylab as pl

# Command line options
if __name__=="__main__":
    parser = optparse.OptionParser("usage: %prog [options] arg1 arg2")
    parser.add_option("-p", "--probsfname", dest="probsfname", default=None,
                      type="string", 
                      help="Filename for output probabilities.")
    parser.add_option("-o", "--output", dest="output", default=1,
                      type="int", 
                      help="Output this report to file.")
    parser.add_option("-b", "--binary", dest="binary", default=1,
                      type="int", 
                      help="Binary input or not (default = 1).")
    parser.add_option("-q", "--qubits", dest="qubits", default=None,
                      type="int", 
                      help="Number of qubits/neurons.")
    (options, args) = parser.parse_args()
    probsfname = options.probsfname
    output = options.output
    binary = options.binary
    qubits = options.qubits

def spins2bitstr(vec):
    """ Return a converted bitstring from @vec, a list of spins. """
    return ''.join([ '0' if k == 1 else '1' for k in vec ])

def hamdist(a,b):
    """ Calculate Hamming distance. """
    return np.sum(abs(np.array(a)-np.array(b))/2.0)

def analysis(filename, binary):
    """
    """
    # Get probability data
    if binary:
        probs = np.load(filename)
    else:
        probs = np.loadtxt(filename)
    # probstr = [ '%1.16E' % e for e in probs.ravel() ]
    bstrs = [ line.rstrip('\n') for line in open('statelabels.txt') ]

    # Get problem outputs
    props = json.load(open('problem_outputs.dat'))
    neurons = props['nQubits']
    lrule = props['learningRule']
    mems = props['memories']
    instate = props['inputState']
    annealtime = props['annealTime']
    answer = props['answer']

    # Calculate average Hamming distance
    avgdist = np.average([ hamdist(m, instate) for m in mems ])
    
    # Calculate success
    sortidx = np.argsort(probs)[::-1]
    sorted_bstrs = np.array(bstrs)[sortidx]
    success = False
    if sorted_bstrs[0] is spins2bitstr(answer):
        success = True
    elif ((sorted_bstrs[1] is spins2bitstr(answer)) and 
          (sorted_bstrs[0] is spins2bitstr(instate))):
        success = True

    return success, avgdist, len(mems)

# Set up plot
pl.xlabel(r'$|P|$')
pl.ylabel(r'$Avg. Hamming distance$')

# Loop through all data directories
for root, dirs, files in os.walk('.'):
    if root == '.':
        continue
    # If we are in a dir with no children..
    if dirs == [] and (int(root.split('n')[1].split('p')[0]) == qubits):
        os.chdir(root)
        print root
        success, dist, kmems = analysis(probsfname, binary)
        pl.scatter(kmems, dist, c=('r' if success is False else 'b'))
        os.chdir('../../')

if output:
    pl.savefig('p_hamming_plot.png')
pl.show()
