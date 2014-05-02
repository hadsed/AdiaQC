'''
File: failure_analysis.py
Author: Hadayat Seddiqi
Date: 4.28.14
Description: Take failed instances from Hopfield experiments #1
             and show some things about it.
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
        try:
            probs = np.load(filename)
        except IOError, e:
            print "IOError: No file "+filename+".npy. \nYou probably forgot "+\
                "to specify the appropriate command-line argument: -b 0."
            sys.exit(0)
    else:
        probs = np.loadtxt(filename)
    # probstr = [ '%1.16E' % e for e in probs.ravel() ]
    bstrs = [ line.rstrip('\n') for line in open('statelabels.txt') ]

    # Get the eigenspectrum
    if binary:
        try:
            eigs = np.load('eigenspectrum.dat.npy')
        except IOError, e:
            print "IOError: No file eigenspectrum.dat.npy. \nYou probably forgot "+\
                "to specify the appropriate command-line argument: -b 0."
            sys.exit(0)
    else:
        eigs = np.loadtxt('eigenspectrum.dat')

    # Organize the eigenstuff
    times = eigs[:,0]/eigs[-1,0]
    eigs = eigs[:,1:]

    # Get problem outputs
    props = json.load(open('problem_outputs.dat'))
    lrule = props['learningRule']
    mems = props['memories']
    instate = props['inputState']
    answer = props['answer']

    # Calculate average Hamming distance
    avgdist = np.average([ hamdist(m, instate) for m in mems ])
    
    # Calculate success
    sortidx = np.argsort(probs)[::-1]
    sorted_bstrs = np.array(bstrs)[sortidx]
    success = False
    prob = 0.0
    if sorted_bstrs[0] == spins2bitstr(answer):
        success = True
        prob = probs[sortidx][0]
    elif ((sorted_bstrs[1] == spins2bitstr(answer)) and 
          (sorted_bstrs[0] == spins2bitstr(instate))):
        success = True
        prob = probs[sortidx][1]

    return success, avgdist, len(mems), prob, lrule, times, eigs

# Initialize some variables we'll want to look at
csuccess = {'hebb': [0]*qubits,
            'stork': [0]*qubits,
            'proj': [0]*qubits}
cfailure = {'hebb': [0]*qubits,
            'stork': [0]*qubits,
            'proj': [0]*qubits}
data = []
faildata = []
# Loop through all data directories
for root, dirs, files in os.walk('.'):
    if root == '.':
        continue
    # If we are in a dir with no children..
    if dirs == [] and (int(root.split('n')[1].split('p')[0]) == qubits):
        os.chdir(root)
        success, dist, kmems, prob, lrule, times, eigs = analysis(probsfname, binary)
        if success:
            continue  # skip successful guys for now
            # plot eigenspectrum
            pl.clf()
            for col in xrange(eigs.shape[1]):
                pl.plot(times, eigs[:,col])
            pl.title('Eigenspectrum in '+root)
            pl.xlabel('Time')
            pl.ylabel('Energy')
            pl.show()
            # csuccess[lrule][kmems-1] += 1
            # data.append([ success, dist, kmems, prob, lrule ])
        else:
            # plot eigenspectrum
            pl.clf()
            for col in xrange(eigs.shape[1]):
                pl.plot(time, eigs[:,col])
            pl.title('Eigenspectrum in '+root)
            pl.xlabel('Time')
            pl.ylabel('Energy')
            pl.show()
            # cfailure[lrule][kmems-1] += 1
            # faildata.append([ success, dist, kmems, prob, lrule ])

        os.chdir('../../')
