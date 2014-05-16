'''
File: exp1_phase.py
Author: Hadayat Seddiqi
Date: 5.3.14
Description: Plot a phase diagram between the answer
             and input state.
'''

import os, optparse, sys
import numpy as np
import matplotlib.pyplot as pl
import itertools
import json


def spins2bitstr(vec):
    """ Return a converted bitstring from @vec, a list of spins. """
    return ''.join([ '0' if k == 1 else '1' for k in vec ])

def analyze(filename, binary):
    """
    Grab the right data when inside a problem instance directory.
    """
    # Get probability and eigspec data
    if binary:
        try:
            probs = np.load(filename)
            eigs = np.load('eigenspectrum.dat.npy')
        except (IOError):
            print ("IOError: No file "+filename+".npy. \nYou probably forgot "+\
                "to specify the appropriate command-line argument: -b 0.")
            sys.exit(0)
    else:
        probs = np.loadtxt(filename)
        eigs = np.loadtxt('eigenspectrum.dat')

    # probstr = [ '%1.16E' % e for e in probs.ravel() ]
    bstrs = [ line.rstrip('\n') for line in open('statelabels.txt') ]

    # Organize the eigenstuff
    times = eigs[:,0]/eigs[-1,0]
    eigs = eigs[:,1:]
    specgap = eigs[:,1]-eigs[:,0]

    # Get state overlap labels
    axis = []
    with open('state_overlap_labels.txt') as f:
        axis = f.readlines()
    axis = [ line.rstrip('\n') for line in axis ]

    # Cycle through state overlap files
    data = []
    ftimes = []
    for fname in os.listdir('.'):
        if fname[:15] == 'state_overlap_T' and fname[-4:] == '.txt':
            prob = np.loadtxt(fname)
            t = fname[15:].split('.txt')[0].split('T')[0]
            ftimes.append(float(t))
            data.append(prob)

    # Fix up the state overlap data
    data = zip(*data)
    sortidx = np.array(ftimes).argsort()
    for icol, col in enumerate(data):
        data[icol] = np.array(col)[sortidx]
    time = np.linspace(0,1,len(data[0]))

    # Figure out what states we want
    props = json.load(open('problem_outputs.dat'))
    answer = props['answer']
    instate = props['inputState']
    answerIdx = axis.index(spins2bitstr(answer))
    inputIdx = axis.index(spins2bitstr(instate))
    print spins2bitstr(answer), spins2bitstr(instate), axis
    print answerIdx, inputIdx
    # Phase diagram tuples
    phase = (data[answerIdx], data[inputIdx])

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
    else:
        pidx = np.where(sorted_bstrs == spins2bitstr(answer))
        prob = probs[sortidx][pidx]

    return success, prob, times, specgap, phase


# Command line options
if __name__=="__main__":
    parser = optparse.OptionParser("usage: %prog [options] arg1 arg2")
    parser.add_option("-q", "--qubits", dest="qubits", default=None,
                      type="int", 
                      help="Number of qubits/neurons.")
    parser.add_option("-b", "--binary", dest="binary", default=1,
                      type="int", 
                      help="Whether data file is in binary (1) or not (0).")
    (options, args) = parser.parse_args()
    binary = options.binary
    qubits = options.qubits

    # Loop through all data directories
    for root, dirs, files in os.walk('.'):
        if root == '.':
            continue
        # If we are in a dir with no children..
        if dirs == [] and (int(root.split('n')[1].split('p')[0]) == qubits):
            # Get in there
            os.chdir(root)
            # Get some data
            success, prob, times, specgap, phase = analyze('probsT15.0.dat', binary)
            # Plot what we see
            color = 'r'
            if success:
                color = 'b'
            pl.plot(data[answerIdx], data[inputIdx], color, alpha=0.5)

            # Head back out
            os.chdir('../../')
    pl.title('Phase diagram of state vector over time', fontweight='bold')
    pl.xlabel('Answer State Probability')
    pl.ylabel('Input State Probability')
    pl.show()

