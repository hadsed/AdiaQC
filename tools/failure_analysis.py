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


def spins2bitstr(vec):
    """ Return a converted bitstring from @vec, a list of spins. """
    return ''.join([ '0' if k == 1 else '1' for k in vec ])

def hamdist(a,b):
    """ Calculate Hamming distance. """
    return np.sum(abs(np.array(a)-np.array(b))/2.0)

def lrule2int(lrule):
    """
    Convert learning rule string to int.
    """
    if lrule == 'hebb':
        return 0
    if lrule == 'stork':
        return 1
    if lrule == 'proj':
        return 2

def analysis(filename, binary):
    """
    Grab the right data when inside a problem instance directory.
    """
    # Get probability and eigspec data
    if binary:
        try:
            probs = np.load(filename)
            eigs = np.load('eigenspectrum.dat.npy')
        except IOError, e:
            print "IOError: No file "+filename+".npy. \nYou probably forgot "+\
                "to specify the appropriate command-line argument: -b 0."
            sys.exit(0)
    else:
        probs = np.loadtxt(filename)
        eigs = np.loadtxt('eigenspectrum.dat')

    # probstr = [ '%1.16E' % e for e in probs.ravel() ]
    bstrs = [ line.rstrip('\n') for line in open('statelabels.txt') ]

    # Organize the eigenstuff
    times = eigs[:,0]/eigs[-1,0]
    eigs = eigs[:,1:]
    mgIdx = np.argmax(eigs[:,1]-eigs[:,0])
    mingap = eigs[mgIdx,1]-eigs[mgIdx,0]
    mgt = times[mgIdx]

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
    else:
        pidx = np.where(sorted_bstrs == spins2bitstr(answer))
        prob = probs[sortidx][pidx]

    return success, avgdist, len(mems), prob, lrule, times, eigs, mingap, mgt

def getdata(qubits):
    """
    Walk through all directories and get relevant data, outputting to files. Calls analysis().
    """

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
            success, dist, kmems, prob, lrule, times, eigs, mingap, mgt = \
                analysis(probsfname, binary)
            if success:
                data.append([ dist, prob, lrule2int(lrule), mingap, mgt ])
                # plot eigenspectrum
                # pl.clf()
                # for col in xrange(eigs.shape[1]):
                #     pl.plot(times, eigs[:,col])
                # pl.title('Eigenspectrum in '+root)
                # pl.xlabel('Time')
                # pl.ylabel('Energy')
                # pl.show()
                # csuccess[lrule][kmems-1] += 1
                # data.append([ success, dist, kmems, prob, lrule ])
            else:
                faildata.append([ dist, prob, lrule2int(lrule), mingap, mgt])
                # plot eigenspectrum
                # pl.clf()
                # for col in xrange(eigs.shape[1]):
                #     pl.plot(times, eigs[:,col])
                # pl.title('Eigenspectrum in '+root)
                # pl.xlabel('Time')
                # pl.ylabel('Energy')
                # pl.show()
                # cfailure[lrule][kmems-1] += 1
                # faildata.append([ success, dist, kmems, prob, lrule ])

            os.chdir('../../')

    np.save('successdata_n'+str(qubits)+'.dat')
    np.save('failuredata_n'+str(qubits)+'.dat')

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

    # Grab the data
    getdata(qubits)

