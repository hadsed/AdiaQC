'''
File: failed_inst_listgen.py
Author: Hadayat Seddiqi
Date: 5.20.14
Description: Find the failed instances and output their paths
             to a file for further processing (e.g. gamma tuning).
'''

import os, optparse, sys
import numpy as np
import json
import pylab as pl


def spins2bitstr(vec):
    """ Return a converted bitstring from @vec, a list of spins. """
    return ''.join([ '0' if k == 1 else '1' for k in vec ])

def success1(bstrs, probs, answer, instate):
    """
    Success defined by looking at the top two states
    and making sure if the top one isn't the answer,
    then it is the input state and the second top
    answer is the true answer, otherwise failure.
    """
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
    return success, prob

def success2(bstrs, probs, answer, instate):
    """
    Success only if the top state is the answer.
    """
    sortidx = np.argsort(probs)[::-1]
    sorted_bstrs = np.array(bstrs)[sortidx]
    success = False
    prob = probs[sortidx][sorted_bstrs == spins2bitstr(answer)]
    if sorted_bstrs[0] == spins2bitstr(answer):
        success = True
    return success, prob

def analysis(filename, binary):
    """
    Grab the right data when inside a problem instance directory.
    """
    # Get probability and eigspec data
    if binary:
        try:
            probs = np.load(filename)
            eigs = np.load('eigenspectrum.dat.npy')
        except (IOError):
            print ("IOError: No file "+filename+". \nYou probably forgot "+\
                "to specify the appropriate command-line argument: -b 0.")
            print ("Directory:")
            print (os.getcwd())
            sys.exit(0)
    else:
        try:
            probs = np.loadtxt(filename)
            eigs = np.loadtxt('eigenspectrum.dat')
        except (IOError):
            print ("IOError: No file "+filename+". \nYou probably forgot "+\
                "to specify the appropriate command-line argument: -b 0.")
            print ("Directory:")
            print (os.getcwd())
            sys.exit(0)

    # probstr = [ '%1.16E' % e for e in probs.ravel() ]
    bstrs = [ line.rstrip('\n') for line in open('statelabels.txt') ]
    # Get problem outputs
    props = json.load(open('problem_outputs.dat'))
    answer = props['answer']
    instate = None
    # Calculate success
    success, sprob = success2(bstrs, probs, answer, instate)
    return success

def getdata(qubits, binary, probsfname):
    """
    Walk through all directories and get relevant data, outputting to files. Calls
    analysis().
    """
    # Create file
    output_fname = 'failed_instances_n'+str(qubits)+'.txt'
    open(output_fname, 'w').close()
    # Loop through all data directories
    for root, dirs, files in os.walk('.'):
        if root == '.':
            continue
        # If we are in a dir with no children..
        if dirs == [] and (int(root.split('n')[1].split('p')[0]) == qubits):
            os.chdir(root)
            success = analysis(probsfname, binary)
            if not success:
                # output to file
                with open('../../'+output_fname, 'a') as outfile:
                    outfile.write(os.getcwd()+'\n')
            os.chdir('../../')


if __name__=="__main__":
    # Command line options
    parser = optparse.OptionParser("usage: %prog [options] arg1 arg2")
    parser.add_option("-f", "--fname", dest="probsfname", default=None,
                      type="string", 
                      help="Filename for output probabilities.")
    parser.add_option("-b", "--binary", dest="binary", default=1,
                      type="int", 
                      help="Binary input or not (default = 1).")
    parser.add_option("-q", "--qubits", dest="qubits", default=None,
                      type="int", 
                      help="Number of qubits/neurons.")
    (options, args) = parser.parse_args()
    probsfname = options.probsfname
    binary = options.binary
    qubits = options.qubits
    # Grab the data
    getdata(qubits, binary, probsfname)
