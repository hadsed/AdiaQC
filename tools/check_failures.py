'''
File: check_failures.py
Author: Hadayat Seddiqi
Date: 8.17.14
Description: Check failed instances to see if their encodings
             would have given the correct answer if not for a
             bad biasing.
'''

import os, optparse, sys
import json
import cPickle as pickle
import numpy as np
import scipy as sp
from scipy import linalg


def spins2bitstr(vec):
    """ Return a converted bitstring from @vec, a list of spins. """
    return ''.join([ '0' if k == 1 else '1' for k in vec ])

def hamdist(a,b):
    """ Calculate Hamming distance. """
    return np.sum(abs(np.array(a)-np.array(b))/2.0)

def fsuccess(bstrs, probs, answer, thresh):
    """
    Success only if the top state is the answer.
    """
    sortidx = np.argsort(probs)[::-1]
    sorted_bstrs = np.array(bstrs)[sortidx]
    success = False
    prob = probs[sortidx][sorted_bstrs == spins2bitstr(answer)]
    # Check if we have degenerate states that are XORs of each other
    t1 = np.array([ int(k) for k in sorted_bstrs[0] ])
    t2 = np.array([ int(k) for k in sorted_bstrs[1] ])
    if (t1^1 == t2).all():
        prob = probs[sortidx][0] + probs[sortidx][1]
        if (sorted_bstrs[0] == spins2bitstr(answer)
            or sorted_bstrs[1] == spins2bitstr(answer)) and (prob > thresh):
           success = True
    # Was it the top state, and was its probability high enough?
    elif sorted_bstrs[0] == spins2bitstr(answer) and prob > thresh:
        success = True
    return success, prob

def analysis(filename, binary, threshold):
    """
    """
    # Get probability data
    if binary:
        probs = np.load(filename+".dat.npy")
    else:
        probs = np.loadtxt(filename+".dat")

    # Get state labels
    bstrs = [ line.rstrip('\n') for line in open('statelabels.txt') ]

    # Get problem outputs
    props = json.load(open('problem_outputs.dat'))
    answer = props['answer']
    gamma = props['gamma']
    mems = props['memories']
    inp = props['inputState']

    # Calculate success
    success, sprob = fsuccess(bstrs, probs, answer, threshold)

    return success, [sprob, gamma, mems, inp, answer]


# Command line options
if __name__=="__main__":
    parser = optparse.OptionParser("usage: %prog [options] arg1 arg2")
    parser.add_option("-q", "--qubits", dest="qubits", default=2,
                      type="int", 
                      help="Number of qubits.")
    parser.add_option("-t", "--threshold", dest="threshold", default=0.66666666,
                      type="float", 
                      help="Threshold probability for determining success "+\
                          "(0 < thresh < 1).")
    parser.add_option("-b", "--binary", dest="binary", default=1,
                      type="int", 
                      help="Whether data file is in binary (1) or not (0).")
    parser.add_option("-f", "--fname", dest="fname", default=None,
                      type="string", 
                      help="Filename for output probabilities.")
    parser.add_option("-g", "--getdata", dest="getdata", default=0,
                      type="int", 
                      help="Do we need to get the success/failure data.")
    parser.add_option("-l", "--lrule", dest="lrule", default=0,
                      type="int", 
                      help="Learning rule to plot (0: all, 1: hebb, "+\
                          "2: stork, 3: proj).")
    parser.add_option("-n", "--nsamples", dest="nsamples", default=250,
                      type="float", 
                      help="Number of samples (to normalize f_x).")
    (options, args) = parser.parse_args()
    qubits = options.qubits
    threshold = options.threshold
    binary = options.binary
    probsfname = options.fname
    getdataopt = options.getdata
    lrule = options.lrule
    nsamples = options.nsamples

    if lrule == 0:
        lrlist = ['hebb', 'proj', 'stork']
    elif lrule == 1:
        lrlist = ['hebb']
    elif lrule == 2:
        lrlist = ['stork']
    elif lrule == 3:
        lrlist = ['proj']
    elif lrule == 4:
        lrlist = ['hebb', 'stork']
    elif lrule == 5:
        lrlist = ['hebb', 'proj']
    elif lrule == 6:
        lrlist = ['proj', 'stork']
    else:
        print("Learning rule number not recognized.")
        sys.exit(0)

    if getdataopt:
        # We're going to output properties of failed instances in the end
        fail_list = []
        pref = 'n'+str(qubits)+'p'
        # Loop over pattern numbers
        for pnum in range(qubits):
            print("Pnum: ", str(pnum))
            # Loop over learning rules
            for lr in lrlist:
                os.chdir(pref+str(pnum+1)+lr)
                # Go through all instances
                for root, dirs, files in os.walk('.'):
                    if root == '.':
                        continue
                    os.chdir(root)
                    # Do what we came here to do
                    # try:
                    success, props = analysis(probsfname, 
                                              binary, 
                                              threshold)
                    if not success:
                        # Add some properties to a list we will
                        # output at the end to a file
                        fail_list.append([root[2:], lr]+props)
                    # except:
                    #     print("Failure at: ", root)
                    #     os.chdir('../')
                    #     continue
                    os.chdir('../')
                os.chdir('../')

        # Save
        pickle.dump(fail_list, open('failure_list_n'+str(qubits)+'.pk', 'w'))
        print("Exiting. Run with -g 0 to do analysis.")
    else:
        def problem_params(lr, gam, memories, inpst, neurons):
            """
            Return the lowest eigenvector of the classical Hamiltonian
            constructed by the learning rule, gamma, memories, and input.
            """
            # Bias Hamiltonian
            alpha = gam*np.array(inpst)

            # Memory Hamiltonian
            beta = np.zeros((qubits, qubits))
            if lr == 'hebb':
                # Hebb rule
                memMat = sp.matrix(memories).T
                beta = sp.triu(memMat*memMat.T)/float(neurons)
            elif lr == 'stork':
                # Storkey rule
                Wm = sp.zeros((neurons,neurons))
                for m, mem in enumerate(memories):
                    Am = sp.outer(mem,mem) - sp.eye(neurons)
                    Wm += (Am - Am*Wm - Wm*Am)/float(neurons)
                beta = sp.triu(Wm)
            elif lr == 'proj':
                # Moore-Penrose pseudoinverse rule
                memMat = sp.matrix(memories).T
                beta = sp.triu(memMat * sp.linalg.pinv(memMat))

            # Find the eigenvectors
            evals, evecs = sp.linalg.eig(np.diag(alpha)+beta)
            idx = evals.argsort()

            return evals[idx], evecs[:,idx], np.diag(alpha), beta

        # Load data
        data = pickle.load(open('failure_list_n'+str(qubits)+'.pk', 'r'))

        #
        # The data is structured like so:
        #    [ instance_num, lrule, sprob, gamma, mems, inp, answer ]
        #

        # Loop through them
        for row in data[::-1]:
            instance_num, lrule, sprob, gamma, mems, inp, answer = row
            neurons = len(inp)
            print "\n\nInstance:"
            print instance_num
            energies, states, alpha, beta  = problem_params(lrule, 
                                                            gamma, 
                                                            mems, 
                                                            inp, 
                                                            neurons)
            print 'energies'
            print energies
            print 'states'
            for col in range(neurons):
                vec0 = np.matrix(np.around(states[:,col], 10)).T
                print vec0.T
                print -np.sum(vec0.T*beta*vec0)
                print -np.sum(vec0.T*beta*vec0) - np.sum(alpha)
                print ""

            for m in mems:
                km = np.matrix(m).T
                print -np.sum(km.T*beta*km), m

            print "\nProperties: \n"
            print lrule
            print gamma
            print 'memories:'
            print mems
            print 'input: '
            print inp
            print 'intended answer: '
            print answer
            # print 1-enc_ans

            break
