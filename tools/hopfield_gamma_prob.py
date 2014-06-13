'''
File: hopfield_gamma_prob.py
Author: Hadayat Seddiqi
Date: 5.10.14
Description: Show the data of a single Hopfield instance
             where the bias term, Gamma, is varied. Make
             a plot of that vs. the answer state probability.
'''

import os, optparse, sys
import json
import cPickle as pickle
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as pl
import matplotlib.colors as mcolors
import matplotlib.cm as cm

def spins2bitstr(vec):
    """ Return a converted bitstring from @vec, a list of spins. """
    return ''.join([ '0' if k == 1 else '1' for k in vec ])

def hamdist(a,b):
    """ Calculate Hamming distance. """
    return np.sum(abs(np.array(a)-np.array(b))/2.0)

def fsuccess(bstrs, probs, answer, thresh):
    """
    Success only if the top state is the answer, or it is
    in a cluster of degenerate ground states.
    """
    sortidx = np.argsort(probs)[::-1]
    sorted_bstrs = np.array(bstrs)[sortidx]
    success = False
    prob = probs[sortidx][sorted_bstrs == spins2bitstr(answer)]
    # Was it the top state, and was its probability high enough? Or do
    # we have a degenerate cluster of states at the top?
    if((sorted_bstrs[0] == spins2bitstr(answer) or 
        np.isclose(probs[sortidx][0], prob)) 
       and prob > thresh):
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
    gamma = props['bias']

    # Calculate success
    success, sprob = fsuccess(bstrs, probs, answer, threshold)

    return success, sprob, gamma


# Command line options
if __name__=="__main__":
    parser = optparse.OptionParser("usage: %prog [options] arg1 arg2")
    parser.add_option("-q", "--qubits", dest="qubits", default=2,
                      type="int", 
                      help="Number of qubits.")
    parser.add_option("-t", "--threshold", dest="threshold", default=0.7,
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
                      help="Learning rule to plot (0: all, 1: hebb, 2: stork, 3: proj).")
    parser.add_option("-n", "--nsamples", dest="nsamples", default=250,
                      type="float", 
                      help="Number of samples (to normalize f_x).")
    parser.add_option("-p", "--patterns", dest="patterns", default=1,
                      type="int", 
                      help="Number of patterns.")
    (options, args) = parser.parse_args()
    qubits = options.qubits
    threshold = options.threshold
    binary = options.binary
    probsfname = options.fname
    getdataopt = options.getdata
    lrule = options.lrule
    nsamples = options.nsamples
    patterns = options.patterns

    # This needs to be hardcoded for now
    gammarng = np.arange(0.0, 1.05, 0.01)
    lengamma = len(gammarng)

    if getdataopt:
        data = {'hebb':  np.zeros(lengamma),
                'stork': np.zeros(lengamma),
                'proj':  np.zeros(lengamma)}
        pref = 'n'+str(qubits)+'p'+str(patterns)
        # Loop over pattern numbers
        for lr in ['hebb', 'stork', 'proj']:
            os.chdir(pref+lr)
            for root, dirs, files in os.walk('.'):
                if root == '.':
                    continue
                os.chdir(root)
                try:
                    success, sprob, gamma = analysis(probsfname, binary, threshold)
                except:
                    print("Failure at: ", root)
                    os.chdir('../')
                    continue
                # if success:
                if True:
                    try:
                        gidx = np.where(np.isclose(gammarng, gamma))[0][0]
                        # sidx = int(root[2:]) % lengamma
                    except:
                        os.chdir('../')
                        continue
                    data[lr][gidx] = sprob
                os.chdir('../')
            os.chdir('../')

        # Save
        pickle.dump(data, open('gamma_prob_n'+str(qubits)+'p'+str(patterns)+'.pk', 'w'))
        print("Exiting. Run with -g 0 to get the plots.")
    else:
        # Load data
        data = pickle.load(open('gamma_prob_n'+str(qubits)+'p'+str(patterns)+'.pk', 'r'))

        pl.rc('text', usetex=True)
        pl.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
        pl.rcParams['xtick.labelsize'] = 24
        pl.rcParams['ytick.labelsize'] = 24
        lwidth = 2.5
        lstyle = '-'
        marker = '.'
        msize = 0
        color = ['green', 'red', 'blue']

        pl.rcParams['figure.figsize'] = 9,7
        # pl.title(r'$\textbf{Fractional success as a function of input bias'+\
        #              ' (N = '+str(qubits)+')}$', fontsize=18)
        pl.xlabel(r'$\boldsymbol{\Gamma}$', fontsize=20)
        pl.ylabel(r'$|\langle \phi_{ans} | \psi \rangle|^2$', fontsize=24)
        pl.ylim([0.0, 1.05])
        pl.xlim([-0.05, 1.05])

        for ilr, lr in enumerate(['proj', 'hebb', 'stork']):
            if lr == 'stork':
                lstyle = '--'
            else:
                lstyle = '-'
            pl.plot(gammarng, data[lr], c=color[ilr], label=lr, marker=marker, 
                    markersize=msize, linewidth=lwidth, linestyle=lstyle)
        pl.legend(prop={'size':12}, loc=1)
        # pl.show()
        pl.savefig('gamma_prob_n'+str(qubits)+'p'+str(patterns)+'.png')
