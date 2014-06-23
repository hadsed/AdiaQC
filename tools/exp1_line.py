'''
File: exp1_line.py
Author: Hadayat Seddiqi
Date: 5.23.14
Description: Capacity lineplots, |P| vs. capacity (fractional success f_x)
'''

import os, optparse, sys
import numpy as np
import json
import pylab as pl
import cPickle as pickle

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
        probs = np.load(filename)
    else:
        probs = np.loadtxt(filename)
    # probstr = [ '%1.16E' % e for e in probs.ravel() ]
    bstrs = [ line.rstrip('\n') for line in open('statelabels.txt') ]

    # Get problem outputs
    props = json.load(open('problem_outputs.dat'))
    lrule = props['learningRule']
    mems = props['memories']
    answer = props['answer']

    # Success/failure?
    success, sprob = fsuccess(bstrs, probs, answer, threshold)

    return success, len(mems), sprob, lrule

if __name__=="__main__":
    # Command line options
    parser = optparse.OptionParser("usage: %prog [options] arg1 arg2")
    parser.add_option("-p", "--probsfname", dest="probsfname", default=None,
                      type="string", 
                      help="Filename for output probabilities.")
    parser.add_option("-t", "--threshold", dest="threshold", default=0.7,
                      type="float", 
                      help="Threshold probability for determining success "+\
                          "(0 < thresh < 1).")
    parser.add_option("-o", "--output", dest="output", default=1,
                      type="int", 
                      help="Output this report to file.")
    parser.add_option("-b", "--binary", dest="binary", default=1,
                      type="int", 
                      help="Binary input or not (default = 1).")
    parser.add_option("-q", "--qubits", dest="qubits", default=None,
                      type="int", 
                      help="Number of qubits/neurons.")
    parser.add_option("-g", "--getdata", dest="getdata", default=0,
                      type="int", 
                      help="Do we need to get the success/failure data.")
    parser.add_option("-n", "--nsamples", dest="nsamples", default=1000,
                      type="float", 
                      help="Number of samples (to normalize f_x).")
    (options, args) = parser.parse_args()
    probsfname = options.probsfname
    threshold = options.threshold
    output = options.output
    binary = options.binary
    qubits = options.qubits
    getdataopt = options.getdata
    nsamples = options.nsamples

    if getdataopt:
        # Initialize some variables we'll want to look at
        data = {'hebb': [0]*qubits,
                'stork': [0]*qubits,
                'proj': [0]*qubits}

        # Loop through all data directories
        for root, dirs, files in os.walk('.'):
            print root
            if root == '.':
                continue
            # If we are in a dir with no children..
            if (dirs == [] and 
                (int(root.split('n')[1].split('p')[0]) == qubits) and
                (probsfname in files)):
                os.chdir(root)
                try:
                    success, kmems, prob, lrule = analysis(probsfname, binary, threshold)
                except IOError, e:
                    print "\nFailure on: " + root
                    print e
                    continue
                if success:
                    data[lrule][kmems-1] += 1
                os.chdir('../../')
        # Save
        pickle.dump(data, open('success_line_n'+str(qubits)+'.pk', 'w'))
        print("Exiting. Run with -g 0 to get the plots.")
    else:
        # Load data
        data = pickle.load(open('success_line_n'+str(qubits)+'.pk', 'r'))
        prng = np.arange(qubits)+1
        # Plot settings
        lwidth = 2.5
        lstyle = '-'
        marker = 'o'
        msize = 3
        pl.rcParams['xtick.major.pad'] = 8
        pl.rcParams['ytick.major.pad'] = 8
        pl.rcParams['xtick.labelsize'] = 24
        pl.rcParams['ytick.labelsize'] = 24
        pl.rc('text', usetex=True)
        pl.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
        colors = ['black', 'blue', 'green', 'orange', 'red']*2
        # Create figure
        fig, ax = pl.subplots(figsize=(10,8))
        fig.subplots_adjust(wspace=0.05, bottom=0.18)
        # Hebb
        ax.set_ylabel(r'$\langle f_x \rangle$', fontsize=24, rotation='horizontal')
        ax.yaxis.set_label_coords(-0.1, 0.45)
        ax.set_ylim([0.0, 1.05])
        ax.set_xticks(range(1,qubits+1))
        ax.set_xlabel(r'$\textbf{P}$', fontsize=24)
        ax.plot(prng, np.array(data['hebb'])/nsamples, lstyle, 
                c='r', linewidth=lwidth, label='Hebb',
                marker=marker, markersize=msize)
        # Stork
        ax.plot(prng, np.array(data['stork'])/nsamples, lstyle, 
                c='b', linewidth=lwidth, label='Storkey',
                marker=marker, markersize=msize)
        # Proj
        ax.plot(prng, np.array(data['proj'])/nsamples, lstyle, 
                c='g', linewidth=lwidth, label='Projection',
                marker=marker, markersize=msize)
        # ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., prop={'size':12})
        # ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
        #           ncol=4, mode="expand", borderaxespad=0.)
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
                  fancybox=True, shadow=True, ncol=3)
        # pl.show()
        pl.savefig('fx_vs_pnum_n'+str(qubits)+'.png')
