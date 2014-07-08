'''
File: exp4_line.py
Author: Hadayat Seddiqi
Date: 6.23.14
Description: Show plot of varying Gamma over counts of
             incorrect retrieval of the input state when
             it's not included in the memory set.
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
    inp = props['inputState']
    gamma = props['gamma']

    # Calculate success
    # success, sprob = fsuccess(bstrs, probs, answer, threshold)
    success, sprob = fsuccess(bstrs, probs, inp, threshold)

    return success, sprob, gamma


# Command line options
if __name__=="__main__":
    parser = optparse.OptionParser("usage: %prog [options] arg1 arg2")
    parser.add_option("-q", "--qubits", dest="qubits", default=2,
                      type="int", 
                      help="Number of qubits.")
    parser.add_option("-t", "--threshold", dest="threshold", default=0.0,
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
    parser.add_option("-d", "--hdist", dest="hdist", default=1,
                      type="int", 
                      help="Only affects filenaming (which Hamming distance did we use).")
    parser.add_option("-l", "--lrule", dest="lrule", default=0,
                      type="int", 
                      help="Learning rule to plot (0: all, 1: hebb, 2: stork, 3: proj).")
    parser.add_option("-n", "--nsamples", dest="nsamples", default=250,
                      type="float", 
                      help="Number of samples (to normalize f_x).")
    (options, args) = parser.parse_args()
    qubits = options.qubits
    threshold = options.threshold
    binary = options.binary
    probsfname = options.fname
    getdataopt = options.getdata
    hdist = options.hdist
    lrule = options.lrule
    nsamples = options.nsamples

    # This needs to be hardcoded for now
    # gammarng = [ 0.01, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0 ]
    gammarng = np.arange(0.0, 1.0, 0.05)
    lengamma = len(gammarng)

    if getdataopt:
        # data = {'hebb': [ [0]*lengamma for k in range(qubits) ],
        #         'stork': [ [0]*lengamma for k in range(qubits) ],
        #         'proj': [ [0]*lengamma for k in range(qubits) ]}
        data = {'hebb':  [ np.zeros((lengamma, nsamples)) for k in range(qubits) ],
                'stork': [ np.zeros((lengamma, nsamples)) for k in range(qubits) ],
                'proj':  [ np.zeros((lengamma, nsamples)) for k in range(qubits) ]}
        pref = 'n'+str(qubits)+'p'
        # Loop over pattern numbers
        for pnum in range(qubits):
            print("Pnum: ", str(pnum))
            for lr in ['hebb', 'stork', 'proj']:
                os.chdir(pref+str(pnum+1)+lr)
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
                    if success:
                        try:
                            gidx = np.where(np.isclose(gammarng, gamma))[0][0]
                            sidx = int(root[2:]) % lengamma
                        except:
                            os.chdir('../')
                            continue
                        data[lr][pnum][gidx, sidx] += 1
                    os.chdir('../')
                os.chdir('../')

        # Save
        pickle.dump(data, open('exp4_line_data_n'+str(qubits)+'_hd'+str(hdist)+'.pk', 'w'))
        print("Exiting. Run with -g 0 to get the plots.")
    else:
        # Load data
        data = pickle.load(open('exp4_line_data_n'+str(qubits)+'_hd'+str(hdist)+'.pk', 'r'))
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
        fig, (ax1, ax2, ax3) = pl.subplots(1, 3, sharey=True, figsize=(24,11))
        fig.subplots_adjust(wspace=0.05, left=0.05, right=0.95, bottom=0.15)
        # Hebb
        ax1.set_title('Hebb', fontsize=24)
        ax1.set_ylabel(r'$\textbf{Failure Probability}$', fontsize=24)
        ax1.yaxis.set_label_coords(-0.1, 0.45)
        ax1.set_ylim([0.0, 1.05])
        ax1.set_xlim([-0.05, 1.0])
        for pnum in range(qubits):
            ax1.plot(gammarng, np.sum(data['hebb'][pnum], axis=1)/nsamples, lstyle, 
                     c=colors[pnum], linewidth=lwidth,
                     marker=marker, markersize=msize)
        # Stork
        ax2.set_title('Storkey', fontsize=24)
        ax2.set_xlim([-0.05, 1.0])
        ax2.set_xlabel(r'$\mathbf{\Gamma}$', fontsize=24)
        for pnum in range(qubits):
            ax2.plot(gammarng, np.sum(data['stork'][pnum], axis=1)/nsamples, lstyle, 
                     c=colors[pnum], linewidth=lwidth, 
                     label=str(pnum+1)+" patterns",
                     marker=marker, markersize=msize)
        # Proj
        ax3.set_title('Projection', fontsize=24)
        ax3.set_xlim([-0.05, 1.0])
        for pnum in range(qubits):
            ax3.plot(gammarng, np.sum(data['proj'][pnum], axis=1)/nsamples, lstyle, 
                     c=colors[pnum], linewidth=lwidth, 
                     marker=marker, markersize=msize)
        # ax2.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., prop={'size':12})
        # ax2.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
        #           ncol=4, mode="expand", borderaxespad=0.)
        ax2.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1),
                   fancybox=True, shadow=True, ncol=10)
        # pl.show()
        pl.savefig('exp4_line_n'+str(qubits)+'_hd'+str(hdist)+'.png')
