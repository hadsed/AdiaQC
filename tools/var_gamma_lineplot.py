'''
File: var_gamma_lineplot.py
Author: Hadayat Seddiqi
Date: 5.22.14
Description: Show the "fractional success" (capacity) over Gamma
             for various problem types (n,p) for different learning
             rules.
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
    gamma = props['gamma']

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
    (options, args) = parser.parse_args()
    qubits = options.qubits
    threshold = options.threshold
    binary = options.binary
    probsfname = options.fname
    getdataopt = options.getdata
    lrule = options.lrule
    nsamples = options.nsamples

    # This needs to be hardcoded for now
    # gammarng = [ 0.01, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0 ]
    gammarng = np.arange(0.0, 2.0, 0.05)
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
        pickle.dump(data, open('var_gamma_n'+str(qubits)+'.pk', 'w'))
        print("Exiting. Run with -g 0 to get the plots.")
    else:
        # Load data
        data = pickle.load(open('var_gamma_n'+str(qubits)+'.pk', 'r'))

        pl.rc('text', usetex=True)
        pl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
        lwidth = 1.8
        lstyle = '-'
        marker = '.'
        msize = 6
        color = ['black', 'blue', 'green', 'orange', 'red']

        if lrule == 0:
            # Create figure
            fig, (ax1, ax2, ax3) = pl.subplots(1, 3, sharey=True, figsize=(20,9))
            fig.subplots_adjust(wspace=0.05, left=0.05, right=0.95, bottom=0.15)
            fig.suptitle(r'$\textbf{Fractional success as a function of input bias'+\
                         ' (N = '+str(qubits)+')}$', fontsize=18)
            pl.ylim([0.0,1.05])
            # Hebb
            ax1.set_title('Hebb')
            ax1.set_ylabel(r'$\langle f_x \rangle$', fontsize=24, rotation='horizontal')
            ax1.set_xlim([-0.05, 2.05])
            for pnum in range(qubits):
                ax1.plot(gammarng, np.array(data['hebb'][pnum])/nsamples, lstyle, 
                         c=color[pnum], linewidth=lwidth, label=str(pnum+1)+' patterns', 
                         marker=marker, markersize=msize)
            # Stork
            ax2.set_title('Stork')
            ax2.set_xlabel(r'$\boldsymbol{\Gamma}$', fontsize=20)
            ax2.set_xlim([-0.05, 2.05])
            for pnum in range(qubits):
                ax2.plot(gammarng, np.array(data['stork'][pnum])/nsamples, lstyle, 
                         c=color[pnum], linewidth=lwidth, label=str(pnum+1)+' patterns', 
                         marker=marker, markersize=msize)
            # Proj
            ax3.set_title('Proj')
            ax3.set_xlim([-0.05, 2.05])
            for pnum in range(qubits):
                ax3.plot(gammarng, np.array(data['proj'][pnum])/nsamples, lstyle, 
                         c=color[pnum], linewidth=lwidth, label=str(pnum+1)+' patterns', 
                         marker=marker, markersize=msize)
            pl.legend(prop={'size':12})
            # pl.show()
            pl.savefig('var_gamma_n'+str(qubits)+'.png')
        else:
            if lrule == 1:
                lrule = 'hebb'
            elif lrule == 2:
                lrule = 'stork'
            elif lrule == 3:
                lrule = 'proj'
            else:
                print("Learning rule number not recognized (0: all,"+
                      "1: hebb, 2: stork, 3: proj).")
                sys.exit(0)
            pl.rcParams['figure.figsize'] = 9,6
            pl.title(r'$\textbf{Fractional success as a function of input bias'+\
                         ' (N = '+str(qubits)+'), '+lrule+'}$', fontsize=18)
            pl.xlabel(r'$\boldsymbol{\Gamma}$', fontsize=20)
            pl.ylabel(r'$\langle f_x \rangle$', fontsize=24)#, rotation='horizontal')
            pl.ylim([0.0, 1.05])
            pl.xlim([-0.05, 2.05])
            for pnum in range(qubits):
            # for pnum in [3]:
                fx = data[lrule][pnum]
                fxavg = np.average(fx, axis=1)
                yerr = np.average(fx**2, axis=1) - fxavg**2
                # print yerr == fx.std(1)
                print gammarng[np.argmax(fxavg)]
                # Do the error margins
                # pl.errorbar(gammarng, fxavg, marker=marker, markersize=msize,
                #             linewidth=lwidth, yerr=yerr, color='grey')
                # Plot with some color
                pl.plot(gammarng, fxavg, c=color[pnum],
                        label=str(pnum+1)+' patterns', marker=marker, markersize=msize,
                        linewidth=lwidth)
            pl.legend(prop={'size':12})
            pl.show()
            pl.savefig('var_gamma_n'+str(qubits)+'_'+lrule+'.png')
