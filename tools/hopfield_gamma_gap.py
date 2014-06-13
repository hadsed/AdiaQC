'''
File: hopfield_gamma_gap.py
Author: Hadayat Seddiqi
Date: 5.10.14
Description: Show the data of a single Hopfield instance
             where the bias term, Gamma, is varied. Generate
             a plot of the energy gap vs. Gamma.
'''

import os, optparse, sys
import json
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
    
def analysis(filename, binary, thresh):
    """
    """
    # Get probability data
    if binary:
        probs = np.load(filename+".dat.npy")
        eigs = np.load('eigenspectrum.dat.npy')
    else:
        probs = np.loadtxt(filename+".dat")
        eigs = np.loadtxt('eigenspectrum.dat')
    # probstr = [ '%1.16E' % e for e in probs.ravel() ]
    bstrs = [ line.rstrip('\n') for line in open('statelabels.txt') ]

    # Fix up eigenvals
    times = eigs[:,0]/eigs[-1,0]
    eigs = eigs[:,1:]

    # Get problem outputs
    props = json.load(open('problem_outputs.dat'))
    lrule = props['learningRule']
    # mems = props['memories']
    instate = props['inputState']
    answer = props['answer']
    gamma = props['bias']

    # # Calculate average Hamming distance
    # avgdist = np.average([ hamdist(m, instate) for m in mems ])
    
    # Calculate success
    success, sprob = fsuccess(bstrs, probs, answer, thresh)

    return success, sprob, lrule, times, eigs, gamma

def make_colormap(seq):
    """
    Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)


# Command line options
if __name__=="__main__":
    parser = optparse.OptionParser("usage: %prog [options] arg1 arg2")
    parser.add_option("-q", "--qubits", dest="qubits", default=2,
                      type="int", 
                      help="Number of qubits.")
    parser.add_option("-p", "--patterns", dest="patterns", default=1,
                      type="int", 
                      help="Number of patterns.")
    parser.add_option("-b", "--binary", dest="binary", default=1,
                      type="int", 
                      help="Whether data file is in binary (1) or not (0).")
    parser.add_option("-f", "--fname", dest="fname", default=None,
                      type="string", 
                      help="Filename for output probabilities.")
    parser.add_option("-t", "--threshold", dest="threshold", default=0.7,
                      type="float", 
                      help="Threshold probability for determining success "+\
                          "(0 < thresh < 1).")
    (options, args) = parser.parse_args()
    qubits = options.qubits
    threshold = options.threshold
    patterns = options.patterns
    binary = options.binary
    probsfname = options.fname
    if not probsfname:
        print "You must specify the -f flag (filename for output probabilities)."
        sys.exit(0)

    # Make some new colormaps
    c = mcolors.ColorConverter().to_rgb
    cg = make_colormap([c('darkgrey'), c('gray'), 0.33,
                        c('gray'), c('lightgray'), 0.66,
                        c('lightgray')])
    c1 = make_colormap([c('darkred'), c('red'), 0.33, 
                        c('red'), c('salmon'), 0.66, 
                        c('salmon')])
    c2 = make_colormap([c('navy'), c('blue'), 0.33, 
                        c('blue'), c('skyblue'), 0.66, 
                        c('skyblue')])
    c3 = make_colormap([c('darkgreen'), c('green'), 0.33, 
                        c('green'), c('lime'), 0.66, 
                        c('lime')])

    # Loop through all data directories
    data = []
    pref = 'n'+str(qubits)+'p'+str(patterns)
    lwidth = 1
    separately = False
    pl.rcParams['xtick.major.pad'] = 8
    pl.rcParams['ytick.major.pad'] = 8
    pl.rcParams['xtick.labelsize'] = 24
    pl.rcParams['ytick.labelsize'] = 24

    fig, (ax1, ax2, ax3) = pl.subplots(1, 3, sharey=True, figsize=(23,11))
    fig.subplots_adjust(wspace=0.12, left=0.05, right=0.95, bottom=0.2)
    # fig.suptitle('Energy gap as a function of Gamma',
    #              fontsize=16, fontweight='bold')

    # Do the Hebb rule
    os.chdir(pref+'hebb')
    for root, dirs, files in os.walk('.'):
        if root == '.':
            continue
        os.chdir(root)
        success, sprob, lrule, times, eigs, gamma = analysis(probsfname, binary, threshold)
        # data.append([ success, dist, kmems, prob, lrule ])
        color = cg(gamma) #cm.gray(gamma)
        lstyle = '-'
        if success:
            color = c1(gamma)
        ax1.text(0.5, 0.93,'Hebb', horizontalalignment='center', 
                 transform=ax1.transAxes, fontsize=24)
        ax1.set_ylabel('Energy Gap', fontweight='bold', fontsize=24)
        ax1.plot(times, eigs[:,1]-eigs[:,0], lstyle, c=color, linewidth=lwidth)
        ax1_cb = fig.add_axes([0.048, 0.08, 0.277, 0.015])
        ax1_cb_b = mpl.colorbar.ColorbarBase(ax1_cb, cmap=c1, orientation='horizontal')
        ax1_cb_b.set_ticks(np.arange(0.0,1.15,0.2))
        os.chdir('../')
    if separately: 
        pl.show()
        pl.clf()
    os.chdir('../')
    # Do the Storkey rule
    os.chdir(pref+'stork')
    for root, dirs, files in os.walk('.'):
        if root == '.':
            continue
        os.chdir(root)
        success, sprob, lrule, times, eigs, gamma = analysis(probsfname, binary, threshold)
        # data.append([ success, dist, kmems, prob, lrule ])
        color = cg(gamma) #'k'
        if success:
            color = c2(gamma)
        ax2.text(0.5, 0.93,'Storkey', horizontalalignment='center', 
                 transform=ax2.transAxes, fontsize=24)
        ax2.set_xlabel('t/T', fontweight='bold', fontsize=24)
        ax2.plot(times, eigs[:,1]-eigs[:,0], c=color, linewidth=lwidth)
        ax2_cb = fig.add_axes([0.361, 0.08, 0.277, 0.015])
        ax2_cb_b = mpl.colorbar.ColorbarBase(ax2_cb, cmap=c2, orientation='horizontal')
        ax2_cb_b.set_label('Gamma (gray indicates failed instance)', 
                           labelpad=8, fontweight='bold', fontsize=24)
        ax2_cb_b.set_ticks(np.arange(0.0,1.15,0.2))
        os.chdir('../')
    os.chdir('../')
    # Do the projection rule
    os.chdir(pref+'proj')
    for root, dirs, files in os.walk('.'):
        if root == '.':
            continue
        os.chdir(root)
        success, sprob, lrule, times, eigs, gamma = analysis(probsfname, binary, threshold)
        # data.append([ success, dist, kmems, prob, lrule ])
        color = cg(gamma)
        if success:
            color = c3(gamma)
        ax3.text(0.5, 0.93,'Projection', horizontalalignment='center', 
                 transform=ax3.transAxes, fontsize=24)
        ax3.plot(times, eigs[:,1]-eigs[:,0], c=color, linewidth=lwidth)
        ax3_cb = fig.add_axes([0.675, 0.08, 0.277, 0.015])
        ax3_cb_b = mpl.colorbar.ColorbarBase(ax3_cb, cmap=c3, orientation='horizontal')
        ax3_cb_b.set_ticks(np.arange(0.0,1.15,0.2))
        os.chdir('../')
    os.chdir('../')
    # pl.show()
    pl.savefig('gamma_gap_'+pref+'.png')
