'''
File: hopfield_gamma_analysis.py
Author: Hadayat Seddiqi
Date: 5.10.14
Description: Show the data of a single Hopfield instance
             where the bias term, Gamma, is varied. Analyze
             the eigenspectrum and probabilities.
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
    success, sprob = success2(bstrs, probs, answer, instate)

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
    parser.add_option("-n", "--nqubits", dest="nqubits", default=2,
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
    (options, args) = parser.parse_args()
    nqubits = options.nqubits
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
    pref = 'n'+str(nqubits)+'p'+str(patterns)
    lwidth = 1
    separately = False
    fig, (ax1, ax2, ax3) = pl.subplots(1, 3, sharey=True, figsize=(20,9))
    fig.subplots_adjust(wspace=0.05, left=0.05, right=0.95, bottom=0.15)
    fig.suptitle('Energy gap as a function of Gamma',
                 fontsize=16, fontweight='bold')

    # Do the Hebb rule
    os.chdir(pref+'hebb')
    for root, dirs, files in os.walk('.'):
        if root == '.':
            continue
        os.chdir(root)
        success, sprob, lrule, times, eigs, gamma = analysis(probsfname, binary)
        # data.append([ success, dist, kmems, prob, lrule ])
        color = cg(gamma) #cm.gray(gamma)
        lstyle = '-'
        # pl.subplot(1,3,1)
        if success:
            color = c1(gamma)
        ax1.set_title('Hebb')
        ax1.set_ylabel('Energy Gap', fontweight='bold')
        ax1.plot(times, eigs[:,1]-eigs[:,0], lstyle, c=color, linewidth=lwidth)
        ax1_cb = fig.add_axes([0.05, 0.08, 0.29, 0.01])
        ax1_cb_b = mpl.colorbar.ColorbarBase(ax1_cb, cmap=c1, orientation='horizontal')
        # print success, sprob, gamma
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
        success, sprob, lrule, times, eigs, gamma = analysis(probsfname, binary)
        # data.append([ success, dist, kmems, prob, lrule ])
        color = cg(gamma) #'k'
        # pl.subplot(1,3,2)
        if success:
            color = c2(gamma)
        ax2.set_title('Storkey')
        ax2.set_xlabel('Annealing Time', fontweight='bold')
        ax2.plot(times, eigs[:,1]-eigs[:,0], c=color, linewidth=lwidth)
        ax2_cb = fig.add_axes([0.355, 0.08, 0.29, 0.01])
        ax2_cb_b = mpl.colorbar.ColorbarBase(ax2_cb, cmap=c2, orientation='horizontal')
        ax2_cb_b.set_label('Gamma (gray indicates failed instance)', 
                           labelpad=8, fontweight='bold')
        # print success, sprob, gamma
        os.chdir('../')
    if separately: 
        pl.show()
        pl.clf()
    os.chdir('../')
    # Do the projection rule
    os.chdir(pref+'proj')
    for root, dirs, files in os.walk('.'):
        if root == '.':
            continue
        os.chdir(root)
        success, sprob, lrule, times, eigs, gamma = analysis(probsfname, binary)
        # data.append([ success, dist, kmems, prob, lrule ])
        color = cg(gamma) #'k'
        # pl.subplot(1,3,3)
        if success:
            color = c3(gamma)
        ax3.set_title('Projection')
        ax3.plot(times, eigs[:,1]-eigs[:,0], c=color, linewidth=lwidth)
        ax3_cb = fig.add_axes([0.66, 0.08, 0.29, 0.01])
        ax3_cb_b = mpl.colorbar.ColorbarBase(ax3_cb, cmap=c3, orientation='horizontal')
        # print success, sprob, gamma
        os.chdir('../')
    os.chdir('../')
    pl.savefig(pref+'_gamma_analysis.png')

# # Plot success count bar graph
# width = 0.3
# fontsize = 16
# ymax = 0.0
# for key in ['hebb', 'stork', 'proj']:
#     if np.amax(csuccess[key]) > ymax:
#         ymax = np.amax(csuccess[key])
#     if np.amax(cfailure[key]) > ymax:
#         ymax = np.amax(cfailure[key])

# fig = pl.figure(figsize=(8,9))
# fig.suptitle('Success vs. # of patterns in '+str(qubits)+
#              '-qubit network', fontsize=fontsize, fontweight='bold')
# # Loop over rules
# for irule, rule in enumerate(['hebb', 'stork', 'proj']):
#     # Plot properties
#     pl.subplot(3,1,irule)
#     pl.title(rule)
#     if irule == 2:
#         pl.ylabel('Success count', fontweight='bold', fontsize=fontsize)
#     if irule == 0:
#         pl.xlabel('Number of patterns', fontweight='bold', fontsize=fontsize)
#     pl.xlim([0.5,qubits+0.5])
#     pl.ylim([0,1.3*ymax])
#     # Loop over memory count
#     for kmems in range(1,qubits):
#         r1 = pl.bar(kmems-width, csuccess[rule][kmems], width, color='b', alpha=0.5)
#         r2 = pl.bar(kmems, cfailure[rule][kmems], width, color='r', alpha=0.5)
#     pl.grid()
#     leg = pl.legend([r1,r2], ["Success", "Failure"], prop={'size':10})
#     leg.get_frame().set_alpha(0.5)

# # Output to file
# if output:
#     pl.savefig('success_bargraph_n'+str(qubits)+'.png')
