'''
File: exp1_analysis.py
Author: Hadayat Seddiqi
Date: 4.20.14
Description: Plots and facts on hopfield experiment #1.
'''

import os, optparse, sys
import numpy as np
import json
import pylab as pl

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

def spins2bitstr(vec):
    """ Return a converted bitstring from @vec, a list of spins. """
    return ''.join([ '0' if k == 1 else '1' for k in vec ])

def hamdist(a,b):
    """ Calculate Hamming distance. """
    return np.sum(abs(np.array(a)-np.array(b))/2.0)

def analysis(filename, binary):
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

    return success, avgdist, len(mems), prob, lrule

# Initialize some variables we'll want to look at
csuccess = {'hebb': [0]*qubits,
            'stork': [0]*qubits,
            'proj': [0]*qubits}
cfailure = {'hebb': [0]*qubits,
            'stork': [0]*qubits,
            'proj': [0]*qubits}
data = []

# Loop through all data directories
for root, dirs, files in os.walk('.'):
    if root == '.':
        continue
    # If we are in a dir with no children..
    if dirs == [] and (int(root.split('n')[1].split('p')[0]) == qubits):
        os.chdir(root)
        success, dist, kmems, prob, lrule = analysis(probsfname, binary)
        data.append([ success, dist, kmems, prob, lrule ])
        # pl.scatter(kmems, dist, c=('r' if success is False else 'b'), marker='x')
        # pl.scatter(kmems, prob, c=('r' if success is False else 'b'), marker='x')
        if success:
            csuccess[lrule][kmems-1] += 1
        else:
            cfailure[lrule][kmems-1] += 1
        os.chdir('../../')

# print '\n'
# print '=============='
# print '   Analysis   '
# print '==============\n'
# print 'Successes: ', str(csuccess)
# print 'Failures: ', str(cfailure)
# print '\n'

# Plot success count bar graph
width = 0.3
fontsize = 16
ymax = 0.0
for key in ['hebb', 'stork', 'proj']:
    if np.amax(csuccess[key]) > ymax:
        ymax = np.amax(csuccess[key])
    if np.amax(cfailure[key]) > ymax:
        ymax = np.amax(cfailure[key])

fig = pl.figure(figsize=(8,9))
fig.suptitle('Success vs. # of patterns in '+str(qubits)+
             '-qubit network', fontsize=fontsize, fontweight='bold')
# Loop over rules
for irule, rule in enumerate(['hebb', 'stork', 'proj']):
    # Plot properties
    pl.subplot(3,1,irule)
    pl.title(rule)
    if irule == 2:
        pl.ylabel('Success count', fontweight='bold', fontsize=fontsize)
    if irule == 0:
        pl.xlabel('Number of patterns', fontweight='bold', fontsize=fontsize)
    pl.xlim([0.5,qubits+0.5])
    pl.ylim([0,1.3*ymax])
    # Loop over memory count
    for kmems in range(1,qubits):
        r1 = pl.bar(kmems-width, csuccess[rule][kmems], width, color='b', alpha=0.5)
        r2 = pl.bar(kmems, cfailure[rule][kmems], width, color='r', alpha=0.5)
    pl.grid()
    leg = pl.legend([r1,r2], ["Success", "Failure"], prop={'size':10})
    leg.get_frame().set_alpha(0.5)

# Output to file
if output:
    pl.savefig('success_bargraph_n'+str(qubits)+'.png')

