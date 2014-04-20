'''
File: hopfield_analysis.py
Author: Hadayat Seddiqi
Date: 4.18.14
Description: Some useful facts about Hopfield simulations.
'''

import os, optparse, sys
import numpy as np
import json

# Command line options
if __name__=="__main__":
    parser = optparse.OptionParser("usage: %prog [options] arg1 arg2")
    parser.add_option("-p", "--probsfname", dest="probsfname", default=None,
                      type="string", 
                      help="Filename for output probabilities.")
    parser.add_option("-o", "--output", dest="output", default=0,
                      type="int", 
                      help="Output this report to file.")
    parser.add_option("-b", "--binary", dest="binary", default=1,
                      type="int", 
                      help="Binary input or not (default = 1).")
    parser.add_option("-k", "--kmems", dest="kmems", default=2,
                      type="int", 
                      help="How many states to show (multiple of memories).")
    (options, args) = parser.parse_args()
    filename = options.probsfname
    output = options.output
    binary = options.binary
    kmems = options.kmems

# Get probability data
if binary:
    probs = np.load(filename)
else:
    probs = np.loadtxt(filename)
probstr = [ '%1.16E' % e for e in probs.ravel() ]
bstrs = [ line.rstrip('\n') for line in open('statelabels.txt') ]

# Get problem outputs
props = json.load(open('problem_outputs.dat'))
neurons = props['nQubits']
lrule = props['learningRule']
mems = props['memories']
instate = props['inputState']
annealtime = props['annealTime']

def spins2bitstr(vec):
    """
    Return a converted bitstring from @vec, a list of spins.
    """
    return ''.join([ '0' if k == 1 else '1' for k in vec ])

# Convert list of spins to bitstring
distances = []
for im, m in enumerate(mems):
    distances.append(np.sum(abs(np.array(instate)-np.array(m))/2))
    mems[im] = spins2bitstr(m)
instate = spins2bitstr(instate)

# Build the output string
outstr = ('\n====================================\n' +
          '  Hopfield network analysis report  \n' +
          '====================================\n')
outstr += '\nNeurons:'+ ' '*10 + str(neurons)
outstr += '\nMemories:'+ ' '*9 + str(len(mems))
outstr += '\nRatio (alpha):'+ ' '*4 + str(len(mems)/float(neurons))
outstr += '\nLearning rule:' + ' '*4 + str(lrule)
outstr += '\nAnnealing time:' + ' '*3 + str(annealtime)
outstr += '\n\nMean distance:' + ' '*4 + str(np.mean(distances))
outstr += '\nMedian distance:' + ' '*2 + str(np.median(distances))

# Memories
outstr += '\n\nPatterns : associated probabilities : rank : Hamming distance\n'
sortidx = np.argsort(probs)[::-1]
bstrs_sorted = np.array(bstrs)[sortidx]
for im, m in enumerate(mems):
    outstr += m + '\t    ' + probstr[bstrs.index(m)] + '     ' + \
        str(np.where(bstrs_sorted == m)[0][0]) + '\t' + str(distances[im]) + '\n'

# Input state
outstr += '\n\nInput state:\n'
outstr += instate + '\t    ' + probstr[bstrs.index(instate)] + '     ' +\
    str(np.where(bstrs_sorted == instate)[0][0]) + '\n'

# Most likely states
k = min(kmems*len(mems), 2**neurons)
outstr += '\n\nTop ' + str(k) + ' most likely states (sorted):\n'
sortidx = np.argsort(probs)[::-1][:k]
ptop = np.array(probstr)[sortidx]
btop = np.array(bstrs)[sortidx]
for p,b in zip(ptop, btop):
    outstr += b + '\t' + str(p) + '\n'

# Print it all out
outstr += '\n'
print outstr

# Output
if output:
    with open('hopfield_report.txt', 'w') as f:
        f.write(outstr)
