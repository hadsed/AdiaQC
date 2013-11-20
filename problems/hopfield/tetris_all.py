'''

File: tetris_hebbrule.py
Author: Hadayat Seddiqi
Date: 4.5.13
Description: Recognize tetris pieces with a quantum Hopfield network
             using the Hebbian learning rule.

'''

import os
import scipy as sp
import itertools
import random
import collections

import params

learningRule = params.learningRules[params.rule]

nQubits = params.numQubits
T = params.annealTime
dt = 0.01*T

# Output parameters
output = 1 # Turn on/off all output except final probabilities
eigspecdat = 1 # Output data for eigspec
eigspecplot = 0 # Plot eigspec
eigspecnum = params.numEigvals # Number of eigenvalues
fidelplot = 1 # Plot fidelity
fideldat = 1 # Output fidelity data
fidelnumstates = 2**nQubits # Check fidelity with this number of eigenstates
overlapdat = 0 # Output overlap data
overlapplot = 0 # Plot overlap

# Output directory stuff
probdir = 'data/hopfield/tetris/n'+str(nQubits)+'p'+params.simCase+learningRule
if isinstance(T, collections.Iterable):
    probdir += 'MultiT'
if os.path.isdir(probdir):
    outlist = sorted([ int(name) for name in os.listdir(probdir) if name.isdigit() ])
else:
    outlist = []
outnum = outlist[-1] + 1 if outlist else 0
outputdir = probdir + '/' + str(outnum) + '/'

probout = 0 # Print final state probabilities
mingap = 1 # Output minimum spectral gap
outdat = 1 # Output probabilities and mingap

errchk = 0 # Error-checking on/off
eps = 0.01 # Numerical error in normalization condition (1 - norm < eps)

# Specify a QUBO (convert to Ising = True), or alpha, beta directly 
# (convert = False), and also specify the signs on the Ising Hamiltonian 
# terms (you can specify coefficients too for some problems if needed)
isingConvert = 0
isingSigns = {'hx': -1, 'hz': -1, 'hzz': -1}

neurons = nQubits
memories = []

# All possible 2^n patterns
patternSet = ["".join(seq) for seq in itertools.product("01", repeat=neurons)]

# Pick P number of random memories
for k in range(params.numMemories):
    bitstring = random.choice(patternSet)
    spins = [ 1 if k == '1' else -1 for k in bitstring ]
    # Make sure we have a unique set
    while spins in memories:
        bitstring = random.choice(patternSet)
        spins = [ 1 if k == '1' else -1 for k in bitstring ]
    memories.append(spins)

# Make the input the last memory recorded
inputstate = params.inputState

# Add in the input state
if params.includeInput and inputstate not in memories:
    memories[0] = inputstate
    
# This is gamma, the appropriate weighting on the input vector
isingSigns['hz'] *= 1 - (len(inputstate) - inputstate.count(0))/(2*neurons)

alpha = sp.array(inputstate)
beta = sp.zeros((neurons,neurons))
delta = sp.array([])

if learningRule == 'hebb':
    # Construct pattern matrix according to Hebb's rule
    for i in range(neurons):
        for j in range(neurons):
            for p in range(len(memories)):
                beta[i,j] += ( memories[p][i]*memories[p][j] -
                               len(memories)*(i == j) )

    beta = sp.triu(beta)/float(neurons)
elif learningRule == 'stork':
    # Construct the memory matrix according to the Storkey learning rule
    memMat = sp.zeros((neurons,neurons))
    for m, mem in enumerate(memories):
        for i in range(neurons):
            for j in range(neurons):
                hij = sp.sum([ memMat[i,k]*mem[k] for k in range(neurons) ])
                hji = sp.sum([ memMat[j,k]*mem[k] for k in range(neurons) ])
                # Don't forget to make the normalization a float!
                memMat[i,j] += 1./neurons*(mem[i]*mem[j] - mem[i]*hji - hij*mem[j])

    beta = sp.triu(memMat)
elif learningRule == 'proj':
    # Construct the memory matrix according to the Moore-Penrose pseudoinverse rule
    memMat = sp.matrix(memories).T
    beta = sp.triu(memMat * sp.linalg.pinv(memMat))

# Calculate Hamming distance between input state and each memory
hammingDistance = []
for mem in memories:
    dist = sp.sum(abs(sp.array(inputstate)-sp.array(mem))/2)
    hammingDistance.append(dist)

hamMean = sp.average(hammingDistance)
hamMed = sp.median(hammingDistance)

# Some outputs
outputs = {
    'input': inputstate,
    'memories': memories,
    'hammingDistance': {'dist': hammingDistance,
                        'mean': hamMean,
                        'median': hamMed },
    'annealTime': list(T)
           }
