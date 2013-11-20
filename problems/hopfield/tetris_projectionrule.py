'''

File: tetris.py
Author: Hadayat Seddiqi
Date: 9.8.13
Description: Recognize tetris pieces with a quantum Hopfield network but
             using the projection (pseudo-inverse) learning rule.

'''

import scipy as sp
import itertools
import random

import params

caseNumber = 16

nQubits = 4
T = params.annealTime
#T = sp.arange(0.1,20,0.1) # Output a sequence of anneal times
dt = 0.01*T

# Output parameters
output = 1 # Turn on/off all output except final probabilities
eigspecdat = 1 # Output data for eigspec
eigspecplot = 0 # Plot eigspec
eigspecnum = params.numEigvals # Number of eigenvalues
fidelplot = 0 # Plot fidelity
fideldat = 0 # Output fidelity data
fidelnumstates = 2**nQubits # Check fidelity with this number of eigenstates
overlapdat = 0 # Output overlap data
overlapplot = 0 # Plot overlap
outputdir = 'data/hopfield/tetris/n4p'+params.simCase+'proj' # In relation to run.py
probout = 0 # Calculate final state probabilities
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
for k in range(params.numMemories):
    bitstring = random.choice(patternSet)
    spins = [ 1 if k == '1' else -1 for k in bitstring ]
    memories.append(spins)

# Case number determines how many patterns we input to the memory matrix
#for p in range(caseNumber):
#    bitstring = list(patternSet[p])
#    spins = [ 1 if k == '1' else -1 for k in bitstring ]
#    memories.append(spins)
    
# Make the input the last memory recorded
#inputstate = memories[-1]
inputstate = params.inputState

# This is gamma, the appropriate weighting on the input vector
isingSigns['hz'] *= (1 - (len(inputstate) - inputstate.count(0))/(2*neurons))

alpha = sp.array(inputstate)
beta = sp.zeros((neurons,neurons))
delta = sp.array([])

# Construct the memory matrix according to the Moore-Penrose pseudoinverse rule
memMat = sp.matrix(memories).T
beta = sp.triu(memMat * sp.linalg.pinv(memMat))
print memMat
print beta

# Calculate Hamming distance between input
# state and each memory
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
    'annealTime': T
           }
