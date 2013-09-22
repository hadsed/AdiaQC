'''

File: tetris.py
Author: Hadayat Seddiqi
Date: 9.14.13
Description: Recognize tetris pieces with a quantum Hopfield network but
             using the Storkey learning rule.

'''

import scipy as sp
import itertools

nQubits = 4
#T = 10.0
T = sp.arange(0.1,20,0.1) # Output a sequence of anneal times
dt = 0.01*T

# Output parameters
output = 1 # Turn on/off all output except final probabilities
eigspecdat = 0 # Output data for eigspec
eigspecplot = 0 # Plot eigspec
eigspecnum = 2**nQubits # Number of eigenvalues
fidelplot = 0 # Plot fidelity
fideldat = 0 # Output fidelity data
fidelnumstates = 2**nQubits # Check fidelity with this number of eigenstates
overlapdat = 0 # Output overlap data
overlapplot = 0 # Plot overlap
outputdir = 'data/hopfield/tetris/n4p1stork' # In relation to run.py
probout = 0 # Calculate final state probabilities
mingap = 0 # Output minimum spectral gap
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
# Case number determines how many patterns we input to the memory matrix
caseNum = 1
for p in range(caseNum):
    bitstring = list(patternSet[p])
    spins = [ 1 if k == '1' else -1 for k in bitstring ]
    memories.append(spins)
    
# Make the input the last memory recorded
inputstate = memories[-1]

# This is gamma, the appropriate weighting on the input vector
isingSigns['hz'] *= (1 - (len(inputstate) - inputstate.count(0))/(2*neurons))

alpha = sp.array(inputstate)
beta = sp.zeros((neurons,neurons))
delta = sp.array([])

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
print beta
