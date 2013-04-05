'''

File: hopfield.py
Author: Hadayat Seddiqi
Date: 4.5.13
Description: Parameters for a Hopfield neural network.

'''

import scipy as sp

nQubits = 4
T = 100.0
#T = sp.arange(2,23,4.0) # Output a sequence of anneal times
dt = 0.01

# Specify Ising coefficients (Q, a), or alpha, beta directly
ising = 1
Q = sp.empty((nQubits,nQubits))
#Q = sp.ones((nQubits,nQubits)) - sp.identity(nQubits)
a = sp.empty(nQubits)

neurons = nQubits
memories = [ [1,1,1,1], [0,0,0,0] ]
inputstate = [1,0,1,1]

# Construct pattern matrix
for p in range(len(memories)):
    for i in range(neurons):
        for j in range(neurons):
            Q[i,j] += memories[p][i]*memories[p][j]

# No self-connections
sp.fill_diagonal(Q, 0)

a = sp.array(inputstate)

# Output parameters
output = 1 # Turn on/off all output except final probabilities
eigspecdat = 1 # Output data for eigspec
eigspecplot = 1 # Plot eigspec
eigspecnum = 2**nQubits # Number of eigenvalues
fidelplot = 1 # Plot fidelity
fideldat = 1 # Output fidelity data
fidelnumstates = 2**nQubits # Check fidelity with this number of eigenstates
overlapdat = 1 # Output overlap data
overlapplot = 1 # Plot overlap
outputdir = 'data/' # In relation to run.py
probout = 1 # Calculate final state probabilities

errchk = 1 # Error-checking on/off
eps = 0.01 # Numerical error in normalization condition (1 - norm < eps)

