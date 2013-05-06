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

neurons = nQubits
G = 10 # Scale factor for input memory 'a'
memories = [ [1,1,1,1] ]
inputstate = [1,1,1,1]

# Construct pattern matrix
for p in range(len(memories)):
    for i in range(neurons):
        for j in range(neurons):
            Q[i,j] += memories[p][i]*memories[p][j]

# No self-connections, encode input state
sp.fill_diagonal(Q, G*inputstate)


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
outputdir = 'data/hopfield/' # In relation to run.py
probout = 1 # Calculate final state probabilities

errchk = 0 # Error-checking on/off
eps = 0.01 # Numerical error in normalization condition (1 - norm < eps)

