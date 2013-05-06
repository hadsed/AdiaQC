'''

File: threevarqubo.py
Author: Hadayat Seddiqi
Date: 5.1.13
Description: Test problem with 3 variables.

'''

import scipy as sp

nQubits = 3
T = 30.0
#T = sp.arange(2,103,10) # Output a sequence of anneal times
dt = 0.1

# Output parameters
output = 1 # Turn on/off all output except final probabilities
eigspecdat = 0 # Output data for eigspec
eigspecplot = 1 # Plot eigspec
eigspecnum = 2**nQubits # Number of eigenvalues
fidelplot = 1 # Plot fidelity
fideldat = 0 # Output fidelity data
fidelnumstates = 2**nQubits # Check fidelity with this number of eigenstates
overlapdat = 0 # Output overlap data
overlapplot = 1 # Plot overlap
outputdir = 'threevarqubodata/' # In relation to run.py
probout = 1 # Calculate final state probabilities

errchk = 0 # Error-checking on/off
eps = 0.01 # Numerical error in normalization condition (1 - norm < eps)

# Specify QUBO matrix Q or alpha, beta directly
ising = 1
Q = sp.triu(sp.ones((nQubits,nQubits)))

Q[0,0] = -6
Q[1,1] = 3
Q[2,2] = 7

Q[0,1] = 2
Q[0,2] = 5
Q[1,2] = -12
