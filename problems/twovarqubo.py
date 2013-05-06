'''

File: twovarqubo.py
Author: Hadayat Seddiqi
Date: 5.1.13
Description: Test problem with 2 variables.

'''

import scipy as sp

nQubits = 2
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
outputdir = 'data/twovarqubodata/' # In relation to run.py
probout = 1 # Calculate final state probabilities

errchk = 0 # Error-checking on/off
eps = 0.01 # Numerical error in normalization condition (1 - norm < eps)

# Specify QUBO matrix Q or alpha, beta directly
ising = 1
Q = 0.5*sp.triu(sp.ones((nQubits,nQubits))) + -3*sp.identity(nQubits)
