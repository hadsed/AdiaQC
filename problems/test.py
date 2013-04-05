'''

File: test.py
Author: Hadayat Seddiqi
Date: 3.18.13
Description: The test problem.

'''

import scipy as sp

nQubits = 2
T = 10.0
#T = sp.arange(2,23,4.0) # Output a sequence of anneal times
dt = 0.01

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

# Specify Ising coefficients (Q, a), or alpha, beta directly
ising = 0
Q = sp.zeros((nQubits,nQubits))
#Q = sp.ones((nQubits,nQubits)) - sp.identity(nQubits)
a = sp.ones(nQubits)

# Only if ising = 0; set all to empty SciPy arrays for default coefficients
alpha = sp.ones(nQubits)
#delta = sp.zeros(nQubits)
#beta = sp.array([])
delta = sp.array([])
#alpha = beta = delta = sp.array([])
beta = sp.zeros((nQubits,nQubits))
beta[0,1] = 1
beta[1,0] = beta[0,1]
