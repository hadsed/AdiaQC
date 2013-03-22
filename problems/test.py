'''

File: test.py
Author: Hadayat Seddiqi
Date: 3.18.13
Description: The test problem.

'''

import scipy as sp

nQubits = 2
#T = 20.0
T = range(1, 110, 9) # Output a sequence of anneal times
dt = 0.1

# Output parameters
output = 1 # Turn on/off all output except final probabilities
eigspecdat = 1 # Output data for eigspec
eigspecplot = 1 # Plot eigspec
eigspecnum = 2**nQubits # Number of eigenvalues
fidelplot = 1 # Plot fidelity
fideldat = 1 # Output fidelity data
fidelnumstates = 1 #2**nQubits # Check fidelity with this number of eigenstates
outputdir = 'data/' # In relation to run.py

errchk = 1 # Error-checking on/off
eps = 0.01 # Numerical error in normalization condition (1 - norm < eps)

# Specify Ising coefficients (Q, a), or alpha, beta directly
ising = 1
Q = sp.zeros((nQubits,nQubits))
#Q = sp.ones((nQubits,nQubits)) - sp.identity(nQubits)
a = sp.ones(nQubits)

# Only if ising = 0; set all to empty SciPy arrays for default coefficients
#alpha = beta = delta = sp.array([])
alpha = sp.ones(2**nQubits)
beta = sp.zeros((2**nQubits,2**nQubits))
delta = sp.array([])
