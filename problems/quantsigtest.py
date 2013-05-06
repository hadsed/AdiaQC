'''

File: quantsigtest.py
Author: Hadayat Seddiqi
Date: 4.9.13
Description: This is a test problem derived from Boxio et al. from
             "Experimental signature of programmable quantum annealing,"
             (2012) where we try to show that the quantum annealing
             really is quantum. We should get a 17-fold degenerate
             ground state.

'''

import scipy as sp

nQubits = 8
T = 20.0
#T = sp.arange(1,10,0.5) # Output a sequence of anneal times
dt = 0.01

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
outputdir = 'data/' # In relation to run.py
probout = 1 # Calculate final state probabilities

errchk = 0 # Error-checking on/off
eps = 0.01 # Numerical error in normalization condition (1 - norm < eps)

# Specify Ising coefficients (Q, a), or alpha, beta directly
ising = 0
Q = sp.zeros((nQubits,nQubits))
a = sp.ones(nQubits)

# Only if ising = 0; set all to empty SciPy arrays for default coefficients
alpha = sp.ones(nQubits)
delta = sp.array([])
beta = sp.zeros((nQubits,nQubits))
alpha[0:4] = 1
alpha[4:] = -1
beta[0,1] = beta[1,2] = beta[2,3] = beta[3,0] = beta[0,4] = beta[1,5] = beta[2,6] = beta[3,7] = 1
beta[1,0] = beta[2,1] = beta[3,2] = beta[0,3] = beta[4,0] = beta[5,1] = beta[6,2] = beta[7,3] = 1
