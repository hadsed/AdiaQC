'''

File: sixvarqubo.py
Author: Hadayat Seddiqi
Date: 5.1.13
Description: A six-variable QUBO test problem. Should
             have minimum vectors of the following:
             (1,0,0,1,0,1)   and
             (1,1,0,1,0,1)

'''

import scipy as sp

nQubits = 6
T = 30.0
#T = sp.arange(2,103,10) # Output a sequence of anneal times
dt = 0.1

# Output parameters
output = 1 # Turn on/off all output except final probabilities
eigspecdat = 0 # Output data for eigspec
eigspecplot = 1 # Plot eigspec
eigspecnum = 2**nQubits # Number of eigenvalues
fidelplot = 0 # Plot fidelity
fideldat = 0 # Output fidelity data
fidelnumstates = 2**nQubits # Check fidelity with this number of eigenstates
overlapdat = 0 # Output overlap data
overlapplot = 1 # Plot overlap
outputdir = 'data/sixvarqubodata/' # In relation to run.py
probout = 1 # Calculate final state probabilities
mingap = 1 # Output the minimum spectral gap

errchk = 0 # Error-checking on/off
eps = 0.01 # Numerical error in normalization condition (1 - norm < eps)

# Specify Ising coefficients (Q, a), or alpha, beta directly
isingConvert = 1
isingSigns = {'hx': 1, 'hz': 1, 'hzz': -1}

Q = sp.zeros((nQubits,nQubits))
a = sp.ones(nQubits)

a[0] = 2
a[1] = 1
a[2] = -2
a[3] = -1
a[4] = 1
a[5] = -1

Q[0,1] = Q[1,0] = -1
Q[0,2] = Q[2,0] = 2
Q[0,3] = Q[3,0] = -2
Q[0,4] = Q[4,0] = 2
Q[0,5] = Q[5,0] = -1

Q[1,2] = Q[2,1] = 1
Q[1,3] = Q[3,1] = -1
Q[1,4] = Q[4,1] = -1
Q[1,5] = Q[5,1] = 1

Q[2,3] = Q[3,2] = 2
Q[2,4] = Q[4,2] = -2
Q[2,5] = Q[5,2] = 1

Q[3,4] = Q[4,3] = 2
Q[3,5] = Q[5,3] = -1

Q[4,5] = Q[5,4] = 2

Q = sp.triu(Q) + a*sp.identity(6)

print Q
