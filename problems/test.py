'''

File: test.py
Author: Hadayat Seddiqi
Date: 3.18.13
Description: The test problem.

'''

import scipy as sp

# Output parameters
output = 1 # Turn on/off all output except final probabilities
eigspecdat = 1 # Output data for eigspec
eigspecplot = 1 # Plot eigspec
eigspecnum = 4 # Number of eigenvalues
fidelplot = 1 # Plot fidelity
fideldat = 1 # Output fidelity data
fidelnumstates = 1 # Check fidelity with this number of eigenstates
outputdir = 'data/' # In relation to run.py

#T = 20.0
T = range(1, 110, 9)
dt = 0.1

errchk = 1 # Error-checking on/off
eps = 0.01 # Numerical error in normalization condition (1 - norm < eps)

nQubits = 2
#Q = sp.ones((nQubits,nQubits)) - sp.identity(nQubits)
Q = sp.zeros((nQubits,nQubits))
a = sp.ones(nQubits)
