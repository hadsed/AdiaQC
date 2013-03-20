'''

File: test.py
Author: Hadayat Seddiqi
Date: 3.18.13
Description: The test problem.

'''

import scipy as sp

# Output parameters
eigspecdat = 1 # Output data
eigspecplot = 1 # Output a plot
eigspecnum = 2 # Number of eigenvalues

outputdir = 'data/' # In relation to run.py

T = 10.0
dt = 0.1

errchk = 1 # Error-checking on/off
eps = 0.01 # Numerical error in normalization condition (1 - norm < eps)

nQubits = 2
#Q = sp.ones((nQubits,nQubits)) - sp.identity(nQubits)
Q = sp.zeros((nQubits,nQubits))
a = sp.ones(nQubits)

print (Q, a)
