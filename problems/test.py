'''

File: test.py
Author: Hadayat Seddiqi
Date: 3.18.13
Description: The test problem.

'''

import scipy as sp

# Output parameters
eigspecdat = 1 # Data for eigenspectrum (or not)
eigspecplot = 1 # Output a plot
eigspecnum = 2

outputdir = 'data/'

T = 10.0
dt = 0.1

nQubits = 1
#Q = sp.ones((nQubits,nQubits)) - sp.identity(nQubits)
Q = sp.zeros((nQubits,nQubits))
a = sp.ones(nQubits)

print (Q, a)
