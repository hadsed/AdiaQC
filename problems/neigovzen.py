'''

File: neigovzen.py
Author: Hadayat Seddiqi
Date: 7.2.13
Description: Running the experimental results from the paper
             'Quantum pattern recognition with liquid-state 
             nuclear magnetic resonance', Neigovzen et al.


'''

import scipy as sp

nQubits = 2
T = 10.0
#T = sp.arange(2,23,4.0) # Output a sequence of anneal times
dt = 0.01

plots = 0

# Output parameters
output = 1 # Turn on/off all output except final probabilities
eigspecdat = 0 # Output data for eigspec
eigspecplot = plots # Plot eigspec
eigspecnum = 2**nQubits # Number of eigenvalues
fidelplot = 0 # Plot fidelity
fideldat = 0 # Output fidelity data
fidelnumstates = 2**nQubits # Check fidelity with this number of eigenstates
overlapdat = 0 # Output overlap data
overlapplot = plots # Plot overlap
outputdir = 'data/neigovzen/' # In relation to run.py
probout = 1 # Calculate final state probabilities

errchk = 0 # Error-checking on/off
eps = 0.01 # Numerical error in normalization condition (1 - norm < eps)

# Specify Ising coefficients (Q, a), or alpha, beta directly
ising = 1

G = 0.5
eta1 = 0
eta2 = 0
w = -1

# Construct Q matrix
Q = sp.matrix([[G*eta1, w], [0, G*eta2]])
#Q = sp.matrix([[G*eta1, w], [w, G*eta2]])

print Q

