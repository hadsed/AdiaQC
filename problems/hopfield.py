'''

File: hopfield.py
Author: Hadayat Seddiqi
Date: 4.5.13
Description: Parameters for a Hopfield neural network.

'''

import scipy as sp

nQubits = 4
T = 10.0
#T = sp.arange(2,23,4.0) # Output a sequence of anneal times
dt = 0.01

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
outputdir = 'data/hopfield/' # In relation to run.py
probout = 1 # Calculate final state probabilities

errchk = 0 # Error-checking on/off
eps = 0.01 # Numerical error in normalization condition (1 - norm < eps)

# Specify Ising coefficients (Q, a), or alpha, beta directly
ising = 1
Q = sp.zeros((nQubits,nQubits))
#Q = sp.ones((nQubits,nQubits)) - sp.identity(nQubits)

neurons = nQubits
G = -0.25 # Scale factor for input memory 'a'
memories = [sp.array([1,1,-1,-1])*-1]
#memories = [ [-1,-1,-1,-1,-1,-1,-1], [1,1,1,1,1,1,1] ]
#memories = [ [1,1,1,-1,-1,-1,-1] ]
inputstate = sp.array([-1,-1,-1,-1])*1
#inputstate = [-1,-1,-1,-1,-1,-1,-1]
#inputstate = [0,0,0,0,0,0,0]

# Construct pattern matrix
for p in range(len(memories)):
    for i in range(neurons):
        for j in range(neurons):
            # Encode input state on diagonal
            if (i == j) : Q[i,j] = G*inputstate[i]
            # Encode memories using Hebb's rule
            else: Q[i,j] += memories[p][i]*memories[p][j]


# Q must be triangular, normalize the weightings too
Q = sp.triu(Q, 1)/neurons + Q.diagonal()*sp.identity(neurons)
#Q = sp.triu(Q)
print(Q)

