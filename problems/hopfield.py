'''

File: hopfield.py
Author: Hadayat Seddiqi
Date: 3.15.13
Description: Parameters for a Hopfield neural network.

'''
import scipy as sp

# Output parameters
eigspecflag = 1 # Plot eigenspectrum (or not)
outputdir = 'data/'

# Simulation parameters
T = 100 # Total time
dt = 0.1 # Timestep
neurons = 4 # qubits
memories = [ [1,0,0,1] ]
inputstate = [1,0,1,1]

# Generate the QUBO
Q = sp.zeros((neurons,neurons))
a = sp.array(inputstate)
nQubits = neurons

# Construct pattern matrix
for p in range(len(memories)):
    for i in range(neurons):
        for j in range(neurons):
            Q[i,j] += memories[p][i]*memories[p][j]

# No self-connections
sp.fill_diagonal(Q, 0)
