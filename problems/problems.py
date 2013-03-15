'''

File: problems.py
Author: Hadayat Seddiqi
Date: 3.15.13
Description: Return a QUBO for some problem.

'''

import scipy as sp
from scipy import linalg

def HopfieldNetwork(n, memories, inpattern):
    " Construct a QUBO corresponding to a Hopfield neural network problem. \
      Returns n (neurons/qubits), Q (memories), and a (input vector). "

    Q = sp.zeros((n,n))
    a = sp.array(inpattern)
    
    for p in range(len(memories)):
        for i in range(n):
            for j in range(n):
                if (i != j): Q[i,j] += memories[p][i]*memories[p][j]
                else: Q[i,j] = 0

    return [n, Q, a]

def CorrelationClustering(n):
    " Construct a QUBO for a correlation clustering problem. "

    return n
