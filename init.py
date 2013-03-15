'''

File: init.py
Author: Hadayat Seddiqi
Date: 3.15.13
Description: Initialize stuff (state, Hamiltonian, etc.)

'''

import scipy as sp
from scipy import linalg

def InitialState(delta):
    " Set initial state as ground state of transverse field Hamiltonian. " 

    # Get ground state of H
    eigval, eigvec = sp.linalg.eig(delta)
    sortperm = eigval.argsort()
    eigval = eigval[sortperm]
    eigvec = eigvec[:, sortperm]

    # Ugly casting hack
    Psi = sp.transpose(sp.matrix(sp.transpose(eigvec))[0])
    return Psi

def QUBO2Ising(n, Q, a):
    " Convert QUBO problem to Ising model, construct and output \
      the Ising Hamiltonian. "

    ##### Convert to Ising ######
    s = 2*a - 1
    J = sp.array(0.25*Q)
    h = 0.5*s + J.sum(axis=1)
    g = J.sum() + 0.5*s.sum()

    ### Construct Hamiltonian ###
    Z = sp.matrix([[1,0],[0,-1]])
    X = sp.matrix([[0,1],[1,0]])
    I = sp.identity(2)
    alpha = sp.zeros((2**n,2**n))
    beta = sp.zeros((2**n,2**n))
    delta = sp.zeros((2**n,2**n))
    matrices = []

    # Calculate alpha
    for i in range(0,n):
        for m in range(0,n-1): matrices.append(I)
        matrices.insert(i, Z)
        temp = matrices[0]
        matrices.pop(0)

        while (len(matrices) != 0):
            temp = sp.kron(temp, matrices[0])
            matrices.pop(0)

        alpha = (alpha + temp)*h[i]

    temp = 0

    # Calculate beta
    for i in range(0,n):
        for j in range(0,n):
            if (i != j):
                for m in range(0,n-2): matrices.append(I)
                matrices.insert(i, Z)
                matrices.insert(j, Z)
                temp = matrices[0]
                matrices.pop(0)

                while (len(matrices) != 0):
                    temp = sp.kron(temp, matrices[0])
                    matrices.pop(0)

                beta = (beta + temp)*J[i,j]

    temp = 0

    # Calculate delta                                                                          
    for i in range(0,n) :
        for m in range(0,n-1) : matrices.append(I)
        matrices.insert(i, X)
        temp = matrices[0]
        matrices.pop(0)

        while (len(matrices) != 0) :
            temp = sp.kron(temp, matrices[0])
            matrices.pop(0)

        delta += temp

    return [alpha, beta, delta]

    # Save matrices to output
    sp.savez('isingcoeffs.npz', alpha=alpha, beta=beta, delta=delta)
