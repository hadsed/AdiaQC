'''

File: initialize.py
Author: Hadayat Seddiqi
Date: 3.15.13
Description: Initialize stuff (state, Hamiltonian (Ising coefficients), etc.)

'''

import scipy as sp
from scipy import linalg

def InitialState(delta):
    " Set initial state as ground state of transverse field Hamiltonian. " 

    # Get ground state of H
    eigval, eigvec = sp.linalg.eigh(delta)
    sortperm = eigval.argsort()
    eigval = eigval[sortperm]
    eigvec = eigvec[:, sortperm]

    # Ugly casting hack
    Psi = sp.transpose(sp.matrix(sp.transpose(eigvec))[0])

    return Psi

def QUBO2Ising(Q):
    " Convert QUBO problem to Ising model, construct and output \
      the Ising Hamiltonian. "

    # Make sure we've got floats
    Q = Q.astype(sp.float64)

    ##### Convert to Ising ######
    a = Q.diagonal()
    nq = len(a)
    Q = Q - sp.array(a)*sp.identity(nq)
    J = 0.25*Q
    h = 0.5*a + (J + J.T).sum(axis=0)
    g = J.sum() + 0.5*a.sum()

    # Quick hack to make 'h' an actual array
    h = sp.squeeze(sp.asarray(h))

    return [h, J, g]

def IsingHamiltonian(n, h, J, g):
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

        alpha = alpha + temp*h[i]

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

                beta = beta + temp*J[i,j]

    beta = beta + g*sp.identity(2**n)
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

def HamiltonianGen(n, a, b, d):
    " Generate default Hamiltonian coefficients. "

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

        alpha = alpha + temp*a[i]

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

                beta = beta + temp*b[i,j]

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

        delta += temp*d[i]

    return [alpha, beta, delta]

def AlphaBetaCoeffs(n, a, b):
    " Construct the alpha and beta coefficient matrices. "

    Z = sp.matrix([[1,0],[0,-1]])
    I = sp.identity(2)
    alpha = sp.zeros((2**n,2**n))
    beta = sp.zeros((2**n,2**n))
    m1 = []
    m2 = []

    for i in range(0,n):
        for m in range(0,n-1): m1.append(I)
        m1.insert(i, Z)
        temp1 = m1[0]
        m1.pop(0)

        while (len(m1) != 0):
            temp1 = sp.kron(temp1, m1[0])
            m1.pop(0)

        alpha += temp1*a[i]

        for j in range(i+1, n):
            for m in range(0, n-2): m2.append(I)

            m2.insert(i, Z)
            m2.insert(j, Z)
            temp2 = m2[0]
            m2.pop(0)

            while (len(m2) != 0):
                temp2 = sp.kron(temp2, m2[0])
                m2.pop(0)

            beta += (temp2)*b[i,j]

    return [alpha, beta]

def AlphaCoeffs(n, a):
    " Construct the alpha coefficient matrix. "
    Z = sp.matrix([[1,0],[0,-1]])
    I = sp.identity(2)
    alpha = sp.zeros((2**n,2**n))
    matrices = []

    for i in range(0,n):
        for m in range(0,n-1): matrices.append(I)
        matrices.insert(i, Z)
        temp = matrices[0]
        matrices.pop(0)

        while (len(matrices) != 0):
            temp = sp.kron(temp, matrices[0])
            matrices.pop(0)

        alpha += temp*a[i]

    return alpha

def BetaCoeffs(n, b):
    " Calculate beta coefficients. "

    Z = sp.matrix([[1,0],[0,-1]])
    I = sp.identity(2)
    beta = sp.zeros((2**n,2**n))
    matrices = []

    for i in range(0,n):
        for j in range(i+1, n):
            for m in range(0, n-2): matrices.append(I)
            matrices.insert(i, Z)
            matrices.insert(j, Z)
            temp = matrices[0]
            matrices.pop(0)

            while (len(matrices) != 0):
                temp = sp.kron(temp, matrices[0])
                matrices.pop(0)

            beta += (temp)*b[i,j]

    return beta

def DeltaCoeffs(n, d):
    " Calculate delta coefficients. "
    X = sp.matrix([[0,1],[1,0]])
    I = sp.identity(2)
    delta = sp.zeros((2**n,2**n))
    matrices = []

    for i in range(0,n) :
        for m in range(0,n-1) : matrices.append(I)
        matrices.insert(i, X)
        temp = matrices[0]
        matrices.pop(0)

        while (len(matrices) != 0) :
            temp = sp.kron(temp, matrices[0])
            matrices.pop(0)

        delta += temp*d[i]

    return delta
