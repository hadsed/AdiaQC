'''

File: initialize.py
Author: Hadayat Seddiqi
Date: 3.15.13
Description: Initialize stuff (state, Hamiltonian (Ising coefficients), etc.)

'''

import scipy as sp
from scipy import linalg
from scipy import sparse as sps
from scipy.sparse import linalg as spla

def InitialState(delta):
    " Set initial state as ground state of transverse field Hamiltonian. " 

    # Get ground state of H
    eigval, eigvec = spla.eigsh(delta)
    # eigval, eigvec = sp.linalg.eigh(delta.todense())
    sortperm = eigval.argsort()
    eigval = eigval[sortperm]
    eigvec = eigvec[:, sortperm]
    # Ugly casting hack
    Psi = sp.matrix(eigvec[:,0]).T
    return Psi

def QUBO2Ising(Q):
    " Convert QUBO problem to Ising model, construct and output \
      the Ising Hamiltonian. "

    # Make sure we've got floats
    Q = Q.astype(sp.float64)
    # Do the conversion
    a = Q.diagonal()
    nq = len(a)
    Q = Q - sp.array(a)*sp.identity(nq)
    J = 0.25*Q
    h = 0.5*a + (J + J.T).sum(axis=0)
    g = J.sum() + 0.5*a.sum()
    # Quick hack to make 'h' an actual array
    h = sp.squeeze(sp.asarray(h))
    return [h, J, g]

def _formZi(nq, i, z):
    """
    Helper function to form 2**(@nq) matrix Zi:

    Zi = I \kron I ... \kron @z \kron I ..

    i.e. where @z is in the @i-th position.
    """
    if i == 0:
        return sps.kron(z, sps.eye(2**(nq-1)))
    elif i == (nq-1):
        return sps.kron(sps.eye(2**(nq-1)), z)
    else:
        k1 = i
        k2 = nq - i - 1
        return sps.kron(sps.eye(2**k1),
                        sps.kron(z, sps.eye(2**k2)))

def _formZiZj(nq, i, j, z):
    """
    Helper function to form 2**(@nq) matrix Zi:

    ZiZj = I \kron I ... \kron @z \kron I .. \kron @z ..

    i.e. where two @z's are in @i-th and @j-th positions.
    """
    imin = min(i,j)
    jmax = max(i,j)

    # Special case of diagonal
    if i == j:
        if i == 0:
            return sps.kron(z*z, sps.eye(2**(nq-1)))
        elif i == (nq-1):
            return sps.kron(sps.eye(2**(nq-1)), z*z)
        else:
            # k1 = i - 1
            # k2 = nq - i
            k1 = i
            k2 = nq - i - 1
            return sps.kron(sps.eye(2**k1), 
                            sps.kron(z*z, sps.eye(2**k2)))
    # Construct first part
    if imin == 0:
        t1 = z
    else:
        t1 = sps.kron(sps.eye(2**(imin)), z)
    # Construct third part
    if jmax == (nq-1):
        t3 = z
    else:
        t3 = sps.kron(z, sps.eye(2**(nq-jmax-1)))
    # Put them all together
    if jmax == (imin+1):
        return sps.kron(t1,t3)
    else:
        k1 = (jmax-1) - imin
        t2 = sps.eye(2**k1)
        return sps.kron(t1, sps.kron(t2,t3))

def _formX(nq, x, i):
    """
    Helper function to form the X part of the Hamiltonian.
    """
    if nq == 1:
        return x
    elif nq == 2:
        return sps.kron(x,i) + sps.kron(i,x)
    elif nq == 3:
        return sps.kron(x, sps.kron(i,i)) + \
            sps.kron(i, sps.kron(x,i)) + \
            sps.kron(i, sps.kron(i,x))
    elif nq == 4:
        return sps.kron(x, sps.kron(i, sps.kron(i,i))) + \
            sps.kron(i, sps.kron(x, sps.kron(i,i))) + \
            sps.kron(i, sps.kron(i, sps.kron(x,i))) + \
            sps.kron(i, sps.kron(i, sps.kron(i,x)))
    else:
        h = _formX(nq-1, x, i)
        return sps.kron(x, sps.eye(h.shape[0])) + sps.kron(i,h)

def IsingHamiltonian(n, h, J, g):
    # Construct Hamiltonian
    Z = sps.csr_matrix([[1,0],[0,-1]])
    X = sps.csr_matrix([[0,1],[1,0]])
    I = sps.eye(2)
    alpha = sps.csr_matrix((2**n, 2**n))
    beta = sps.csr_matrix((2**n, 2**n)) 
    delta = sps.csr_matrix((2**n, 2**n))
    # Form alpha and beta operators
    for q1 in xrange(n):
        alpha = alpha + h[q1]*_formZi(n,q1,Z)
        for q2 in xrange(n):
            beta = beta + J[q1,q2]*_formZiZj(n,q1,q2,Z)
    # Add energy shift
    beta.setdiag(beta.diagonal() + g)
    # Form delta operator
    delta = _formX(n, X, I)
    return alpha, beta, delta

def HamiltonianGen(n, a, b, d):
    # Construct Hamiltonian
    Z = sps.csr_matrix([[1,0],[0,-1]])
    X = sps.csr_matrix([[0,1],[1,0]])
    I = sps.eye(2)
    alpha = sps.csr_matrix((2**n, 2**n))
    beta = sps.csr_matrix((2**n, 2**n)) 
    delta = sps.csr_matrix((2**n, 2**n))
    # Form alpha and beta operators
    for q1 in xrange(n):
        alpha = alpha + a[q1]*_formZi(n,q1,Z)
        for q2 in xrange(q1+1,n):
            beta = beta + b[q1,q2]*_formZiZj(n,q1,q2,Z)
    # Form delta operator
    delta = _formX(n, X, I)
    return alpha, beta, delta

def AlphaBetaCoeffs(n, a, b):
    # Construct alpha and beta operators
    Z = sps.csr_matrix([[1,0],[0,-1]])
    I = sps.eye(2)
    alpha = sps.csr_matrix((2**n, 2**n))
    beta = sps.csr_matrix((2**n, 2**n)) 
    # Form alpha and beta operators
    for q1 in xrange(n):
        alpha = alpha + a[q1]*_formZi(n,q1,Z)
        for q2 in xrange(n):
            beta = beta + b[q1,q2]*_formZiZj(n,q1,q2,Z)
    return alpha, beta

def AlphaCoeffs(n, a):
    # Construct alpha operator
    Z = sps.csr_matrix([[1,0],[0,-1]])
    I = sps.eye(2)
    alpha = sps.eye(2**n)
    # Form alpha and beta operators
    for q1 in xrange(n):
        alpha = alpha + a[q1]*_formZi(n,q1,Z)
    return alpha

def BetaCoeffs(n, b):
    # Construct beta operator
    Z = sps.csr_matrix([[1,0],[0,-1]])
    I = sps.eye(2)
    beta = sps.csr_matrix((2**n, 2**n)) 
    # Form alpha and beta operators
    for q1 in xrange(n):
        for q2 in xrange(n):
            beta = beta + b[q1,q2]*_formZiZj(n,q1,q2,Z)
    return beta

def DeltaCoeffs(n, d):
    # Construct X/delta operator
    X = sps.csr_matrix([[0,1],[1,0]])
    I = sps.eye(2)
    delta = sps.csr_matrix((2**n, 2**n))
    # Form delta operator
    delta = _formX(n, X, I)
    return delta


#############################################################
##############     Old stuff.     ###########################
#############################################################


def IsingHamiltonian_old(n, h, J, g):
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
        for m in range(0,n-1):
            matrices.append(I)
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
                for m in range(0,n-2):
                    matrices.append(I)
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
        for m in range(0,n-1):
            matrices.append(I)
        matrices.insert(i, X)
        temp = matrices[0]
        matrices.pop(0)
        while (len(matrices) != 0):
            temp = sp.kron(temp, matrices[0])
            matrices.pop(0)
        delta += temp
    return [alpha, beta, delta]

def HamiltonianGen_old(n, a, b, d):
    """ Generate default Hamiltonian coefficients. """
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
        for j in range(i+1,n):
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
    for i in range(0,n):
        for m in range(0,n-1):
            matrices.append(I)
        matrices.insert(i, X)
        temp = matrices[0]
        matrices.pop(0)
        while (len(matrices) != 0):
            temp = sp.kron(temp, matrices[0])
            matrices.pop(0)
        delta += temp*d[i]
    return [alpha, beta, delta]

def AlphaCoeffs_old(n, a):
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

def BetaCoeffs_old(n, b):
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

def DeltaCoeffs_old(n, d):
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

def AlphaBetaCoeffs_old(n, a, b):
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
