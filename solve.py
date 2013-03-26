'''

File: solve.py
Author: Hadayat Seddiqi
Date: 3.7.13
Description: Simulate quantum annealing with some technique.

'''

import os
import scipy as sp
from scipy import linalg

import output

def DiagHam(alpha, beta, delta, t, T):
    " Get exact eigen states/energies from H. "

    H = 1/T*(t*(alpha + beta) + (T - t)*delta) # Exact Hamiltonian
    #if (~(sp.matrix.getH(H) == H).any()): print ("H is not Hermitian at time "+str(t))
    return sp.linalg.eigh(H)

def CheckNorm(t, nQubits, Psi, Hvecs, eps):
    " Check for numerical error with normalization. "
    norm = 0.0

    for i in range(0, 2**nQubits): 
        norm += abs(sp.vdot(Psi, Hvecs[i]))**2

    if ( (1.0 - norm) > eps ): 
        print (str((1.0 - norm)) + " (norm error) > " + str(eps) + " (eps) @ t = " + str(t))

def ExpPert(nQubits, alpha, beta, delta, Psi, T, dt, errchk, eps, outinfo):
    " Solve using exponential perturbation theory (i.e. Magnus expansion). "

    if (outinfo['eigdat']): eigspec = []
    if (outinfo['overlapdat']): overlap = []

    N = T/dt # steps

    # Loop over time
    for i in range(0, int(sp.floor(N)) + 1):
        t = i*dt
        t0 = (i-1)*dt

        # Approximate Hamiltonian to first term in Magnus expansion (OPTIMIZE)
        Hexp = 1/(2*T)*((t**2 - t0**2)*(alpha + beta) + (2*T*(t - t0) + t0**2 - t**2)*delta)

        A = linalg.expm(-1j*Hexp)
        Psi = A*Psi

        # Get eigendecomposition of true Hamiltonian if necessary
        if (errchk | outinfo['eigdat'] | outinfo['fiddat']): 
            Hvals, Hvecs = DiagHam(alpha, beta, delta, t, T)

            # Sort by eigenvalues
            idx = Hvals.argsort()
            Hvals = Hvals[idx]
            Hvecs = Hvecs[:,idx]
            Hvecs = sp.transpose(Hvecs) # So we can grab them as vectors
        
        # Check for numerical error
        if (errchk):
            CheckNorm(t, nQubits, Psi, Hvecs, eps)

        # Construct eigenspectrum datapoint = [t, eigval 1, ... , eigval n]
        if (outinfo['eigdat']):
            eigspec.append(output.ConstructEigData(t, Hvals, outinfo['eignum']))

        if (outinfo['overlapdat']):
            overlap.append(output.ConstructOverlapData(t, Psi, Hvecs[0]))

    if (outinfo['eigdat']): output.RecordEigSpec(eigspec, outinfo['outdir'], T)
    if (outinfo['eigplot']): output.PlotEigSpec(eigspec, outinfo['outdir'], T)

    if (outinfo['overlapdat']): output.RecordOverlap(overlap, outinfo['outdir'], T)
    if (outinfo['overlapplot']): output.PlotOverlap(overlap, outinfo['outdir'], T)

    return Psi
