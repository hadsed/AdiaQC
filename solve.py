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

def DiagHam(hz, hzz, hx, t, T, n):
    """ 
    Get exact eigen states/energies from H. 
    """

    H = 1/T*(t*(hz + hzz) + (T - t)*hx)

    return sp.linalg.eigh(H)

def CheckNorm(t, nQubits, Psi, Hvecs, eps):
    """ 
    Check for numerical error with normalization. 
    """
    norm = 0.0

    for i in range(0, 2**nQubits): 
        norm += abs(sp.vdot(Psi, Hvecs[i]))**2

    if ( (1.0 - norm) > eps ): 
        print (str((1.0 - norm)) + " (norm error) > " + str(eps) + 
               " (eps) @ t = " + str(t))

def ExpPert(nQubits, hz, hzz, hx, Psi, T, dt, errchk, eps, outinfo):
    """ 
    Solve using exponential perturbation theory (i.e. Magnus expansion).
    """

    if outinfo['eigdat'] or outinfo['eigplot']:
        eigspec = []
    if outinfo['overlapdat'] or outinfo['overlapplot']:
        overlap = []

    N = T/dt # steps
    mingap = None

    # Loop over time
    for i in range(0, int(sp.floor(N))):
        t = i*dt
        t0 = (i-1)*dt

        # Approximate Hamiltonian to first term in Magnus expansion (OPTIMIZE)

        # The sign on alpha is flipped so that the output probabilities are
        # interpreted as 1 being up and 0 being down, which is more intuitive.
#        Hexp = -1/(2*T)*((t**2 - t0**2)*(alpha - beta) + \
#                        (2*T*(t - t0) + t0**2 - t**2)*delta)

        # This Hamiltonian should have user-specified signs already
        Hexp = 1/(2*T)*((t**2 - t0**2)*(hz + hzz) + \
                        (2*T*(t - t0) + t0**2 - t**2)*hx)

        A = linalg.expm(-1j*Hexp)
        Psi = A*Psi

        # Get eigendecomposition of true Hamiltonian if necessary
        if (errchk or outinfo['mingap'] or outinfo['outdat']
            or outinfo['eigdat'] or outinfo['eigplot']
            or outinfo['fiddat'] or outinfo['fidplot']):
            Hvals, Hvecs = DiagHam(hz, hzz, hx, t, T, nQubits)

            # Sort by eigenvalues
            idx = Hvals.argsort()
            Hvals = Hvals[idx]
            Hvecs = Hvecs[:,idx]
            Hvecs = sp.transpose(Hvecs) # So we can grab them as vectors

            if mingap is None:
                mingap = sp.absolute(Hvals[1] - Hvals[0])
            elif mingap > sp.absolute(Hvals[1] - Hvals[0]):
                mingap = sp.absolute(Hvals[1] - Hvals[0])

        # Check for numerical error
        if (errchk):
            CheckNorm(t, nQubits, Psi, Hvecs, eps)

        # Construct eigenspectrum datapoint = [t, eigval 1, ... , eigval n]
        if (outinfo['eigdat'] or outinfo['eigplot']):
            eigspec.append(output.ConstructEigData(t, Hvals, outinfo['eignum']))

        if (outinfo['overlapdat'] or outinfo['overlapplot']):
            overlap.append(output.ConstructOverlapData(t, Psi, Hvecs[0]))

        # Output our progress, if specified
        if outinfo['progressout']:
            output.ProgressOutput(t, T, outinfo['outdir'])

    # Output stuff as needed
    if (outinfo['eigdat']): 
        output.RecordEigSpec(eigspec, outinfo['outdir'], outinfo['binary'])
    if (outinfo['eigplot']):
        output.PlotEigSpec(eigspec, outinfo['outdir'], T)
    if (outinfo['overlapdat']): 
        output.RecordOverlap(overlap, outinfo['outdir'], T, outinfo['binary'])
    if (outinfo['overlapplot']): 
        output.PlotOverlap(overlap, outinfo['outdir'], T)

    return Psi, mingap
