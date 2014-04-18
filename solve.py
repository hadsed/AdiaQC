'''

File: solve.py
Author: Hadayat Seddiqi
Date: 3.7.13
Description: Simulate quantum annealing with some technique.

'''

import os
import scipy as sp
import numpy as np
from scipy.sparse import linalg as sla

import output

def CheckNorm(t, nQubits, Psi, Hvecs, eps):
    """ 
    Check for numerical error with normalization. 
    """
    norm = 0.0

    for i in range(0, 2**nQubits): 
        norm += abs(sp.vdot(Psi, Hvecs[:,i]))**2

    if (1.0 - norm) > eps:
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

        # Approximate Hamiltonian to first term in Magnus expansion
        cz = (t**2 - t0**2)/(2*T)
        cx = (2*T*(t - t0) + t0**2 - t**2)/(2*T)
        Psi = sla.expm_multiply(-1j*(cx*hx + cz*(hz + hzz)), Psi)
        # This is a HUGE performance loss -- requires sparse to dense
        # A = sla.expm(-1j*(cx*hx + cz*(hz + hzz)))
        # Psi = A*Psi

        # Get eigendecomposition of true Hamiltonian if necessary
        if (errchk or outinfo['mingap']
            or outinfo['eigdat'] or outinfo['eigplot']
            or outinfo['fiddat'] or outinfo['fidplot']):
            # Unfortunately we cannot compute all eigenpairs
            Hvals, Hvecs = sla.eigsh(cx*hx + cz*(hz + hzz), 
                                     k=2**nQubits-1,
                                     which='SM')
            # Sort by eigenvalues
            idx = Hvals.argsort()
            Hvals = Hvals[idx]
            Hvecs = Hvecs[:,idx]

            if mingap is None:
                mingap = [sp.absolute(Hvals[1] - Hvals[0]), t/T]
            elif mingap[0] > sp.absolute(Hvals[1] - Hvals[0]):
                mingap = [sp.absolute(Hvals[1] - Hvals[0]), t/T]

        # Check for numerical error
        if (errchk):
            CheckNorm(t, nQubits, Psi, Hvecs, eps)

        # Construct eigenspectrum datapoint = [t, eigval 1, ... , eigval n]
        if (outinfo['eigdat'] or outinfo['eigplot']):
            eigspec.append(output.ConstructEigData(t, Hvals, outinfo['eignum']))

        if (outinfo['overlapdat'] or outinfo['overlapplot']):
            overlap.append(output.ConstructOverlapData(t, Psi, Hvecs[:,0]))

        # Output our progress, if specified
        if outinfo['progressout']:
            output.ProgressOutput(t, T, outinfo['outdir'])

        # Output the overlap with pattern vectors
        if outinfo['stateoverlap'] is not None:
            output.StateOverlapOutput(t, outinfo, Psi)

    # Output stuff as needed
    if (outinfo['eigdat']): 
        output.RecordEigSpec(eigspec, outinfo['outdir'], outinfo['binary'])
    if (outinfo['eigplot']):
        output.PlotEigSpec(eigspec, outinfo['outdir'], T)
    if (outinfo['overlapdat']): 
        output.RecordOverlap(overlap, outinfo['outdir'], T, outinfo['binary'])
    if (outinfo['overlapplot']): 
        output.PlotOverlap(overlap, outinfo['outdir'], T)
    if outinfo['stateoverlap'] is not None:
        output.StateOverlapLabelsOutput(t, outinfo)

    return Psi, mingap

# def edsolver(nQubits, hz, hzz, hx, Psi, T, dt, errchk, eps, outinfo):
#     """ 
#     Solve using edsolver.m.
#     """

#     from oct2py import octave
#     N = sp.floor(T/dt+1)
#     probs, eigk, Psi = octave.call('edsolver/solve1b.m', 
#                                    float(nQubits), hz, hzz, 
#                                    float(N), 1., dt, float(2**nQubits-2))

#     dgap = eigk[1,:] - eigk[0,:]
#     mingap = [ np.min(dgap), np.argmin(dgap)/float(dgap.size) ]

#     times = np.atleast_2d(np.arange(0,N+1))
#     eigspec = np.hstack((times.T, eigk.T))
#     output.PlotEigSpec(eigspec, outinfo['outdir'], T)

#     return Psi, mingap
