'''

File: solver.py
Author: Hadayat Seddiqi
Date: 3.7.13
Description: Simulate quantum annealing with some technique.

'''

import os
import scipy as sp
from scipy import linalg

def GetEigSpec(t, A, n):
    datapoint = [t]

    E = sp.linalg.eigvalsh(A)
    for i in range(0, n): datapoint.append(E[i].real)

    return datapoint # [t, eigval 1, eigval 2, ... , eigval n]

def RecordEigSpec(eigspec, outputdir):
    eigpath = os.path.dirname(os.path.realpath(__file__)) + "/" + outputdir + "/eigenspectrum.dat"
    sp.savetxt(eigpath, eigspec)

def ExpEvolve(alpha, beta, delta, Psi, T, dt, eigspecparams, outputdir):
    " Evolve in time using sequential matrix exponential. "

    if (eigspecparams[0]): eigspec = []

    N = T/dt # steps

    # Loop over time
    for i in range(1, int(sp.floor(N)) + 1):
        t = i*dt
        t0 = (i-1)*dt

        # Approximate Hamiltonian to first term in Magnus expansion (OPTIMIZE)
        H = 1/(2*T)*((t**2 - t0**2)*(alpha + beta) + (2*T*(t - t0) + t0**2 - t**2)*delta)

        A = linalg.expm(-1j*H)
        Psi = A*Psi

        # Record eigenvalues
        if (eigspecparams[0]): eigspec.append(GetEigSpec(t, H, eigspecparams[2]))


    if (eigspecparams[0]): RecordEigSpec(eigspec, outputdir)

    return Psi
