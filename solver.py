'''

File: solver.py
Author: Hadayat Seddiqi
Date: 3.7.13
Description: Simulate quantum annealing with some technique.

'''

import os
import scipy as sp
from scipy import linalg

def GetEigSpec(t, A):
    E = sp.linalg.eigvals(A)
    return [t, E[0], E[1]]

def RecordEigSpec(eigspec):
    eigpath = os.path.dirname(os.path.realpath(__file__)) + "/data/eigenspectrum.dat"
    sp.savetxt(eigpath, eigspec)

def ExpEvolve(alpha, beta, delta, Psi, T, dt, eigspecflag):
    " Evolve in time using sequential matrix exponential. "
    i = 0
    t = 0
    N = T/dt

    if (eigspecflag): eigspec = []

    # Loop over time
    while (i <= N) :
        if (i != 1) : t0 = t - dt
        else : t0 = 0.0

        # Approximated Hamiltonian
        H = -1j*(((t**2 - t0**2)/(2*T))*(alpha + beta) + (t - t0 + (t0**2 - t**2)/(2*T))*delta)

        A = linalg.expm(H)
        Psi = A*Psi

        # Record eigenvalues
        if (eigspecflag): eigspec.append(GetEigSpec(t, H))

        i += 1
        t += dt

    if (eigspecflag): RecordEigSpec(eigspec)

    return Psi
