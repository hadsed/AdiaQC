'''

File: solver.py
Author: Hadayat Seddiqi
Date: 3.7.13
Description: Simulate quantum annealing with some technique.

'''

import os
import scipy as sp
from scipy import linalg

def RecordEigSpec(eigspec, outputdir):
    eigpath = os.path.dirname(os.path.realpath(__file__)) + "/" + outputdir + "/eigenspectrum.dat"
    sp.savetxt(eigpath, eigspec)

def PlotEigSpec(eigspec, outputdir):
    import pylab as pl

    eigpath = os.path.dirname(os.path.realpath(__file__)) + "/" + outputdir + "/eigenspectrum.png"

    # Get columns of eigspec to plot
    t = [ row[0] for row in eigspec ]
    for i in range(1,len(eigspec[0])): 
        pl.plot(t, [ row[i] for row in eigspec ])

    pl.xlabel('Time')
    pl.ylabel('Energy')
    pl.savefig(eigpath)

def ExpEvolve(nQubits, alpha, beta, delta, Psi, T, dt, eigspecparams, outputdir, errchk, eps):
    " Evolve in time using sequential matrix exponential. "

    if (eigspecparams['dat']): eigspec = []

    N = T/dt # steps

    # Loop over time
    for i in range(0, int(sp.floor(N)) + 1):
        t = i*dt
        t0 = (i-1)*dt

        # Approximate Hamiltonian to first term in Magnus expansion (OPTIMIZE)
        Hexp = 1/(2*T)*((t**2 - t0**2)*(alpha + beta) + (2*T*(t - t0) + t0**2 - t**2)*delta)

        A = linalg.expm(-1j*Hexp)
        Psi = A*Psi

        if (errchk | eigspecparams['dat']):
            # Get exact eigen states/energies from H
            H = 1/T*(t*(alpha + beta) + (T - t)*delta)
            Hvals, Hvecs = sp.linalg.eigh(H)
        
        if (errchk):
            # Normalization condition
            norm = 0.0
            for i in range(0, 2**nQubits): 
                braket = sp.vdot(Psi, Hvecs[i])
                norm += sp.vdot(braket, braket).real

            if ( (1.0 - norm) > eps ): 
                print ("Error is greater than eps = " + str(eps) + ": " + str((1.0 - norm)))

        # Construct eigenspectrum datapoint = [t, eigval 1, ... , eigval n]
        if (eigspecparams['dat']): 
            datapoint = [t]
            for i in range(0, eigspecparams['num']): datapoint.append(Hvals[i].real)

            eigspec.append(datapoint)


    if (eigspecparams['dat']): RecordEigSpec(eigspec, outputdir)
    if (eigspecparams['plot']): PlotEigSpec(eigspec, outputdir)

    return Psi
