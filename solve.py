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
from scipy import sparse as sps

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
            if outinfo['eignum'] == 2**nQubits:
                # This is very expensive!!
                Hvals, Hvecs = sp.linalg.eigh((cx*hx + cz*(hz + hzz)).todense())
            else:
                Hvals, Hvecs = sla.eigsh(cx*hx + cz*(hz + hzz), 
                                         k=outinfo['eignum'],
                                         which='SA')
            # Sort by eigenvalues
            idx = Hvals.argsort()
            Hvals = Hvals[idx]/dt
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

def SuzTrot(nQubits, hz, hzz, hx, Psi, T, dt, errchk, eps, outinfo):
    """ 
    Solve using classical Suzuki-Trotter of 2nd order.

    Given H(t) = tmid * Hf + (1-timd)*Hi

    approximate  expm(cc*H(t))*psi as:

        expm( cc/2 * (tmid*HF) ) * ...
        expm( cc * ((1-tmid)*HI) ) * ...
        expm( cc/2 * (tmid*HF) ) * psi

    """

    if outinfo['eigdat'] or outinfo['eigplot']:
        eigspec = []
    if outinfo['overlapdat'] or outinfo['overlapplot']:
        overlap = []

    N = T/dt # steps
    mingap = None
    expm = sla.expm

    # Hamiltonian coefficients
    # acoef = np.array([0, 1, 0, 0])
    # bcoef = np.array([0.5, 0.5, 0, 0])

    # Holder operators (requires more memory, but it's faster)
    Az = (hz + hzz).copy()
    Ax = hx.copy().todense()

    # Hadamard matrix
    h = sp.linalg.hadamard(2**nQubits)/(2.**(nQubits/2.))

    # Loop over time
    for i in xrange(0, int(sp.floor(N))):
        t = i*dt
        t0 = (i-1)*dt
        
        # Time-dependent coefficients
        cz = (-1j*dt)*t/T
        cx = (-1j*dt)*(1 - t/T)

        # Construct Z and X operators
        Az = 0.5*cz*(hz + hzz)
        Az.setdiag(np.exp(Az.diagonal()))
        Ax = h*(cx*hx.todense())*h
        Ax = h*np.matrix(np.diag(np.exp(np.diag(Ax))))*h

        # Apply
        Psi = (Az*Ax*Az)*Psi

        # Get eigendecomposition of true Hamiltonian if necessary
        if (errchk or outinfo['mingap']
            or outinfo['eigdat'] or outinfo['eigplot']
            or outinfo['fiddat'] or outinfo['fidplot']):
            # Unfortunately we cannot compute all eigenpairs
            if outinfo['eignum'] == 2**nQubits:
                # This is very expensive!!
                Hvals, Hvecs = sp.linalg.eigh((cx*hx + cz*(hz + hzz)).todense())
            else:
                Hvals, Hvecs = sla.eigsh(cx*hx + cz*(hz + hzz), 
                                         k=outinfo['eignum'],
                                         which='SA')
            # Sort by eigenvalues
            idx = Hvals.argsort()
            Hvals = Hvals[idx]/dt
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

def ForRuth(nQubits, hz, hzz, hx, Psi, T, dt, errchk, eps, outinfo):
    """ 
    Forest and Ruth
    4th order, Physica D 43, 105 (1990)

    psi_new = C1 * C2 * C3 * C4 * psi

    C1 =  expm( cc*(1-tmid)*a1 * HI ) * ...
          expm( (cc*tmid)*b1 * HF )

    C2 = expm( (cc*(1-tmid))*a2 * HI ) * ...
         expm( (cc*tmid)*b2 * HF )

    C3 = expm( (cc*(1-tmid)* a3 *HI ) * ...
         expm( (cc*timid)* b3 *HF )

    C4 = expm( (cc*(1-tmid)) a4 *HI ) * ...
         expm( (cc*tmid)* b4 * HF )

    where theta = 1/(2 - 2^(1/3) )
    and cc = -i*dt
    """

    if outinfo['eigdat'] or outinfo['eigplot']:
        eigspec = []
    if outinfo['overlapdat'] or outinfo['overlapplot']:
        overlap = []

    N = T/dt # steps
    mingap = None
    expm = sla.expm

    # Hamiltonian coefficients
    theta = -1./(2-2.**(1/3.))  # requires a negative sign in front?
    acoef = np.array([0, theta, 1-2.*theta, theta])
    bcoef = np.array([theta/2., (1-theta)/2., (1-theta)/2., theta/2.])

    # Holder operators (requires more memory, but it's faster)
    Az = (hz + hzz).copy()
    Ax = hx.copy().todense()
    h = sp.linalg.hadamard(2**nQubits)/(2.**(nQubits/2.))

    # Loop over time
    for i in range(0, int(sp.floor(N))):
        t = i*dt
        t0 = (i-1)*dt
        
        # Time-dependent coefficients
        cz = (-1j*dt)*t/T
        cx = (-1j*dt)*(1 - t/T)
        
        # Run through coefficient indices backwards
        for idx in range(len(acoef))[::-1]:
            if bcoef[idx] != 0:
                # Psi = expm(bcoef[idx]*(hz + hzz))*Psi
                # A = bcoef[idx]*cz*(hz + hzz).copy()
                # A.setdiag(np.exp(A.diagonal()))
                Az = bcoef[idx]*cz*(hz + hzz)
                Az.setdiag(np.exp(Az.diagonal()))
                Psi = Az*Psi
            if acoef[idx] != 0:
                # Psi = expm(acoef[idx]*hx)*Psi
                Ax = acoef[idx]*cx*hx.todense()
                Ax = h*Ax*h
                Ax = np.matrix(np.diag(np.exp(np.diag(Ax))))
                Psi = (h*Ax*h)*Psi

        # Get eigendecomposition of true Hamiltonian if necessary
        if (errchk or outinfo['mingap']
            or outinfo['eigdat'] or outinfo['eigplot']
            or outinfo['fiddat'] or outinfo['fidplot']):
            # Unfortunately we cannot compute all eigenpairs
            if outinfo['eignum'] == 2**nQubits:
                # This is very expensive!!
                Hvals, Hvecs = sp.linalg.eigh((cx*hx + cz*(hz + hzz)).todense())
            else:
                Hvals, Hvecs = sla.eigsh(cx*hx + cz*(hz + hzz), 
                                         k=outinfo['eignum'],
                                         which='SA')
            # Sort by eigenvalues
            idx = Hvals.argsort()
            Hvals = Hvals[idx]/dt
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

def BCM(nQubits, hz, hzz, hx, Psi, T, dt, errchk, eps, outinfo):
    """ 

    Blanes-Casas-Nurua, P_38 2 formula (2nd order with large stability)

    Symplectic splitting operator methods for the time-dependent 
    Schrodinger equation by Sergio Blanes, Fernando Casas, and 
    Ander Nurua

    Citation: J. Chem. Phys. 124, 234105 (2006)

    """

    if outinfo['eigdat'] or outinfo['eigplot']:
        eigspec = []
    if outinfo['overlapdat'] or outinfo['overlapplot']:
        overlap = []

    N = T/dt # steps
    mingap = None
    expm = sla.expm

    # Hamiltonian coefficients
    acoef = np.zeros(40)
    bcoef = np.zeros(40)
    aclist = [0.0215672851797585075705350295278,
              0.0431726343853101639735369714998,
              0.0431324297795690599949127838602,
              0.0427852961505675320118200419401,
              0.0449747930772476869948630891275,
              0.521477840977180737598212898081,
              -0.460297865581209561666776462059,
              0.0476657723717784446737564703982,
              -0.299809415632442402707251772031,
              0.360890555491738732398154005651,
              0.0355310860247975525993505717327,
              0.0451459109591929143698396854787,
              0.151663982419594313475358779605,
              -0.122723981192628473398202625228,
              -0.0342003644722802255132523920962,
              0.0514702802470565594888643277103,
              -0.00346916149683374374401491713903,
              0.0201046430669616823814202845610,
              -0.0245251277750599926319683675996]
    aclist.append(1 - 2*(np.sum(aclist)))
    acoef[0:19] = aclist
    acoef[39:20:-1] = aclist[0:19]

    bclist = [0.0431461454881085359990876258277,
              0.0431853234593364152087490292063,
              0.0429704744650982147539363885468,
              0.0430364300871454499243887883740,
              0.0532805678508921227350798781968,
              -0.0000741632590652008982349604299511,
              0.0549252685049280768846009673282,
              0.0572922318289063436814214008313,
              -0.000216083699929765754852184048464,
              0.0429262827299850710231689679598,
              0.0509590583382259625517957082533,
              0.0125876466303119396367352929903,
              -0.00110143601875055751217588524309,
              0.0589864485893508739845735668507,
              -0.00393919091210338198661577774009,
              0.0909189791588641823686791563103,
              -0.107654717879545729464023522278,
              0.0254278113893309936197644680648]
    bclist.append(0.5 - np.sum(bclist))
    bcoef[0:19] = bclist
    bcoef[38:19:-1] = bclist[0:19]

    # Holder operators (requires more memory, but it's faster)
    Az = (hz + hzz).copy()
    Ax = hx.copy().todense()
    h = sp.linalg.hadamard(2**nQubits)/(2.**(nQubits/2.))

    # Loop over time
    for i in range(0, int(sp.floor(N))):
        t = i*dt
        t0 = (i-1)*dt
        
        # Time-dependent coefficients
        cz = (-1j*dt)*t/T
        cx = (-1j*dt)*(1 - t/T)
        
        # Run through coefficient indices backwards
        for idx in range(len(acoef))[::-1]:
            if bcoef[idx] != 0:
                # Psi = expm(bcoef[idx]*(hz + hzz))*Psi
                # A = bcoef[idx]*cz*(hz + hzz).copy()
                # A.setdiag(np.exp(A.diagonal()))
                Az = bcoef[idx]*cz*(hz + hzz)
                Az.setdiag(np.exp(Az.diagonal()))
                Psi = Az*Psi
            if acoef[idx] != 0:
                # Psi = expm(acoef[idx]*hx)*Psi
                Ax = acoef[idx]*cx*hx.todense()
                Ax = h*Ax*h
                Ax = np.matrix(np.diag(np.exp(np.diag(Ax))))
                Psi = (h*Ax*h)*Psi

        # Get eigendecomposition of true Hamiltonian if necessary
        if (errchk or outinfo['mingap']
            or outinfo['eigdat'] or outinfo['eigplot']
            or outinfo['fiddat'] or outinfo['fidplot']):
            # Unfortunately we cannot compute all eigenpairs
            if outinfo['eignum'] == 2**nQubits:
                # This is very expensive!!
                Hvals, Hvecs = sp.linalg.eigh((cx*hx + cz*(hz + hzz)).todense())
            else:
                Hvals, Hvecs = sla.eigsh(cx*hx + cz*(hz + hzz), 
                                         k=outinfo['eignum'],
                                         which='SA')
            # Sort by eigenvalues
            idx = Hvals.argsort()
            Hvals = Hvals[idx]/dt
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
