'''

File: run.py
Author: Hadayat Seddiqi
Date: 3.15.13
Description: Runs everything.

'''

import os
import optparse
import collections
import scipy as sp
from scipy import linalg

import initialize
import solve
import statelabels

# Command line options
if __name__=="__main__":
    parser = optparse.OptionParser("usage: %prog [options] arg1 arg2")
    parser.add_option("-p", "--problem", dest="problem", default="",
                      type="string", help="The problem that you want to run.")
    (options, args) = parser.parse_args()
    problem = options.problem

# Import proper problem file
if problem.endswith('.py'): problem = problem[:-3]

try:
    params = __import__("problems." + problem, fromlist=[problem])
except ImportError:
    print ("Unable to import config file for '%s'." % problem)
    raise SystemExit

# Create data directory
pathpref =  os.path.dirname(os.path.realpath(__file__)) + "/"
#os.makedirs(pathpref + params.outputdir, exist_ok=True) # The Python 3 way
try:
    os.makedirs(pathpref + params.outputdir)
except OSError:
    if not os.path.isdir(pathpref + params.outputdir):
        raise

# Get parameters from problem file
nQubits = params.nQubits
Q = params.Q
a = params.a
T = params.T
dt = params.dt
ising = params.ising
errchk = params.errchk
eps = params.eps

# Get user-specified coefficients
if (ising):
    alpha = beta = delta = 0
else:
    alpha = params.alpha
    beta = params.beta
    delta = params.delta

# Construct output parameters dictionary
outinfo = { 'eigdat': params.eigspecdat, 
            'eigplot': params.eigspecplot, 
            'eignum': params.eigspecnum, 
            'fiddat': params.fideldat, 
            'fidplot': params.fidelplot, 
            'fidnumstates': params.fidelnumstates,
            'overlapdat': params.overlapdat,
            'overlapplot': params.overlapplot,
            'outdir': params.outputdir,
            'probout': params.probout }

# Turn off all outputs (potentially)
if (params.output == 0):
    for param in outinfo: outinfo[param] = 0

# Get our initial Hamiltonian coefficients
if (ising):
    # Get Ising coefficients
    h, J = initialize.QUBO2Ising(Q, a)
    alpha, beta, delta = initialize.IsingHamiltonian(nQubits, h, J)
elif (alpha.size == 0 & beta.size == 0 & delta.size == 0):
    # Get default generated coefficients
    alpha, beta, delta = initialize.HamiltonianGen(nQubits, alpha, beta, delta)
else:
    # Check if we need to generate individually
    if (alpha.size == 0): alpha = sp.ones(nQubits)
    if (beta.size == 0): beta = sp.ones((nQubits, nQubits))
    if (delta.size == 0): delta = sp.ones(nQubits)

    alpha = initialize.AlphaCoeffs(nQubits, alpha)
    beta = initialize.BetaCoeffs(nQubits, beta)
    delta = initialize.DeltaCoeffs(nQubits)

# Initial state
Psi0 = initialize.InitialState(delta)
Psi = sp.empty(2**nQubits)

if outinfo['probout']:
    print ("Initial state:")
    print (Psi0)

# Determine if we're doing multiple simulations over T
if isinstance(T, collections.Iterable):
    if outinfo['fiddat']: fidelitydata = []

    for i in range(0, len(T)): # Go through all the T's
        Psi = solve.ExpPert(nQubits, alpha, beta, delta, Psi0, T[i], dt, \
                             errchk, eps, outinfo)

        # Do fidelity stuff
        if ( outinfo['fiddat'] | outinfo['fidplot'] ):
            from solve import output

            Hvals, Hvecs = sp.linalg.eigh(alpha + beta)
            
            # Sort by eigenvalues
            idx = Hvals.argsort()
            Hvals = Hvals[idx]
            Hvecs = Hvecs[:,idx]
            Hvecs = sp.transpose(Hvecs) # So we can grab them as vectors

            # Construct fidelity data
            if outinfo['fiddat']:
                d = solve.output.ConstructFidelityData(Psi, Hvecs[0:outinfo['fidnumstates']], 
                                                       T[i], outinfo['outdir'])

                for i in range(0, outinfo['fidnumstates']): fidelitydata.append(d[i])

    # Sort fidelity data
    if (outinfo['fiddat'] | outinfo['fidplot']): 
        fidelitydata, fidelitydataplot = solve.output.SortFidelity(outinfo['fidnumstates'], fidelitydata)

    # Write out fidelity data
    if outinfo['fiddat']: solve.output.RecordFidelity(fidelitydata, outinfo['outdir'])

    # Plot fidelity(T)
    if outinfo['fidplot']: solve.output.PlotFidelity(fidelitydataplot, outinfo['outdir'],
                                                     outinfo['fidnumstates'])

    if outinfo['probout']:
        # Get state labelings, sort them in descending order
        bitstring = statelabels.GenerateLabels(nQubits)
        bitstring, density = statelabels.SortStateProbabilities(nQubits, Psi, bitstring)

        print ("\nProbability (T = "+str(list(T)[-1])+"):")
        for i in range(2**nQubits):
            outstr = bitstring[i] + '\t' + '%.8E' % density[i]
            print (outstr)

else:
    Psi = solve.ExpPert(nQubits, alpha, beta, delta, Psi0, T, dt, 
                        errchk, eps, outinfo)
    if outinfo['probout']:
        sp.set_printoptions(precision=16)
        print (Psi)

        # Get state labelings, sort them in descending order
        bitstring = statelabels.GenerateLabels(nQubits)
        bitstring, density = statelabels.SortStateProbabilities(nQubits, Psi, bitstring)

        print ("Probability (T = "+str(T)+"):")
        for i in range(2**nQubits):
            outstr = bitstring[i] + '\t' + '%.8E' % density[i]
            print (outstr)
