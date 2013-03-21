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

errchk = params.errchk
eps = params.eps

# Construct output parameters dictionary
outinfo = { 'eigdat': params.eigspecdat, 
            'eigplot': params.eigspecplot, 
            'eignum': params.eigspecnum, 
            'fiddat': params.fideldat, 
            'fidplot': params.fidelplot, 
            'fidnumstates': params.fidelnumstates,
            'outdir': params.outputdir }

# Turn off all outputs (potentially)
if (params.output == 0):
    for param in outinfo: outinfo[param] = 0

# Get Ising coefficients
alpha, beta, delta = initialize.QUBO2Ising(nQubits, Q, a)

# Initial state
Psi = initialize.InitialState(delta)

print ("Initial state:")
print (Psi)

# Determine if we're doing multiple simulations over T
if isinstance(T, collections.Iterable):
    if outinfo['fiddat']: fidelitydata = []

    for i in range(0, len(T)): # Go through all the T's
        Psi = solve.ExpPert(nQubits, alpha, beta, delta, Psi, T[i], dt, \
                             errchk, eps, outinfo)

        # Do fidelity stuff
        if ( outinfo['fiddat'] | outinfo['fidplot'] ):
            from solve import output

            Hvals, Hvecs = sp.linalg.eigh(alpha + beta)

            if outinfo['fiddat']: # Output fidelity data
                d = solve.output.RecordFidelity(Psi, Hvecs[0:outinfo['fidnumstates']], 
                                                T[i], outinfo['outdir'])

                for i in range(0, outinfo['fidnumstates']): fidelitydata.append(d[i])

    # Write out fidelity data
    if outinfo['fiddat']:
        path = os.path.dirname(os.path.realpath(__file__)) + "/" + outinfo['outdir'] + "/fidelity.dat"
        sp.savetxt(path, fidelitydata)

    # Plot fidelity(T)
    if outinfo['fidplot']:
        solve.output.PlotFidelity(fidelitydata, outinfo['outdir'])

else:
    Psi = solver.ExpPert(nQubits, alpha, beta, delta, Psi, T[i], dt, 
                         errchk, eps, outinfo)

# Get state labelings, sort them in descending order
bitstring = statelabels.GenerateLabels(nQubits)
bitstring, density = statelabels.SortStateProbabilities(nQubits, Psi, bitstring)

print ("Probability:")
for i in range(2**nQubits):
    outstr = bitstring[i] + '\t' + '%.8E' % density[i]
    print (outstr)
