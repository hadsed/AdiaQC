'''

File: run.py
Author: Hadayat Seddiqi
Date: 3.15.13
Description: Runs everything.

'''

import os, shutil, optparse
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

# Parse the dirs
problemClean = problem.replace('/', '.')
problemPath = ''
while problem.rfind('.') > 0:
    idx = problem.rfind('.') + 1
    problemPath = problem[0:idx]
    problem = problem[idx:]

# Import proper problem file
if problemClean.endswith('.py'): 
    problemClean = problem[:-3]

try:
    params = __import__("problems." + problemPath + problemClean, fromlist=[problem])
except ImportError:
    print ("Unable to import config file for '%s'." % problem)
    raise SystemExit

# Create data directory
pathpref =  os.path.dirname(os.path.realpath(__file__)) + "/"

try:
    os.makedirs(pathpref + params.outputdir)
except OSError:
    if not os.path.isdir(pathpref + params.outputdir):
        raise

# Get parameters from problem file
nQubits = params.nQubits
T = params.T
dt = params.dt
errchk = params.errchk
eps = params.eps
isingConvert = params.isingConvert
isingSigns = params.isingSigns

# Get user-specified coefficients
if (isingConvert):
    alpha = beta = delta = gamma = 0
    Q = params.Q
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
            'probout': params.probout,
            'mingap': params.mingap,
            'outdat': params.outdat }

# Copy the input file to the output dir
shutil.copyfile('problems/'+problem+'.py', 
                outinfo['outdir']+'/'+problemClean+'.out')

# Turn off all outputs (potentially)
if (params.output == 0):
    for param in outinfo: outinfo[param] = 0

# Get our initial Hamiltonian coefficients
if (isingConvert):
    # Get Ising coefficients
    h, J, g = initialize.QUBO2Ising(Q)
    hz, hzz, hx = initialize.IsingHamiltonian(nQubits, h, J, g)
elif (alpha.size == 0 & beta.size == 0 & delta.size == 0):
    # Get default generated coefficients
    hz, hzz, hx = initialize.HamiltonianGen(nQubits, alpha, beta, delta)
else:
    # Check if we need to generate individually
    if (alpha.size == 0): alpha = sp.ones(nQubits)
    if (beta.size == 0): beta = sp.ones((nQubits, nQubits))
    if (delta.size == 0): delta = sp.ones(nQubits)
    if (alpha.size & beta.size == 1): 
        hz, hzz = initialize.AlphaBetaCoeffs(nQubits, alpha, beta)
        hx = initialize.DeltaCoeffs(nQubits, delta)
    else:
        hz = initialize.AlphaCoeffs(nQubits, alpha)
        hzz = initialize.BetaCoeffs(nQubits, beta)
        hx = initialize.DeltaCoeffs(nQubits, delta)

# Initial state
Psi0 = initialize.InitialState(-hx)
Psi = sp.empty(2**nQubits)

# Apply signs to our operators
hz *= isingSigns['hz']
hzz *= isingSigns['hzz']
hx *= isingSigns['hx']

# Output a string to file
def RecordStr(string, fname):
    with open(outinfo['outdir']+'/'+fname, "w") as file:
        file.write(string)

# Output (append) the minimum gap to file
def RecordMingap(time, gap, fname, it):
    filepath = outinfo['outdir'] + '/' + fname
    # Kill ghost data files
    if (it is None or it == 0) and os.path.isfile(filepath):
        os.remove(filepath)
    # Write to file
    with open(filepath, "a") as file:
        file.write(str(time)+'\t'+str(gap)+'\n')

if outinfo['probout']:
    print ("Initial state:")
    print (Psi0)

# Determine if we're doing multiple simulations over T
if isinstance(T, collections.Iterable):
    if (outinfo['fiddat'] or outinfo['fidplot']):
        fidelitydata = []

    for i in range(0, len(T)): # Go through all the T's
        Psi, mingap = solve.ExpPert(nQubits, hz, hzz, hx, Psi0, T[i], dt[i],
                                    errchk, eps, outinfo)
        # Do fidelity stuff
        if outinfo['fiddat'] or outinfo['fidplot']:
            # Yeah, yeah, I know this is bad practice, god
            from solve import output

            Hvals, Hvecs = sp.linalg.eigh(hz + hzz)
            
            # Sort by eigenvalues
            idx = Hvals.argsort()
            Hvals = Hvals[idx]
            Hvecs = Hvecs[:,idx]
            Hvecs = sp.transpose(Hvecs) # So we can grab them as vectors

            # Construct fidelity data
            if (outinfo['fiddat'] or outinfo['fidplot']):
                d = solve.output.ConstructFidelityData(Psi, 
                                                       Hvecs[0:outinfo['fidnumstates']], 
                                                       T[i], 
                                                       outinfo['outdir'])

                for j in range(0, outinfo['fidnumstates']):
                    fidelitydata.append(d[j])
        # Record the mingap and probabilities
        if outinfo['outdat']:
            # Record the minimum spectral gap
            RecordMingap(T[i], mingap, 'mingap.dat', i)

            # Get state labelings, sort them in descending order
            bitstring = statelabels.GenerateLabels(nQubits)
            bitstring, density = statelabels.SortStateProbabilities(nQubits, Psi, bitstring)
            finalOutputStr = ''

            for j in range(2**nQubits):
                outstr = bitstring[j] + '\t' + '%.8E' % density[j]
                finalOutputStr += outstr + '\n'
            if outinfo['probout']:
                # Print out the probabilities
                print "\nProbability (T = "+str(T[i])+"):"
                print finalOutputStr

            RecordStr(finalOutputStr, 'probsT'+str(T[i])+'.dat')

    # Sort fidelity data
    if (outinfo['fiddat'] or outinfo['fidplot']): 
        fidelitydata, fidelitydataplot = solve.output.SortFidelity(outinfo['fidnumstates'], 
                                                                   fidelitydata)
    # Write out fidelity data
    if outinfo['fiddat']: 
        solve.output.RecordFidelity(fidelitydata, outinfo['outdir'])
    # Plot fidelity(T)
    if outinfo['fidplot']: 
        solve.output.PlotFidelity(fidelitydataplot, outinfo['outdir'],
                                  outinfo['fidnumstates'])
else:
    Psi, mingap = solve.ExpPert(nQubits, hz, hzz, hx, Psi0, T, dt, 
                                errchk, eps, outinfo)
    # Output the minimal spectral gap
    if outinfo['mingap']:
        print ("\nMinimum spectral gap: "+str(mingap))
        if outinfo['outdat']:
            RecordMingap(T, mingap, 'mingap.dat', None)

    # Output the probabilities
    if outinfo['probout']:
        sp.set_printoptions(precision=16)
        # Get state labelings, sort them in descending order
        bitstring = statelabels.GenerateLabels(nQubits)
        bitstring, density = statelabels.SortStateProbabilities(nQubits, Psi, bitstring)
        finalOutputStr = ''

        print ("Probability (T = "+str(T)+"):\n")
        for i in range(2**nQubits):
            outstr = bitstring[i] + '\t' + '%.8E' % density[i]
            finalOutputStr += outstr + '\n'
        print finalOutputStr

        if outinfo['outdat']:
            RecordStr(finalOutputStr, 'probsT'+str(T)+'.dat')
