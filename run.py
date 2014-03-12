'''

File: run.py
Author: Hadayat Seddiqi
Date: 3.11.14
Description: Runs everything.

'''

import os, shutil, optparse
import json
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
                      type="string", 
                      help="The problem that you want to run.")
    parser.add_option("-r", "--relpath", dest="relpath", default="./",
                      type="string", 
                      help="Relative path of run.py from current dir.")
    parser.add_option("-i", "--instance", dest="instance", default=None,
                      type="int", 
                      help="Instance number for large scale simulations.")
    (options, args) = parser.parse_args()
    problem = options.problem
    relpath = options.relpath
    instance = options.instance

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

# Now import it
try:
    fparams = __import__("problems."+problemPath+problemClean, fromlist=[problem])
except ImportError:
    print ("Unable to import config file for '%s'." % problem)
    raise SystemExit

# Get parameter dict from problem file
cmdargs = {'problem': problem, 
           'relpath': relpath, 
           'instance': instance}
params = fparams.parameters(cmdargs)

# Create data directory
pathpref =  os.path.dirname(os.path.realpath(__file__)) + "/"
try:
    os.makedirs(pathpref + params['outputdir'])
except OSError:
    if not os.path.isdir(pathpref + params['outputdir']):
        raise

# Get the other variables from the problem file
nQubits = params['nQubits']
T = params['T']
dt = params['dt']
errchk = params['errchk']
eps = params['eps']
isingConvert = params['isingConvert']
isingSigns = params['isingSigns']

# Output some params to a file
try:
    params['outputs']
except NameError:
    pass
else:
    with open(pathpref + params['outputdir'] +
              '/networkProperties.dat', 'w') as handle:
        json.dump(params['outputs'], handle)

# Get user-specified coefficients
if (isingConvert):
    alpha = beta = delta = gamma = 0
    Q = params['Q']
else:
    alpha = params['alpha']
    beta = params['beta']
    delta = params['delta']

# Construct output parameters dictionary
outinfo = { 'eigdat': params['eigdat'], 
            'eigplot': params['eigplot'], 
            'eignum': params['eignum'], 
            'fiddat': params['fiddat'], 
            'fidplot': params['fidplot'], 
            'fidnumstates': params['fidnumstates'],
            'overlapdat': params['overlapdat'],
            'overlapplot': params['overlapplot'],
            'outdir': params['outdir'],
            'probout': params['probout'],
            'mingap': params['mingap'],
            'outdat': params['outdat'] }

# Copy the input file to the output dir
shutil.copyfile(relpath+'/problems/'+problem+'.py', 
                relpath+outinfo['outdir']+'/'+problemClean+'.out')

# Turn off all outputs (potentially)
if (params['outputs'] == 0):
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
def RecordStr(string, fname, rpath):
    with open(rpath+outinfo['outdir']+'/'+fname, "w") as file:
        file.write(string)

# Output (append) the minimum gap to file
def RecordMingap(time, gap, fname, it, rpath):
    filepath = rpath+outinfo['outdir'] + '/' + fname
    # Kill ghost data files
    if (it is None or it == 0) and os.path.isfile(filepath):
        os.remove(filepath)
    # Write to file
    with open(filepath, "w") as file:
        file.write(str(gap))

if outinfo['probout']:
    print ("Initial state:")
    print (Psi0)

# Determine if we're doing multiple simulations over T
if isinstance(T, collections.Iterable):
    # Keep the fidelity data somewhere
    if (outinfo['fiddat'] or outinfo['fidplot']):
        fidelitydata = []
    # Go through all the T's
    for i in range(0, len(T)): 
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
                d = solve.output.ConstructFidelityData(
                                    Psi, 
                                    Hvecs[0:outinfo['fidnumstates']], 
                                    T[i], 
                                    outinfo['outdir'])

                for j in range(0, outinfo['fidnumstates']):
                    fidelitydata.append(d[j])
        # Record the mingap and probabilities
        if outinfo['outdat']:
            # Record the minimum spectral gap
            RecordMingap(T[i], mingap, 'mingap.dat', i, relpath)

            # Get state labelings, sort them in descending order
            bitstring = statelabels.GenerateLabels(nQubits)
            bitstring, density = statelabels.SortStateProbabilities(nQubits, Psi, 
                                                                    bitstring)
            finalOutputStr = ''

            for j in range(2**nQubits):
                outstr = bitstring[j] + '\t' + '%.8E' % density[j]
                finalOutputStr += outstr + '\n'
            if outinfo['probout']:
                # Print out the probabilities
                print "\nProbability (T = "+str(T[i])+"):"
                print finalOutputStr

            # Record the final output probabilities
            RecordStr(finalOutputStr, 'probsT'+str(T[i])+'.dat', relpath)

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
            RecordMingap(T, mingap, 'mingap.dat', None, relpath)

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
            RecordStr(finalOutputStr, 'probsT'+str(T)+'.dat', relpath)
