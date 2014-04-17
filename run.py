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
from scipy import linalg as sla
from scipy import sparse as sps

import initialize
import solve
import statelabels
from solve import output

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
    parser.add_option("-k", "--simtype", dest="simtype", default=None,
                      type="string", 
                      help="Further specification of problem type, if needed.")
    parser.add_option("-q", "--qubits", dest="qubits", default=None,
                      type="int", 
                      help="Number of qubits, if needed for batch scripts.")
    parser.add_option("-f", "--farg", dest="farg", default=None,
                      type="string", 
                      help="Just another parameter, if needed.")
    (options, args) = parser.parse_args()
    problem = options.problem
    relpath = options.relpath
    instance = options.instance
    simtype = options.simtype
    qubits = options.qubits
    farg = options.farg

# Clean up the problem path
problemClean = problem.replace('/', '.')
problemPath = ''

# Separate problem path and problem itself
if problemClean.endswith('.py'): 
    problemClean = problemClean[:-3]
while problemClean.rfind('.') > 0:
    idx = problemClean.rfind('.') + 1
    problemPath = problemClean[0:idx]
    problemClean = problemClean[idx:]

# Now import it
try:
    fparams = __import__("problems."+problemPath+problemClean, 
                         fromlist=[problem])
except ImportError:
    print ("Unable to import config file for '%s'." % problem)
    raise SystemExit

# Get parameter dict from problem file
cmdargs = {'problem': problem, 
           'relpath': relpath, 
           'instance': instance,
           'simtype': simtype,
           'qubits': qubits,
           'farg': farg}
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
              '/problem_outputs.dat', 'w') as handle:
        json.dump(params['outputs'], handle)

# Construct output parameters dictionary
outinfo = { 
    'eigdat': params['eigdat'], 
    'eigplot': params['eigplot'], 
    'eignum': params['eignum'], 
    'fiddat': params['fiddat'], 
    'fidplot': params['fidplot'], 
    'fidnumstates': params['fidnumstates'],
    'overlapdat': params['overlapdat'],
    'overlapplot': params['overlapplot'],
    'outdir': params['outdir'],
    'probshow': params['probshow'],
    'probout': params['probout'],
    'mingap': params['mingap'],
    'binary': params['binary'],
    'progressout': params['progressout']
    }

# Copy the input file to the output dir
shutil.copyfile(relpath+'/problems/'+problemClean+'.py', 
                relpath+outinfo['outdir']+'/'+problemClean+'.out')

# Get our initial Hamiltonian coefficients
if (isingConvert):
    # Get Ising coefficients
    h, J, g = initialize.QUBO2Ising(params['Q'])
    # Construct operators
    hz, hzz, hx = initialize.IsingHamiltonian(nQubits, h, J, g)
    # Free memory
    del h, J, g, params['Q']
elif (params['alpha'].size == 0 and
      params['beta'].size == 0 and
      params['delta'].size == 0):
    # Get default generated coefficients
    hz, hzz, hx = initialize.HamiltonianGen(nQubits, 
                                            params['alpha'], 
                                            params['beta'], 
                                            params['delta'])
else:
    alpha = beta = delta = 0
    # Check if we need to generate individually
    if (params['alpha'].size == 0): 
        alpha = sp.ones(nQubits)
    else:
        alpha = params['alpha']
    if (params['beta'].size == 0):
        beta = sp.ones((nQubits, nQubits))
    else:
        beta = params['beta']
    if (params['delta'].size == 0):
        delta = sp.ones(nQubits)
    else:
        delta = params['delta']
    if (params['alpha'].size == 1 and params['beta'].size == 1):
        hz, hzz = initialize.AlphaBetaCoeffs(nQubits, alpha, beta)
        hx = initialize.DeltaCoeffs(nQubits, delta)
        # Free up some memory
        del params['alpha'], params['beta'], params['delta']
        del alpha, beta, delta
    else:
        hz, hzz, hx = initialize.HamiltonianGen(nQubits, alpha, beta, delta)
        del params['alpha'], params['beta'], params['delta']
        del alpha, beta, delta

# Initial state
Psi0 = initialize.InitialState(-hx)
Psi = sp.empty(2**nQubits)

# Apply signs to our operators
hz *= isingSigns['hz']
hzz *= isingSigns['hzz']
hx *= isingSigns['hx']

if outinfo['probshow']:
    print ("Initial state:")
    print (Psi0)

# Determine if we're doing multiple simulations over T
if isinstance(T, collections.Iterable):
    # Keep the fidelity data somewhere
    if (outinfo['fiddat'] or outinfo['fidplot']):
        fidelitydata = []
    # Keep the user-specified values for eigspec stuff
    ueigdat = outinfo['eigdat']
    ueigplot = outinfo['eigplot']
    # Go through all the T's
    for i in range(0, len(T)): 
        # If user wants eigenspectrum data/plots and is also doing multiple T's,
        # make sure we output only one file since spectrum is independent of T.
        if ueigdat:
            if i == (len(T) - 1):
                outinfo['eigdat'] = 1
            else:
                outinfo['eigdat'] = 0
        if ueigplot:
            if i == (len(T) - 1):
                outinfo['eigplot'] = 1
            else:
                outinfo['eigplot'] = 0
        # Solve the Schrodinger equation, get back the final state and mingap
        Psi, mingap = solve.ExpPert(nQubits, hz, hzz, hx, Psi0, T[i], dt[i],
                                    errchk, eps, outinfo)
        # Do fidelity stuff
        if outinfo['fiddat'] or outinfo['fidplot']:
            # Get the eigenpairs
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
        # Record the mingap
        if outinfo['mingap']:
            solve.output.RecordMingap(T[i], mingap, 'mingap.dat', i, 
                                      relpath, outinfo)
        # Record probabilities
        if outinfo['probout']:
            # Get state labelings and probabilities
            bitstring = statelabels.GenerateLabels(nQubits)
            density = statelabels.GetProbabilities(nQubits, Psi)

            # Record the final output probabilities
            solve.output.RecordProbs(bitstring, density, 
                                     'probsT'+str(T[i])+'.dat', 
                                     relpath, outinfo)

            # Check if the bitstrings were recorded, if not, record them
            if not os.path.isfile(relpath+outinfo['outdir']+'/statelabels.txt'):
                # Build a nice string for output to file
                finalstr = ''
                for j in range(2**nQubits):
                    finalstr += bitstring[j] + '\n'
                with open(relpath+outinfo['outdir']+'/statelabels.txt', 
                          "w") as file:
                    file.write(finalstr)

            # Incase the user wants this printed to screen
            if outinfo['probshow']:
                finalOutputStr = ''
                # Sort by probability
                bitstring, density = statelabels.SortStates(nQubits, 
                                                            Psi, 
                                                            bitstring,
                                                            density)
                # Construct a nice-looking string
                for j in range(2**nQubits):
                    outstr = bitstring[j] + '\t' + '%.8E' % density[j]
                    finalOutputStr += outstr + '\n'
                # Print out the probabilities
                print "\nProbability (T = "+str(T[i])+"):"
                print finalOutputStr

    # Sort fidelity data (so the plots come out correctly)
    if (outinfo['fiddat'] or outinfo['fidplot']): 
        fidelitydata, fidelitydataplot = \
            solve.output.SortFidelity(outinfo['fidnumstates'], fidelitydata)
    # Write out fidelity data
    if outinfo['fiddat']: 
        solve.output.RecordFidelity(fidelitydata, outinfo['outdir'], 
                                    outinfo['binary'])
    # Plot fidelity(T)
    if outinfo['fidplot']: 
        solve.output.PlotFidelity(fidelitydataplot, outinfo['outdir'],
                                  outinfo['fidnumstates'])
else:
    Psi, mingap = solve.ExpPert(nQubits, hz, hzz, hx, Psi0, T, dt, 
                                errchk, eps, outinfo)
    # Psi, mingap = solve.edsolver(nQubits, hz, hzz, hx, Psi0, T, dt, 
    #                              errchk, eps, outinfo)

    # Output the minimal spectral gap
    if outinfo['mingap']:
        solve.output.RecordMingap(T, mingap, 'mingap.dat', None, 
                                  relpath, outinfo)
    # Output the probabilities
    if outinfo['probout']:
        sp.set_printoptions(precision=16)
        # Get state labelings, sort them in descending order
        bitstring = statelabels.GenerateLabels(nQubits)
        # Get probability densities
        density = statelabels.GetProbabilities(nQubits, Psi)
        # Output to file
        solve.output.RecordProbs(bitstring, density, 
                                 'probsT'+str(T)+'.dat', 
                                 relpath, outinfo)
        # Check if the bitstrings were recorded, if not, record them
        if not os.path.isfile(relpath+outinfo['outdir']+'/statelabels.txt'):
            # Build a nice string for output to file
            finalstr = ''
            for j in range(2**nQubits):
                finalstr += bitstring[j] + '\n'
            with open(relpath+outinfo['outdir']+'/statelabels.txt', 
                      "w") as file:
                file.write(finalstr)

        # Output probabilities to screen if user wants it
        if outinfo['probshow']:
            finalOutputStr = ''
            bitstring, density = statelabels.SortStates(nQubits, 
                                                        bitstring,
                                                        density)
            print ("Probability (T = "+str(T)+"):\n")
            for i in range(2**nQubits):
                outstr = bitstring[i] + '\t' + '%.8E' % density[i]
                finalOutputStr += outstr + '\n'
            print finalOutputStr
