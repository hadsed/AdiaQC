'''

File: run.py
Author: Hadayat Seddiqi
Date: 3.15.13
Description: Runs everything.

'''

import os
import optparse
import scipy as sp
from scipy import linalg

import initialize
import solver
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
Q = params.Q
a = params.a
T = params.T
dt = params.dt
nQubits = params.nQubits
eigspecflag = params.eigspecflag

# Get Ising coefficients
alpha, beta, delta = initialize.QUBO2Ising(nQubits, Q, a)

# Initial state
Psi = initialize.InitialState(delta)

print ("Initial state:")
print (Psi)

# Evolve in time
Psi = solver.ExpEvolve(alpha, beta, delta, Psi, T, dt, eigspecflag)

# Get state labelings, sort them in descending order
bitstring = statelabels.GenerateLabels(nQubits)
bitstring, density = statelabels.SortStateProbabilities(nQubits, Psi, bitstring)

print ("Probability:")
for i in range(2**nQubits):
    outstr = bitstring[i] + '\t' + '%.8E' % density[i]
    print (outstr)
