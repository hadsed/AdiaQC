'''

File: random.py
Author: Hadayat Seddiqi
Date: 3.25.14
Description: Run a random case of n-qubits.

'''

import os
import scipy as sp
import itertools
import random
import collections

def parameters(cmdargs):
    """
    """

    # Basic simulation params
    nQubits = int(cmdargs['simtype'])
    T = 10.0  # sp.arange(0.1, 15, 0.5)
    dt = 0.01*T

    # Output parameters
    binary = 1 # Save as binary Numpy
    progressout = 0 # Output simulation progress over anneal timesteps

    eigspecdat = 0 # Output data for eigspec
    eigspecplot = 0 # Plot eigspec
    eigspecnum = 2**nQubits # Number of eigenvalues
    fidelplot = 0 # Plot fidelity
    fideldat = 0 # Output fidelity data
    fidelnumstates = 2**nQubits # Check fidelity with this number of eigenstates
    overlapdat = 0 # Output overlap data
    overlapplot = 0 # Plot overlap

    # Output directory stuff
    probdir = 'data/random_n'+str(nQubits)
    if isinstance(T, collections.Iterable):
        probdir += 'MultiT'
    if os.path.isdir(probdir):
        outlist = sorted([ int(name) for name in os.listdir(probdir) 
                           if name.isdigit() ])
    else:
        outlist = []
    outnum = outlist[-1] + 1 if outlist else 0
    outputdir = probdir + '/' + str(outnum) + '/'

    probshow = 0 # Print final state probabilities to screen
    probout = 0 # Output probabilities to file
    mingap = 0 # Record the minimum spectral gap

    errchk = 0 # Error-checking on/off (for simulation accuracy)
    eps = 0.01 # Numerical error in normalization condition (1 - norm < eps)

    # Specify a QUBO (convert to Ising = True), or alpha, beta directly 
    # (convert = False), and also specify the signs on the Ising Hamiltonian 
    # terms (you can specify coefficients too for some problems if needed)
    isingConvert = 0
    isingSigns = {'hx': -1, 'hz': -1, 'hzz': -1}

    # Generate random matrix of coefficents between [-1,1]
    beta = sp.triu(2*sp.random.ranf((nQubits,nQubits)) - 1)
    alpha = 2*sp.random.ranf(nQubits) - 1
    delta = sp.array([])

    # Some outputs
    outputs = None

    ############################################################################
    ######## All variables must be specified here, do NOT change the keys ######
    ############################################################################

    return {
        'nQubits': nQubits,
        'Q': None,
        'T': T,
        'dt': dt,
        'outputdir': outputdir,
        'errchk': errchk,
        'eps': eps,
        'isingConvert': isingConvert,
        'isingSigns': isingSigns,
        'outputs': outputs,
        'alpha': alpha,
        'beta': beta,
        'delta': delta,
        'eigdat': eigspecdat,
        'eigplot': eigspecplot,
        'eignum': eigspecnum,
        'fiddat': fideldat,
        'fidplot': fidelplot,
        'fidnumstates': fidelnumstates,
        'overlapdat': overlapdat,
        'overlapplot': overlapplot,
        'outdir': outputdir,
        'binary': binary,
        'progressout': progressout,
        'probshow': probshow,
        'probout': probout,
        'mingap': mingap,
        }
