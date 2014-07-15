'''

File: hopfield_exp2.py
Author: Hadayat Seddiqi
Date: 5.3.14
Description: Hopfield experiment #2. This experiment is just
             like experiment #1, except that a particular
             problem instance is defined by its memory set, 
             and three sub-instances are defined with differing
             input states that correspond to Hamming distance
             equal to 1 from each pattern, showing that all
             patterns can be retrieved (or not).

'''

import os
import scipy as sp
import scipy.random.random_integers as sp_rint
import itertools
import random
import collections
import cPickle as pickle

def parameters(cmdargs):
    """
    """

    # The Hopfield parameters
    hparams = {
        'numNeurons': cmdargs['qubits'],
        'inputState': [ 2*sp.random.random_integers(0,1)-1 
                        for k in xrange(cmdargs['qubits']) ],
        'learningRule': cmdargs['simtype'],
        'numMemories': int(cmdargs['farg'])
        }

    # Instances will be actual instances x # of patterns
    instIdx = cmdargs['instance'] / hparams['numMemories']
    # Subinstance will be given by the modulus
    subinst = cmdargs['instances'] % hparams['numMemories']

    # Get memories
    memdat = pickle.load(
        open('exp2_memset_n'+str(hparams['numNeurons'])+'.dat', 'rb'))
    memories = memdat[instIdx]
    del memdat

    # Define the input state as a pattern with Hamming dist. 1 from
    # the pattern with index = subinst
    hparams['inputState'] = memories[subinst]
    hparams['inputState'][sp_rint(0,hparams['numNeurons']-1)] *= -1

    # Simulation variables
    nQubits = hparams['numNeurons']
    T = 15.0
    dt = 0.01*T

    # Define states for which to track probabilities in time
    import statelabels
    label_list = statelabels.GenerateLabels(nQubits)
    stateoverlap = []
    for mem in memories:
        # Convert spins to bits
        bitstr = ''.join([ '0' if k == 1 else '1' for k in mem ])
        # Get the index of the current (converted) memory and add it to list
        stateoverlap.append([ label_list.index(bitstr), bitstr ])

    # Output parameters
    binary = 1 # Save as binary Numpy
    progressout = 0 # Output simulation progress over anneal timesteps

    eigspecdat = 1 # Output data for eigspec
    eigspecplot = 0 # Plot eigspec
    eigspecnum = 2**nQubits # Number of eigenvalues
    fidelplot = 0 # Plot fidelity
    fideldat = 0 # Output fidelity data
    fidelnumstates = 2**nQubits # Check fidelity with this number of eigenstates
    overlapdat = 0 # Output overlap data
    overlapplot = 0 # Plot overlap
    solveMethod = 'ExpPert' # 'ExpPert', 'SuzTrot', 'ForRuth', 'BCM'

    # Output directory stuff
    probdir = 'data/hopfield_exp2/n'+str(nQubits)+'p'+str(pnum)+learningRule
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
    probout = 1 # Output probabilities to file
    mingap = 1 # Record the minimum spectral gap

    errchk = 0 # Error-checking on/off (for simulation accuracy)
    eps = 0.01 # Numerical error in normalization condition (1 - norm < eps)

    # Specify a QUBO (convert to Ising = True), or alpha, beta directly 
    # (convert = False), and also specify the signs on the Ising Hamiltonian 
    # terms (you can specify coefficients too for some problems if needed)
    isingConvert = 0
    isingSigns = {'hx': -1, 'hz': -1, 'hzz': -1}

    # This is gamma, the appropriate weighting on the input vector
    isingSigns['hz'] *= 1.0/(5*nQubits)

    neurons = nQubits

    # Initialize Hamiltonian parameters
    alpha = sp.array(hparams['inputState'])
    beta = sp.zeros((neurons,neurons))
    delta = sp.array([])

    # Construct the memory matrix according to a learning rule
    if hparams['learningRule'] == 'hebb':
        # Construct pattern matrix according to Hebb's rule
        for i in range(neurons):
            for j in range(neurons):
                for p in range(len(memories)):
                    beta[i,j] += ( memories[p][i]*memories[p][j] -
                                   len(memories)*(i == j) )
        beta = sp.triu(beta)/float(neurons)
    elif hparams['learningRule'] == 'stork':
        # Construct the memory matrix according to the Storkey learning rule
        memMat = sp.zeros((neurons,neurons))
        for m, mem in enumerate(memories):
            for i in range(neurons):
                for j in range(neurons):
                    hij = sp.sum([ memMat[i,k]*mem[k] for k in range(neurons) ])
                    hji = sp.sum([ memMat[j,k]*mem[k] for k in range(neurons) ])
                    # Don't forget to make the normalization a float!
                    memMat[i,j] += 1./neurons*(mem[i]*mem[j] - mem[i]*hji - 
                                               hij*mem[j])
        beta = sp.triu(memMat)
    elif hparams['learningRule'] == 'proj':
        # Construct memory matrix according to the Moore-Penrose pseudoinverse rule
        memMat = sp.matrix(memories).T
        beta = sp.triu(memMat * sp.linalg.pinv(memMat))

    # Calculate Hamming distance between input state and each memory
    hammingDistance = []
    for mem in memories:
        dist = sp.sum(abs(sp.array(hparams['inputState'])-sp.array(mem))/2)
        hammingDistance.append(dist)

    hamMean = sp.average(hammingDistance)
    hamMed = sp.median(hammingDistance)

    # Some outputs
    outputs = {
        'nQubits': nQubits,
        'learningRule': hparams['learningRule'],
        'outdir': probdir,
        'hparams['inputState']': hparams['inputState'],
        'memories': memories,
        'hammingDistance': {'dist': hammingDistance,
                            'mean': hamMean,
                            'median': hamMed },
        'annealTime': list(T) if isinstance(T, collections.Iterable) else T
               }

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
        'stateoverlap': stateoverlap,
        'hzscale': None,
        'hzzscale': None,
        'hxscale': None,
        'solveMethod': solveMethod
        }
