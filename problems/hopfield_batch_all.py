'''

File: lr_all.py
Author: Hadayat Seddiqi
Date: 3.11.14
Description: Run Hopfield net.

'''

import os
import scipy as sp
import itertools
import random
import collections
import cPickle as pickle

def parameters(cmdargs):
    """
    """

    # import problems.hopfield.params as params

    # learningRule = params.learningRules[params.rule]
    learningRule = cmdargs['simtype']
    # nQubits = params.numQubits
    nQubits = 5
    # T = params.annealTime
    T = sp.arange(0.1, 15, 0.5)
    T = 10.0
    dt = 0.01*T
    inputstate = [1,1,-1,-1,1]

    # Output parameters
    output = 1 # Turn on/off all output except final probabilities
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

    # Output directory stuff
    pnum = 0
    if 32 > cmdargs['instance']:
        pnum = 1
    elif (496-1) < cmdargs['instance']:
        pnum = 3
    else:
        pnum = 2
    # probdir = 'data/all_n'+str(nQubits)+'p'+str(pnum)+learningRule
    probdir = 'problems/all_n'+str(nQubits)+'p'+str(pnum)+learningRule
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

    neurons = nQubits
    memories = []

    # Load the list of memories
    patternSet = pickle.load(open('problems/n5p'+str(pnum)+'mems.dat', 'rb'))
    # Define the right index for each pnum case
    psetIdx = cmdargs['instance']
    if 31 < cmdargs['instance'] < 496:
        psetIdx -= 32
    elif cmdargs['instance'] >= 496:
        psetIdx -= 496
    # Iterate through the set that corresponds to the [corrected] instance number
    for bitstring in patternSet[psetIdx]:
        # Convert binary to Ising spins
        spins = [ 1 if k == '1' else -1 for k in bitstring ]
        memories.append(spins)

    # Make the input the last memory recorded
    # inputstate = params.inputState
    # Add in the input state
    # if params.includeInput and inputstate not in memories:
    #     memories[0] = inputstate
    if inputstate not in memories:
        memories[0] = inputstate

    # This is gamma, the appropriate weighting on the input vector
    isingSigns['hz'] *= 1 - (len(inputstate) - inputstate.count(0))/(2*neurons)

    # Initialize Hamiltonian parameters
    alpha = sp.array(inputstate)
    beta = sp.zeros((neurons,neurons))
    delta = sp.array([])

    # Construct the memory matrix according to a learning rule
    if learningRule == 'hebb':
        # Construct pattern matrix according to Hebb's rule
        for i in range(neurons):
            for j in range(neurons):
                for p in range(len(memories)):
                    beta[i,j] += ( memories[p][i]*memories[p][j] -
                                   len(memories)*(i == j) )
        beta = sp.triu(beta)/float(neurons)
    elif learningRule == 'stork':
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
    elif learningRule == 'proj':
        # Construct memory matrix according to the Moore-Penrose pseudoinverse rule
        memMat = sp.matrix(memories).T
        beta = sp.triu(memMat * sp.linalg.pinv(memMat))

    # Calculate Hamming distance between input state and each memory
    hammingDistance = []
    for mem in memories:
        dist = sp.sum(abs(sp.array(inputstate)-sp.array(mem))/2)
        hammingDistance.append(dist)

    hamMean = sp.average(hammingDistance)
    hamMed = sp.median(hammingDistance)

    # Some outputs
    outputs = {
        'nQubits': nQubits,
        'learningRule': learningRule,
        'outdir': probdir,
        'inputState': inputstate,
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
        }
