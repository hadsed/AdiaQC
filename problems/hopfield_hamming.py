'''

File: hopfield_hamming.py
Author: Hadayat Seddiqi
Date: 8.05.14
Description: This problem file is meant to test the
             invariance of choosing "clusters" of
             memories with the same Hamming distance.

'''

import os
import scipy as sp
import itertools
import random
import collections

def parameters(cmdargs):
    """
    cmdargs:
             -q, qubits 
             -k, lrule
             -f, nmems
             -g, hdist
    """

    # The Hopfield parameters
    hparams = {
        'numNeurons': cmdargs['qubits'],
        'inputState': [],
        'learningRule': cmdargs['simtype'],
        'numMemories': int(cmdargs['farg'])
        }

    # Construct memories
    memories = [ [ 2*sp.random.random_integers(0,1)-1 
                   for k in xrange(hparams['numNeurons']) ] ]
    # Add a few memories that are some Hamming distance away from each other
    def hamdist(a,b):
        """ Calculate Hamming distance. """
        return sp.sum(abs(sp.array(a)-sp.array(b))/2.0)
    for m in range(hparams['numMemories']-1):
        # pick a random memory
        newmem = memories[sp.random.random_integers(0,len(memories)-1)][:]
        # flip a random bit
        rndbit = sp.random.random_integers(0,hparams['numNeurons']-1)
        newmem[rndbit] *= -1
        # make sure no duplicates
        while newmem in memories:
            rndbit = sp.random.random_integers(0,hparams['numNeurons']-1)
            newmem[rndbit] *= -1
        # add it
        memories.append(newmem)

    # Choose a random memory for the input state
    hparams['inputState'] = memories[
        sp.random.random_integers(0,len(memories)-1)]
        
    # Basic simulation params
    nQubits = hparams['numNeurons']
    T = 10.0 # sp.arange(0.1, 15, 0.5)
    # T = sp.array([10.0, 20.0, 50.0, 100.0])
    dt = 0.01*T

    # Define states for which to track probabilities in time
    # import statelabels
    # label_list = statelabels.GenerateLabels(nQubits)
    # stateoverlap = []
    # for mem in memories:
    #     # Convert spins to bits
    #     bitstr = ''.join([ '0' if k == 1 else '1' for k in mem ])
    #     # Get the index of the current (converted) memory and add it to list
    #     stateoverlap.append([ label_list.index(bitstr), bitstr ])
    stateoverlap = None

    # Output parameters
    binary = 1 # Save output files as binary Numpy format
    progressout = 0 # Output simulation progress over anneal timesteps

    eigspecdat = 1 # Output data for eigspec
    eigspecplot = 0 # Plot eigspec
    eigspecnum = 16 # Number of eigenvalues to output
    fidelplot = 0 # Plot fidelity
    fideldat = 0 # Output fidelity data
    fidelnumstates = 2**nQubits # Check fidelity with this number of eigenstates
    overlapdat = 0 # Output overlap data
    overlapplot = 0 # Plot overlap
    solveMethod = 'ExpPert' # 'ExpPert', 'SuzTrot', 'ForRuth', 'BCM'

    # Output directory stuff
    probdir = 'data/hopfield_hamming/n'+str(nQubits)+'p'+\
        str(hparams['numMemories'])+hparams['learningRule']
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
    mingap = 0 # Record the minimum spectral gap

    errchk = 0 # Error-checking on/off (for simulation accuracy)
    eps = 0.01 # Numerical error in normalization condition (1 - norm < eps)

    # Specify a QUBO (convert to Ising = True), or alpha, beta directly 
    # (convert = False), and also specify the signs on the Ising Hamiltonian 
    # terms (you can specify coefficients too for some problems if needed)
    isingConvert = 0
    isingSigns = {'hx': -1, 'hz': -1, 'hzz': -1}

    # Construct network Ising parameters
    neurons = nQubits

    # This is gamma, the appropriate weighting on the input vector
    # isingSigns['hz'] *= 1 - (len(hparams['inputState']) - 
    #                          hparams['inputState'].count(0))/(2*neurons)
    # isingSigns['hz'] *= 1.0/(5*neurons)

    alpha = sp.array(hparams['inputState'])
    beta = sp.zeros((neurons,neurons))
    delta = sp.array([])

    # Construct the memory matrix according to a learning rule
    if hparams['learningRule'] == 'hebb':
        # Hebb rule
        isingSigns['hz'] *= 0.7
        memMat = sp.matrix(memories).T
        beta = sp.triu(memMat*memMat.T)/float(neurons)
    elif hparams['learningRule'] == 'stork':
        # Storkey rule
        isingSigns['hz'] *= 0.15
        Wm = sp.zeros((neurons,neurons))
        for m, mem in enumerate(memories):
            Am = sp.outer(mem,mem) - sp.eye(neurons)
            Wm += (Am - Am*Wm - Wm*Am)/float(neurons)
        beta = sp.triu(Wm)
    elif hparams['learningRule'] == 'proj':
        isingSigns['hz'] *= 0.15
        # Moore-Penrose pseudoinverse rule
        memMat = sp.matrix(memories).T
        beta = sp.triu(memMat * sp.linalg.pinv(memMat))

    # Some outputs
    outputs = {
        'nQubits': nQubits,
        'learningRule': hparams['learningRule'],
        'outdir': probdir,
        'inputState': hparams['inputState'],
        'memories': memories,
        'answer': memories[0],
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
