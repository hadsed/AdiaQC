'''

File: hopfield_opt_gamma.py
Author: Hadayat Seddiqi
Date: 8.20.14
Description: This Hopfield experiment seeks to find the
             optimal Gamma bias weight using classical
             methods. This should indicate whether all
             problems are solvable even in principle.

'''

import os
import numpy as np
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
    """

    # The Hopfield parameters
    hparams = {
        'numNeurons': cmdargs['qubits'],
        'inputState': [ 2*sp.random.random_integers(0,1)-1 
                        for k in xrange(cmdargs['qubits']) ],
        'learningRule': cmdargs['simtype'],
        'numMemories': int(cmdargs['farg'])
        }

    # Construct memories
    memories = [ [ 2*sp.random.random_integers(0,1)-1 
                   for k in xrange(hparams['numNeurons']) ]
                 for j in xrange(hparams['numMemories']) ]

    # At least one pattern must be one Hamming unit away from the input
    memories[0] = list(hparams['inputState'])
    memories[0][sp.random.random_integers(0,hparams['numNeurons']-1)] *= -1

    # Make sure all other patterns have Hamming distance > 1
    def hamdist(a,b):
        """ Calculate Hamming distance. """
        return sp.sum(abs(sp.array(a)-sp.array(b))/2.0)
    # Loop over additional memories, if there are any
    for imem, mem in enumerate(memories[1:]):
        # Avoid duplication and ensure we're far away enough from the input
        while (memories.count(mem) > 1 or 
               hamdist(mem, hparams['inputState']) < 2.0):
            # Flip a random spin
            rndbit = sp.random.random_integers(0,hparams['numNeurons']-1)
            memories[imem+1][rndbit] *= -1
        
    # Basic simulation params
    nQubits = hparams['numNeurons']
    T = 1000.0 # sp.arange(0.1, 15, 0.5)
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
    eigspecnum = 8 # Number of eigenvalues to output
    fidelplot = 0 # Plot fidelity
    fideldat = 0 # Output fidelity data
    fidelnumstates = 2**nQubits # Check fidelity with this number of eigenstates
    overlapdat = 0 # Output overlap data
    overlapplot = 0 # Plot overlap
    solveMethod = 'ExpPert' # 'ExpPert', 'SuzTrot', 'ForRuth', 'BCM'

    # Output directory stuff
    probdir = 'data/hopfield_opt_gamma/n'+str(nQubits)+'p'+\
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

    alpha = sp.array(hparams['inputState'])
    beta = sp.zeros((neurons,neurons))
    delta = sp.array([])

    def kbits(n, k):
        " kbits(n: length of bitstring, k: number of ones) "
        result = []
        for bits in itertools.combinations(range(n), k):
            s = ['0'] * n
            for bit in bits:
                s[bit] = '1'
            result.append(''.join(s))
        return result

    def bitstr2spins(vec):
        """ Return converted list of spins from bitstring @vec. """
        return [ 1 if k == '1' else -1 for k in vec ]

    def calculate_gamma(Jmat, inp, mems):
        """
        Return a threshold+eps for Gamma by analyzing energies
        of the given memory matrix.
        """
        step = 0.001
        eps = 0.1*step
        Glist = np.arange(0,1.+step, step)
        Gthresh = []
        Minp = np.matrix(inp).T
        mem0 = np.matrix(mems[0]).T
        # The input state's energy
        Einp = -Glist*np.sum(Minp.T*Minp) + -np.sum(Minp.T*Jmat*Minp)
        # Loop backwards through Gamma (go high to low)
        for i in range(len(Glist))[::-1]:
            # Calculate new bias matrix
            bias = np.sum(np.array(inp)*np.array(mems[0]))
            # Compute the memory state energies
            lower = -np.sum(mem0.T*Jmat*mem0) - bias
            # Generate bitstrings
            bitstring = []
            for i in range(0,neurons+1): 
                bitstring.append(kbits(neurons, i))
            # Flatten, convert to spins
            spinstr = [ bitstr2spins(item) for item in 
                        list(itertools.chain.from_iterable(bitstring)) ]
            alllist = [ [-np.sum(np.matrix(v)*Jmat*np.matrix(v).T)-bias, v]
                        for v in spinstr ]
            # Sort by energy
            alllist = sorted(alllist, key=lambda x: x[0])
            indices = []
            upper = 0
            # Find the next energy level
            for irec, rec in enumerate(alllist):
                if len(indices) == 2*len(mems):
                    upper = rec[0]
                    break
                if rec[1] in mems:
                    indices.append(irec)
            # If this Gamma value gives the right input state energy,
            # keep it for later analysis
            if lower < Einp[i] and Einp[i] <= upper:
                Gthresh.append(Glist[i])
        Gthresh = list(set(Gthresh))
        # Now pick a Gamma in the middle somewhere
        if len(Gthresh) == 0:
            return 0
        elif len(Gthresh) == 1:
            return Gthresh[0]
        else:
            return Gthresh[(len(Gthresh)-1)/2]

    # Construct the memory matrix according to a learning rule
    if hparams['learningRule'] == 'hebb':
        # Hebb rule
        memMat = sp.matrix(memories).T
        beta = sp.triu(memMat*memMat.T)/float(neurons)
        isingSigns['hz'] *= calculate_gamma(beta, hparams['inputState'], memories)
    elif hparams['learningRule'] == 'stork':
        # Storkey rule
        Wm = sp.zeros((neurons,neurons))
        for m, mem in enumerate(memories):
            Am = sp.outer(mem,mem) - sp.eye(neurons)
            Wm += (Am - Am*Wm - Wm*Am)/float(neurons)
        beta = sp.triu(Wm)
        isingSigns['hz'] *= calculate_gamma(beta, hparams['inputState'], memories)
    elif hparams['learningRule'] == 'proj':
        # Moore-Penrose pseudoinverse rule
        memMat = sp.matrix(memories).T
        beta = sp.triu(memMat * sp.linalg.pinv(memMat))
        isingSigns['hz'] *= calculate_gamma(beta, hparams['inputState'], memories)

    # Some outputs
    outputs = {
        'nQubits': nQubits,
        'learningRule': hparams['learningRule'],
        'outdir': probdir,
        'inputState': hparams['inputState'],
        'memories': memories,
        'bias': -isingSigns['hz'],
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
