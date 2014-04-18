'''

File: sixvarqubo.py
Author: Hadayat Seddiqi
Date: 3.15.14
Description: A six-variable QUBO test problem. Should
             have the following minimum vectors:
             (1,0,0,1,0,1) and (1,1,0,1,0,1).

'''

import scipy as sp

def parameters(cmdargs):
    """
    """
    nQubits = 6
    T = 30.0
    #T = sp.arange(2,103,10) # Output a sequence of anneal times
    dt = 0.1

    # Output parameters
    binary = 1 # Save output files as binary Numpy format
    progressout = 0 # Output simulation progress over anneal timesteps

    outputdir = 'data/sixvarqubo_overlap/' # In relation to run.py
    eigspecdat = 0 # Output data for eigspec
    eigspecplot = 0 # Plot eigspec
    eigspecnum = 2**nQubits # Number of eigenvalues to output
    fidelplot = 0 # Plot fidelity
    fideldat = 0 # Output fidelity data
    fidelnumstates = 2**nQubits # Check fidelity with this number of eigenstates
    overlapdat = 0 # Output overlap data
    overlapplot = 0 # Plot overlap

    probshow = 1 # Print final state probabilities to screen
#### BUG ##### probshow requires probout
    probout = 1 # Output probabilities to file

    mingap = 0 # Record the minimum spectral gap

    errchk = 0 # Error-checking on/off (for simulation accuracy)
    eps = 0.01 # Numerical error in normalization condition (1 - norm < eps)

    # Specify a QUBO (convert to Ising = True), or alpha, beta directly 
    # (convert = False), and also specify the signs on the Ising Hamiltonian 
    # terms (you can specify coefficients too for some problems if needed)
    isingConvert = 1
    isingSigns = {'hx': 1, 'hz': 1, 'hzz': -1}

    # Define pattern list for overlap thing
    import statelabels
    label_list = statelabels.GenerateLabels(nQubits)
    stateoverlap = []
    for bitstr in ['100101', '110101']:
        # Get the index of the current (converted) memory and add it to list
        stateoverlap.append([ label_list.index(bitstr), bitstr ])

    # Define the QUBO and its diagonal
    Q = sp.zeros((nQubits,nQubits))
    a = sp.zeros(nQubits)

    a[0] = 2
    a[1] = 1
    a[2] = -2
    a[3] = -1
    a[4] = 1
    a[5] = -1

    Q[0,1] = Q[1,0] = -1
    Q[0,2] = Q[2,0] = 2
    Q[0,3] = Q[3,0] = -2
    Q[0,4] = Q[4,0] = 2
    Q[0,5] = Q[5,0] = -1

    Q[1,2] = Q[2,1] = 1
    Q[1,3] = Q[3,1] = -1
    Q[1,4] = Q[4,1] = -1
    Q[1,5] = Q[5,1] = 1

    Q[2,3] = Q[3,2] = 2
    Q[2,4] = Q[4,2] = -2
    Q[2,5] = Q[5,2] = 1

    Q[3,4] = Q[4,3] = 2
    Q[3,5] = Q[5,3] = -1

    Q[4,5] = Q[5,4] = 2

    Q = sp.triu(Q) + a*sp.identity(6)

    # Usually we specify outputs that may be of interest in the form of a dict, 
    # but we don't need any for this problem
    outputs = None

    # We could specify these below, but we want to avoid modifying that part
    alpha = beta = delta = None

    ############################################################################
    ######## All variables must be specified here, do NOT change the keys ######
    ############################################################################

    return {
        'nQubits': nQubits,
        'Q': Q,
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
        'stateoverlap': stateoverlap
        }
