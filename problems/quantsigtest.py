'''

File: quantsigtest.py
Author: Hadayat Seddiqi
Date: 3.15.14
Description: This is a test problem derived from Boxio et al. from
             "Experimental signature of programmable quantum annealing,"
             (2012) where we try to show that the quantum annealing
             really is quantum. We should get a 17-fold degenerate
             ground state.

'''

import scipy as sp

def parameters(cmdargs):
    """
    """
    # Some basic simulation params
    nQubits = 8
    T = 10.0
    dt = 0.1

    # Output parameters
    binary = 1 # Save output files as binary Numpy format
    progressout = 0 # Output simulation progress over anneal timesteps

    outputdir = 'data/quantsigtest_new/' # In relation to run.py
    eigspecdat = 0 # Output data for eigspec
    eigspecplot = 0 # Plot eigspec
    eigspecnum = 2**nQubits # Number of eigenvalues to output
    fidelplot = 0 # Plot fidelity
    fideldat = 0 # Output fidelity data
    fidelnumstates = 2**nQubits # Check fidelity with this number of eigenstates
    overlapdat = 0 # Output overlap data
    overlapplot = 0 # Plot overlap
    solveMethod = 'ExpPert' # 'ExpPert', 'SuzTrot', 'ForRuth', 'BCM'

    probshow = 1 # Print final state probabilities to screen
    probout = 1 # Output probabilities to file
    mingap = 0 # Record the minimum spectral gap

    errchk = 0 # Error-checking on/off (for simulation accuracy)
    eps = 0.01 # Numerical error in normalization condition (1 - norm < eps)

    # Specify a QUBO (convert to Ising = True), or alpha, beta directly 
    # (convert = False), and also specify the signs on the Ising Hamiltonian 
    # terms (you can specify coefficients too for some problems if needed)
    isingConvert = 0
    isingSigns = {'hx': -1, 'hz': 1, 'hzz': 1}

    # Only if ising = 0; set all to empty SciPy arrays for default coefficients
    alpha = sp.ones(nQubits)
    beta = sp.zeros((nQubits,nQubits))
    delta = sp.array([])
    alpha[0:4] = 1
    alpha[4:] = -1
    beta[0,1] = beta[1,2] = beta[2,3] = beta[3,0] = beta[0,4] = \
        beta[1,5] = beta[2,6] = beta[3,7] = 1
    beta[1,0] = beta[2,1] = beta[3,2] = beta[0,3] = beta[4,0] = \
        beta[5,1] = beta[6,2] = beta[7,3] = 1

    # We must do this because of Boixo's definition of the Ising Hamiltonian,
    # which has an overall negative sign
    alpha = -alpha
    beta = -beta

    # Usually we specify outputs that may be of interest in the form of a dict, 
    # but we don't need any for this problem
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
        'stateoverlap': None,
        'hzscale': None,
        'hzzscale': None,
        'hxscale': None,
        'solveMethod': solveMethod
        }
