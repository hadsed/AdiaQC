'''

File: hopfield.py
Author: Hadayat Seddiqi
Date: 4.5.13
Description: Parameters for a Hopfield neural network.

'''

import scipy as sp

def parameters(cmdargs):
    """
    """
    nQubits = 8
    T = 10.0
    #T = sp.arange(2,23,4.0) # Output a sequence of anneal times
    dt = 0.01*T

    # Output parameters
    binary = 0 # Save output files as binary Numpy format
    progressout = 0 # Output simulation progress over anneal timesteps

    outputdir = 'data/hopfield_hamming_inv/' # In relation to run.py
    eigspecdat = 1 # Output data for eigspec
    eigspecplot = 0 # Plot eigspec
    eigspecnum = 16 # Number of eigenvalues to output
    fidelplot = 0 # Plot fidelity
    fideldat = 0 # Output fidelity data
    fidelnumstates = 2**nQubits # Check fidelity with this number of eigenstates
    overlapdat = 0 # Output overlap data
    overlapplot = 0 # Plot overlap
    solveMethod = 'ExpPert' # 'ExpPert', 'SuzTrot', 'ForRuth', 'BCM'
    # solveMethod = 'SuzTrot' # 'ExpPert', 'SuzTrot', 'ForRuth', 'BCM'

    probshow = 1 # Print final state probabilities to screen
    probout = 1 # Output probabilities to file
    mingap = 0 # Record the minimum spectral gap

    errchk = 0 # Error-checking on/off (for simulation accuracy)
    eps = 0.01 # Numerical error in normalization condition (1 - norm < eps)

    # Specify a QUBO (convert to Ising = True), or alpha, beta directly 
    # (convert = False), and also specify the signs on the Ising Hamiltonian 
    # terms (you can specify coefficients too for some problems if needed)
    isingConvert = 0
    isingSigns = {'hx': -1, 'hz': -1, 'hzz': -1}

    def spins2bitstr(vec):
        """ Return a converted bitstring from @vec, a list of spins. """
        return ''.join([ '0' if k == 1 else '1' for k in vec ])

    # Construct network Ising parameters
    neurons = nQubits

    # if nQubits % 2 == 0:
    #     memories = sp.linalg.hadamard(neurons).tolist()
    # inputstate = memories[0]

    # memories = [[ 1,-1,-1, 1, 1,-1,-1, 1],
    #             [-1,-1,-1, 1, 1,-1,-1, 1],
    #             [ 1,-1,-1, 1, 1,-1,-1,-1]]
    memories = [[ 1,-1,-1, 1,-1,-1,-1, 1],
                [ 1,-1,-1, 1,-1,-1,-1, 1],
                [ 1,-1,-1, 1,-1,-1,-1,-1]]
    inputstate = memories[0]

    print("Input:")
    print(spins2bitstr(inputstate))
    print('')
    print("Memories:")
    for mem in memories:
        print(spins2bitstr(mem))
    print('')
          
    # This is gamma, the appropriate weighting on the input vector
    # isingSigns['hz'] *= 1 - (len(inputstate) - inputstate.count(0))/(2*neurons)
    # isingSigns['hz'] *= 0.6772598397
    isingSigns['hz'] *= 1.0
    print("Gamma: ", isingSigns['hz'])
    print('')

    alpha = sp.array(inputstate)
    beta = sp.zeros((neurons,neurons))
    delta = sp.array([])

    # Hebb rule - even better matrix style
    memMat = sp.matrix(memories).T
    beta = (sp.triu(memMat*memMat.T)-len(memories)*sp.eye(nQubits))/float(neurons)
    print beta

    # # Hebb rule - matrix style
    # for m, mem in enumerate(memories):
    #     beta += (sp.outer(mem, mem))-sp.eye(neurons)
    # beta = sp.triu(beta)/float(neurons)

    # beta = sp.zeros((neurons,neurons))

    # # Construct pattern matrix according to the Hebb learning rule
    # for i in range(neurons):
    #     for j in range(neurons):
    #         for p in range(len(memories)):
    #             beta[i,j] += ( memories[p][i]*memories[p][j] -
    #                            len(memories)*(i == j) )
    # beta = sp.triu(beta)/float(neurons)

    # Storkey rule (incorrect, but may be useful later)
    # memMat = sp.zeros((neurons,neurons))
    # for m, mem in enumerate(memories):
    #     for i in range(neurons):
    #         for j in range(neurons):
    #             hij = sp.sum([ memMat[i,k]*mem[k] for k in range(neurons) ])
    #             hji = sp.sum([ memMat[j,k]*mem[k] for k in range(neurons) ])
    #             # Don't forget to make the normalization a float!
    #             memMat[i,j] += (mem[i]*mem[j] - mem[i]*hji - hij*mem[j])/float(neurons)
    # beta = sp.triu(memMat)

    # # Storkey rule
    # Wm = sp.zeros((neurons,neurons))
    # for m, mem in enumerate(memories):
    #     Am = sp.mat((sp.outer(mem,mem) - sp.eye(neurons)))
    #     Wm += (Am - Am*Wm - Wm*Am)/float(neurons)
    # beta = sp.triu(Wm)

    # # Construct memory matrix according to the Moore-Penrose pseudoinverse rule
    # memMat = sp.matrix(memories).T
    # beta = sp.triu(memMat * sp.linalg.pinv(memMat))

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
