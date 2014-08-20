'''
Analyze Gamma thresholds and other energetic properties of
classical Hopfield Hamiltonians.
'''

import numpy as np
import scipy as sp
from scipy import linalg
import itertools

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

def hamdist(a,b):
    """ Calculate Hamming distance. """
    return sp.sum(abs(sp.array(a)-sp.array(b))/2.0)

def hebb(neurons, memories):
    W = []
    for m, mem in enumerate(memories):
        W.append(sp.triu(sp.outer(mem, mem)-sp.eye(neurons)))
    return np.triu(np.sum(W, axis=0))/float(neurons)

def stork(neurons, memories):
    Wm = sp.zeros((neurons,neurons))
    for m, mem in enumerate(memories):
        Am = sp.outer(mem,mem) - sp.eye(neurons)
        Wm += (Am - Am*Wm - Wm*Am)/float(neurons)
    return sp.triu(Wm)

def stork2(neurons, memories):
    memMat = sp.zeros((neurons,neurons))
    for m, mem in enumerate(memories):
        for i in range(neurons):
            for j in range(i, neurons):
                hij = sp.sum([ memMat[i,k]*mem[k] for k in range(neurons) ])
                hji = sp.sum([ memMat[j,k]*mem[k] for k in range(neurons) ])
                # Don't forget to make the normalization a float!
                memMat[i,j] += (mem[i]*mem[j] - mem[i]*hji - hij*mem[j])/float(neurons)
    return sp.triu(memMat)

def proj(neurons, memories):
    memMat = sp.matrix(memories).T
    return sp.triu(memMat * sp.linalg.pinv(memMat))

# Neurons
neurons = 6

# Memories
# hadamard = sp.linalg.hadamard(neurons).tolist()
# memories = np.array(hadamard)[[0,1,2]].tolist()
# memories = np.array(hadamard)[:].tolist()
patterns = 4
rand = False

memories = [[1, -1, -1, 1, -1, -1],
            [-1, -1, -1, 1, -1, -1],
            # [-1, 1, -1, 1, -1, 1],
            # [-1, 1, -1, -1, -1, -1],
            [1, 1, 1, -1, -1, -1],
            [-1, 1, 1, 1, 1, 1]]

if rand is True:
    memories = [ [ 2*sp.random.random_integers(0,1)-1 
                   for k in xrange(neurons) ]
                 for j in xrange(patterns) ]
    # memories = zip(*rmems)
    memories = np.matrix(rmems).T.tolist()
print "Memories:"
print memories

# Ising matrix
# J = hebb(neurons, memories)
J = stork(neurons, memories)
# J = proj(neurons, memories)
print "Ising matrix:"
print -J
print "Ising energy:"
print np.sum(-J)
print "Fixed-point energies:"
for mem in memories:
    m = np.matrix(mem).T
    print -np.sum(m.T*J*m), mem

# Generate bitstrings
bitstring = []
for i in range(0,neurons+1): 
    bitstring.append(kbits(neurons, i))
# Flatten, convert to spins
spinstr = [ bitstr2spins(item) for item in 
            list(itertools.chain.from_iterable(bitstring)) ]
print "Other guys:"
alllist = []
for vec in spinstr:
    v = np.matrix(vec).T
    alllist.append([-np.sum(v.T*J*v), vec])
alllist = sorted(alllist, key=lambda x: x[0])
for vec in alllist:
    print vec

# Set input state
inp = np.matrix(memories[0]).T
# Flip a random spin
# inp[sp.random.random_integers(0,neurons-1)] *= -1
# inp *= -1
# inp = np.matrix(hadamard[-1]).T
# Calculate bias energies
step = 0.001
G = np.arange(0,1.+step,step)
Eb = -G*(inp.T*inp)[0,0]
Ej = -inp.T*J*inp
E = np.ravel(Eb+Ej)
print "Input:"
print Ej[0,0], inp.T
print "Gamma / Total Energy threshold:"
print (G[E >= alllist[0][0]][::-1][0], E[E >= alllist[0][0]][::-1][0])
# print "All G-energy pairs:"
# print [ np.around(k, 5).tolist() for k in zip(G,E) ]
test = np.ravel(inp)
# test[sp.random.random_integers(0,neurons-1)] *= -1
print hamdist(test, memories[0])
print np.inner(test, memories[0])
