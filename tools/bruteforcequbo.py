'''
A brute force QUBO solution finder thing.
'''

import numpy as np
import itertools


# Number of bits
N = 4

# Set the QUBO problem
Q = np.zeros((N,N))
a = np.zeros(N)

a[0] = 1
a[1] = 1
a[2] = 2
a[3] = 2

Q[0,1] = Q[1,0] = -1
Q[0,2] = Q[2,0] = -1
Q[0,3] = Q[3,0] = 5

Q[1,2] = Q[2,1] = 5
Q[1,3] = Q[3,1] = -1

Q[2,3] = Q[3,2] = -1

Q = np.matrix(np.triu(Q) + a*np.identity(N))

# Go through all possible solutions
def kbits(n, k):
    """
    Generate list of n k-length bitstrings.
    """
    result = []
    for bits in itertools.combinations(range(n), k):
        s = ['0'] * n
        for bit in bits:
            s[bit] = '1'
        result.append(''.join(s))
    return result

bitstring = []
for i in range(0,N+1):
    bitstring.append(kbits(N, i))

# Flatten
bitstring = list(itertools.chain.from_iterable(bitstring))

# Compute costs
results = []
for candy in bitstring:
    bitlist = np.matrix([ 0 if b == '0' else 1
                          for b in candy ]).T
    results.append([candy, np.sum(bitlist.T*Q*bitlist)])

results = sorted(results, key=lambda x: x[1])
print "[Bits, Cost]"
for r in results:
    print r
