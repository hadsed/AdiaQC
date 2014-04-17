'''

File: statelabels.py
Author: Hadayat Seddiqi
Date: 3.15.13
Description: Construct proper labeling scheme for states.
             The NumPy kron() does a kronecker product like:
             kron(A,B) = [ [a00*B, a01*B], [a10*B, a11*B] ]
             So we construct our states similarly, with the 
             same kron() routine.

'''

import scipy as sp
import itertools
from scipy import linalg

def kbits(n, k):
    " kbits(n: length of bitstring, k: number of ones) "
    result = []
    for bits in itertools.combinations(range(n), k):
        s = ['0'] * n
        for bit in bits:
            s[bit] = '1'
        result.append(''.join(s))
    return result

def GenerateLabels(n):
    " Get proper labeling for output states. "
    # Generate bitstrings
    bitstring = []
    for i in range(0,n+1): 
        bitstring.append(kbits(n, i))
    # Flatten
    bitstring = list(itertools.chain.from_iterable(bitstring))
    # Generate unit vectors
    statelist = []
    poslist = []
    pos = 0
    unit0 = sp.array([1,0])
    unit1 = sp.array([0,1])
    for i in range(len(bitstring)):
        # Construct unit vector corresponding to bitstring
        state = unit1 if (bitstring[i][0] == '1') else unit0
        for j in range(n-1):
            state = sp.kron(state, 
                            unit1 if (bitstring[i][j+1] == '1') else unit0)
        statelist.append(state)
        # Record orientation of unit vector (using position of 1 value)
        for j in range(2**n):
            if (statelist[-1][j]):
                pos = j
                break
        poslist.append(pos)
    # Sort the states
    sortperm = sp.array(poslist).argsort()
    bitstring = [ bitstring[i] for i in sortperm ]

    return bitstring


def GetProbabilities(n, Psi):
    """
    Get probabilites by doing mod-squared of @Psi.
    """
    density = sp.zeros(2**n)
    for i in range(2**n): 
        density[i] = sp.vdot(Psi[i], Psi[i]).real
    return density


def SortStates(n, bitstring, density):
    """ 
    Sort the probabilities of Psi in descending order with labels. 
    """
    sortperm = density.argsort()[::-1] # Flip order of argsort (descending)
    bitstring = [ bitstring[i] for i in sortperm ]
    density = [ density[i] for i in sortperm ]
    return [bitstring, density]
