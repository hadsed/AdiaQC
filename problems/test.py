'''

File: test.py
Author: Hadayat Seddiqi
Date: 3.18.13
Description: The test problem.

'''

import scipy as sp

eigspecflag = 0
outputdir = 'data/'

T = 10.0
dt = 0.1

nQubits = 2
Q = sp.ones((nQubits,nQubits)) - sp.identity(nQubits)
a = sp.ones(nQubits)

print (Q, a)
