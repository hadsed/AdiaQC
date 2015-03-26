import json
import numpy as np

# Read from textfile directly to be sure
loaded = np.loadtxt('boixo16.txt')
# Construct Ising matrix
J = np.zeros((16,16))
for i,j,val in loaded:
    J[i-1,j-1] = val
# get diagonal, apply bias
alpha = np.diag(J).copy()
alpha[:8] *= 0.44
beta = []
# go through upper triangle, no diags
for row in xrange(16):
    for col in xrange(row+1,16):
        beta.append(-1.0*J[row,col])
# this never changes really
delta = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
                  1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
# output this
jout = {
    "scalars": {
        "nq" : 16,
        "lrgs" : 10,
        "t" : 500.0,
        "dt" : 0.01
    },
    "coefficients" : {
        "alpha" : list(alpha),
        "beta" : list(beta),
        "delta" : list(delta)
    }
}

with open('boixo16.json', 'w') as outfile:
    json.dump(jout, outfile)
