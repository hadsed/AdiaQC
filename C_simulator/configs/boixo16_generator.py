import json
import numpy as np

# Read from textfile directly to be sure
loaded = np.loadtxt('boixo16.txt')
# Construct Ising matrix
J = np.zeros((16,16))
# this simulator has a negative somewhere
for i,j,val in loaded:
    J[i-1,j-1] = -1.0*val
# get diagonal, apply bias
alpha = 0.44*np.diag(J)
beta = []
# go through upper triangle, no diags
for row in xrange(16):
    for col in xrange(row+1,16):
        beta.append(J[row,col])
# this never changes really
delta = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
                  1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
# output this
jout = {
    "scalars": {
        "nq" : 16,
        "lrgs" : 10,
        "t" : 400.0,
        "dt" : 0.1
    },
    "coefficients" : {
        "alpha" : list(alpha),
        "beta" : list(beta),
        "delta" : list(delta)
    }
}

with open('boixo16.json', 'w') as outfile:
    json.dump(jout, outfile)
