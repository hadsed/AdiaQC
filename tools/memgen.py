"""
Generate a list of Hopfield memory matrices to
iterate through for large instances.
"""

import itertools
import cPickle

def kbits(n, k):
    " kbits(n: length of bitstring, k: number of ones) "
    result = []
    for bits in itertools.combinations(range(n), k):
        s = ['0'] * n
        for bit in bits:
            s[bit] = '1'
        result.append(''.join(s))
    return result

# Define some stuff
nqubits = 5

# Every possible bitstring
pbits = []
for i in range(0,nqubits+1):
    pbits += kbits(nqubits, i)

# Now generate all possible combinations of those
# cmb = [ k for k in itertools.chain(itertools.combinations(pbits, 1), 
#                                    itertools.combinations(pbits, 2), 
#                                    itertools.combinations(pbits, 3)) ]
cmb1 = [ k for k in itertools.combinations(pbits, 1) ]
cmb2 = [ k for k in itertools.combinations(pbits, 2) ]
cmb3 = [ k for k in itertools.combinations(pbits, 3) ]

cPickle.dump(cmb1, open('n'+str(nqubits)+'p1mems.dat', 'wb'))
cPickle.dump(cmb2, open('n'+str(nqubits)+'p2mems.dat', 'wb'))
cPickle.dump(cmb3, open('n'+str(nqubits)+'p3mems.dat', 'wb'))

incmb1 = cPickle.load(open('n'+str(nqubits)+'p1mems.dat', 'rb'))
incmb2 = cPickle.load(open('n'+str(nqubits)+'p2mems.dat', 'rb'))
incmb3 = cPickle.load(open('n'+str(nqubits)+'p3mems.dat', 'rb'))

print "Generated.."
print "p1 len: " + str(len(cmb1))
print "p2 len: " + str(len(cmb2))
print "p3 len: " + str(len(cmb3))
print "From file.."
print "p1 len: " + str(len(incmb1))
print "p2 len: " + str(len(incmb2))
print "p3 len: " + str(len(incmb3))
print "Done."
