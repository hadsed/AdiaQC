"""
Generate a file full of random instances. Final output
file will be a pickled nested list structure that looks
like:

[ [ p1, p2, ..., pN ] x num_instances ]

where pk = [ 1, -1, -1, ... ].
"""

import os, optparse, sys
import numpy as np
import cPickle

def hamdist(a,b):
    """ Calculate Hamming distance. """
    return sp.sum(abs(sp.array(a)-sp.array(b))/2.0)

# Command line options
if __name__=="__main__":
    parser = optparse.OptionParser("usage: %prog [options] arg1 arg2")
    parser.add_option("-q", "--qubits", dest="qubits", default=None,
                      type="int", 
                      help="Number of qubits/neurons.")
    parser.add_option("-i", "--instances", dest="instances", default=1000,
                      type="int", 
                      help="Number of instances.")
    (options, args) = parser.parse_args()
    qubits = options.qubits
    instances = options.instances

    # Generate random sets
    for i in xrange(instances):
        # Generate set of random patterns
        memset = [ [ 2*sp.random.random_integers(0,1)-1
                   for k in xrange(qubits) ]
                 for m in xrange(qubits) ]
        # Make sure all patterns are at least Hamming distance = 2 away from
        # all other patterns in the set
        for m in xrange(qubits):
            # Loop over other memories
            othermems = range(qubits)
            del othermems[m]
            for k in othermems:
                while hamdist(memset[m], memset[k]) < 2.0:
                    # Flip a random spin
                    rndbit = sp.random.random_integers(0,qubits-1)
                    memories[k][rndbit] *= -1

    # Output to pickled file
    cPickle.dump(memset, open('exp2_memset_n'+str(qubits)+'.dat', 'wb'))
    # Sometimes we just need some closure
    print "Done."
