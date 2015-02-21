import random as r
import numpy as np
import sys

def print_config( comments, nq, lrgs, t, dt, alpha, beta, delta ):    
    print "/*"
    print comments
    print "*/"
    print '{' 
    print '    "scalars" : {'
    print '        "nq" : ' + str(nq) + "," 
    print '        "lrgs" : ' + str(lrgs) + "," 
    print '        "t" : ' + str(t) + "," 
    print '        "dt" : ' + str(dt)
    print '    },'
    print '    "coefficients" : {'
    print '        "alpha" : ' + str(alpha) + ',' 
    print '        "beta" : ' + str(beta) + ',' 
    print '        "delta" : ' + str(delta)
    print '    }'
    print '}' 


def rand_config( nq, alpha, beta, delta ):
    for i in xrange( nq ):
        alpha.append( round( r.random(), 3 ) ) 
        delta.append( 1.0 )

    for i in xrange( nq ):
        for j in xrange( i+1, nq ):
            beta.append( round( r.random(), 3 ) );

    return "Random config file"


def hopf_config( nq, alpha, beta, delta, p ):
    comment = ''

    for i in xrange(0, nq):
        #alpha.append( 2*r.randint(0, 1) - 1 )
        #alpha.append( 1 if i < nq//2 else -1 )
        alpha.append( -1 )
        delta.append( 1.0 )

    comment += "input: " + str(alpha) + "\n"
    
    mem = []
    for i in xrange( 0, p ):
        temp = []
        for j in xrange( 0, nq ):
            #temp.append( 2*r.randint(0,1) - 1 )
            temp.append( 1 )
        mem.append( temp )
    #mem = [ [1]*nq, [1,-1]*(nq/2), [1,1,-1,-1]*(nq/4), [1]*(nq/2) + [-1]*(nq/2-1)+[1] ]
    # hebb rule

    comment += "memories: " + str( np.array(mem) ) + "\n"
    
    memMat = np.matrix(mem).T
    ising_off = np.triu(memMat*memMat.T, 1)/float(nq)
    for i in xrange( 0, nq ):
        for j in xrange( i+1, nq ):
            beta.append( -ising_off[i,j] );
    #ising_diag = np.array(inputstate)

    return comment


if __name__ == "__main__":
    nq = int(sys.argv[1])
    lrgs = 10
    t = 10.0
    dt = 0.1 
    alpha = []
    beta = []
    delta = []

    r.seed( 0 ) 
    #comments = rand_config( nq, alpha, beta, delta )
    p = int( 0.14*nq )
    #comments = hopf_config( nq, alpha, beta, delta, p )
    comments = hopf_config( nq, alpha, beta, delta, 1 )

    print_config( comments, nq, lrgs, t, dt, alpha, beta, delta )
