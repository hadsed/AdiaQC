import numpy as np


def _fwht( vec, start, end ):
    iters = (end - start)/2;
    if iters >= 1:
        for i in xrange(start, start+iters):
            temp = vec[i]
            vec[i] += vec[i+iters]
            vec[i+iters] = -vec[i+iters] + temp

        _fwht( vec, start, start+iters ) #upper
        _fwht( vec, start+iters, end ) #lower


def FWHTI( vec, n ):
    iter = 1

    while( iter <= n ):
        temp = 2**(n-iter)
        i = 0
        while( i < 2**n ):
        #for i in xrange( 0, 2**n ):
            #if( i//temp % 2 == 1 ):
            #    vec[i] = vec[i - temp] - 2*vec[i]
            #else:
            #    vec[i] += vec[i + temp]

            veci = vec[i]
            vec[i] += vec[i + temp]
            vec[i + temp] = veci - vec[i + temp]
            if( (i+1)//temp % 2 == 1 ):
                i += temp + 1
            else:
                i += 1
        iter += 1

            
def FWHT( vec, n ):
    _fwht( vec, 0, 2**n )


if __name__ == "__main__":
    vec = [1,0,1,0,0,1,1,0]

    FWHTI( vec, 3 )
    FWHTI( vec, 3 )

    print (1/8.0)*np.array(vec)
