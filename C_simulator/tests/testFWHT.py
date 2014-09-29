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


def FWHT( vec, n ):
    _fwht( vec, 0, 2**n )


if __name__ == "__main__":
    vec = [1,0,1,0,0,1,1,0]

    FWHT( vec, 3 )

    print np.array(vec)
