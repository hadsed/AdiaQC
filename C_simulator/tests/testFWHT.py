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

            
def FWHTI2( vec, n ):
    iter = 0

    while( iter < n ):
        temp = 2**(iter)
        i = 0
        while( i < 2**n ):
            veci = vec[i]
            vec[i] += vec[i + temp]
            vec[i + temp] = veci - vec[i + temp]
            if( (i+1)//temp % 2 == 1 ):
                i += temp + 1
            else:
                i += 1
        iter += 1

            
def FWHTIP( vec, n, np, id, ln ):
    iter = ln
    half = np/2;

    FWHTI2( vec, ln ) #perform FWHT on local array

    while( iter < n ):
        #calculate where the needed data resides, src
        temp = 2**(iter-ln)
        src = (id + temp) % np
        if id < src:
            send( id, src, vec )
            recv( src, vec2 )
        else:
            recv( src, vec2 )
            send( id, src, vec )
        #swap( vec, vec2 )

        i = 0
        while( i < 2**ln ):
            if id < src:
                vec[i] += vec2[i]
            else:
                vec[i] = vec2[i] - vec[i]
            i += 1
        iter += 1
        

if __name__ == "__main__":
    vec = np.array([1,0,1,0,0,1,1,0])
    vec2 = np.array([1,0,1,0,0,1,1,0])
    #vecR = np.array([1,0,1,0,0,1,1,0])
    #vec2 = np.array([1,0,1,0,0,1,1,0])
    #vec = [1,0,1,1]
    
    FWHTI2( vec2, 3 )
    print vec2
    #FWHTI( vec, 3 )

    # for n = n1 + n2 := 3 = 1 + 2
    #FWHTI( vec[0:2], 1 ) #proc 0
    #FWHTI( vec[2:4], 1 ) #proc 0
    #FWHTI( vec[4:6], 1 ) #proc 1
    #FWHTI( vec[6:8], 1 ) #proc 1
    '''
    temp = np.copy(vec)
    # L_{2}^{4*2}
    for i in xrange(0, 4):
        for j in xrange(0, 2):
            #temp[j*4 + i] = vec[i*2 + j]
            temp[j*4 + i] = i*2 + j

    #FWHTI( temp[:4], 2 ) #proc 0
    #FWHTI( temp[4:], 2 ) #proc 1
    for i in xrange(0, 2):
        for j in xrange(0, 4):
            #vec[j*2 + i] = temp[i*4 + j]
            vec[j*2 + i] = i*4 + j

    #print np.array([1, 1, 0, 1, 0, 0, 1, 0])
    print "L_{2}^{4*2} = ", temp
    print "L_{4}^{2*4} = ", vec
    #print vec2
    #print vecR
    #print temp
    #print (1/8.0)*np.array(vec)
    #print np.array(vec)
    '''
