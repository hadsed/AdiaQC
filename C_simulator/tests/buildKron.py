import numpy as np
import sys


def f( i, j, k, n ):
    testi = 2**(n - i - 1)
    testj = 2**(n - j - 1)

    dzi = 1 if (k//testi) % 2 == 0 else -1
    dzj = 1 if (k//testj) % 2 == 0 else -1
    
    return dzi * dzj


def makeKron(n, i, j):
    testi = 2**(n - i - 1)
    testj = 2**(n - j - 1)
    dzi = dzj = -1
    for k in xrange(0, 2**n):
        dzi = (-1)*dzi if( k % testi == 0) else dzi

        dzj = (-1)*dzj if( k % testj == 0) else dzj
        print dzi*dzj


def makeHam(n, al, be, de):
    hz = [0 for i in xrange(0, 2**n)]
    hzz = [0 for i in xrange(0, 2**n)]
    hhxh = [0 for i in xrange(0, 2**n)]

    for i in xrange(0, n):
        testi = 2**(n - i - 1)
        #dzi = -1
        #for k in xrange(0, 2**n):
        #    dzi = (-1)*dzi if( k % testi == 0) else dzi
        #    hz[k] += al[i]*dzi
        #    hhxh[k] += de[i]*dzi

            
        flag = True
        for j in xrange(i, n):
            testj = 2**(n - j - 1)
            dzi = dzj = -1
            for k in xrange(0, 2**n):
                dzi = (-1)*dzi if( k % testi == 0) else dzi
                if flag:
                    hz[k] += al[i]*dzi
                    hhxh[k] += de[i]*dzi
                else:
                    dzj = (-1)*dzj if( k % testj == 0) else dzj
                    hzz[k] += be[i][j]*dzi*dzj
            flag = False


    return hz, hzz, hhxh



def main1(n):
    n = int(sys.argv[1])

    al = [1 for i in xrange(n)]
    de = [1 for i in xrange(n)]
    be = [[1 for j in xrange(n)] for i in xrange(n)]

    hz, hzz, hhxh = makeHam(n, al, be, de)

    print "hz =", np.array(hz)
    print
    print "hzz =", np.array(hzz)
    print
    print "hhxh =", np.array(hhxh)
    print


def main2():
    n = int(sys.argv[1])
    i = int(sys.argv[2])
    j = int(sys.argv[3])

    h = [0 for k in xrange(0, 2**n)]

    for k in xrange(0, 2**n):
        h[k] = f(i, j, k, n)

    print np.array(h)

if __name__ == "__main__":
    main2()
