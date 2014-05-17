'''
File: failure_analysis.py
Author: Hadayat Seddiqi
Date: 4.28.14
Description: Take failed instances from Hopfield experiments #1
             and show some things about it.
'''

import os, optparse, sys
import numpy as np
import json
import pylab as pl


def spins2bitstr(vec):
    """ Return a converted bitstring from @vec, a list of spins. """
    return ''.join([ '0' if k == 1 else '1' for k in vec ])

def hamdist(a,b):
    """ Calculate Hamming distance. """
    return np.sum(abs(np.array(a)-np.array(b))/2.0)

def lrule2int(lrule):
    """
    Convert learning rule string to int.
    """
    if lrule == 'hebb':
        return 0
    if lrule == 'stork':
        return 1
    if lrule == 'proj':
        return 2

def success1(bstrs, probs, answer, instate):
    """
    Success defined by looking at the top two states
    and making sure if the top one isn't the answer,
    then it is the input state and the second top
    answer is the true answer, otherwise failure.
    """
    sortidx = np.argsort(probs)[::-1]
    sorted_bstrs = np.array(bstrs)[sortidx]
    success = False
    prob = 0.0
    if sorted_bstrs[0] == spins2bitstr(answer):
        success = True
        prob = probs[sortidx][0]
    elif ((sorted_bstrs[1] == spins2bitstr(answer)) and 
          (sorted_bstrs[0] == spins2bitstr(instate))):
        success = True
        prob = probs[sortidx][1]
    return success, prob

def success2(bstrs, probs, answer, instate):
    """
    Success only if the top state is the answer.
    """
    sortidx = np.argsort(probs)[::-1]
    sorted_bstrs = np.array(bstrs)[sortidx]
    success = False
    prob = probs[sortidx][sorted_bstrs == spins2bitstr(answer)]
    if sorted_bstrs[0] == spins2bitstr(answer):
        success = True
    return success, prob

def analysis(filename, binary):
    """
    Grab the right data when inside a problem instance directory.
    """
    # Get probability and eigspec data
    if binary:
        try:
            probs = np.load(filename)
            eigs = np.load('eigenspectrum.dat.npy')
        except (IOError):
            print ("IOError: No file "+filename+". \nYou probably forgot "+\
                "to specify the appropriate command-line argument: -b 0.")
            print ("Directory:")
            print (os.getcwd())
            sys.exit(0)
    else:
        try:
            probs = np.loadtxt(filename)
            eigs = np.loadtxt('eigenspectrum.dat')
        except (IOError):
            print ("IOError: No file "+filename+". \nYou probably forgot "+\
                "to specify the appropriate command-line argument: -b 0.")
            print ("Directory:")
            print (os.getcwd())
            sys.exit(0)

    # probstr = [ '%1.16E' % e for e in probs.ravel() ]
    bstrs = [ line.rstrip('\n') for line in open('statelabels.txt') ]

    # Organize the eigenstuff
    times = eigs[:,0]/eigs[-1,0]
    eigs = eigs[:,1:]
    mgIdx = np.argmin(eigs[:,1]-eigs[:,0])
    mingap = eigs[mgIdx,1]-eigs[mgIdx,0]
    mgt = times[mgIdx]

    # Get problem outputs
    props = json.load(open('problem_outputs.dat'))
    lrule = props['learningRule']
    mems = props['memories']
    instate = props['inputState']
    answer = props['answer']

    # Calculate average Hamming distance
    avgdist = np.average([ hamdist(m, instate) for m in mems ])
    
    # Calculate success
    success, sprob = success2(bstrs, probs, answer, instate)

    return success, avgdist, len(mems), sprob, lrule, times, eigs, mingap, mgt

def getdata(qubits, binary, probsfname):
    """
    Walk through all directories and get relevant data, outputting to files. Calls
    analysis().
    """

    # Initialize some variables we'll want to look at
    csuccess = {'hebb': [0]*qubits,
                'stork': [0]*qubits,
                'proj': [0]*qubits}
    cfailure = {'hebb': [0]*qubits,
                'stork': [0]*qubits,
                'proj': [0]*qubits}
    data = []
    faildata = []
    # Loop through all data directories
    for root, dirs, files in os.walk('.'):
        if root == '.':
            continue
        # If we are in a dir with no children..
        if dirs == [] and (int(root.split('n')[1].split('p')[0]) == qubits):
            os.chdir(root)
            success, dist, kmems, prob, lrule, times, eigs, mingap, mgt = \
                analysis(probsfname, binary)
            if success:
                data.append([ dist, prob, lrule2int(lrule), mingap, mgt, kmems ])
            else:
                faildata.append([ dist, prob, lrule2int(lrule), mingap, mgt, kmems])

            os.chdir('../../')

    np.save('successdata_n'+str(qubits)+'.dat', data)
    np.save('failuredata_n'+str(qubits)+'.dat', faildata)

# Command line options
if __name__=="__main__":
    parser = optparse.OptionParser("usage: %prog [options] arg1 arg2")
    parser.add_option("-f", "--fname", dest="probsfname", default=None,
                      type="string", 
                      help="Filename for output probabilities.")
    parser.add_option("-g", "--getdata", dest="getdata", default=0,
                      type="int", 
                      help="Do we need to get the success/failure data.")
    parser.add_option("-b", "--binary", dest="binary", default=1,
                      type="int", 
                      help="Binary input or not (default = 1).")
    parser.add_option("-q", "--qubits", dest="qubits", default=None,
                      type="int", 
                      help="Number of qubits/neurons.")
    parser.add_option("-s", "--separate", dest="separate", default=0,
                      type="int", 
                      help="Separate by patterns or learning rules (or both):"+
                      "0: no separation, 1: patterns, 2: learning rules, 3: both")
    parser.add_option("-p", "--patterns", dest="patterns", default=None,
                      type="int", 
                      help=("Number of patterns (default ignores this and "+
                            "does all)."))
    (options, args) = parser.parse_args()
    probsfname = options.probsfname
    getdataopt = options.getdata
    binary = options.binary
    qubits = options.qubits
    patterns = options.patterns
    separate = options.separate

    # Grab the data
    if getdataopt:
        getdata(qubits, binary, probsfname)
        print ("Data collected and outputted to [success,failure]data_n"+
               "[qubits].dat.npy.")
        print ("Exiting. To run plots, run this script again with -g 0.")
        sys.exit(0)
    else:
        successdat = np.load('successdata_n'+str(qubits)+'.dat.npy')
        failuredat = np.load('failuredata_n'+str(qubits)+'.dat.npy')

    # Data is in this format:
    # [ Hamming dist., answer prob., learning rule, mingap, mingap time, kmems]
    # [ d, p, l, m, mt, km ]

    if separate == 1:
        # Build data list based on number of patterns
        sdat = sorted(successdat, key=lambda x: int(x[5]))
        fdat = sorted(failuredat, key=lambda x: int(x[5]))
        zsdat = zip(*sdat)
        zfdat = zip(*fdat)
        sdat = np.matrix(sdat)
        fdat = np.matrix(fdat)
        # Make them into nice paired bounds
        sbounds = [ (zsdat[5].index(k), zsdat[5].index(k+1))
                    for k in range(2,qubits) ]
        fbounds = [ (zfdat[5].index(k), zfdat[5].index(k+1)) 
                    for k in range(2,qubits) ]
        sbounds.append((zsdat[5].index(qubits), len(zsdat[5])))
        fbounds.append((zfdat[5].index(qubits), len(zfdat[5])))
        # Loop over pattern numbers
        for pidx in range(len(sbounds)):
            # Array bounds for particular number of patterns
            starts, ends = sbounds[pidx]
            startf, endf = fbounds[pidx]
            # Bar graphs
            width = 0.1
            htype = 'bar'
            fontsize = 16
            fig = pl.figure(figsize=(8,8))
            fig.suptitle('Failure Analysis: '+str(qubits)+' qubits, '+
                         str(pidx+2)+' patterns',
                         fontsize=fontsize, fontweight='bold')

            # Hamming distances
            ax = fig.add_subplot(4,1,1)
            recIdx = 0
            nbins = 15
            ax.set_title('Average Hamming distance')
            rmax = max(sdat[starts:ends][:,recIdx].max(), 
                       fdat[startf:endf][:,recIdx].max())
            rmin = min(sdat[starts:ends][:,recIdx].min(), 
                       fdat[startf:endf][:,recIdx].min())
            ax.hist(sdat[starts:ends][:,recIdx], nbins, normed=0, alpha=0.25, 
                    facecolor='blue', histtype=htype, label="Successes", 
                    range=(rmin,rmax))
            ax.hist(fdat[startf:endf][:,recIdx], nbins, normed=0, alpha=0.25, 
                    facecolor='red', histtype=htype, label="Failures", 
                    range=(rmin,rmax))
            pl.grid()
            pl.legend(prop={'size':8})

            # Answer probabilities
            ax = fig.add_subplot(4,1,2)
            recIdx = 1
            nbins = 35
            ax.set_title('Answer probability (log)')
            rmax = max(sdat[starts:ends][:,recIdx].max(), 
                       fdat[startf:endf][:,recIdx].max())
            rmin = min(sdat[starts:ends][:,recIdx].min(), 
                       fdat[startf:endf][:,recIdx].min())
            ax.hist(sdat[starts:ends][:,recIdx], nbins, normed=0, alpha=0.25, 
                    facecolor='blue', histtype=htype, label="Successes", 
                    log=True, range=(rmin,rmax))
            ax.hist(fdat[startf:endf][:,recIdx], nbins, normed=0, alpha=0.25, 
                    facecolor='red', histtype=htype, label="Failures", 
                    log=True, range=(rmin,rmax))
            pl.grid()

            # Mingap
            ax = fig.add_subplot(4,1,3)
            recIdx = 3
            nbins = 15
            ax.set_title('Minimum Spectral Gap')
            rmax = max(sdat[starts:ends][:,recIdx].max(), 
                       fdat[startf:endf][:,recIdx].max())
            rmin = min(sdat[starts:ends][:,recIdx].min(), 
                       fdat[startf:endf][:,recIdx].min())
            ax.hist(sdat[starts:ends][:,recIdx], nbins, normed=0, alpha=0.25, 
                    facecolor='blue', histtype=htype, label="Successes", 
                    log=False, range=(rmin,rmax))
            ax.hist(fdat[startf:endf][:,recIdx], nbins, normed=0, alpha=0.25, 
                    facecolor='red', histtype=htype, label="Failures", 
                    log=False, range=(rmin,rmax))
            pl.grid()

            # Mingap time
            ax = fig.add_subplot(4,1,4)
            recIdx = 4
            nbins = 15
            ax.set_title('Time of Minimum Spectral Gap')
            rmax = max(sdat[starts:ends][:,recIdx].max(), 
                       fdat[startf:endf][:,recIdx].max())
            rmin = min(sdat[starts:ends][:,recIdx].min(), 
                       fdat[startf:endf][:,recIdx].min())
            ax.hist(sdat[starts:ends][:,recIdx], nbins, normed=0, alpha=0.25, 
                    facecolor='blue', histtype=htype, label="Successes", 
                    log=False, range=(rmin,rmax))
            ax.hist(fdat[startf:endf][:,recIdx], nbins, normed=0, alpha=0.25, 
                    facecolor='red', histtype=htype, label="Failures", 
                    log=False, range=(rmin,rmax))
            pl.grid()

            pl.subplots_adjust(hspace=0.5)
            pl.savefig('failure_analysis_n'+str(qubits)+'p'+str(pidx+1)+'.png')
    elif separate == 2:
        # Build data list based on learning rule
        successdat = sorted(successdat, key=lambda x: int(x[2]))
        failuredat = sorted(failuredat, key=lambda x: int(x[2]))
        sdat = {'hebb': [], 'stork': [], 'proj': []}
        fdat = {'hebb': [], 'stork': [], 'proj': []}
        for row in successdat:
            if row[2] == 0:
                sdat['hebb'].append(row)
                fdat['hebb'].append(row)
            elif row[2] == 1:
                sdat['stork'].append(row)
                fdat['stork'].append(row)
            else:
                sdat['proj'].append(row)
                fdat['proj'].append(row)
        for k,v in sdat.iteritems():
            sdat[k] = np.matrix(sdat[k])
            fdat[k] = np.matrix(fdat[k])
        for lr in ['hebb', 'stork', 'proj']:
            # Bar graphs
            width = 0.1
            htype = 'bar'
            fontsize = 16
            fig = pl.figure(figsize=(8,8))
            fig.suptitle('Failure Analysis: '+str(qubits)+' qubits, '+lr, 
                         fontsize=fontsize, fontweight='bold')

            # Hamming distances
            ax = fig.add_subplot(4,1,1)
            recIdx = 0
            nbins = 15
            ax.set_title('Average Hamming distance')
            rmax = max(sdat[lr][:,recIdx].max(), fdat[lr][:,recIdx].max())
            rmin = min(sdat[lr][:,recIdx].min(), fdat[lr][:,recIdx].min())
            ax.hist(sdat[lr][:,recIdx], nbins, normed=0, alpha=0.25, 
                    facecolor='blue', histtype=htype, label="Successes", 
                    range=(rmin,rmax))
            ax.hist(fdat[lr][:,recIdx], nbins, normed=0, alpha=0.25, 
                    facecolor='red', histtype=htype, label="Failures", 
                    range=(rmin,rmax))
            pl.grid()
            pl.legend(prop={'size':8})

            # Answer probabilities
            ax = fig.add_subplot(4,1,2)
            recIdx = 1
            nbins = 35
            ax.set_title('Answer probability (log)')
            rmax = max(sdat[lr][:,recIdx].max(), fdat[lr][:,recIdx].max())
            rmin = min(sdat[lr][:,recIdx].min(), fdat[lr][:,recIdx].min())
            ax.hist(sdat[lr][:,recIdx], nbins, normed=0, alpha=0.25, 
                    facecolor='blue', histtype=htype, label="Successes", 
                    log=True, range=(rmin,rmax))
            ax.hist(fdat[lr][:,recIdx], nbins, normed=0, alpha=0.25, 
                    facecolor='red', histtype=htype, label="Failures", 
                    log=True, range=(rmin,rmax))
            pl.grid()

            # Mingap
            ax = fig.add_subplot(4,1,3)
            recIdx = 3
            nbins = 15
            ax.set_title('Minimum Spectral Gap')
            rmax = max(sdat[lr][:,recIdx].max(), fdat[lr][:,recIdx].max())
            rmin = min(sdat[lr][:,recIdx].min(), fdat[lr][:,recIdx].min())
            ax.hist(sdat[lr][:,recIdx], nbins, normed=0, alpha=0.25, 
                    facecolor='blue', histtype=htype, label="Successes", 
                    log=False, range=(rmin,rmax))
            ax.hist(fdat[lr][:,recIdx], nbins, normed=0, alpha=0.25, 
                    facecolor='red', histtype=htype, label="Failures", 
                    log=False, range=(rmin,rmax))
            pl.grid()

            # Mingap time
            ax = fig.add_subplot(4,1,4)
            recIdx = 4
            nbins = 15
            ax.set_title('Time of Minimum Spectral Gap')
            rmax = max(sdat[lr][:,recIdx].max(), fdat[lr][:,recIdx].max())
            rmin = min(sdat[lr][:,recIdx].min(), fdat[lr][:,recIdx].min())
            ax.hist(sdat[lr][:,recIdx], nbins, normed=0, alpha=0.25, 
                    facecolor='blue', histtype=htype, label="Successes", 
                    log=False, range=(rmin,rmax))
            ax.hist(fdat[lr][:,recIdx], nbins, normed=0, alpha=0.25, 
                    facecolor='red', histtype=htype, label="Failures", 
                    log=False, range=(rmin,rmax))
            pl.grid()

            pl.subplots_adjust(hspace=0.5)
            pl.savefig('failure_analysis_n'+str(qubits)+'_'+lr+'.png')
    elif separate == 0:
        # Put everything on one graph

        # Bar graphs
        width = 0.1
        htype = 'bar'
        fontsize = 16
        fig = pl.figure(figsize=(8,8))
        fig.suptitle('Failure Analysis: '+str(qubits)+' qubits', 
                     fontsize=fontsize, fontweight='bold')

        # Hamming distances
        ax = fig.add_subplot(4,1,1)
        recIdx = 0
        nbins = 15
        ax.set_title('Average Hamming distance')
        rmax = max(successdat[:,recIdx].max(), failuredat[:,recIdx].max())
        rmin = min(successdat[:,recIdx].min(), failuredat[:,recIdx].min())
        ax.hist(successdat[:,recIdx], nbins, normed=0, alpha=0.25, facecolor='blue', 
                histtype=htype, label="Successes", range=(rmin,rmax))
        ax.hist(failuredat[:,recIdx], nbins, normed=0, alpha=0.25, facecolor='red', 
                histtype=htype, label="Failures", range=(rmin,rmax))
        pl.grid()
        pl.legend(prop={'size':8})

        # Answer probabilities
        ax = fig.add_subplot(4,1,2)
        recIdx = 1
        nbins = 35
        ax.set_title('Answer probability (log)')
        rmax = max(successdat[:,recIdx].max(), failuredat[:,recIdx].max())
        rmin = min(successdat[:,recIdx].min(), failuredat[:,recIdx].min())
        ax.hist(successdat[:,recIdx], nbins, normed=0, alpha=0.25, facecolor='blue', 
                histtype=htype, label="Successes", log=True, range=(rmin,rmax))
        ax.hist(failuredat[:,recIdx], nbins, normed=0, alpha=0.25, facecolor='red', 
                histtype=htype, label="Failures", log=True, range=(rmin,rmax))
        pl.grid()

        # Mingap
        ax = fig.add_subplot(4,1,3)
        recIdx = 3
        nbins = 15
        ax.set_title('Minimum Spectral Gap')
        rmax = max(successdat[:,recIdx].max(), failuredat[:,recIdx].max())
        rmin = min(successdat[:,recIdx].min(), failuredat[:,recIdx].min())
        ax.hist(successdat[:,recIdx], nbins, normed=0, alpha=0.25, facecolor='blue', 
                histtype=htype, label="Successes", log=False, range=(rmin,rmax))
        ax.hist(failuredat[:,recIdx], nbins, normed=0, alpha=0.25, facecolor='red', 
                histtype=htype, label="Failures", log=False, range=(rmin,rmax))
        pl.grid()

        # Mingap time
        ax = fig.add_subplot(4,1,4)
        recIdx = 4
        nbins = 15
        ax.set_title('Time of Minimum Spectral Gap')
        rmax = max(successdat[:,recIdx].max(), failuredat[:,recIdx].max())
        rmin = min(successdat[:,recIdx].min(), failuredat[:,recIdx].min())
        ax.hist(successdat[:,recIdx], nbins, normed=0, alpha=0.25, facecolor='blue', 
                histtype=htype, label="Successes", log=False, range=(rmin,rmax))
        ax.hist(failuredat[:,recIdx], nbins, normed=0, alpha=0.25, facecolor='red', 
                histtype=htype, label="Failures", log=False, range=(rmin,rmax))
        pl.grid()

        pl.subplots_adjust(hspace=0.5)
        pl.savefig('failure_analysis_n'+str(qubits)+'.png')
