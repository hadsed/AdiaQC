"""

"""

import subprocess, os, time, optparse, shutil
import numpy as np
import json

# Read instance directories from file
# loop over instances
## pass instance number, which equals the line number in the file being read
## (yes it's read twice)
## also pass gamma value
## grab all relevant problem parameters from json file
## run the instance
# vary gamma until success
# start with coarse-grained grid, keep interpolating up to some level
# keep a running count of successes

def spins2bitstr(vec):
    """ Return a converted bitstring from @vec, a list of spins. """
    return ''.join([ '0' if k == 1 else '1' for k in vec ])

def success(bin, ans, thresh, dpath, filename):
    """
    Determine whether a particular instance of a problem was successful.
    Returns True if successful (i.e. if the answer state is most likely),
    False otherwise.
    """
    # Get into the right directory
    os.chdir(dpath)
    # Get probability data
    if binary:
        probs = np.load(filename)
    else:
        probs = np.loadtxt(filename)
    # Read in proper state labeling
    bstrs = [ line.rstrip('\n') for line in open('statelabels.txt') ]
    sortidx = np.argsort(probs)[::-1]
    sorted_bstrs = np.array(bstrs)[sortidx]
    success = False
    prob = probs[sortidx][sorted_bstrs == spins2bitstr(ans)]
    if sorted_bstrs[0] == spins2bitstr(ans) and prob > thresh:
        success = True
    os.chdir('../../')
    return success

if __name__=="__main__":
    # Command line options
    parser = optparse.OptionParser("usage: %prog [options] arg1 arg2")
    parser.add_option("-l", "--listfname", dest="listfname", default='',
                      type="string", help="List of paths to failed instances.")
    parser.add_option("-q", "--qubits", dest="qubits", default=0,
                      type="int", help="Number of qubits/neurons.")
    parser.add_option("-b", "--binary", dest="binary", default=1,
                      type="int", 
                      help="Binary input or not (default = 1).")
    parser.add_option("-t", "--threshold", dest="threshold", default=0.7,
                      type="float", 
                      help="Probability of answer state must be above this threshold"+\
                      " in order to be successful. Must be 0 < threshold < 1.")
    (options, args) = parser.parse_args()
    listfname = options.listfname
    qubits = options.qubits
    binary = options.binary
    threshold = options.threshold
    if listfname is None or qubits == 0:
        print("You must specify both flags: qubit number and list filename.")
        sys.exit(0)
    devnull = open('/dev/null', 'w')
    # Create list of gammas to search over
    gamma_grid = [0.1, 0.9, 0.5, 0.3, 0.7,
                  0.01, 0.05, 0.005]
    # Keep count of how many successes we got
    successes = {
        'hebb': [0]*qubits,
        'stork': [0]*qubits,
        'proj': [0]*qubits
        }
    # Get main path (assuming we're in problems/)
    os.chdir('../')
    maindir = os.getcwd()
    datadir = maindir+'/data/hopfield_gamma_tune'
    os.chdir(datadir)
    # Keep the ones that still failed for further analysis
    failed_fname = 'superfailed_instances_n'+str(qubits)+'.txt'
    open(datadir+'/'+failed_fname, 'w').close()
    # Read the paths
    paths = [ line.strip() for line in open(listfname) ]
    # Loop over paths of failed instances
    for ipath, path in enumerate(paths):
        print("\nPath number: "+str(ipath)+"/"+str(len(paths)))
        shutil.copyfile(path+'/problem_outputs.dat', 'gamma_problem_inputs.json')
        # Get some attributes
        properties = json.load(open('gamma_problem_inputs.json'))
        pnum = len(properties['memories'])
        answer = properties['answer']
        lrule = properties['learningRule']
        annealtime = properties['annealTime']
        # Search for a successful gamma
        for gamma in gamma_grid:
            os.chdir(maindir)  # go to main directory
            print("\t\t Gamma = "+str(gamma))
            subprocess.call('python2 run.py -p hopfield_gamma_tune -f '+
                            str(gamma)+' -q '+str(qubits), shell=True,
                            stdout=devnull, stderr=devnull)
            os.chdir(datadir)
            probfname = 'probsT'+str(annealtime)+'.dat'
            if binary:
                probfname += '.npy'
            if success(binary, answer, threshold, datadir, probfname):
                # increase some counters here
                successes[lrule][pnum-1] += 1
                print("\n##############\n"+
                      "   SUCCESS\n"+
                      "##############\n")
                break
            if gamma == gamma_grid[-1]:
                with open(datadir+'/'+failed_fname, 'a') as outfile:
                    outfile.write(path+'\n')
                
        print(successes)
        # Sum up successes so far
        success_sum = 0
        for k,v in successes.items():
            success_sum += np.sum(v)
        print(str(success_sum)+'/'+str(len(paths)))
        # Get back to where we started
        os.chdir(maindir+'/problems')

    # Output our results
    print(successes)
    with open('gamma_tune_results.txt', 'w') as jfile:
        json.dump(successes, jfile)

