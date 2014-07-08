"""
Run hopfield_gamma.py using subprocesses.
"""

import subprocess, os, time, optparse
import numpy as np

# Command line options
if __name__=="__main__":
    parser = optparse.OptionParser("usage: %prog [options] arg1 arg2")
    # parser.add_option("-i", "--instances", dest="instances", default=1,
    #                   type="int", help="Number of problem instances.")
    # parser.add_option("-q", "--qubits", dest="qubits", default=1,
    #                   type="int", help="Number of neurons/qubits.")
    # parser.add_option("-n", "--nummems", dest="nummems", default=1,
    #                   type="int", help="Number of memories.")
    # Learning rule mode is:
    # 0 : all
    # 1 : hebb
    # 2 : stork
    # 3 : proj
    # 4 : hebb & stork
    # 5 : hebb & proj
    # 6 : stork & proj
    # parser.add_option("-l", "--lmode", dest="lmode", default=0,
    #                   type="int", help="Learning rule mode.")
    parser.add_option("-m", "--maxproc", dest="maxproc", default=1,
                      type="int", help="Number of processes to run at a time.")
    (options, args) = parser.parse_args()
    max_processes = options.maxproc
    # n_instances = options.instances
    # nummems = options.nummems
    # qubits = options.qubits
    # lmode = options.lmode

processes = set()
bias_rng = np.arange(0.0, 1.05, 0.01)

# Get into the right dir (we're in problems)
os.chdir('../')

# Python command
python = "python2"
hop = "hopfield_gamma"

# Run Hebb rule instances
for gamma in bias_rng:
    time.sleep(1)  # wait for last dir to be created
    command = python+" run.py -p "+hop+" -k hebb -f "+str(gamma)
    print "Hopfield Hebb, max(G) = "+str(bias_rng[-1])+", G = "+str(gamma)
    processes.add(subprocess.Popen([command], shell=True))
    if len(processes) >= max_processes:
        os.wait()
        processes.difference_update(
            p for p in processes.copy() if p.poll() is not None)

# Run Storkey rule instances
for gamma in bias_rng:
    time.sleep(1)  # wait for last dir to be created
    command = python+" run.py -p "+hop+" -k stork -f "+str(gamma)
    print "Hopfield Stork, max(G) = "+str(bias_rng[-1])+", G = "+str(gamma)
    processes.add(subprocess.Popen([command], shell=True))
    if len(processes) >= max_processes:
        os.wait()
        processes.difference_update(
            p for p in processes.copy() if p.poll() is not None)

# Run Projection rule instances
for gamma in bias_rng:
    time.sleep(1)  # wait for last dir to be created
    command = python+" run.py -p "+hop+" -k proj -f "+str(gamma)
    print "Hopfield Proj, max(G) = "+str(bias_rng[-1])+", G = "+str(gamma)
    processes.add(subprocess.Popen([command], shell=True))
    if len(processes) >= max_processes:
        os.wait()
        processes.difference_update(
            p for p in processes.copy() if p.poll() is not None)

# Make sure we have no stragglers before we exit
for p in processes:
    if p.poll() is None:
        p.wait()
