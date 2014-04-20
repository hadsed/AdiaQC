"""
Run hopfield_exp1.py using subprocesses.
"""

import subprocess, os, time, optparse

# Command line options
if __name__=="__main__":
    parser = optparse.OptionParser("usage: %prog [options] arg1 arg2")
    parser.add_option("-i", "--instances", dest="instances", default=1,
                      type="int", help="Number of problem instances.")
    parser.add_option("-m", "--maxproc", dest="maxproc", default=1,
                      type="int", help="Number of problem instances.")
    parser.add_option("-q", "--qubits", dest="qubits", default=1,
                      type="int", help="Number of neurons/qubits.")
    parser.add_option("-n", "--nummems", dest="nummems", default=1,
                      type="int", help="Number of memories.")
    # Learning rule mode is:
    # 0 : all
    # 1 : hebb
    # 2 : stork
    # 3 : proj
    # 4 : hebb & stork
    # 5 : hebb & proj
    # 6 : stork & proj
    parser.add_option("-l", "--lmode", dest="lmode", default=0,
                      type="int", help="Learning rule mode.")
    (options, args) = parser.parse_args()
    n_instances = options.instances
    max_processes = options.maxproc
    nummems = options.nummems
    qubits = options.qubits
    lmode = options.lmode

# n_instances = 5488  # 32 + 496 + 4960
processes = set()
nrng = range(n_instances)

# Get into the right dir (we're in problems)
os.chdir('../')

# Python command
python = "python2"
hop = "hopfield_exp1"

# Run Hebb rule instances
if (lmode == 0) or (lmode == 1) or (lmode == 4) or (lmode == 5):
    for instance in nrng:
        time.sleep(1)  # wait for last dir to be created
        command = python+" run.py -p "+hop+" -i "+str(instance)+" -k hebb" + \
            " -q "+str(qubits)+" -f "+str(nummems)
        print "HopfieldHebb" + str(instance)
        processes.add(subprocess.Popen([command], shell=True))
        if len(processes) >= max_processes:
            os.wait()
            processes.difference_update(
                p for p in processes.copy() if p.poll() is not None)

# Run Storkey rule instances
if (lmode == 0) or (lmode == 2) or (lmode == 4) or (lmode == 6):
    for instance in nrng:
        command = python+" run.py -p "+hop+" -i "+str(instance)+" -k stork" + \
            " -q "+str(qubits)+" -f "+str(nummems)
        print "HopfieldStork" + str(instance)
        processes.add(subprocess.Popen([command], shell=True))
        if len(processes) >= max_processes:
            os.wait()
            processes.difference_update(
                p for p in processes.copy() if p.poll() is not None)

# Run Projection rule instances
if (lmode == 0) or (lmode == 3) or (lmode == 5) or (lmode == 6):
    for instance in nrng:
        command = python+" run.py -p "+hop+" -i "+str(instance)+" -k proj" + \
            " -q "+str(qubits)+" -f "+str(nummems)
        print "HopfieldProj" + str(instance)
        processes.add(subprocess.Popen([command], shell=True))
        if len(processes) >= max_processes:
            os.wait()
            processes.difference_update(
                p for p in processes.copy() if p.poll() is not None)

# Make sure we have no stragglers before we exit
for p in processes:
    if p.poll() is None:
        p.wait()
