"""
Run lr_all.py using subprocesses.
"""

import subprocess, os, time

n_instances = 5488  # 32 + 496 + 4960
processes = set()
max_processes = 4

# Get into the right dir (we're in problems)
os.chdir('../')

# Run Hebb rule instances
for instance in range(n_instances):
    command = "python2 run.py -p lr_all -i "+str(instance)+" -k hebb"
    print "HopfieldHebb" + str(instance)
    processes.add(subprocess.Popen([command], shell=True))
    if len(processes) >= max_processes:
        os.wait()
        processes.difference_update(
            p for p in processes.copy() if p.poll() is not None)

# Run Storkey rule instances
for instance in range(n_instances):
    command = "python2 run.py -p lr_all -i "+str(instance)+" -k stork"
    print "HopfieldStork" + str(instance)
    processes.add(subprocess.Popen([command], shell=True))
    if len(processes) >= max_processes:
        os.wait()
        processes.difference_update(
            p for p in processes.copy() if p.poll() is not None)

# Run Projection rule instances
for instance in range(n_instances):
    command = "python2 run.py -p lr_all -i "+str(instance)+" -k proj"
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
