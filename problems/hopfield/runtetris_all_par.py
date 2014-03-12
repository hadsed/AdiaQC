'''
Run all learning rules with params from params.py in
parallel with X number of processors.
'''

import subprocess, os, time

n_instances = 1000
processes = set()
max_processes = 4

# Get into the right directory
os.chdir('../../')

# Run Hebb rule instances
command = "python2 run.py -p hopfield/tetris_hebbrule"
for instance in range(n_instances):
    print "HopfieldHebb" + str(instance)
    processes.add(subprocess.Popen([command], shell=True))
    if len(processes) >= max_processes:
        os.wait()
        processes.difference_update(
            p for p in processes.copy() if p.poll() is not None)

# Run Storkey rule instances
command = "python2 run.py -p hopfield/tetris_storkeyrule"

for instance in range(n_instances):
    print "HopfieldStork" + str(instance)
    processes.add(subprocess.Popen([command], shell=True))
    # processes.add(subprocess.Popen([command, name]))
    if len(processes) >= max_processes:
        os.wait()
        processes.difference_update(
            p for p in processes.copy() if p.poll() is not None)

# Run Projection rule instances
command = "python2 run.py -p hopfield/tetris_projectionrule"

for instance in range(n_instances):
    print "HopfieldProj" + str(instance)
    processes.add(subprocess.Popen([command], shell=True))
    # processes.add(subprocess.Popen([command, name]))
    if len(processes) >= max_processes:
        os.wait()
        processes.difference_update(
            p for p in processes.copy() if p.poll() is not None)

# Make sure we have no stragglers before we exit
for p in processes:
    if p.poll() is None:
        p.wait()
