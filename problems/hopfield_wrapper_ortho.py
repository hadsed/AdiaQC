'''
Experiment with fully orthogonal patterns
'''

import subprocess, os, time

# Subprocess stuff
processes = set()
max_processes = 6

# Get into the right dir (we're in problems)
os.chdir('../')

# Number of (orthogonal) memories
instances = range(1,10)

# Start timer
t1 = time.time()

# hebb
for inst in instances:
    command = "python2 run.py -p hopfield_ortho -i "+str(inst)+" -k hebb"
    processes.add(subprocess.Popen([command], shell=True))
    # block
    if len(processes) >= max_processes:
        os.wait()
        processes.difference_update(
            p for p in processes.copy() if p.poll() is not None)

# stork
for inst in instances:
    command = "python2 run.py -p hopfield_ortho -i "+str(inst)+" -k stork"
    processes.add(subprocess.Popen([command], shell=True))
    # block
    if len(processes) >= max_processes:
        os.wait()
        processes.difference_update(
            p for p in processes.copy() if p.poll() is not None)

# proj
for inst in instances:
    command = "python2 run.py -p hopfield_ortho -i "+str(inst)+" -k proj"
    processes.add(subprocess.Popen([command], shell=True))
    # block
    if len(processes) >= max_processes:
        os.wait()
        processes.difference_update(
            p for p in processes.copy() if p.poll() is not None)

# Make sure we have no stragglers before we exit
for p in processes:
    if p.poll() is None:
        p.wait()

# End timer
t2 = time.time()

# Write time to file
with open('hop_n10_time.txt', 'w') as file:
    file.write('Time elapsed (sec): '+str(t2-t1))

