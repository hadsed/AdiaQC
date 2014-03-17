'''
File: histogram1.py
Author: Hadayat Seddiqi
Date: 3.16.14
Description: Builds success probability vs. # of instances
             histograms.

'''

import os, optparse
import itertools
import numpy as np
import pylab as pl
import json

def getData(fpath, t, lrule):
    """
    Get final state probability data.
    """
    data = []
    for num in range(0,sims):
        try:
            # Load final state probabilities
            probs = np.load(fpath+str(num)+'/probsT'+str(t)+'.dat.npy')
        except Exception as e:
            print "Error in sim instance: " + str(num)
            print "Learning rule: " + str(lrule)
            print "\t", e
            continue
        # Pick out the highest probability--this is almost surely the answer
        data.append(np.amax(probs))
    return data

# Command line options
if __name__=="__main__":
    parser = optparse.OptionParser("usage: %prog [options] arg1 arg2")
    parser.add_option("-n", "--nqubits", dest="nqubits", default=None,
                      type="int", 
                      help="Number of qubits/neurons.")
    parser.add_option("-p", "--patterns", dest="patterns", default=None,
                      type="int", 
                      help="Number of patterns in network.")
    parser.add_option("-t", "--time", dest="time", default=None,
                      type="float", 
                      help="Final annealing time.")
    parser.add_option("-i", "--instances", dest="instances", default=None,
                      type="int", 
                      help="Number of instances.")
    parser.add_option("-a", "--animate", dest="animate", default=1,
                      type="int", 
                      help="Animate or no (1 or 0, only for multiT).")
    parser.add_option("-d", "--destination", dest="destination", default="./",
                      type="string", 
                      help="Path for saving image files.")
    (options, args) = parser.parse_args()

n = options.nqubits
p = options.patterns
T = options.time
sims = options.instances
anim = options.animate
dest = options.destination + '/'  # to be safe

# Construct the right path(s)
pref = 'n'+str(n)+'p'+str(p)
dirs = [ pref+'hebb/', pref+'stork/', pref+'proj/']
dirs = [ prefix[:-1] + 'MultiT/' for prefix in dirs ]

# If annealing time isn't specified, grab it from problem_outputs.dat. This
# is usually done if we did a MultiT simulation.
if T is None:
    props = json.load(open(dirs[0]+'/0/problem_outputs.dat'))
    T = props['annealTime']
else:
    T = list(T)

# Storing filenames for animation
fnames = []

# Loop through all annealing times possible
for time in T:
    print "Starting T = "+str(time)+" for "+pref+".\n"
    hebbData = getData(dirs[0], time, 'hebb')
    storkData = getData(dirs[1], time, 'storkey')
    projData = getData(dirs[2], time, 'projection')

    bins = 20
    xlabel = 'Ground state probability'
    ylabel = 'Population'
    xlim = [0,1]
    ylim = [0,sims]
    labelfsize = 16
    normed = False
    fig = pl.figure(figsize=(8,9))
    fig.suptitle(r'$N_{qubits} = '+str(n)+'$,  $T_{final} = ' + 
                 str(time) + '$,  $P = 3$,  $G = N-n/2N$', 
                 fontsize=14)
    pl.grid(True)

    # Plot the stuff
    pl.subplot(3,1,1)
    pl.title('Hebb rule')
    pl.xlim(xlim)
    pl.ylim(ylim)
    pl.hist(x=hebbData, bins=bins, range=(0,1), normed=normed)

    pl.subplot(3,1,2)
    pl.title('Storkey rule')
    pl.ylabel(ylabel, fontweight='bold', fontsize=labelfsize)
    pl.xlim(xlim)
    pl.ylim(ylim)
    pl.hist(x=storkData, bins=bins, range=(0,1), normed=normed)

    pl.subplot(3,1,3)
    pl.title('Projection rule')
    pl.xlabel(xlabel, fontweight='bold', fontsize=labelfsize)
    pl.xlim(xlim)
    pl.ylim(ylim)
    pl.hist(x=projData, bins=bins, range=(0,1), normed=normed)

    # Save it
    fn = 'hist1_n'+str(n)+'p'+str(p)+'T'+str(time)+'.png'
    fnames.append(fn)
    pl.savefig(dest+fn)

# If multiT, put together animated GIF
if anim and len(T) > 1:
    # Save textfile with filenames
    outfile = open("hist1_list.txt", "wb")
    outfile.writelines([ name+'\n' for name in fnames ])
    outfile.close()
    # Make into animation using imagemagick
    os.system("convert -delay 30 @hist1_list.txt "+
              dest+"hist1_"+pref+"_anim.gif")
