import numpy as np
import pylab as pl
import json
import itertools
import os

n = 5
p = 3
# T = 4.1
T = np.arange(0.1,15,0.5)
sims = 1000
hamMetric = 'mean'
pref = 'n'+str(n)+'p'+str(p)
dirs = [ pref+'hebb/', pref+'stork/', pref+'proj/']
# if we want the multiT dirs
dirs = [ prefix[:-1] + 'MultiT/' for prefix in dirs ]
pthresh = 0.9

hebbData = []
storkData = []
projData = []

#
# Plot the gap vs. T vs. mean or median Hamming distance
#

def getData(fpath, t):
    fpath += '/' # Just incase!
    data = []
    pthreshCount = 0
    for num in range(0,sims):
        time, e0, e1 = np.loadtxt(fpath+str(num)+'/eigenspectrum'+str(t)+'.dat',
                                  usecols=[0,1,2],
                                  delimiter=' ',
                                  unpack=True)
        bitstr, prob = np.loadtxt(fpath+str(num)+'/probsT'+str(t)+'.dat',
                                  usecols=[0,1],
                                  unpack=True,
                                  dtype=str)
        prob = np.array(prob, dtype=float)
        props = json.load(open(fpath+str(num)+'/networkProperties.dat'))

        # Check some conditions:
        # if probability is high enough
        # if, given input state is in mems, did it get right answer
        # if, not given that, then whatever (can ignore this)
#        if prob[0] < pthresh:
#            pthreshCount = pthreshCount + 1
#            continue
#        if not bitstr[0] == props['input']:
#            continue
#        if (props['input'] in props['memories']):
#            continue
#        weight = prob[0]
#        weight = props['hammingDistance'][hamMetric]
        energy = e1 - e0
        # data.append([ time, energy, weight ])
        data.append(energy)
    return data

# Store the filenames for animation
fnames = []

# Loop through all annealing times possible
for time in T:
    print "Starting time " + str(time) +".\n"
    hebbData = getData(dirs[0], time)
    storkData = getData(dirs[1], time)
    projData = getData(dirs[2], time)

    bins = 20
    xlabel = 'Minimum gap (E1 - E0)'
    ylabel = 'Population'
    xlim = [0,10]
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

    #pl.rcParams.update({'font.size': 10})
    #pl.show()
    fn = 'hist2_n'+str(n)+'p'+str(p)+'T'+str(time)+'.png'
    fnames.append(fn)
    pl.savefig(fn)

# Save the filenames for animation
outfile = open("hist2_list.txt", "wb")
outfile.writelines([ name+'\n' for name in fnames ])
outfile.close()

# Make into animation
os.system("convert -delay 30 @hist2_list.txt hist2_"+pref+"_anim.gif")
