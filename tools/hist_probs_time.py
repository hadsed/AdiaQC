'''
File: hist_probs_time.py
Author: Hadayat Seddiqi
Date: 4.27.14
Description: Plot a series of histograms showing the probability
             to be in a particular state (so 2^nqubits bar graphs)
             over the annealing time.
'''

import os, optparse, sys
import numpy as np
import matplotlib.pyplot as pl

# Command line options
if __name__=="__main__":
    parser = optparse.OptionParser("usage: %prog [options] arg1 arg2")
    parser.add_option("-b", "--binary", dest="binary", default=1,
                      type="int", 
                      help="Whether data file is in binary (1) or not (0).")
    (options, args) = parser.parse_args()
    binary = options.binary

# Get state overlap labels
axis = []
with open('state_overlap_labels.txt') as f:
    axis = f.readlines()

# Get probabilities
data = []
ftimes = []
if binary:
    try:
        for fname in os.listdir('.'):
            if fname[:15] == 'state_overlap_T' and fname[-4:] == '.npy':
                prob = np.load(fname)
                t = fname[15:].split('.txt')[0].split('T')[0]
                ftimes.append(float(t))
                data.append(prob)
    except IOError, e:
        print "IOError: No file eigenspectrum.dat.npy. \nYou probably forgot "+\
            "to specify the appropriate command-line argument: -b 0."
        sys.exit(0)
else:
    for fname in os.listdir('.'):
        if fname[:15] == 'state_overlap_T' and fname[-4:] == '.txt':
            prob = np.loadtxt(fname)
            t = fname[15:].split('.txt')[0].split('T')[0]
            ftimes.append(float(t))
            data.append(prob)

# Fix up the data
data = zip(*data)
sortidx = np.array(ftimes).argsort()
data = [ np.array(rec)[sortidx] for rec in data ]
time = np.linspace(0,1,len(data[0]))

# Visualize as a timeseries plot
pl.title('Six-var QUBO State Probabilities vs. Time')
pl.xlabel('Time')
pl.ylabel('Probability')
pl.ylim([0,1])
pl.xlim([0,1])
for idat in range(len(data)):
    pl.plot(time, data[idat], label=(axis[idat] if idat <= 5 else None))
pl.legend(loc=2)
pl.savefig('stateoverlap_timeseries.png')
pl.clf()

# Now do a histogram
xlabels = []
fnames = []
width = 0.1
# Loop through time
for itime in xrange(time.size):
    ax = pl.figure().add_axes([0.1, 0.1, 0.8, 0.8])
    ax.set_ylim([0,1])
    ax.set_title('Probability to be in a pattern state, t = '+\
                     str(np.around(time[itime], 3)))
    # Loop through pattern data
    for idat in range(len(data)):
        ax.bar(float(idat)/len(data), data[idat][itime], width, 
               color=None, alpha=0.5, align='center')
        xlabels.append(axis[idat])
    # Set up bar labelings
    ax.set_xticks(np.arange(0,1,1./len(data)))
    ax.set_xticklabels(xlabels)
    fn = 'hist_probs_time_t'+str(itime)+'.png'
    fnames.append(fn)
    pl.savefig(fn)
    pl.clf()

# Save textfile with filenames
outfile = open("hist_probs_list.txt", "wb")
outfile.writelines([ name+'\n' for name in fnames ])
outfile.close()
# Make into animation using imagemagick
os.system("convert -delay 10 @hist_probs_list.txt "+\
          "hist_probs_anim.gif")
