'''
File: 2dplot.py
Author: Hadayat Seddiqi
Date: 3.16.14
Description: Builds 2D plots of energy spectra vs. annealing time
             of the three learning rules with colors, which represent
             success probabilities.

'''

import os, optparse
import itertools
import numpy as np
import pylab as pl
import json
import matplotlib.colors as mcolors

def getData(fpath, t):
    """
    Get data for the spectral gap and probabilities.
    """
    data = []
    for num in range(0,sims):
        try:
            # Load final state probabilities
            eigmat = np.load(fpath+str(num)+'/eigenspectrum.dat.npy')
            probs = np.load(fpath+str(num)+'/probsT'+str(t)+'.dat.npy')
            # props = json.load(open(fpath+str(num)+'/problem_outputs.dat'))
        except Exception as e:
            print "Error in sim instance: " + str(num)
            print "Learning rule: " + str(lrule)
            print "\t", e
            continue
        # Grab the relevant data
        weight = np.amax(probs)
#        weight = props['hammingDistance'][hamMetric]
        gap = eigmat[:,2] - eigmat[:,1]
        time = eigmat[:,0] / eigmat[:,0][-1]
        data.append([ time, gap, weight ])
    return data

def normalize(data):
    """
    Normalize the weights (could be probabilities or some dist. metric)
    """
    zData = zip(*data)
    vec = zData[2]
    sortedList = sorted(vec)
    hi, lo = sortedList[-1], sortedList[0]
    vecNormed = [ (w-lo)/(hi-lo) for w in vec ]
    del zData[-1]
    zData.append(vecNormed)
    return zip(*zData)

def make_colormap(seq):
    """
    Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)

def resetFig():
    """
    Clears the last figure and redefines the axes correctly.
    """
    pl.clf()
    pl.figure(figsize=(8,6))
    pl.rcParams.update({'font.size': 9.5})
    pl.subplots_adjust(left=0.1, right=1.0, top=0.9, bottom=0.1)
    pl.tick_params(axis='both', which='major', labelsize=9)
    pl.tick_params(axis='both', which='minor', labelsize=9)
    sm1 = pl.cm.ScalarMappable(cmap=c1,norm=pl.normalize(vmin=0.0, vmax=1.0))
    sm1._A = []
    sm2 = pl.cm.ScalarMappable(cmap=c2,norm=pl.normalize(vmin=0.0, vmax=1.0))
    sm2._A = []
    sm3 = pl.cm.ScalarMappable(cmap=c3,norm=pl.normalize(vmin=0.0, vmax=1.0))
    sm3._A = []
    return sm1,sm2,sm3

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

# Hamming distance metric (for weightings)
hamMetric = 'mean'

# Problem prefix
pref = 'n'+str(n)+'p'+str(p)

# Build the paths
dirs = []
if T is None:
    postf = 'MultiT/'
    dirs = [ pref+'hebb'+postf, pref+'stork'+postf, pref+'proj'+postf]
    props = json.load(open(dirs[0]+'/0/problem_outputs.dat'))
    T = props['annealTime']
else:
    dirs = [ pref+'hebb/', pref+'stork/', pref+'proj/']
    T = list(T)

# Filename prefix
fname_pref = '2dplot_'+pref+'_'

# Record filenames for animation
fnames = []

# Loop over times
for time in T:
    hebbData = getData(dirs[0], time)
    storkData = getData(dirs[1], time)
    projData = getData(dirs[2], time)

    hebbData = normalize(hebbData)
    storkData = normalize(storkData)
    projData = normalize(projData)

    # Make the new colormaps
    c = mcolors.ColorConverter().to_rgb
    c1 = make_colormap([c('darkred'), c('red'), 0.33, 
                        c('red'), c('salmon'), 0.66, 
                        c('salmon')])
    c2 = make_colormap([c('navy'), c('blue'), 0.33, 
                        c('blue'), c('skyblue'), 0.66, 
                        c('skyblue')])
    c3 = make_colormap([c('darkgreen'), c('green'), 0.33, 
                        c('green'), c('lime'), 0.66, 
                        c('lime')])
    # c1 = pl.get_cmap('Greens')
    # c2 = pl.get_cmap('Oranges')
    # c3 = pl.get_cmap('Blues')
    dataPairs = [ (hebbData, c1), (storkData, c2), (projData, c3) ]

    fname_pref = '2dplot_'+pref+'_t'+str(time)+'_'
    pmsg = '(N = '+str(n)+', P = '+str(p)+', T = '+str(time)+')'
    marker = ''
    msize = 5
    xlim = [0,1]
    ylim = [0,3]
    lwidth = 1.5

    print "Plotting " + pmsg

    # Format the figure correctly
    sm1,sm2,sm3 = resetFig()
    cbarArgs = {'pad':0.01,'aspect':40, 'fraction':0.09, 'ticks':[0.0,0.5,1.0]}

    # Plot Hebb
    for t, e, w in hebbData:
        pl.plot(t,e,color=c1(w),marker=marker,markersize=msize,linewidth=lwidth)
    pl.xlim(xlim)
    pl.xlabel(u'$Time$')
    pl.ylabel(u'$E_1 - E_0$')
    pl.title('Hebb Rule '+pmsg)
    pl.colorbar(sm1, **cbarArgs).set_label('Hebb probability')
    pl.savefig(fname_pref + 'hebb.png')
    sm1,sm2,sm3 = resetFig()

    # Plot Storkey
    for t, e, w in storkData:
        pl.plot(t,e,color=c2(w), marker=marker,markersize=msize,linewidth=lwidth)
    pl.xlim(xlim)
    pl.xlabel(u'$Time$')
    pl.ylabel(u'$E_1 - E_0$')
    pl.title('Storkey Rule '+pmsg)
    pl.colorbar(sm2, **cbarArgs).set_label('Storkey probability')
    pl.savefig(fname_pref + 'stork.png')
    sm1,sm2,sm3 = resetFig()

    # Plot Proj
    for t, e, w in projData:
        pl.plot(t,e,color=c3(w), marker=marker,markersize=msize,linewidth=lwidth)
    pl.xlim(xlim)
    pl.xlabel(u'$Time$')
    pl.ylabel(u'$E_1 - E_0$')
    pl.title('Projection Rule '+pmsg)
    pl.colorbar(sm3, **cbarArgs).set_label('Projection probability')
    pl.savefig(fname_pref + 'proj.png')
    sm1,sm2,sm3 = resetFig()

    # Plot Hebb and Stork
    for data, cmap in dataPairs[0:2]:
        for t, e, w in data:
            pl.plot(t,e,color=cmap(w), marker=marker, markersize=msize,
                    linewidth=lwidth)
    pl.xlim(xlim)
    pl.ylim(ylim)
    pl.xlabel(u'$Time$')
    pl.ylabel(u'$E_1 - E_0$')
    pl.title('Hebb and Storkey '+pmsg)
    pl.colorbar(sm1, **cbarArgs).set_label('Hebb probability')
    pl.colorbar(sm2, **cbarArgs).set_label('Storkey probability')
    pl.savefig(fname_pref + 'hebbstork.png')
    sm1,sm2,sm3 = resetFig()

    # Plot Hebb and Proj
    for data, cmap in [dataPairs[0],dataPairs[2]]:
        for t, e, w in data:
            pl.plot(t,e,color=cmap(w), marker=marker, markersize=msize,
                    linewidth=lwidth)
    pl.xlim(xlim)
    pl.ylim(ylim)
    pl.xlabel(u'$Time$')
    pl.ylabel(u'$E_1 - E_0$')
    pl.title('Hebb and Projection '+pmsg)
    pl.colorbar(sm1, **cbarArgs).set_label('Hebb probability')
    pl.colorbar(sm3, **cbarArgs).set_label('Projection probability')
    pl.savefig(fname_pref + 'hebbproj.png')
    sm1,sm2,sm3 = resetFig()

    # Plot Stork and Proj
    for data, cmap in dataPairs[1:3]:
        for t, e, w in data:
            pl.plot(t,e,color=cmap(w), marker=marker, markersize=msize,
                    linewidth=lwidth)
    pl.xlim(xlim)
    pl.ylim(ylim)
    pl.xlabel(u'$Time$')
    pl.ylabel(u'$E_1 - E_0$')
    pl.title('Storkey and Projection '+pmsg)
    pl.colorbar(sm2, **cbarArgs).set_label('Storkey probability')
    pl.colorbar(sm3, **cbarArgs).set_label('Projection probability')
    pl.savefig(fname_pref + 'storkproj.png')
    sm1,sm2,sm3 = resetFig()

    # Plot all three
    for data, cmap in dataPairs:
        for t, e, w in data:
            pl.plot(t,e,color=cmap(w), marker=marker, markersize=msize,
                    linewidth=lwidth)
    pl.xlim(xlim)
    pl.ylim(ylim)
    pl.xlabel(u'$Time$')
    pl.ylabel(u'$E_1 - E_0$')
    pl.title('Hebb, Storkey, Projection '+pmsg)
    pl.colorbar(sm1, **cbarArgs).set_label('Hebb probability')
    pl.colorbar(sm2, **cbarArgs).set_label('Storkey probability')
    pl.colorbar(sm3, **cbarArgs).set_label('Projection probability')
    pl.savefig(fname_pref + 'all.png')
    sm1,sm2,sm3 = resetFig()

    # Record filenames
    fnames.append(dest+fname_pref + 'all.png')

# If multiT, put together animated GIF
if anim and len(T) > 1:
    # Save textfile with filenames
    outfile = open("2dplot_list.txt", "wb")
    outfile.writelines([ name+'\n' for name in fnames ])
    outfile.close()
    # Make into animation using imagemagick
    os.system("convert -delay 30 @2dplot_list.txt "+
              dest+"2dplot_"+pref+"_anim.gif")
