import numpy as np
import pylab as pl
import json
import itertools
import os

n = 5
p = 5
sims = 990
T = np.arange(0.1, 15., 0.5)
hamMetric = 'mean'
pref = 'n'+str(n)+'p'+str(p)
multiT = True
dirs = []
if multiT:
    postf = 'MultiT/'
    dirs = [ pref+'hebb'+postf, pref+'stork'+postf, pref+'proj'+postf]
else:
    dirs = [ pref+'hebb/', pref+'stork/', pref+'proj/']
    
pthresh = 0.9
fname_pref = '2dplot_'+pref+'_'

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
        weight = prob[0]
#        weight = props['hammingDistance'][hamMetric]
        energy = e1 - e0
        data.append([ time, energy, weight ])
    return data

# Normalize weights
def normalize(data):
    zData = zip(*data)
    vec = zData[2]
    sortedList = sorted(vec)
    hi, lo = sortedList[-1], sortedList[0]
    vecNormed = [ (w-lo)/(hi-lo) for w in vec ]
    del zData[-1]
    zData.append(vecNormed)
    return zip(*zData)

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

    c1 = pl.get_cmap('Greens')
    c2 = pl.get_cmap('Oranges')
    c3 = pl.get_cmap('Blues')
    dataPairs = [ (hebbData, c1), (storkData, c2), (projData, c3) ]

    fname_pref = '2dplot_'+pref+'_t'+str(time)+'_'
    pmsg = '(N = '+str(n)+', P = '+str(p)+', T = '+str(time)+')'
    marker = ''
    msize = 5
    xlim = [0,time]
    ylim = [0,3]
    lwidth = 2.0

    print "Plotting " + pmsg

    def resetFig():
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
            pl.plot(t,e,color=cmap(w), marker=marker,markersize=msize,linewidth=lwidth)
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
            pl.plot(t,e,color=cmap(w), marker=marker,markersize=msize,linewidth=lwidth)
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
            pl.plot(t,e,color=cmap(w), marker=marker,markersize=msize,linewidth=lwidth)
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
            pl.plot(t,e,color=cmap(w), marker=marker,markersize=msize,linewidth=lwidth)
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
    fnames.append(fname_pref + 'all.png')

# Save the filenames for animation
outfile = open("2dplot_list.txt", "wb")
outfile.writelines([ name+'\n' for name in fnames ])
outfile.close()

# Make into animation
os.system("convert -delay 30 @2dplot_list.txt 2dplot_"+pref+"_anim.gif")
