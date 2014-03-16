from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
import json

n = 4
p = 3
T = np.arange(0.1,15,0.5)
sims = 10
hamMetric = 'mean'
pref = 'n'+str(n)+'p'+str(p)
suf = 'MultiT/'
dirs = [ pref+'hebb'+suf, pref+'stork'+suf, pref+'proj'+suf]

hebbData = []
storkData = []
projData = []

def getData(fpath):
    fpath += '/' # Just incase!
    data = []
    for num in range(0,sims):
        for t in T:
#        for t in [4.1,6.6,8.6,10.1,12.1,13.6]:
            time, e0, e1 = np.loadtxt(fpath+str(num)+'/eigenspectrum'+str(t)+'.dat',
                                      usecols=[0,1,2],
                                      delimiter=' ',
                                      unpack=True)
            bitstr, prob = np.loadtxt(fpath+str(num)+'/probsT'+str(t)+'.dat',
                                      usecols=[0,1],
                                      unpack=True,
                                      dtype=str)
            prob = np.array(prob, dtype=float)
            weight = prob[0]
            energy = e1 - e0
            data.append([ list(energy[0:99]), list(time[0:99]), t, weight ])
    return data

#hebbData = getData(dirs[0])
#storkData = getData(dirs[1])
projData = getData(dirs[2])

# Normalize weights
def normalize(data):
    zData = zip(*data)
    vec = zData[-1]
    sortedList = sorted(vec)
    hi, lo = sortedList[-1], sortedList[0]
    vecNormed = [ (w-lo)/(hi-lo) for w in vec ]
    del zData[-1]
    zData.append(vecNormed)
    return zip(*zData)

#hebbData = normalize(hebbData)
#storkData = normalize(storkData)
projData = normalize(projData)
#hebbData = storkData
hebbData = projData
# Get everything out properly
zHD = zip(*hebbData)
energies = reduce(lambda x,y: x+y, zHD[0])
times = reduce(lambda x,y: x+y, zHD[1])
atimes = list(np.tile(np.array(zHD[2]), len(energies)/len(zHD[2])))
weights = list(np.tile(np.array(zHD[3]), len(energies)/len(zHD[3])))

# x = mingap
# y = time
# z = T
# c = prob.

def ndmesh(*args):
   args = map(np.asarray,args)
   return np.broadcast_arrays(*
       [ x[(slice(None),)+(None,)*i] for i, x in enumerate(args) ])


#X,Y = ndmesh(energies, times)
fig = plt.figure()
ax = fig.gca(projection='3d')
x = np.arange(len(energies))
y = np.arange(len(energies))
scat = ax.scatter(energies, times, atimes, c=weights, s=30, cmap='gist_heat')
fig.colorbar(scat, aspect=40).set_label('Final state probability')
ax.set_xlabel('Min. Gap')
ax.set_ylabel('Simulation Time')
ax.set_zlabel('Final Annealing Time')


#X, Y = np.meshgrid(X, Y)
#R = np.sqrt(x**2 + y**2)
#Z = np.sin(R)
#surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
#        linewidth=0, antialiased=False)
#ax.set_zlim(-1.01, 1.01)

#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

#fig.colorbar(surf, shrink=0.5, aspect=5)

#plt.savefig('3d.png')
plt.show()
