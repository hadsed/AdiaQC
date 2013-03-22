'''

File: output.py
Author: Hadayat Seddiqi
Date: 3.21.13
Description: Output methods (e.g. fidelity, eigenspectrum, etc.)

'''

import os
import scipy as sp
from scipy import linalg

def ConstructEigData(t, vals, n):
    " Construct proper datapoint for eigenspectrum. "
    datapoint = [t]
    for i in range(0, n): datapoint.append(vals[i].real)

    return datapoint

def RecordEigSpec(eigspec, outdir, T):
    " Output eigenspectrum to data file. "
    eigpath = os.path.dirname(os.path.realpath(__file__)) + "/" + \
        outdir + "/eigenspectrum" + str(T) + ".dat"
    sp.savetxt(eigpath, eigspec)

def PlotEigSpec(eigspec, outdir, T):
    " Plot eigenspectrum. "
    import pylab as pl

    eigpath = os.path.dirname(os.path.realpath(__file__)) + "/" + \
        outdir + "/eigenspectrum" + str(T) + ".png"

    # Clear figure
    pl.clf()

    # Get columns of eigspec to plot
    t = [ row[0] for row in eigspec ]
    for i in range(1,len(eigspec[0])): 
        pl.plot(t, [ row[i] for row in eigspec ])

    pl.xlim([0, T])
    pl.xlabel(r'$Time$')
    pl.ylabel(r'$Energy$')
    pl.savefig(eigpath)

def ConstructOverlapData(t, Psi, vec):
    " Construct proper datapoint for overlap. "
    datapoint = [t]

    overlap = sp.vdot(Psi, vec)
    overlapr = sp.vdot(overlap, overlap).real

    datapoint.append(overlapr)

    return datapoint

def RecordOverlap(overlap, outdir, T):
    " Output overlap (in time) to data file. "
    path = os.path.dirname(os.path.realpath(__file__)) + "/" + \
        outdir + "/overlap" + str(T) + ".dat"
    sp.savetxt(path, overlap)

def PlotOverlap(overlap, outdir, T):
    " Plot the overlap. "
    import pylab as pl

    path = os.path.dirname(os.path.realpath(__file__)) + "/" + \
        outdir + "/overlap" + str(T) + ".png"

    # Clear figure
    pl.clf()

    # Get columns of eigspec to plot
    t = [ row[0] for row in overlap ]
    for i in range(1,len(overlap[0])): 
        pl.plot(t, [ row[i] for row in overlap ], marker='.', linestyle='')

    pl.xlim([0, T])
    pl.xlabel(r'$Time$')
    pl.ylabel(r'Overlap $\|\langle \psi \| \phi_0\rangle\|^2$')
    pl.savefig(path)

def RecordFidelity(Psi, vecs, T, outdir):
    " Output the fidelity for multi-T simulations by calculating \
      the overlap between however many eigenstates we want. "

    data = []

    for i in range(0, len(vecs)):
        overlap = sp.vdot(Psi, vecs[i])
        overlapr = sp.vdot(overlap, overlap).real
        data.append([T, overlapr])

    return data

def PlotFidelity(data, outdir):
    " Plot the fidelity. "
    import pylab as pl

    pl.clf()

    path = os.path.dirname(os.path.realpath(__file__)) + "/" + outdir + "/fidelity.png"

    T = [ row[0] for row in data ]
    for i in range(1, len(data[0])):
        pl.plot(T, [ row[i] for row in data ], marker='o', linestyle='')

    pl.xlabel(r'T (anneal time)')
    pl.ylabel(r'Fidelity $\|\langle \psi \| \phi_n\rangle\|^2$')
    pl.savefig(path)
