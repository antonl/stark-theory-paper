import dask
import dask.array as da
import h5py
import numpy as np
import click
from pathlib import Path
from collections import namedtuple
from concurrent.futures import ProcessPoolExecutor

StarkData = namedtuple('StarkData', ['fieldon', 'fieldoff'])

import matplotlib
matplotlib.use('Agg')

from matplotlib.pyplot import *
from matplotlib.colors import SymLogNorm
from matplotlib.ticker import SymmetricalLogLocator
from cycler import cycler

from pyqcfp.delayed import SimulationCheckpoint
from pyqcfp.runqcfp import QcfpConfig

from pyqcfp.plotting import plot_result, plot_absorption

params = {
    'savefig.dpi': 300,
    'figure.figsize': (12,10),
    'font.size': 12,
    #'xtick.major.size': 3,
    #'xtick.minor.size': 3,
    #'xtick.major.width': 1,
    #'xtick.minor.width': 1,
    'axes.linewidth': 1
}

'''
params.update({
    'ytick.major.size': params['xtick.major.size'],
    'ytick.minor.size': params['xtick.minor.size'],
    'ytick.major.width': params['xtick.major.width'],
    'ytick.minor.width': params['xtick.minor.width'],
    'axes.labelsize': params['font.size'],
})
'''

def normalize(data):
    return data / np.abs(data).max()

def plot_evecs(energies, evecs2, which, path, axlim=None):
    pts_used = energies.shape[0]
    nvecs = energies.shape[1]
    fig = figure()
    ax = fig.add_subplot(111, polar=True)
    cyc = cycler(color=['b', 'r', 'g'])

    #for which in range(nvecs):
    for i,prop in zip(range(nvecs), cyc):
        ax.scatter(3*np.pi/2*evecs2[:pts_used//2, i, which],
               energies[:pts_used//2, which], alpha=0.8, **prop)
    if axlim is None:
        ax.set_rlim(14200, 15100)
    else:
        ax.set_rlim(*axlim)

    ax.set_theta_offset(-np.pi/4)
    ax.set_rlabel_position(np.degrees(3*np.pi/2))
    ax.set_xticklabels(['{:3.0f}%'.format(x) for x in np.linspace(0, 100, 7)] + [''])
    fig.savefig(path, bbox_inches='tight')

def plot_2d(w1, w3, signal, path, invert_w1=False, scale=None, axlim=None):
    signal2 = -1*signal.copy()

    # plot in 1000 of wn
    w1 = w1.copy()/1e3
    w3 = w3.copy()/1e3

    if invert_w1:
        w1 = -w1

    if scale is None: # calculate scale, return it in meta
        scale = np.max(np.abs(signal2))
        signal2 /= scale

    # fiddle parameters to make the plots look good
    linthresh = 0.01
    linscale = 0.1
    nlevels = 100
    norm = SymLogNorm(linthresh=linthresh, linscale=linscale, vmin=-1, vmax=1)
    logloc = SymmetricalLogLocator(linthresh=linthresh, base=2)

    levels = logloc.tick_values(-1, 1)

    fig = figure()
    ax = fig.add_subplot(111, aspect='equal', adjustable='box-forced')

    qset = ax.contourf(w1, w3, signal2, nlevels, norm=norm, cmap='RdBu_r')
    c = ax.contour(w1, w3, signal2, levels=levels, colors='k', alpha=0.4)

    loc = matplotlib.ticker.MaxNLocator(11)
    fmt = matplotlib.ticker.ScalarFormatter(useOffset=False, useMathText=True)

    cb = fig.colorbar(qset, ax=ax, ticks=loc, format=fmt)
    cb.add_lines(c)

    if axlim:
        ypts = xpts = np.array(sorted(axlim))/1e3
        print('Using limits ', xpts, ypts)
        ax.set_xlim(xpts)
        ax.set_ylim(ypts)
    else:
        ypts = xpts = sorted([np.min(w3), np.max(w3)])
        ax.set_xlim(xpts)
        ax.set_ylim(ypts)

    #ax.grid(True)
    ax.add_artist(Line2D(xpts, ypts, linewidth=2, color='k', alpha=0.5,
                              transform=ax.transData))

    mainpath = Path(path)

    fig.savefig(str(mainpath))

def plot_linear(w3, signal, path=None, scale=None, axlim=None, ax=None,
                eigenstates=None):

    w3 = w3.copy()/1e3
    signal2 = signal.copy()

    if ax is None:
        fig = figure()
        ax = fig.add_subplot(111)
    else:
        fig = ax.get_figure()

    if scale is None: # calculate scale, return it in meta
        scale = np.max(np.abs(signal2))
    # scale signal to the passed scale
    signal2 /= scale

    if axlim is None:
        axlim = 15200, 16200

    ax.plot(w3, signal, linewidth=2)
    ax.set_xlim(*axlim)

    # add eigenstate positions
    if eigenstates:
        ax.vlines(eigenstates, -1, 1, linewidth=2, alpha=0.7)

    s = str(path)
    fig.savefig(s)
    return ax, scale

@click.command()
@click.argument('path', type=click.Path(file_okay=False, exists=True))
@click.option('--limits', default=(None, None), type=(float, float))
@click.option('-c','--ncores', default=6)
@click.option('--fudge-factor', default=0.)
def make_figures(path, limits, ncores, fudge_factor):
    # look for pump-probe data file
    path = Path(path)
    pool = ProcessPoolExecutor(max_workers=ncores)
    #rcParams.update(params)

    try:
        ddfile = h5py.File(str(path/'pump-probe.h5'), 'r')
        absfile = h5py.File(str(path/'absorption.h5'), 'r')
    except FileNotFoundError as e:
        print('Datafiles not found in dir {!s}'.format(path))
        return

    # load ref
    ddref = np.array(ddfile['reference']).imag
    shape = (100, *ddref.shape)

    # calculate average pump-probe
    tmp = da.from_array(ddfile['00000/data'], chunks=shape)
    pts_used = tmp.shape[0]
    rdataon = tmp[:pts_used//2].imag.mean(axis=0)
    rdataoff = tmp[pts_used//2:].imag.mean(axis=0)
    dd = StarkData(*dask.compute(rdataon, rdataoff))
    w3, w1 = np.array(ddfile['w3']), np.array(ddfile['w1'])

    # load evecs
    energies = np.array(ddfile['meta/one band energies'])
    nstates = energies.shape[0]
    evecs2 = np.array(ddfile['meta/ge eigenvectors'])**2
    reorgs = np.array(ddfile['meta/reorganization energy matrix'])
    sbcouplingdiag = np.diag(ddfile['meta/sb coupling diagonal'])
    sbcouplingoffdiag = np.diag(ddfile['meta/sb coupling off-diagonal'])
    redfield = np.array(ddfile['meta/redfield relaxation matrix'])
    reorgs = np.diag(reorgs)[1:nstates+1]
    redfield = np.diag(redfield)[1:nstates+1]
    dephasingmat = np.array(ddfile['meta/lifetime dephasing matrix'])[:, ::2] + \
                   1j*np.array(ddfile['meta/lifetime dephasing matrix'])[:, 1::2]
    imagdeph = dephasingmat[1:nstates+1,0].imag

    fixed_energies = energies - reorgs
    fixed_energies2 = energies - reorgs + imagdeph + fudge_factor
    # fudge_factor is the calibrated correction factor that comes from Stokes
    # shift, which moves the location of the monomer

    reorgs_trace = np.diagonal(ddfile['00000/meta/reorganization energy matrix'],
                               axis1=1, axis2=2)[:, 1:nstates+1]
    evecs2_trace = np.array(ddfile['00000/meta/ge eigenvectors'])**2

    energies_trace = np.array(ddfile['00000/meta/one band energies'])
    dephasingmat_trace = np.array(ddfile['00000/meta/lifetime dephasing '
                                         'matrix'])[:,:,::2] + \
                   1j*np.array(ddfile['00000/meta/lifetime dephasing '
                                      'matrix'])[:,:,1::2]
    imagdeph_trace = dephasingmat_trace[:, 1:nstates+1,0].imag
    corr_energies = energies_trace - reorgs_trace + imagdeph_trace + fudge_factor

    # prepare folder for writing things
    figpath = (path / 'figures')
    figpath.mkdir(exist_ok=True)

    with (figpath / 'eigen-energies.info').open('w') as f:
        print('Eigen-energies:', energies, file=f)
        print('GE reorganization energies:', reorgs, file=f)
        print('Reorg\'ed energies:', fixed_energies, file=f)
        print('Reorg\'ed energies + deph + fudge:', fixed_energies2, file=f)
        print(file=f)
        for i in range(evecs2.shape[0]):
            print('Localization of eigenstate {:d}:'.format(i), file=f)
            print(evecs2[i, :], file=f)
            print(file=f)
        print('S-B diagonal couplings:', sbcouplingdiag, file=f)
        print('S-B off-diagonal couplings:', sbcouplingoffdiag, file=f)
        print(file=f)

    # make diagnostic plots to make sure rotational averaging matches analytic
    s = str(figpath / '2d-reference.png')
    pool.submit(plot_2d, w1=w1, w3=w3, signal=ddref, path=s, axlim=limits)

    #s = str(figpath / '2d-reference-old.png')
    #pool.submit(plot_result, w1=w1, w3=w3, signal=ddref, path=s,
    #            show=False)

    s = str(figpath / '2d-fieldon.png')
    pool.submit(plot_2d, w1=w1, w3=w3, signal=dd.fieldon, path=s, axlim=limits)

    s = str(figpath / '2d-fieldoff.png')
    pool.submit(plot_2d, w1=w1, w3=w3, signal=dd.fieldoff, path=s, axlim=limits)

    s = str(figpath / '2d-stark.png')
    pool.submit(plot_2d, w1=w1, w3=w3, signal=dd.fieldon-dd.fieldoff, path=s, axlim=limits)

    for i in range(nstates):
        s = str(figpath / '2d-evecs{:d}.png'.format(i))
        pool.submit(plot_evecs, corr_energies, evecs2_trace, i, s)

    dd_projection = -(ddref).sum(axis=1)
    ddess_projection = -(dd.fieldon - dd.fieldoff).sum(axis=1)

    # do the same for absorption
    absref = np.array(absfile['reference'])
    shape = (100, *absref.shape)

    tmp = da.from_array(absfile['00000/data'], chunks=shape)
    pts_used = tmp.shape[0]
    rdataon = tmp[:pts_used//2].mean(axis=0)
    rdataoff = tmp[pts_used//2:].mean(axis=0)
    abs = StarkData(*dask.compute(rdataon, rdataoff))
    w3 = np.array(absfile['w3'])

    eigenenergies = fixed_energies2/1e3
    s = str(figpath / 'linear-fieldoff.png')
    pool.submit(plot_linear, w3=w3, signal=abs.fieldoff, path=s,
                eigenenergies=eigenenergies)

    s = str(figpath / 'linear-fieldon.png')
    pool.submit(plot_linear, w3=w3, signal=abs.fieldon, path=s,
                eigenenergies=eigenenergies)

    s = str(figpath / 'linear-stark.png')
    pool.submit(plot_linear, w3=w3, signal=abs.fieldon - abs.fieldoff, path=s,
                eigenenergies=eigenenergies)

    s = str(figpath / 'linear-projections.png')
    ax, scale = plot_linear(w3=w3, signal=abs.fieldoff, path=s)
    plot_linear(w3=w3, signal=dd_projection, path=s, ax=ax)

    s = str(figpath / 'linear-stark-projections.png')
    ax, scale = plot_linear(w3=w3, signal=abs.fieldon - abs.fieldoff, path=s)
    plot_linear(w3=w3, signal=ddess_projection, path=s, ax=ax,
                eigenenergies=eigenenergies)

if __name__ == '__main__':
    make_figures()
