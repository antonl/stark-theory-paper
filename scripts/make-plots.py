import dask
import dask.array as da
import h5py
import numpy as np
import click
from pathlib import Path
from collections import namedtuple
from concurrent.futures import ProcessPoolExecutor
from pyqcfp import QcfpConfig

StarkData = namedtuple('StarkData', ['fieldon', 'fieldoff'])

import matplotlib
matplotlib.use('Agg')

from matplotlib.pyplot import *
from matplotlib.colors import SymLogNorm
from matplotlib.ticker import SymmetricalLogLocator
from cycler import cycler

from plotutils import *

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

@click.command()
@click.argument('path', type=click.Path(file_okay=False, exists=True))
@click.option('--limits', default=(None,None), type=(float, float))
@click.option('-c','--ncores', default=6)
@click.option('--fudge-factor', default=0.)
@click.option('--scale', default=-1)
def make_figures(path, limits, ncores, fudge_factor, scale):
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

    if scale < 0:
        scale = None

    # load ref
    ddref = np.array(ddfile['reference']).imag
    shape = (100, *ddref.shape)

    # calculate average pump-probe
    tmp = da.from_array(ddfile['00000/data'], chunks=shape)
    pts_used = tmp.shape[0] - 1 # assume we did Stark averaging
    rdataon = tmp[:pts_used].imag.mean(axis=0)
    rdataoff = tmp[pts_used].imag # last one is the field-off data
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
    mu2_trace = np.linalg.norm(np.array(ddfile['00000/meta/ge dipoles'])[..., 2:], axis=-1)**2

    cfg = QcfpConfig.from_yaml(str(np.array(ddfile['cfg'])))
    if not cfg.include_complex_lifetimes:
        dephasingmat = dephasingmat.real

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

    if not cfg.include_complex_lifetimes:
        dephasingmat_trace = dephasingmat_trace.real

    imagdeph_trace = dephasingmat_trace[:, 1:nstates+1,0].imag
    corr_energies = energies_trace - reorgs_trace + imagdeph_trace + fudge_factor

    # prepare folder for writing things
    figpath = (path / 'figures')
    figpath.mkdir(exist_ok=True)

    with (figpath / 'eigen-energies.info').open('w') as f:
        print('Eigen-energies:', energies, file=f)
        print('GE reorganization energies:', reorgs, file=f)
        print('Reorg\'ed energies:', fixed_energies, file=f)
        print('Dephasing: ', imagdeph, file=f)
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
    pool.submit(plot_2d, w1=w1, w3=w3, signal=ddref, path=s, axlim=limits,
            scale=scale)
    #plot_2d(w1=w1, w3=w3, signal=ddref, path=s, axlim=limits)

    #s = str(figpath / '2d-reference-old.png')
    #pool.submit(plot_result, w1=w1, w3=w3, signal=ddref, path=s,
    #            show=False)

    s = str(figpath / '2d-fieldon.png')
    pool.submit(plot_2d, w1=w1, w3=w3, signal=dd.fieldon, path=s, axlim=limits,
            scale=scale)

    s = str(figpath / '2d-fieldoff.png')
    pool.submit(plot_2d, w1=w1, w3=w3, signal=dd.fieldoff, path=s,
            axlim=limits,
            scale=scale)

    s = str(figpath / '2d-stark.png')
    pool.submit(plot_2d, w1=w1, w3=w3, signal=dd.fieldon-dd.fieldoff, path=s,
            axlim=limits,
            scale=scale)

    for i in range(nstates):
        s = str(figpath / '2d-evecs{:d}.png'.format(i))
        pool.submit(plot_evecs, corr_energies, evecs2_trace, i, s, axlim=limits)

    dd_projection = -(ddref).sum(axis=1)
    ddess_projection = -(dd.fieldon - dd.fieldoff).sum(axis=1)

    # do the same for absorption
    absref = np.array(absfile['reference'])
    shape = (100, *absref.shape)

    tmp = da.from_array(absfile['00000/data'], chunks=shape)
    pts_used = tmp.shape[0] - 1
    rdataon = tmp[:pts_used].mean(axis=0)
    rdataoff = tmp[pts_used]
    abs = StarkData(*dask.compute(rdataon, rdataoff))
    w3 = np.array(absfile['w3'])

    eigenenergies = {'with dephasing': fixed_energies2/1e3,
                     'without dephasing': fixed_energies/1e3}

    # add the localization plot
    s = str(figpath / 'linear-localization.png')
    fig, (ax1, ax2) = subplots(2, 1, sharex=True)
    for i in range(0, nstates):
        weights_trace = mu2_trace*evecs2_trace[:, i, :]
        heights, bins = np.histogram(energies_trace.reshape(-1)/1e3,
                                     bins=80,
                                     weights=weights_trace.reshape(-1),
                                     density=False)
        widths = np.diff(bins)
        ax1.bar(bins[:-1], heights/heights.max(), widths, alpha=0.8,
                label='site {:d}'.format(i+1))

    ax1.legend()

    ax2.plot(w3/1e3, abs.fieldoff, label='field off')
    ax2.plot(w3/1e3, abs.fieldon, label='field on')
    ax2.plot(w3/1e3, abs.fieldon - abs.fieldoff, label='stark')
    ax2.set_xlabel(r'$\omega_\tau$ ($\times 10^3\ \mathrm{cm}^{-1}$)')
    ax2.set_xlim(*limits)

    ax2.legend()
    fig.savefig(str(s))

    s = str(figpath / 'linear-reference.png')
    pool.submit(plot_linear, w3=w3, signal=absref, path=s,
                axlim=limits, eigenenergies=eigenenergies,
                scale=scale)

    s = str(figpath / 'linear-fieldoff.png')
    pool.submit(plot_linear, w3=w3, signal=abs.fieldoff, path=s,
                axlim=limits, eigenenergies=eigenenergies,
                scale=scale)

    s = str(figpath / 'linear-fieldon.png')
    pool.submit(plot_linear, w3=w3, signal=abs.fieldon, path=s,
                axlim=limits, eigenenergies=eigenenergies,
                scale=scale)

    s = str(figpath / 'linear-stark.png')
    pool.submit(plot_linear, w3=w3, signal=abs.fieldon - abs.fieldoff, path=s,
                axlim=limits, eigenenergies=eigenenergies, scale=scale)

    s = str(figpath / 'linear-projections.png')
    ax, scale2 = plot_linear(w3=w3, signal=abs.fieldoff, path=s, axlim=limits)
    plot_linear(w3=w3, signal=dd_projection, path=s, ax=ax, axlim=limits,
            eigenenergies=eigenenergies, scale=scale)

    s = str(figpath / 'linear-stark-projections.png')
    ax, scale2 = plot_linear(w3=w3, signal=abs.fieldon - abs.fieldoff, path=s,
            axlim=limits, scale=scale)
    plot_linear(w3=w3, signal=ddess_projection, path=s, ax=ax,
                axlim=limits, eigenenergies=eigenenergies, scale=scale)
    print('submitted some figures')

    pool.shutdown(wait=True)
if __name__ == '__main__':
    make_figures()
