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

@click.command()
@click.argument('file', type=click.Path(file_okay=True,
                                        dir_okay=False,
                                        exists=True))
@click.option('--limits', default=(None,None), type=(float, float))
@click.option('-c','--ncores', default=6)
@click.option('--fudge-factor', default=0.)
@click.option('--scale', default=-1)
def doit(file, limits, ncores, fudge_factor, scale):
    # look for pump-probe data file
    path = Path(file)
    pool = ProcessPoolExecutor(max_workers=ncores)
    #rcParams.update(params)

    try:
        ddfile = h5py.File(str(path), 'r')
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
    p = Path(path)
    figpath = (p.parent / 'figures' / p.with_suffix('').name)
    #print('Figure path: ', str(figpath))
    figpath.mkdir(exist_ok=True, parents=True)

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

if __name__ == '__main__':
    doit()
