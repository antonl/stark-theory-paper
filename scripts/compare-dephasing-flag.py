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

from .plotutils import *

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
        ddfile = h5py.File(str(path/'pump-probe-no-dephasing.h5'), 'r')
        absfile = h5py.File(str(path/'absorption.h5'), 'r')
    except FileNotFoundError as e:
        print('Datafiles not found in dir {!s}'.format(path))
        return

    if scale < 0:
        scale = None

    # load ref
    ddref = np.array(ddfile['reference']).imag
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
    pool.submit(plot_2d, w1=w1, w3=w3, signal=ddref, path=s, axlim=limits,
                scale=scale)


if __name__ == '__main__':
    make_figures()