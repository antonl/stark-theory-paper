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

from plotutils import *

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
        ddfile_real = h5py.File(str(path/'pump-probe-real.h5'), 'r')
        ddfile_complex = h5py.File(str(path/'pump-probe-complex.h5'), 'r')
        absfile_real = h5py.File(str(path/'absorption-real.h5'), 'r')
        absfile_complex = h5py.File(str(path/'absorption-complex.h5'), 'r')

    except FileNotFoundError as e:
        print('Datafiles not found in dir {!s}'.format(path))
        return

    if scale < 0:
        scale = None

    # load ref
    ddref = np.array(ddfile_real['reference']).imag
    w3, w1 = np.array(ddfile_real['w3']), np.array(ddfile_real['w1'])

    # load evecs
    energies = np.array(ddfile_real['meta/one band energies'])
    nstates = energies.shape[0]
    evecs2 = np.array(ddfile_real['meta/ge eigenvectors'])**2
    reorgs = np.array(ddfile_real['meta/reorganization energy matrix'])
    sbcouplingdiag = np.diag(ddfile_real['meta/sb coupling diagonal'])
    sbcouplingoffdiag = np.diag(ddfile_real['meta/sb coupling off-diagonal'])
    redfield = np.array(ddfile_real['meta/redfield relaxation matrix'])
    reorgs = np.diag(reorgs)[1:nstates+1]
    redfield = np.diag(redfield)[1:nstates+1]
    dephasingmat = np.array(ddfile_real['meta/lifetime dephasing matrix'])[:,
                   ::2] + \
                   1j*np.array(ddfile_real['meta/lifetime dephasing '
                                           'matrix'])[:, 1::2]
    imagdeph = dephasingmat[1:nstates+1,0].imag

    fixed_energies = energies - reorgs
    fixed_energies2 = energies - reorgs + imagdeph + fudge_factor

    # prepare folder for writing things
    figpath = (path / 'figures')
    figpath.mkdir(exist_ok=True)

    with (figpath / 'eigen-energies.info').open('w') as f:
        print('Without complex lifetimes: ', file=f)
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
    s = str(figpath / '2d-real.png')
    pool.submit(plot_2d, w1=w1, w3=w3, signal=ddref, path=s, axlim=limits,
                scale=scale)

    # load ref
    ddref = np.array(ddfile_complex['reference']).imag
    w3, w1 = np.array(ddfile_complex['w3']), np.array(ddfile_complex['w1'])

    # load evecs
    energies = np.array(ddfile_complex['meta/one band energies'])
    nstates = energies.shape[0]
    evecs2 = np.array(ddfile_complex['meta/ge eigenvectors'])**2
    reorgs = np.array(ddfile_complex['meta/reorganization energy matrix'])
    sbcouplingdiag = np.diag(ddfile_complex['meta/sb coupling diagonal'])
    sbcouplingoffdiag = np.diag(ddfile_complex['meta/sb coupling off-diagonal'])
    redfield = np.array(ddfile_complex['meta/redfield relaxation matrix'])
    reorgs = np.diag(reorgs)[1:nstates+1]
    redfield = np.diag(redfield)[1:nstates+1]
    dephasingmat = np.array(ddfile_complex['meta/lifetime dephasing matrix'])[:,
                   ::2] + \
                   1j*np.array(ddfile_complex['meta/lifetime dephasing '
                                           'matrix'])[:, 1::2]
    imagdeph = dephasingmat[1:nstates+1,0].imag

    fixed_energies = energies - reorgs
    fixed_energies2 = energies - reorgs + imagdeph + fudge_factor

    with (figpath / 'eigen-energies.info').open('a+') as f:
        print('With complex lifetimes: ', file=f)
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
    s = str(figpath / '2d-complex.png')
    pool.submit(plot_2d, w1=w1, w3=w3, signal=ddref, path=s, axlim=limits,
                scale=scale)


    # do the same for absorption
    absref = np.array(absfile_real['reference'])
    w3 = np.array(absfile_real['w3'])

    eigenenergies = {'with dephasing': fixed_energies2/1e3,
                     'without dephasing': fixed_energies/1e3}

    s = str(figpath / 'linear-real.png')
    ax, scale2 = plot_linear(w3=w3, signal=absref, path=s,
                        axlim=limits, eigenenergies=eigenenergies,
                        scale=scale)

    absref = np.array(absfile_complex['reference'])
    w3 = np.array(absfile_complex['w3'])

    s = str(figpath / 'linear-complex.png')
    plot_linear(w3=w3, signal=absref, path=s,
                axlim=limits, scale=scale, ax=ax)

if __name__ == '__main__':
    make_figures()
