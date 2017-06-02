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

    if scale < 0:
        scale = None

    try:
        absfile = h5py.File(str(path), 'r')
    except FileNotFoundError as e:
        print('Datafiles not found in dir {!s}'.format(path))
        return

    absref = np.array(absfile['reference'])
    shape = (100, *absref.shape)

    tmp = da.from_array(absfile['00000/data'], chunks=shape)
    pts_used = tmp.shape[0] - 1
    rdataon = tmp[:pts_used].mean(axis=0)
    rdataoff = tmp[pts_used]
    abs = StarkData(*dask.compute(rdataon, rdataoff))
    w3 = np.array(absfile['w3'])

    # prepare folder for writing things
    p = Path(path)
    figpath = (p.parent / 'figures' / p.with_suffix('').name)
    #print('Figure path: ', str(figpath))
    figpath.mkdir(exist_ok=True, parents=True)

    s = str(figpath / 'linear-reference.png')
    pool.submit(plot_linear, w3=w3, signal=absref, path=s,
                axlim=limits, scale=scale)

    s = str(figpath / 'linear-fieldoff.png')
    pool.submit(plot_linear, w3=w3, signal=abs.fieldoff, path=s,
                axlim=limits, scale=scale)

    s = str(figpath / 'linear-fieldon.png')
    pool.submit(plot_linear, w3=w3, signal=abs.fieldon, path=s,
                axlim=limits, scale=scale)

    s = str(figpath / 'linear-stark.png')
    pool.submit(plot_linear, w3=w3, signal=abs.fieldon - abs.fieldoff, path=s,
                axlim=limits, scale=scale)

if __name__ == '__main__':
    doit()
