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
@click.option('--scale', default=-1)
def doit(file, limits, scale):
    # look for pump-probe data file
    path = Path(file)

    try:
        ddfile = h5py.File(str(path), 'r')
    except FileNotFoundError as e:
        print('Datafiles not found in dir {!s}'.format(path))
        return

    if scale < 0:
        scale = None

    # load ref
    ddref = np.array(ddfile['reference']).imag

    # calculate average pump-probe
    w3, w1 = np.array(ddfile['w3']), np.array(ddfile['w1'])

    # prepare folder for writing things
    p = Path(path)
    figpath = (p.parent / 'figures' / p.with_suffix('').name)
    #print('Figure path: ', str(figpath))
    figpath.mkdir(exist_ok=True, parents=True)


    # make diagnostic plots to make sure rotational averaging matches analytic
    s = str(figpath / '2d-reference.png')
    plot_2d(w1=w1, w3=w3, signal=ddref, path=s, axlim=limits, scale=scale)

if __name__ == '__main__':
    doit()

