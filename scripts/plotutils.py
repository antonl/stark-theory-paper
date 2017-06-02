from pathlib import Path
from collections import namedtuple

import matplotlib
matplotlib.use('Agg')

from matplotlib.pyplot import *
from matplotlib.colors import SymLogNorm
from matplotlib.ticker import SymmetricalLogLocator
from cycler import cycler

def normalize(data):
    return data / np.abs(data).max()

def latex_float(f, precision=2):
    float_str = "{0:." + str(int(precision)) + "e}"
    float_str = float_str.format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"{0} \times 10^{{{1}}}".format(base, int(exponent))
    else:
        return float_str

def plot_evecs(energies, evecs2, which, path, axlim=(None, None)):
    pts_used = energies.shape[0]
    nvecs = energies.shape[1]
    fig = figure()
    ax = fig.add_subplot(111, polar=True)
    cyc = cycler(color=['b', 'r', 'g'])

    energies2 = energies/1e3
    #for which in range(nvecs):
    for i,prop in zip(range(nvecs), cyc):
        ax.scatter(3*np.pi/2*evecs2[:pts_used//2, i, which],
                   energies2[:pts_used//2, which], alpha=0.8, **prop)

    ax.set_rlim(*axlim)

    ax.set_theta_offset(-np.pi/4)
    ax.set_rlabel_position(np.degrees(3*np.pi/2))
    ax.set_xticklabels(['{:3.0f}%'.format(x) for x in np.linspace(0, 100, 7)] + [''])
    fig.savefig(path, bbox_inches='tight')

def plot_2d(w1, w3, signal, path, invert_w1=False, scale=None,
            axlim=(None,None)):
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
    ax.set_xlabel(r'$\omega_\tau$ ($\times 10^3\ \mathrm{cm}^{-1}$)')
    ax.set_ylabel(r"$\omega_t$ ($\times 10^3\ \mathrm{cm}^{-1}$)")
    ax.text(0.99, 0.01, r'$\mathrm{{scale:}} {!s}$'.format(latex_float(scale)),
            transform=ax.transAxes,
            horizontalalignment='right',
            verticalalignment='bottom')

    loc = matplotlib.ticker.MaxNLocator(11)
    fmt = matplotlib.ticker.ScalarFormatter(useOffset=False, useMathText=True)

    cb = fig.colorbar(qset, ax=ax, ticks=loc, format=fmt)
    cb.add_lines(c)

    ax.set_xlim(*axlim)
    ax.set_ylim(*axlim)

    #ax.grid(True)
    ax.add_artist(Line2D((0, 1), (0, 1), linewidth=2, color='k', alpha=0.5,
                         transform=ax.transAxes))

    #ax.relim()
    #ax.autoscale_view()
    fig.savefig(str(path))

def plot_linear(w3, signal, path, scale=None, axlim=(None, None), ax=None,
                eigenenergies=None):

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

    ax.plot(w3, signal2, linewidth=2)
    ax.set_xlabel(r'$\omega_\tau$ ($\times 10^3\ \mathrm{cm}^{-1}$)')
    ax.set_ylabel(r'norm. abs.')
    ax.text(0.99, 0.01, r'$\mathrm{{scale:}} {!s}$'.format(latex_float(scale)),
            transform=ax.transAxes,
            horizontalalignment='right',
            verticalalignment='bottom')
    ax.set_xlim(*axlim)
    ax.set_ylim(-1.05, 1.05)
    ax.add_line(Line2D(axlim, [0, 0], alpha=0.8, linestyle='--', linewidth=0.5))

    # add eigenstate positions
    if eigenenergies is not None:
        ax.vlines(eigenenergies['without dephasing'],
                  -1, 1,
                  linewidth=2,
                  alpha=0.7, linestyle='--')
        ax.vlines(eigenenergies['with dephasing'],
                  -1, 1,
                  linewidth=2,
                  alpha=0.7)

    #ax.relim()
    #ax.autoscale_view()
    fig.savefig(str(path))
    return ax, scale
