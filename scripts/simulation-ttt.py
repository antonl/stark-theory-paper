import click
import pathlib
import subprocess
from copy import deepcopy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from  plotutils import *
import pyqcfp
from pyqcfp.runqcfp import render_template

from pkg_resources import Requirement, resource_filename
import numpy as np
from scipy.signal import get_window

BINDIR = resource_filename(Requirement.parse("pyqcfp"), "bin")

bin_path = str(pathlib.Path(BINDIR) / 'qcfp.calculator_2dttt_excitons')

@click.command()
@click.argument('template_yaml', type=click.Path(file_okay=True,
                                         dir_okay=False,
                                         readable=True))
@click.option('--limits', default=(None, None), type=(float, float))
@click.option('--window/--no-window', default=False)
def doit(template_yaml, limits, window):
    path = pathlib.Path(template_yaml)
    cfg = pyqcfp.QcfpConfig.from_yaml_file(str(path))
    cfg.stark_perturbation = False
    cfg.analytic_orientational_averaging = True

    cf1 = 0.5*(cfg.w1_min + cfg.w1_max)
    cf3 = 0.5*(cfg.w3_min + cfg.w3_max)

    cfg_rephasing = cfg
    cfg_rephasing.simulation_type = 'rephasing-ttt'
    cfg_nonrephasing = deepcopy(cfg)
    cfg_nonrephasing.simulation_type = 'nonrephasing-ttt'

    simdir = path.parent

    inpr = simdir / 'rephasing.inp'
    with open(str(inpr), 'w') as f:
        f.write(render_template(cfg_rephasing))

    inpnr = simdir / 'nonrephasing.inp'
    with open(str(inpnr), 'w') as f:
        f.write(render_template(cfg_nonrephasing))

    if window:
        window = get_window(('general_gaussian', 2, 100), cfg.nfreqs,
                fftbins=True)[cfg.nfreqs//2:]
    else:
        window = np.ones((cfg.nfreqs//2,))

    # always scale the ends
    window[0] *= 0.5
    window[-1] *= 0.5

    outputs = []
    for sim in [inpr, inpnr]:
        input = str(sim)
        output = str(sim.with_suffix('.outp'))
        with open(str(sim.with_suffix('.text')), 'w') as f:
            subprocess.check_call([bin_path, input, output],
                                  stdout=f, stderr=subprocess.STDOUT)
        rawdata = np.genfromtxt(output,
                             comments='#',
                             delimiter='\t',
                             names=['t3', 't1', 'real', 'imag'])
        d = int(np.sqrt(rawdata['t1'].shape[0]))
        t1 = 2*rawdata['t1'].reshape(d, d)
        t3 = 2*rawdata['t3'].reshape(d, d)

        # try generating the same time axis:
        freF3 = cfg.w3_max - cfg.w3_min
        freS3 = freF3/cfg.nfreqs
        timeS3 = 2*np.pi/freF3

        freF1 = cfg.w1_max - cfg.w1_min
        freS1 = freF1/cfg.nfreqs
        timeS1 = 2*np.pi/freF1

        #t1 = timeS1*np.arange(0, cfg.nfreqs//2)
        #t3 = timeS3*np.arange(0, cfg.nfreqs//2)
        #t3, t1 = np.meshgrid(t3, t1)

        dt1 = abs(t1[1, 0] - t1[0, 0])
        dt3 = abs(t3[0, 1] - t3[0, 0])
        
        print('F3 range: {:f}, F3 step: {:f}, T3 step: {:f}'.format(freF3,
            freS3, timeS3))
        print('Extracted dt3: {:f}'.format(dt3))
        print('F1 range: {:f}, F1 step: {:f}, T1 step: {:f}'.format(freF1,
            freS1, timeS1))
        print('Extracted dt1: {:f}'.format(dt1))
        
        data = (rawdata['real'] + 1j*rawdata['imag']).reshape(d, d)
        unwindowed = np.real(np.diagonal(data).copy())
        unwindowed /= np.max(np.abs(unwindowed))
        data = np.einsum('ij,j->ij', data, window)
        data = np.einsum('ij,i->ij', data, window)
        windowed = np.real(np.diagonal(data).copy())
        windowed /= np.max(np.abs(windowed))

        # plot windowed trace
        plt.figure()
        plt.plot(np.diag(t3), unwindowed)
        plt.plot(np.diag(t3), windowed)
        plt.plot(np.diag(t3), window)
        plt.savefig(str(sim.with_suffix('')) + '-windowed.png')

        if sim == inpr:
            mod = -1
        else:
            mod = 1

        data = data*np.exp(1j*(mod*cf1*t1 + cf3*t3)) # shift freq
        fdata = np.fft.ifft2(data, s=(2*d, 2*d))
        fdata = np.fft.fftshift(fdata)

        f1 = 2*np.pi*np.fft.fftfreq(2*d, dt1)
        f3 = 2*np.pi*np.fft.fftfreq(2*d, dt3)
        F1 = np.fft.fftshift(f1)
        F3 = np.fft.fftshift(f3)
        F1, F3 = np.meshgrid(F1, F3)
        #F1 += cfg.w1_min
        #F3 += cfg.w3_min
        F1 += cf1
        F3 += cf3

        outputs.append(fdata.T)
        plot_2d(w3=F3, w1=F1, signal=fdata.imag.T, 
                path=str(sim.with_suffix('.png')), axlim=limits)

    # make absorptive
    R, NR = outputs[0], outputs[1]
    ABS = (R[:, ::-1]+ NR)*(2*d)**2

    plot_2d(w3=F3, w1=F1, signal=ABS.imag, path=str(simdir / 'absorptive.png'),
            axlim=limits)

if __name__ == '__main__':
    doit()

