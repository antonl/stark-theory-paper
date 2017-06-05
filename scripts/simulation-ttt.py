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

BINDIR = resource_filename(Requirement.parse("pyqcfp"), "bin")

bin_path = str(pathlib.Path(BINDIR) / 'qcfp.calculator_2dttt_excitons')

@click.command()
@click.argument('template_yaml', type=click.Path(file_okay=True,
                                         dir_okay=False,
                                         readable=True))
def doit(template_yaml):
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
        t3 = rawdata['t3'].reshape(d, d)
        t1 = rawdata['t1'].reshape(d, d)
        data = (rawdata['real'] + 1j*rawdata['imag']).reshape(d, d)
        data[(0, -1), (0, -1)] *= 0.5

        if sim == inpr:
            mod = -1
        else:
            mod = 1
        data = data*np.exp(1j*(mod*cf1*t1 + cf3*t3)) # shift freq
        fdata = np.fft.ifft2(data, s=(2*d, 2*d))
        fdata = np.fft.fftshift(fdata)
        f1 = np.fft.fftfreq(d, abs(t1[1, 0] - t1[0, 0]))
        f3 = np.fft.fftfreq(d, abs(t3[0, 1] - t3[0, 0]))
        F1 = np.fft.fftshift(f1)
        F3 = np.fft.fftshift(f3)
        F1, F3 = np.meshgrid(F1, F3)

        outputs.append(fdata)
        plot_2d(w3=F3, w1=F1, signal=fdata.imag, path=str(sim.with_suffix('.png')))

    # make absorptive
    R, NR = outputs[0], outputs[1]
    ABS = 0.5*(R[::-1] + NR)

    plot_2d(w3=F3, w1=F1, signal=ABS.imag, path=str(simdir / 'absorptive.png'))

if __name__ == '__main__':
    doit()

