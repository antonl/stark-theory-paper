import click
import pathlib
import subprocess
from copy import deepcopy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
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

    cfg_rephasing = cfg
    cfg_rephasing.simulation_type = 'rephasing'
    cfg_nonrephasing = deepcopy(cfg)
    cfg_nonrephasing.simulation_type = 'nonrephasing'

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
        d = int(np.sqrt(rawdata.t1.shape[0]))
        t3 = rawdata.t3.reshape(d, d)
        t1 = rawdata.t1.reshape(d, d)
        data = (rawdata.real + 1j*rawdata.imag).reshape(d, d)
        outputs.append(data)
        data[(0, -1), (0, -1)] *= 0.5

        fdata = np.fft.fft2(data, s=(2*d, 2*d))
        f1 = np.fft.fftfreq(2*d, abs(t1[1, 0] - t1[0, 0]))
        f3 = np.fft.fftfreq(2*d, abs(t3[0, 1] - t3[0, 0]))
        F1, F3 = np.meshgrid(f1, f3)
        fdata = np.fft.fftshift(fdata)
        F1 = np.fft.fftshift(F1)
        F3 = np.fft.fftshift(F3)

        plt.contour(F1, F3, -1*fdata.imag, 50)
        plt.savefig(str(sim.with_suffix('.png')))

if __name__ == '__main__':
    doit()

