import click
import pathlib
import subprocess
from copy import deepcopy
import pyqcfp
from pyqcfp.runqcfp import render_template

from pkg_resources import Requirement, resource_filename
BINDIR = resource_filename(Requirement.parse("pyqcfp"), "bin")

bin_path = str(pathlib.Path(BINDIR) / 'pyqcfp.calculator_2dttt_excitons')

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

    inpr = simdir / 'inp-rephasing.inp'
    with open(str(inpr), 'w') as f:
        f.write(render_template(cfg_rephasing))

    inpnr = simdir / 'inp-nonrephasing.inp'
    with open(str(inpnr), 'w') as f:
        f.write(render_template(cfg_nonrephasing))

    for sim in [inpr, inpnr]:
        input = str(sim)
        output = str(sim.with_suffix('.outp'))
        with open(str(sim.with_suffix('.text')), 'w') as f:
            subprocess.check_call([bin_path, input, output],
                                  stdout=f, stderr=subprocess.STDOUT)

if __name__ == '__main__':
    doit()

