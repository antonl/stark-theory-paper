import click
import pathlib
import straitlets
from pyqcfp.runqcfp import QcfpConfig
import numpy as np

class SitePositionConfig(straitlets.Serializable):
    site_index = straitlets.Integer(default_value=0)
    site_energies = straitlets.List(trait=straitlets.Float, default_value=[0.])
    range = straitlets.Tuple(type=straitlets.Float, default_value=(-400, 400))
    count = straitlets.Integer(default_value=3)

@click.command()
@click.argument('cfgpath', type=click.Path(file_okay=True,
                                           dir_okay=False,
                                           exists=True))
@click.option('--index', '-i', type=int, default=0)
@click.option('--range', '-r', type=(float, float), default=(0.01, 1.1))
@click.option('--count', '-c', type=int, default=3)
def doit(cfgpath, index, range, count):
    cfg = QcfpConfig.from_yaml_file(str(cfgpath))

    if cfg.nsites < index:
        raise ValueError('index too large for system hamiltonian size')

    scfg = SitePositionConfig()
    scfg.range = range
    scfg.count = count
    site_energy = cfg.system_hamiltonian[index][index]
    scfg.site_energies = (site_energy + np.linspace(scfg.range[0], scfg.range[1], scfg.count)).tolist()

    with open('site_energy_config.yaml', 'w') as f:
        scfg.to_yaml(stream=f)

    path = str(pathlib.Path(cfgpath).with_suffix(''))
    for i, m in enumerate(scfg.site_energies):
        p = path + '-{:03d}.yaml'.format(i)
        with open(p, 'w') as f:
            cfg.system_hamiltonian[i][i] = m
            cfg.to_yaml(stream=f)

if __name__ == '__main__':
    doit()
