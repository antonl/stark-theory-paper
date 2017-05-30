import click
import pathlib
import straitlets
from pyqcfp.runqcfp import QcfpConfig
import numpy as np

class VoltageConfig(straitlets.Serializable):
    stark_field_magnitudes = straitlets.List(trait=straitlets.Float, default_value=[1.1])
    range = straitlets.Tuple(type=straitlets.Float, default_value=(0.01, 1.1))
    count = straitlets.Integer(default_value=3)

@click.command()
@click.argument('cfgpath', type=click.Path(file_okay=True,
                                           dir_okay=False,
                                           exists=True))
@click.option('--range', '-r', type=(float, float), default=(0.01, 1.1))
@click.option('--count', '-c', type=int, default=3)
def doit(cfgpath, range, count):
    cfg = QcfpConfig.from_yaml_file(str(cfgpath))

    vcfg = VoltageConfig()
    vcfg.range = range
    vcfg.count = count
    vcfg.stark_field_magnitudes = np.linspace(vcfg.range[0], 
                                              vcfg.range[1],
                                              vcfg.count).tolist()

    with open('voltagecfg.yaml', 'w') as f:
        vcfg.to_yaml(stream=f)

    path = str(pathlib.Path(cfgpath).with_suffix(''))
    for i, m in enumerate(vcfg.stark_field_magnitudes):
        p = path + '-{:03d}.yaml'.format(i)
        with open(p, 'w') as f:
            vec = np.array(cfg.stark_field_vector)
            vec /= np.linalg.norm(vec)
            cfg.stark_field_vector = (m*vec).tolist()
            cfg.to_yaml(stream=f)

if __name__ == '__main__':
    doit()
