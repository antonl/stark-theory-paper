from pyqcfp.runqcfp import QcfpConfig
import pathlib
import sys

path = pathlib.Path(sys.argv[1])

cfg = QcfpConfig.from_yaml_file(str(path))
cfg.simulation_type = 'absorption'
cfg.simulation_code = 'excitons-abs'

linearpath = str(path.with_suffix('').name + '-linear.yaml')

with open(linearpath, 'w') as f:
    cfg.to_yaml(stream=f)
