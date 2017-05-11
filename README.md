# 2DESS Theory Paper

This project collects the resources for running the simulations in the 2DESS
Theory paper, as well as the paper itself. By using the ```Snakemake``` file,
it should be straight-forward to reproduce the results in the future.

## Requirements

The simulation code assumes that everything is being run on a Linux machine.
Although it is possible to run these simulations on Windows, it is more
challenging to setup and may require additional tweaking. 

- Linux machine with gcc compiler
- git
- Anaconda python distribution

## Setting up the environment

To setup the environment, clone the project into some directory and then create
the conda environment which is then used for the simulations. Conda should
obtain all of the required dependencies.

```
# clone workflow into working directory
git clone https://bitbucket.org/user/myworkflow.git path/to/workdir
cd path/to/workdir

# edit config and workflow as needed
vim config.yaml

# install dependencies into isolated environment
conda env create -n simulations --file environment.yaml

# activate environment
source activate myworkflow

# execute workflow
snakemake -n
```

## Running the simulations

To be written
