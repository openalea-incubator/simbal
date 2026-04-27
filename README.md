[![Docs](https://readthedocs.org/projects/simbal/badge/?version=latest)](https://simbal.readthedocs.io/en/latest/?badge=latest)
[![Build Status](https://github.com/openalea-incubator/simbal/actions/workflows/openalea_ci.yml/badge.svg)](https://github.com/openalea-incubator/simbal/actions/workflows/openalea_ci.yml)
[![Python Version](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12%20%7C%203.13-blue)](https://www.python.org/downloads/)
[![Anaconda-Server Badge](https://anaconda.org/openalea3/openalea.simbal/badges/version.svg)](https://anaconda.org/openalea3/openalea.simbal)
[![License](https://img.shields.io/badge/License--CeCILL-C-blue)](https://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html)

# Simbal
FSPM model for crop mixture with root and shoot, integrating C & N cycle

# Final aim 

See the dvt of OSR from seed to seed production

Have the least number of parameters to set

Parameters only about phenology and morphology (which can be considered genetically determined)

See the influence of the environment, in particular light on C, N acquisition and allocation among the organs

=======
# SIMBAL

## Dependencies

Python 3
numpy 
pandas 
scipy
openpyxl
pyvista
openalea.mtg
openalea.plantgl
openalea.caribu
openalea.rhizodep

### Installation

Conda must be installed, please follow the instructions and recommendations [here](https://openalea.readthedocs.io/en/latest/install.html).

#### for user
Creating a new conda environment with simbal and its dependencies installed
```bash
mamba create -n simbal -c openalea3/label/dev -c openalea3 -c conda-forge openalea.simbal
```

activate the environment
```bash
mamba activate simbal
```

#### for developer
```bash
mamba env create -f ./conda/environment.yml
```

activate the environment
```bash
mamba activate simbal_dev
```

This will create a conda environment with dependencies installed and install simbal in editable state.

## Test
- To run a simulation
    * Add a directory simulation/outputs
        mkdir -p simulation/output
    * launch ipython
        ipython
        %gui qt
        %run src/openalea/simbal/simbal_orchestrator.py 

