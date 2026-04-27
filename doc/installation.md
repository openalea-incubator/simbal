# Installation

You must use conda environment : <https://docs.conda.io/en/latest/index.html>

## Users

### Create a new environment with simbal installed in there

```bash

mamba create -n simbal -c openalea3 -c conda-forge  openalea.simbal
mamba activate simbal
```

Install simbal in a existing environment

```bash
mamba install -c openalea3 -c conda-forge openalea.simbal
```

### (Optional) Test your installation

```bash
mamba install -c conda-forge pytest
git clone https://github.com/openalea/simbal.git
cd simbal/test; pytest
```

## Developers

### Install From source

```bash
# Install dependency with conda
mamba env create -n phm -f conda/environment.yml
mamba activate simbal

# Clone simbal and install
git clone https://github.com/openalea/simbal.git
cd simbal
pip install .

# (Optional) Test your installation
cd test; pytest
```
