# SIMBAL

## Dependencies

Python 3
numpy 
pandas 
scipy
pyvista
openalea.mtg
openalea.plantgl
alinea.caribu
rhizodep

## Install

- First create an environment

mamba create -n simbal -c openalea3 -c conda-forge pyvista openalea.mtg openalea.plantgl alinea.caribu numpy scipy pandas

- Activate the env

mamba activate simbal

- Install rhizodep

git clone https://github.com/openalea/rhizodep

cd rhizodep
python setup.py develop
