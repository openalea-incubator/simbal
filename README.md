<<<<<<< HEAD
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

## Test
- To run a simulation
    * Add a directory simulation/outputs
        mkdir -p simulation/output
    * launch ipython
        ipython
        %gui qt
        %run simbal/simbal_orchestrator.py 
>>>>>>> frederic
