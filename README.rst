========================
simbal
========================

.. {# pkglts, doc

.. image:: https://travis-ci.org/openalea-incubator/caracol.svg?branch=master
    :alt: Travis build status
    :target: https://travis-ci.org/openalea-incubator/caracol

.. #}

Simbal is a FSPM model for agroecological simulation

--------------------------
Model layout
--------------------------

**environment.py** : an above-ground environment module 
- spatial light model distribution
- temperature 
- spatialized light acquisition
- degree day cumulation 

**soilN.py** (or soil.py) a below-ground envt / plant interaction module 
- spatialized N acquisition 
- N spatial distribution and disponibility

**balanceCN.py**: a C,N balance module : 
- check available N : stored, acquired, remobilized
- distribute N among organs according to rules and hypothesis 
- check available C : stored, acquired, remobilized 
- distribute C among organs according to rules and hypothesis

**aptos.py** : a  module based on optimization and turnover theory that :
- integrate C-N availability, C-N potential acquisition, C-N re-allocation, processes carbon cost
- modify organs metabolic properties : C-N demand, C-N deficit, C-N content
- modify organs phenotypic properties : color, shape, area, size, growth
- modify organs relations : 
    * nitrogen source-sink status
- modify plant properties : number of organs, structure-function competition status

**main.py**:  run simulations and visualisation of the plant development

