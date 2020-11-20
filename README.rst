========================
caracol
========================

.. {# pkglts, doc

.. image:: https://travis-ci.org/openalea-incubator/caracol.svg?branch=master
    :alt: Travis build status
    :target: https://travis-ci.org/openalea-incubator/caracol

.. #}

Caracol is a FSPM model for agroecological simulation

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
- check C and N stocks
- check C and N acquired
- remobilize C,N from senescing organs 
- distribute C,N among organs of the plant according to rules and hypothesis

**simulation.py** : a  module that
- integrate C,N dispo and thermal time 
- modify organs properties : size, growth rate, senescence state
- modify number of organs : apical growth 
- modify organs relations : 
    no leaf after floral transition but flower 
    flower give fruits after fruit transition

**plot.py**:  visualisation of the plant development

