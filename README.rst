========================
caracol
========================

.. {# pkglts, doc

.. image:: https://travis-ci.org/openalea-incubator/caracol.svg?branch=master
    :alt: Travis build status
    :target: https://travis-ci.org/openalea-incubator/caracol

.. #}

Caracol is a FSPM model for agroecological simulation

# TO DO 
an above-ground environment module : <ag_environment.py> 
- spatial light model distribution
- temperature 

an above-ground envt / plant interaction module : <ag_interaction.py> 
- spatialized light acquisition
- degree day cumulation 

a below-ground environment module : <bg_environment.py> 
- N spatial distribution and disponibility

a below-ground envt / plant interaction module : <bg_interaction.py> 
- spatialized N acquisition 

a C,N balance module : <balanceCN.py>
- check C and N stocks
- check C and N acquired
- remobilize C,N from senescing organs 
- distribute C,N among organs of the plant according to rules and hypothesis

a developing module : <simulation.py>
- integrate C,N dispo and thermal time 
- modify organs properties : size, growth rate, senescence state
- modify number of organs : apical growth 
- modify organs relations : 
    no leaf after floral transition but flower 
    flower give fruits after fruit transition

a visualisation of the plant development

