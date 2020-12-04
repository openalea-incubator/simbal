""" Integration of C,N disponibility and
thermal time
* modification of organs properties
    (1) size
    (2) growth rate
    (3) senescence state
* modification of the number of organs 
    (1) apical growth
    (2) branching ??
* modification of organs relations : 
    (1) rosette -> 
    (2) bolting
    (3) flourishing
    (4) fructification
"""
from openalea import *
from openalea.lpy import *
from openalea.plantgl.all import *
#### Phenological parameters ####
transition1 = 15 #rosette -> bolting (days)
transition2 = 30 #bolting -> flower (days)
transition3 = 45 #flower -> fruits (days)

#### 1st stage of the plantule ####
l = Lsystem('caracolpy\colza_shape_v1.lpy')
lstring1, tt = l.axiom, 1
l.plot(lstring1)

# Development (derivation) if resources and thermaltime #
"""derive(lstring, iteration number, nb of iterations to apply) """
## rosette ###
variables = {'dl': 0.005}
l = Lsystem('caracolpy\colza_shape_v1.lpy', variables)
while tt < transition1 : 
    lstring1 = l.derive(lstring1,0,1)    
    # graphical interpretation
    tt +=1
l.plot(lstring1)

## transition to bolting ###
if tt == transition1 : 
    variables = {'dl':0.1, 'RGRa':0.01}
    l = Lsystem('caracolpy\colza_shape_v1.lpy', variables)

lstring2 = lstring1
while tt<transition2 : 
    lstring2 = l.derive(lstring2,0,1)    
    # graphical interpretation
    tt +=1
l.plot(lstring2)
l.plot(lstring1+lstring2)