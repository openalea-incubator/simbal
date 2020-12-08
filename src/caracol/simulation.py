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
from morpho import *
from soilN import *
from environment import *
#### Phenological parameters ####
transition1 = 50 #rosette -> bolting (degree days) Jullien et al. 2011
transition2 = 300#bolting -> flower (days)
transition3 = 345 #flower -> fruits (days)
phyllochrone1 = 10
#### 1st stage of the plantule ####
parameter = {'Lmax_Internode':0.5}
l = Lsystem('plant.lpy', parameter)


# Development (derivation) if resources and thermaltime #
"""derive(lstring, corresponding iteration number, nb of iterations to apply) """
## rosette ###
lstring1, tt, i = l.axiom, 1, 0
while tt < transition1 : 
    scene = l.sceneInterpretation(lstring1)
    cscene, raw, agg = light_plant(scene)
    _, _=cscene.plot(raw, maxval = 1.01)
    # total BM obtained
    irr = leaf_irradiance(agg, lstring1, leaf_name = 'Leaf')
    # S variation : 
    dS_tot = N_acquisition()/sfqn() #cm2
    nb_tot = len(irr)
    dS_i = dS_tot / nb_tot
    # Compare if Ccost : 
    BM =  irr.apply(sum)['Ei'] * 0.7  # RUE W.cm-2 *
    #apply : 
    l = Lsystem('plant.lpy', {'GRL':dS_i})
    lstring1 = l.derive(lstring1,i,1) 
    #modify lstring
    tt+=1
    i+=1  
l.plot(lstring1)

## transition to bolting ###
if tt == transition1 : 
    variables = {'tip_Leaf':0.1, 'RGRa':0.01}
    l = Lsystem('caracolpy\colza_shape_v1.lpy', variables)

lstring2 = lstring1
while tt<transition2 : 
    lstring2 = l.derive(lstring2,0,1)    
    # graphical interpretation
    tt +=1
l.plot(lstring2)
l.plot(lstring1+lstring2)

## transition to flourishing ###
if tt == transition3 : 
    variables = {'dl':0.1, 'RGRa':0.01}
    l = Lsystem('caracolpy\colza_shape_v1.lpy', variables)

lstring2 = lstring1
while tt<transition2 : 
    lstring2 = l.derive(lstring2,0,1)    
    # graphical interpretation
    tt +=1
l.plot(lstring2)
l.plot(lstring1+lstring2)

### simulation : light on the plant ###
lsys, lstring, lscene = run_lsystem('plant.lpy', {'derivation length'})
lsys.plot(lstring)
cscene, raw, agg = illuminate(lscene, scene_unit= 'cm')
_, _=cscene.plot(raw, maxval = 1.01)

# total BM obtained
irr = leaf_irradiance(agg, lstring, leaf_name = 'Leaf')
C_acquired =  irr.apply(sum)['Ei'] * RUE #W.cm-2 * 

def organs(lstring):
    for mod in enumerate(lstring):
        if len(mod[1].name)>1: # select only the organs
            #do
            return()

def leaves(lstring):
    for mod in enumerate(lstring):
        if len(mod[1].name)=='Leaf': # select only the organs
            #do
            return()




for i in enumerate(l.axiom):
    a = i[1]
    print(type(a), openalea.lpy.__lpy_kernel__.ParamModule.setParameter(a, length, 4))
pass