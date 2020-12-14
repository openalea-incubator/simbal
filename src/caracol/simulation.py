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
from balanceCN import *

# Development (derivation) if resources and thermaltime #
"""derive(lstring, corresponding iteration number, nb of iterations to apply) """

###################################################################################################
""" Define first plant"""
#### Phenological parameters ####
transition1 = 50 #rosette -> bolting (degree days) Jullien et al. 2011
transition2 = 300#bolting -> flower (days)
transition3 = 345 #flower -> fruits (days)
phyllochrone1 = 10
l = Lsystem('interpret.lpy')
# Define Axiom :
pI = ParameterSet(area = 0.3)
pP = ParameterSet(length = 0.01)
pL = ParameterSet(area = 1)
pA = ParameterSet(age = 0)
lstring = Lsystem('plant3.lpy', {'pI':pI, 'pP':pP, 'pL':pL}).axiom +Lsystem('Bud.lpy', {'pA':pA}).axiom
tt,i = 1,0
plant = balance(lstring, 0)
plant = plant.loc[:, ['area','rank', 'Organ']]   
plant.loc[:,'Ntot'] = 0.13

###################################################################################################
biblio = []
init = False
for k in range(5):
    # balance N
    N_disp_aer = 0.13
    step = balance(lstring, N_disp_aer)
    plant = pandas.concat([plant, step.loc[:,['N_attr']] ], axis = 1)
    plant.loc[:,'Ntot'] =plant['Ntot'] + plant['N_attr']
    plant = plant.loc[:, ['area', 'rank', 'Organ', 'Ntot']]

    # surface growth
    plant.loc[:,'dS'] = (plant['Ntot'] * NgNt())/sfqn()

    # leaf initiation
    N4bud = plant[plant['Organ']=='Apex'].iloc[0]['Ntot']
    if N4bud >= BT() : 
        # create new lines in df with increase
        newI = dict(zip(plant.columns,[0, max(plant['rank'])+1, 'Internode', IoL()*N4bud]))
        newL = dict(zip(plant.columns,[0, max(plant['rank'])+1, 'Leaf',(1-IoL())*N4bud]))
        plant = plant.append(newI, ignore_index = True)
        plant = plant.append(newL, ignore_index = True)
        plant.loc[plant['Organ']=='Apex','Ntot']=0
        init = True
    # lstring actualisation
    p=[]
    j= l.axiom #empty AxialTree 
    for i in range(max(plant['rank'])): 
        for org in plant['Organ']:
            a = plant[plant['Organ']==org].iloc[i]['dS'] + plant[plant['Organ']==org].iloc[i]['area'] *1e4
            p += [ ParameterSet(area = a) ]
        pI, pL, pA = p
        j += Lsystem('plant3.lpy', {'pI':pI, 'pP':pP, 'pL':pL}).axiom
    if init : 
        j+= Lsystem('Bud.lpy', {'pA':pA}).axiom
    
    lstring = j
    biblio +=[j]
