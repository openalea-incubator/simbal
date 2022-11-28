# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 10:44:53 2022

@author: FARRA
"""


#MTG
from openalea.mtg import MTG

g = MTG()

# Scale Plant
plant1 = g.add_component(g.root)
# Scale axis or GU or Phytomer
root1 = g.add_component(plant1)

v1 = g.add_child(root1, edge_type='<', carbon=12, label='Phytomer')
v2 = g.add_child(v1, edge_type='+')
v3 = g.add_child(v2, edge_type='<')
v4 = g.add_child(v3, edge_type= '+')

assert g.nb_vertices(scale=1) == 1
assert g.nb_vertices(scale=2) == 5
assert g.complex(v4) == plant1

#FONCTIONS

def CarbonBudget(
    ci: float, 
    cm: float, 
    cd: float, 
    ce: float):
    """ Allocate leaf carbon based on available carbon.

    Parameters
    ----------
    Ci : float
        Carbon import (input)
    Cm : float
        Metabolical C of the leaf = fonctional C (state)
    Cd : float
        Already available C on the leaf = structural C (state)
    Ce : float
        Exported C

    Returns
    -------
    Cm : float
        Metabolical C of the leaf = fonctional C (state)
    Cd : float
        Already available C on the leaf = structural C (state)
    mature: bool
        Is the leaf mature?
    """
    ci: int = 0.02
    cm: int = 0.1
    cd: int = 0.14
    ce: int = 0.01
    mature = False
    if ci > 0 : 
        cd = cm + ci - ce 
    else:
        mature = True

    return cd, mature
        
    
def YoungLeaf_OrganGrowth(lsph: float, 
                          x: float):    
    """ A juvenil leaf can still import = source/sink organ.
    
    Parameters
    ----------
    lsph : float
        Photosynthetic leaf surface (state)
        
    x : float
     Created surface with the available carbon (input)

    Returns
    -------
    lsph : float
        Photosynthetic leaf surface (state)
    """
    lsph :int= 5    # Leaf Surface photosynthetic
    x :int= 0.5 # How much Cc is necessary to construct a leaf surface unit ?
    lsph = lsph + x
    return lsph
    
        
def LeafMature(ce_m: float,
               ee: int, 
               nfonc: float,
               ne: float,
               cd_m: float):   
    """ A mature leaf cannot import = source organ
    
    Parameters
    ----------
    Ce : float
        Exported C
    Ee : int
        exported carbon threshold 
    Nd : float
        Available labil N on the leaf
    Nfonc : float
        Already available fonctionnal N on the leaf = used N for the photosynthetic apparel (state)
    Ne : float
        Exported N
    Cd : float
        Already available C on the leaf (state)

    Returns
    -------
   Nfonc : float
       Already available fonctionnal N on the leaf (state)
   Cd : float
       Already available C on the leaf (state)
   senesce: bool
       Is the leaf senesces ?
    """
    is_senesce= False
    if ce_m > ee: 
        nd = nfonc – ne	
        cd_m = cd-ce_m
    else:    
        senesce = True
    return nd, cd_m

        
        
def Senescence(nstruct: float,
               es: int,
               deltan: float,
               deltal: float,
               lsph: float,
               lsn: float,
               nfonc: float,
               cd: float) :
    """ Senesce is a development state of a leaf = source organ
    
    Parameters
    ----------
    Nstruct : float
        Structural N (state)
    Es : int
        structural N degradation threshold (input)
    deltaN : float
        degradated portion from structural N (input)
    deltaL : float
        dagradated surface from photosynthetic leaf surface (input)
    lsph : float
        photosynthetic leaf surface (state)
    lsn : float
        necrotic leaf surface (state)
    Nfonc : float
        Already available fonctionnal N on the leaf (state)
    Cd : float
        Already available C on the leaf (state)

    Returns
    -------
    degradation: bool
            Is the leaf can still degradate structure ressources ?
    exportation: bool
            Is the leaf can still export ressources ?
    leaffall: bool
            Is the leaf is full necrotic and then can fall ?
    Nstruct : float
         structural N on the leaf (state)
    Nd : float
        Available labil N in the leaf (state)
    lsph : float
        photosynthetic leaf surface (state)
    lsn : float
        necrotic leaf surface (state)
    Cd : float
        Already available C in the leaf (state)
    Nfall : float
        a falling leaf is a non biological state (state)
    
    """
    
    nstruc :int= 0.01
    nd :int=0.015
    es :int= 0.005 
    deltan = 0.0003
    deltal = 0.001
    lsph :int=5
    lsn :int= 0
    if nstruct > Es and lsph > 0:
        degradation = True
        exportation = True
        leaffall = False
        nstruct = nstruct – deltan
        nd = (nd + deltan)- ne
        lsph = lsph – deltal
        lsn = lsn + deltal
        cd = cd - ce
    return nstruct, nd, lsph, lsn, cd
    else:
        degradation = False
        exportation = False
        leaffall = True
    return nfall = nd   


