# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 10:44:53 2022

@author: FARRA
"""
# Coding rules
# Objects : class CarbonBudget
# Functions : def carbon_budget
# Variables : lower case

def light_interception(
    projected_leaf_area: float,
    k: float,
    e_amax: float):
    """ Calculation of PAR intercepted by the plant with Beer Lambert law
    
    Parameters
    ----------
    projected_leaf_area : float
        projected leaf area at the beginning of the time step in cm2
    k : float
        extinction coefficient
    e_amax: float
        maximum intercepted radiation
   
    Returns
    -------
    intercepted_PAR : float
        in MJ.cm-2
    """
    intercepted_PAR = e_amax * (1 - exp(k * projected_leaf_area)) 
    return intercepted_PAR


def photosynthesis(
    intercepted_PAR: float, 
    RUE: float):
    """ Transform radiation intercepted into biomass

    Parameters
    ----------

    intercepted_PAR: float
        intercepted PAR by the plant
    RUE: float
        Radiation use efficiency
    
    Returns
    -------
    producted_biomass : float
        ??????
     
    """
    producted_biomass = RUE * intercepted_PAR
    return producted_biomass




def carbon_budget(
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
    ci: float = 0.02
    cm: float = 0.1
    cd: float = 0.14
    ce: float = 0.01

    mature = False
    if ci > 0 : 
        cd = cm + ci - ce 
    else:
        mature = True

    return cd, mature
        
    
def young_leaf_organ_growth(lsph: float, 
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
    
        
def leaf_mature(ce_m: float,
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

        
        
def senescence(nstruct: float,
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
    
    """
    
    nstruc :int= 0.01
    nd :int=0.015
    es :int= 0.005 
    deltan = 0.0003
    deltal = 0.001
    lsph :int=5
    lsn :int= 0

    degradation = True
    exportation = True
    leaffall = False

    if nstruct > Es and lsph > 0:
        nstruct = nstruct – deltan
        nd = (nd + deltan)- ne
        lsph = lsph – deltal
        lsn = lsn + deltal
        cd = cd - ce
    else:
        degradation = False
        exportation = False
        leaffall = True

    return nstruct, nd, lsph, lsn, cd, degradation, exportation, leaffall
    
    
    
    