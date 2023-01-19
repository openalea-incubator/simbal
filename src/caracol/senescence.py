# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 10:44:53 2022

@author: FARRA
"""

from math import (exp, pi)


#Nitrogen Budget
def potential_N_root_uptake(
    _root_surface: float, 
    Nitrate: float, 
    time_step: int, 
    SNU: float):
    """
    Calculation of N potential uptake by roots

    Parameters
    ----------
    _root_surface: float
        root object
    Nitrate: float
        N quantity in soil in gN
    time_step: float
        time step of the simulation in day
    SNU:


    Returns
    -------
    _potential_N_root_uptake: float
        potential N uptake by roots
    
    """
    _potential_N_root_uptake = SNU * Nitrate * time_step * _root_surface
    print("N root uptake is", _potential_N_root_uptake, "gN")
    return _potential_N_root_uptake


def N_potential_offer(
    _potential_N_root_uptake: float,
    respiration_cost: float):
    """
    Calculation of N available to the plant (in gN)
    
    Parameters
    ----------
    _potential_N_root_uptake: float
        N uptake by roots (in g)   
    respiration_cost: float
        cost of respiration

    Returns
    -------
    _N_potential_offer: float
        N potential offer to the plant in gN
    """
    _N_potential_offer = _potential_N_root_uptake - respiration_cost
    print("N available to the plant is", _N_potential_offer, "gN")
    return _N_potential_offer


def N_potential_root_demand(
    _delta_volume: float, 
    Ndensity: int):
    """
    Calculation of N potential root demand (in gN) in function of the root volume and root N density

    Parameters
    ----------
    delta_volume: float
        potential increase of root volume in m3
    Ndensity: float
        N root density in gN.m-3

    Returns
    -------
    _N_potential_root_demand: float
        N root potential demand in gN

    """
    _N_potential_root_demand = _delta_volume * Ndensity
    print("N potential root demand is", _N_potential_root_demand, "gN")
    return _N_potential_root_demand


def N_potential_canopy_offer(
    _N_potential_offer: float, 
    _N_potential_root_demand: float):
    """
    Calculation of the potential shoot demand (gN) in function of N available to the plant (gN) and potential root N demand (gN)

    Parameters
    ----------
    _N_potential_offer: float
        N offer to the plant in gN
    _N_potential_root_demand: float
        potential root N demand in gN
    
    Returns
    -------
    _N_potential_canopy_offer: float
        N shoot potential demand in gN
    """
    _N_potential_canopy_offer = _N_potential_offer - _N_potential_root_demand
    print("N potential shoot demand is", _N_potential_canopy_offer, "gN")
    return _N_potential_canopy_offer



#Light Budget
def canopy_light_interception(
    projected_canopy_area: float,
    k: float,
    e_amax: float):
    """ Calculation of PAR intercepted by the canopy with Beer Lambert law
    
    Parameters
    ----------
    projected_canopy_area : float
        projected canopy area at the beginning of the time step in m2
    k : float
        extinction coefficient (dimensionless)
    e_amax: float
        maximum intercepted radiation
   
    Returns
    -------
    canopy_intercepted_PAR : float
        in MJ.m-2
    """
    canopy_intercepted_PAR = e_amax * (1 - exp(k * projected_canopy_area)) 
    return canopy_intercepted_PAR


def canopy_photosynthesis(
    canopy_intercepted_PAR: float, 
    RUE: float):
    """ Transform radiation intercepted into biomass by the canopy

    Parameters
    ----------

    canopy_intercepted_PAR: float
        intercepted PAR by the plant
    RUE: float
        Radiation use efficiency
    
    Returns
    -------
    producted_biomass: float
        producted C
    rejacted_C: float
        CO2 rejacted by photosynthesis
     
    """
    producted_biomass = RUE * canopy_intercepted_PAR
    return producted_biomass



#Carbon Budget
def C_potential_root_demand(
    _delta_volume: float, 
    Cdensity: int):
    """
    Calculation of C potential root demand (in gC) in function of the root segment volume and root C density

    Parameters
    ----------
    delta_volume: float
        potential increase of root volume in cm3
    Cdensity: float
        C root density in gC.m-3
    
    Returns
    -------
    _C_potential_root_demand: float
        C root potential demand in gC

    """
    _C_potential_root_demand = _delta_volume * Cdensity
    print("C potential root demand is", _C_potential_root_demand, "gC")
    return _C_potential_root_demand


def potential_leaf_area(
    _N_potential_leaf_demand: float, 
    coef_N_leaf_surface: int):
    """
    Calculation of potential leaf area (cm2) in function of potential N leaf demand (gN) and the coefficient of the relationship of leaf area - leaf N quantity (cm2.g-1)

    Parameters
    ----------
    N_leaf_potential_demand: float
        leaf potential N demand in gN
    coef_N_leaf_area: int
        coefficient between leaf area and leaf N quantity in m2.g-1

    Returns
    -------
    _potential_leaf_surface: float
        potential leaf area in m2
    
    """
    _potential_leaf_surface = {}
    for lid in _N_potential_leaf_demand:
        _potential_leaf_surface = coef_N_leaf_surface * _N_potential_leaf_demand[lid] 
    print("potential leaf surface is", _potential_leaf_surface, "cm2")
    return _potential_leaf_surface


def potential_photosynthetic_leaf_area(
    _potential_leaf_surface: float, 
    necrotic_surface: float):
    """
    Calculation of the photosynthetic leaf area (m2) 

    Parameters
    ----------
    _potential_leaf_surface: float
        potential leaf area in m2
    necrotic_surface: float
        necrose area in m2

    Returns
    -------
    _potential_photosynthetic_leaf_surface: float
        potential leaf area in m2
    
    """
    _potential_photosynthetic_leaf_area =  _potential_leaf_surface - necrotic_surface
    print("potential photosynthetic leaf surface is", _potential_photosynthetic_leaf_area, "m2")
    return _potential_photosynthetic_leaf_area


def C_potential_leaf_demand(
    potential_photosynthetic_leaf_area: float,
    Cdensity: int):
    """
    Calculation of C potential leaf demand (in gC) in function of the potential photosynthetic area and leaf N density
    
    Parameters
    ----------
    potential_photosynthetic_leaf_area: float
        potential photosynthetic area in m2
    Cdensity: float
        C shoot density in gC.cm-2
    
    Returns
    -------
    _C_potential_leaf_demand: float
        C leaves potential demand in gC
    
    """
    _C_potential_leaf_demand = potential_photosynthetic_leaf_area * Cdensity
    print("C potential leaf demand is", _C_potential_leaf_demand, "gC")
    return _C_potential_leaf_demand


def C_potential_offer(
    _producted_biomass: float,
    _C_respiration: float):
    """
    Calculation of C available to the plant (in gC), in function of C uptake (in gC) 
    
    Parameters
    ----------
    _producted_biomass: float
        quantity of producted carbon biomass 
        
    _C_respiration: 
       C respiration cost 

    Returns
    -------
    _C_potential_offer:
        Available carbon in the whole plant

    """
    _C_potential_offer = _producted_biomass - _C_respiration
    print("C available to the plant is", _C_potential_offer, "gC")
    return _C_potential_offer


def C_actual_uptake(
    _C_potential_offer: float,
    _C_potential_leaf_demand: float):
    """
    Calculation of C necessary to the plant development (in gC), in function of C offer (in g) and potential demand
    
    Parameters
    ----------
    _C_potential_offer: float 
        Available carbon in the whole plant
    _C_potential_leaf_demand: float
        C leaves potential demand in gC
    
    Returns
    -------
    _C_actual_uptake:
        Actual necessary carbon for the plant development
    """
    _C_actual_uptake = _C_potential_offer - _C_potential_leaf_demand
    return _C_actual_uptake


   
#Leaf photosynthetic system statut
def leaf_surface(
    above):
    """
    Calculation of the mtg leaf surface 

    Parameters
    ----------
    above :
    length:
    radius:
    
    Returns
    -------
    dS : float
        leaf surface
    _leaf_surface : float
        total leaves surface
    """
    g = above
    max_scale = g.max_scale()
    length = g.property('length')
    radius = g.property('radius')
    def leaf_surface(vid):
        dS = 2 * pi * length.get(vid,0)  * radius.get(vid,0)
        return dS
    
    leaves = [vid for vid in g.vertices(scale=max_scale) if g.label(vid)=='Leaf']
    _leaf_surface = sum(leaf_surface(vid) for vid in leaves)
    return _leaf_surface


def potential_leaf_light_interception(
    dS: float,
    k: int,
    Light: int):
    """
    Intercepted light by the leaf

    Parameters
    ----------
    k: int
        extinction coefficient
    dS : float
        leaf surface
    Light : int
        global incident light
    
    Returns
    -------
    leaf_light_interception: float
        MJ.m-2
    leaf_PNUE: float
        PNUE = Photosynthetic Nitrogen Use efficiency
    """
    leaf_light_interception = k*Light*exp(-k*dS)#(Bertheloot et al. 2008)
    return  leaf_light_interception


def leaf_photosynthesis(
    leaf_light_interception: float, 
    RUE: float):
    """ Transform radiation intercepted into biomass by the canopy

    Parameters
    ----------

    leaf_light_interception: float
        intercepted PAR by the leaf
    RUE: float
        Radiation use efficiency
    
    Returns
    -------
    producted_biomass: float
        producted C
    rejacted_C: float
        CO2 rejacted by photosynthesis
     
    """
    leaf_producted_biomass = RUE * leaf_light_interception
    return leaf_producted_biomass


def leaf_PNUE(
    N_potential_canopy_offer: float,
    Light: float,
    potential_leaf_light_interception: float,
    above):
    """
    Calculation of leaves N potential demand (in gN) in function of the N by light unit necessary

    Parameters
    ----------
    _N_potential_canopy_offer: float
        N shoot potential offer in gN
    Light : int
        global incident light
    potential_leaf_light_interception: float
        leaf incident light
    
    Returns
    -------
    leaf_N_potential_offer: float
        N leaf potential offer in gN
    """
    g = above
    leaves = [vid for vid in g.vertices(scale=4) if g.label(vid)=='Leaf']
    canopy_N_by_lightunit = N_potential_canopy_offer / Light
    leaf_N_potential_offer = dict((lid, canopy_N_by_lightunit*potential_leaf_light_interception[lid]) for lid in leaves)
    
    g.add_property('leaf_N_by_lightunit')
    Nps = g.property('leaf_N_by_lightunit')
    Nps.update(leaf_N_potential_offer)

    return leaf_N_potential_offer, canopy_N_by_lightunit


def Allirand_law(
    _N_potential_leaf_offer,
    leaf_light_interception):
    """
    Allirand law ; relationship between the light intercepted by the leaf and the nitrogen content

    Parameters
    ----------
    _N_potential_leaf_offer: float
        N leaf potential offer in gN
    leaf_light_interception : float
        in MJ.cm-2
    
    Returns
    -------
    Allirand_cst: float
        relationship between the light intercepted by the leaf and the nitrogen content
    """
    Allirand_cst = leaf_light_interception*_N_potential_leaf_offer
    return Allirand_cst







#Leaf Growth statut
def actual_growth(
    _C_actual_uptake: float):
    """
    Calculation of the actual growth of the leaf
    
    Parameters
    ----------
    _C_actual_uptake:
        Actual necessary carbon for the plant development

    Returns
    -------
    _C_actual_uptake:
        Actual necessary carbon for the plant development
    """
    
    if _C_actual_uptake == 0:
        C_actual_uptake = _C_actual_uptake
        print("The leaf is mature. It doesn't import carbon anymore")
    if _C_actual_uptake < 0:
        C_actual_uptake = _C_actual_uptake
        print("The growth is slowed by a carbon deficit")
    else :
        C_actual_uptake = _C_actual_uptake
        print("The leaf is growing and is storaging starch")
    
    return C_actual_uptake

    

#Turnover Budget

if _C_actual_uptake < 0:

def C_import(
    C_potential_leaf_demand,
    C_potential_offer,
    C_source_available,
    C_actual_available,
    leaf_N_potential_offer,
    dC):

dC = C_potential_leaf_demand - C_potential_offer 
C_actual_available = C_source_available - dC

if _C_actual_available == 0:
        C_actual_available = _C_actual_available
        print("The C deficit is compensated thanks to importation but the C storage from the source organ is now empty")
        CN_ratio = C_potential_leaf_demand / leaf_N_potential_offer
if _C_actual_available < 0:
        C_actual_available = _C_actual_available
        print("The importation isn't enougth to face the C deficit")
        CN_ratio = (C_potential_offer + C_source_available) / leaf_N_potential_offer
        print("The actual C/N ratio initiate senescence")
else :
    C_actual_available = _C_actual_available
    print("The C deficit is compensated thanks to importation. The source organ can still storage C")
    CN_ratio = C_potential_leaf_demand / leaf_N_potential_offer

    return C_actual_available, CN_ratio




"""
if _C_actual_available < 0:

def senescence()

"""









