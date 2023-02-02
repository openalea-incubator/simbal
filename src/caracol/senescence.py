# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 10:44:53 2022

@author: FARRA
"""

from math import (exp, pi)


#Nitrogen Budget
def potential_N_root_uptake(
    root_surface: float, 
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
    _potential_N_root_uptake = SNU * Nitrate * time_step * root_surface
    print("N root uptake is", _potential_N_root_uptake, "gN")
    return _potential_N_root_uptake


def N_potential_offer(
    _potential_N_root_uptake: float):
    """
    Calculation of N available to the whole plant (in gN)
    
    Parameters
    ----------
    _potential_N_root_uptake: float
        N uptake by roots (in g)   

    Returns
    -------
    _N_potential_offer: float
        N potential offer to the plant in gN
    """
    _N_potential_offer = _potential_N_root_uptake 
    print("N available to the whole plant is", _N_potential_offer, "gN")
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


def N_potential_plant_offer(
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
    _N_potential_plant_offer: float
        N shoot potential demand in gN
    """
    _N_potential_plant_offer = _N_potential_offer - _N_potential_root_demand
    print("N potential shoot demand is", _N_potential_plant_offer, "gN")
    return _N_potential_plant_offer



#Light Budget
def leaf_actual_surface(
    above,
    vid):
    """
    Calculation of the mtg leaf surface 

    Parameters
    ----------
    above :
    length:
    radius:
    
    Returns
    -------
    _leaf_actual_surface : float
        leaf surface
    """
    g = above
    max_scale = g.max_scale()
    length = g.property('length')
    radius = g.property('radius')

    leaves = [vid for vid in g.vertices(scale=max_scale) if g.label(vid)=='Leaf']
    for vid in leaves:
        _leaf_actual_surface = 2 * pi * length.get(vid,0) * radius.get(vid,0)
        return _leaf_actual_surface
    



def photosynthetic_leaf_surface(
    leaf_actual_surface: float, 
    leaf_necrotic_surface: float):
    """
    Calculation of the photosynthetic leaf area (m2) 

    Parameters
    ----------
    leaf_actual_surface: float
        leaf area in m2
    leaf_necrotic_surface: float
        necrose area in m2

    Returns
    -------
    _potential_photosynthetic_leaf_surface: float
        potential leaf area in m2
    
    """
    _photosynthetic_leaf_surface =  leaf_actual_surface - leaf_necrotic_surface
    print("potential photosynthetic leaf surface is", _photosynthetic_leaf_surface, "m2")
    return _photosynthetic_leaf_surface


def potential_leaf_light_interception(
    _photosynthetic_leaf_surface: float,
    k: int,
    Light: int):
    """
    Intercepted light by the leaf

    Parameters
    ----------
    k: int
        extinction coefficient
    _photosynthetic_leaf_surface : float
        leaf surface
    Light : int
        global incident light
    
    Returns
    -------
    _potential_leaf_light_interception: float
        MJ.m-2
    leaf_PNUE: float
        PNUE = Photosynthetic Nitrogen Use efficiency
    """
    _potential_leaf_light_interception = k*Light*exp(-k*_photosynthetic_leaf_surface)#(Bertheloot et al. 2008)
    return  _potential_leaf_light_interception


def leaf_photosynthesis(
    _leaf_intercepted_PAR: float, 
    RUE: float):
    """ Transform radiation intercepted into biomass by the leaf

    Parameters
    ----------

    leaf_intercepted_PAR: float
        intercepted PAR by the plant
    RUE: float
        Radiation use efficiency
    
    Returns
    -------
    leaf_producted_biomass: float
        producted C
    rejacted_C: float
        CO2 rejacted by photosynthesis
     
    """
    _leaf_producted_biomass = RUE * _leaf_intercepted_PAR
    return _leaf_producted_biomass


def leaf_N_potential_offer(
    _N_potential_plant_offer: float,
    Light: float,
    _potential_leaf_light_interception: float,
    above):
    """
    Calculation of leaves N potential demand (in gN) in function of the N by light unit necessary

    Parameters
    ----------
    _N_potential_plant_offer: float
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
    _plant_N_by_lightunit = _N_potential_plant_offer / Light
    _leaf_N_potential_offer = dict((lid, _plant_N_by_lightunit*_potential_leaf_light_interception[lid]) for lid in leaves)
    
    g.add_property('leaf_N_by_lightunit')
    Nps = g.property('leaf_N_by_lightunit')
    Nps.update(_leaf_N_potential_offer)

    return _leaf_N_potential_offer, _plant_N_by_lightunit


def potential_leaf_surface(
    _leaf_N_potential_offer: float, 
    coef_N_leaf_surface: int):
    """
    Calculation of potential formed leaf area  (cm2) in function of potential N leaf offer (gN) and the coefficient of the relationship of leaf area - leaf N quantity (cm2.g-1)

    Parameters
    ----------
    leaf_N_potential_offer: float
        leaf potential N demand in gN
    coef_N_leaf_surface: int
        coefficient between leaf area and leaf N quantity in m2.g-1

    Returns
    -------
    _potential_leaf_surface: float
        potential leaf area in m2
    
    """
    _potential_leaf_surface = {}
    for lid in _leaf_N_potential_offer:
        _potential_leaf_surface = coef_N_leaf_surface * leaf_N_potential_offer[lid] 
    print("potential leaf surface is", _potential_leaf_surface, "cm2")
    return _potential_leaf_surface


def FSPM_Allirand_law(
    _leaf_N_potential_offer,
    _potential_leaf_light_interception):
    """
    Allirand law ; relationship between the light intercepted by the leaf and the nitrogen content

    Parameters
    ----------
    _N_potential_leaf_offer: float
        N leaf potential offer in gN
    _potential_leaf_light_interception : float
        in MJ.cm-2
    
    Returns
    -------
    _FSPM_Allirand_cst: float
        relationship between the light intercepted by the leaf and the nitrogen content
    """
    _FSPM_Allirand_cst = _potential_leaf_light_interception*_leaf_N_potential_offer
    return _FSPM_Allirand_cst



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

def C_potential_leaf_demand(
    _potential_leaf_surface: float,
    Cdensity: int):
    """
    Calculation of C potential leaf demand (in gC) in function of the potential photosynthetic area and leaf N density
    
    Parameters
    ----------
    _potential_leaf_surface: float
        potential formed area in m2
    Cdensity: float
        C shoot density in gC.cm-2
    
    Returns
    -------
    _C_potential_leaf_demand: float
        C leaves potential demand in gC
    
    """
    _C_potential_leaf_demand = _potential_leaf_surface * Cdensity
    print("C potential leaf demand is", _C_potential_leaf_demand, "gC")
    return _C_potential_leaf_demand


def C_potential_leaf_offer(
    _leaf_producted_biomass: float,
    _C_respiration: float):
    """
    Calculation of C available to the leaf (in gC), in function of C uptake (in gC) 
    
    Parameters
    ----------
    leaf_producted_biomass: float
        quantity of producted carbon biomass 
        
    _C_respiration: 
       C respiration cost 

    Returns
    -------
    _C_potential_leaf_offer:
        Available carbon in the whole plant

    """
    _C_potential_leaf_offer = _leaf_producted_biomass - _C_respiration
    print("C available to the leaf is", _C_potential_leaf_offer, "gC")
    return _C_potential_leaf_offer



#Turnover Budget
def C_actual_uptake(
    g, 
    vid,
     _C_potential_leaf_offer: float,
     _C_potential_leaf_demand: float): 

    """
    Calculation of the actual growth by stress of the leaf
    
    Parameters
    ----------
    _C_potential_leaf_demand:
        Actual necessary carbon for the plant development
     _C_actual_uptake:

    g :
        MTG

    vid :
        Vertices id

    Returns
    -------
    dC:
        Actual necessary carbon for the plant development
     _C_actual_uptake:

    storage:

    growth:

    """  
    dC = _C_potential_leaf_offer - _C_potential_leaf_demand 
    if dC < 0 :
        _C_actual_uptake = dC 
        storage = 0

        def C_actual_uptake_import(
            _C_source_offer,
            dC):
    
            dA = dC + _C_source_offer
            if dA < 0:
             storage = 0
             growth = False

             def leaf_senescence(
                _C_potential_leaf_offer,
                _leaf_N_potential_offer,
                _C_sink_demand,
                _N_sink_demand,
                w):

                _C_source_offer = w*_C_potential_leaf_offer*_C_sink_demand
                _N_source_offer = w*_leaf_N_potential_offer*_N_sink_demand
                if _C_source_offer <= 0 :
                    leaf_fall = True

                if _N_source_offer <= 0 :
                    leaf_fall = True

                else :
                    leaf_fall = False
                
                return _C_source_offer, _N_source_offer, leaf_fall

            else :
                storage = dA
                _C_actual_uptake = _C_potential_leaf_demand
                growth = True
                
            return  _C_actual_uptake, growth
        
    else :
        _C_actual_uptake = _C_potential_leaf_demand
        storage = dC
        growth = True

    return _C_actual_uptake, growth
    
    








