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
    root_surface: float
        root object
    Nitrate: float
        N quantity in soil in gN
    time_step: float
        time step of the simulation in day
    SNU:


    Returns
    -------
    potential_N_root_uptake: float
        potential N uptake by roots
    """
    _potential_N_root_uptake = SNU * Nitrate * time_step * root_surface
    print("N root uptake is", _potential_N_root_uptake, "gN")
    return _potential_N_root_uptake


def N_potential_offer(
    N_potential_root_uptake: float):
    """
    Calculation of N available to the whole plant (in gN)
    
    Parameters
    ----------
    _potential_N_root_uptake: float
        N uptake by roots (in g)   

    Returns
    -------
    N_potential_offer: float
        N potential offer to the plant in gN
    """
    _N_potential_offer = N_potential_root_uptake 
    print("N available to the whole plant is", _N_potential_offer, "gN")
    return _N_potential_offer


def N_potential_root_demand(
    delta_volume: float, 
    N_density: int):
    """
    Calculation of N potential root demand (in gN) in function of the root volume and root N density

    Parameters
    ----------
    delta_volume: float
        potential increase of root volume in m3
    N_density: float
        N root density in gN.m-3

    Returns
    -------
    N_potential_root_demand: float
        N root potential demand in gN
    """
    _N_potential_root_demand = delta_volume * N_density
    print("N potential root demand is", _N_potential_root_demand, "gN")
    return _N_potential_root_demand


def N_potential_plant_offer(
    N_potential_shoot_offer: float, 
    N_potential_root_demand: float):
    """
    Calculation of the potential shoot demand (gN) in function of N available to the plant (gN) and potential root N demand (gN)

    Parameters
    ----------
    N_potential_shoot_offer: float
        N offer to the plant in gN
    N_potential_root_demand: float
        N potential root N demand in gN
    
    Returns
    -------
    N_potential_plant_offer: float
        N shoot potential demand in gN
    """
    _N_potential_plant_offer = N_potential_shoot_offer - N_potential_root_demand
    print("N potential shoot demand is", _N_potential_plant_offer, "gN")
    return _N_potential_plant_offer



#Light Budget
def leaf_actual_surface(
    above, vid):
    """
    Calculation of the mtg leaf surface 

    Parameters
    ----------
    above :
    length:
    radius:
    
    Returns
    -------
    _plant_actual_surface : float
       total shoot surface
    """
    g = above
    max_scale = g.max_scale()
    length = g.property('length')
    radius = g.property('radius')

    def _leaf_actual_surface(vid):
        dS = 2 * pi * length.get(vid,0)  * radius.get(vid,0)
        return dS
    
    leaves = [vid for vid in g.vertices(scale=max_scale) if g.label(vid)=='Leaf']
    _plant_actual_surface = sum(_leaf_actual_surface(vid) for vid in leaves)
    return _plant_actual_surface

def photosynthetic_leaf_surface(
    leaf_surface: float, 
    leaf_necrotic_surface: float):
    """
    Calculation of the photosynthetic leaf area (m2) 

    Parameters
    ----------
    leaf_surface: float
        leaf area in m2
    leaf_necrotic_surface: float
        necrose area in m2

    Returns
    -------
    potential_photosynthetic_leaf_surface: float
        potential leaf area in m2
    """
    _photosynthetic_leaf_surface =  leaf_surface - leaf_necrotic_surface
    print("potential photosynthetic leaf surface is", _photosynthetic_leaf_surface, "m2")
    return _photosynthetic_leaf_surface


def potential_leaf_light_interception(
    photosynthetic_leaf_surface: float,
    k: int,
    Light: int):
    """
    Intercepted light by the leaf

    Parameters
    ----------
    photosynthetic_leaf_surface : float
        leaf surface
    k: int
        extinction coefficient
    Light : int
        global incident light
    
    Returns
    -------
    potential_leaf_light_interception: float
        MJ.m-2
    
    References
    ----------
    Bertheloot et al. 2008
    """
    _potential_leaf_light_interception = k * Light * exp(-k * photosynthetic_leaf_surface)
    
    return  _potential_leaf_light_interception


def leaf_photosynthesis(
    potential_leaf_light_interception: float, 
    RUE: float):
    """ Transform radiation intercepted into biomass by the leaf

    Parameters
    ----------

    potential_leaf_light_interception: float
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
    _leaf_producted_biomass = RUE * potential_leaf_light_interception
    return _leaf_producted_biomass


def leaf_N_potential_offer(
    g,
    N_potential_plant_offer: float,
    Light: float,
    potential_leaf_light_interception: float
    ):
    """
    Calculation of leaves N potential demand (in gN) in function of the N by light unit necessary

    Parameters
    ----------
    g : MTG
        Shoot/above architecture
    N_potential_plant_offer: float
        N shoot potential offer in gN
    Light : float
        global incident light
    potential_leaf_light_interception: float
        leaf incident light
    
    Returns
    -------
    N_leaf_potential_offer: float
        N leaf potential offer in gN
    """
    max_scale = g.max_scale()
    leaves = [vid for vid in g.vertices(scale=max_scale) if g.label(vid)=='Leaf']
    _plant_N_by_lightunit = N_potential_plant_offer / Light
    _leaf_N_potential_offer = dict(
        (lid, _plant_N_by_lightunit * potential_leaf_light_interception[lid]) 
        for lid in leaves
    )
    
    # TODO : define the initial properties once in an init function?
    g.properties().setdefault('leaf_N_by_lightunit', {}).update(_leaf_N_potential_offer)
    #g.add_property('leaf_N_by_lightunit')
    #Nps = g.property('leaf_N_by_lightunit')
    #Nps.update(_leaf_N_potential_offer)

    return _leaf_N_potential_offer, _plant_N_by_lightunit


def potential_leaf_surface(
    N_leaf_potential_offer: dict, 
    coef_N_leaf_surface: float):
    """
    Calculation of potential formed leaf area  (cm2) in function of potential N leaf offer (gN) and the coefficient of the relationship of leaf area - leaf N quantity (cm2.g-1)

    Parameters
    ----------
    N_leaf_potential_offer: dict
        leaf potential N demand in gN
    coef_N_leaf_surface: float
        coefficient between leaf area and leaf N quantity in m2.g-1

    Returns
    -------
    _potential_leaf_surface: dict
        potential leaf area in m2
    """
    _potential_leaf_surface = {}
    for lid in N_leaf_potential_offer:
        _potential_leaf_surface[lid] = coef_N_leaf_surface * N_leaf_potential_offer[lid] 
    print("potential leaf surface is", _potential_leaf_surface, "cm2")
    return _potential_leaf_surface


def FSPM_Allirand_law(
    N_leaf_potential_offer,
    potential_leaf_light_interception):
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
    FSPM_Allirand_cst: float
        relationship between the light intercepted by the leaf and the nitrogen content
    """
    _FSPM_Allirand_cst = potential_leaf_light_interception * N_leaf_potential_offer
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
    Cdensity: float):
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



##Turnover Budget
# internal functions
def leaf_statut(
    C_potential_leaf_offer: float,
    C_potential_leaf_demand: float):
    """
    Evaluation of the leaf state to determinate his value (weigth) in the pool 
    
    Parameters
    ----------
    C_potential_leaf_offer: float
        Available carbon in the whole plant    
    C_potential_leaf_demand: 
       C leaves potential demand in gC 

    Returns
    -------
    growth : bool
        growth state
    sink : bool
        organe state
    source : bool
        organe state
    storage : bool
        C storage value
    dC : float
        Carbon deficit in the leaf
    """
    dC = C_potential_leaf_offer - C_potential_leaf_demand 

    if dC < 0:
        storage = False
        growth = False
        sink = True
        source = False
    
    elif dC == 0:
        storage = False
        growth = True
        sink = False
        source = False
    
    else:
        storage = True
        growth = True
        sink = False
        source = True

    return growth, sink, source, storage, dC

 
def C_pool_budget(
    C_potential_leaf_offer,
    w):
    """
    Available C in the plant
    
    Parameters
    ----------
    C_potential_leaf_offer : float
        Available carbon in the leaf  
    w : float
        leaf ponderation (weight) in the pool for export

    Returns
    -------
    C_pool_offer : float
        Total C available in the plant
    C_leaf_export : float
        C exported to other organs    
    empty : bool
        pool state
    dA : float
        Carbon deficit after pool budget
    """
    #source = leaf_statut()
    dC = leaf_statut()

    #if source = True
    
    #w =#pondération 
    _C_pool_offer = C_potential_leaf_offer #weighted leaf sources based on relative PAR
    _C_leaf_export = w * _C_pool_offer #update

    if _C_pool_offer == 0 :
        empty = True

    else :
        empty = False
        dA = dC + _C_pool_offer
    
    return _C_pool_offer, _C_leaf_export, empty, dA


def C_exportation_budget(
    C_potential_leaf_offer,
    C_leaf_export):

    _C_leaf_actual_content = C_potential_leaf_offer - C_leaf_export#update feuilles sources

    return _C_leaf_actual_content


def C_importation_budget(
    C_potential_leaf_demand,
    C_potential_leaf_offer,
    N_potential_leaf_offer,
    dC, w):
    """
    Evaluation of the stress (deficit) after importation 
    
    Parameters
    ----------
    C_potential_leaf_demand: float
        C leaves potential demand in gC
    C_potential_leaf_offer : float
        Available carbon in the whole plant    
    N_potential_leaf_offer : float
        leaf potential N demand in gN
    dC : float
         Carbon deficit in the leaf
    w : float
        leaf weigth in the pool according to his PAR and statut

    Returns
    -------
    stress : bool
        stress (nutrient deficit) state
    C_actual_deficit : float
        Final deficit after importation
    C_actual_uptake : float
        Final imported carbon quantity
    storage : float
        C storage value
    growth : bool
        growth state
    senescence_initiation : bool
        senescence initiation indicator
    """
    dA = C_pool_budget(
        C_potential_leaf_offer,
        N_potential_leaf_offer,
        w, dC)

    if dA < 0 :
        stress = True
        _C_actual_deficit = dA 
        storage = 0
        growth = False
        senescence_initiation = True

    else :
        stress = False
        _C_actual_uptake = C_potential_leaf_demand
        storage = dA
        growth = True

    return stress, _C_actual_deficit, _C_actual_uptake, storage, growth, senescence_initiation


def leaf_senescence(g, vid,
    N_potential_plant_offer,
    plant_actual_surface,
    C_potential_leaf_offer,
    N_potential_leaf_offer,
    dS, w):

    #sink = leaf_statut()
    #if not sink : initatiation de la feuille avec le PAR le moins élevé
    #update g.property senescence
    #if feuille senescente
    senescence_stage = True
    N_remobilisation = True

    _C_leaf_export = C_pool_budget(
        C_potential_leaf_offer,
        w)

    # TODO : _N_leaf_actual_content is not defined
    _N_senescence_threshold = _N_leaf_actual_content
    _C_senescence_threshold = C_potential_leaf_offer - _C_leaf_export

    #g = above
    max_scale = g.max_scale()
    
    leaves = [vid for vid in g.vertices(scale=max_scale) if g.label(vid)=='Leaf']

    leaf_area = g.property('leaf_area')
    _N_leaf_actual_content = dict(
        (vid, (N_potential_plant_offer / plant_actual_surface)*leaf_area[vid])
        for vid in leaves
    )  

    _N_pool_offer = sum(_N_leaf_actual_content[vid] for vid in leaves)# des feuilles en sénescence

    Np_leaf_offer = g.property('N_potential_leaf_offer')
    _N_plant_demand = sum(_N_leaf_actual_content[vid] - Np_leaf_offer[vid] for vid in leaves)# des feuilles puits

    _N_plant_hunger = _N_pool_offer - _N_plant_demand

    if _N_plant_hunger > 0 :
       
        _N_remobilized = _N_plant_demand
       
        _N_leaf_remobilized = (w * _N_remobilized)# feuille sénescence w le plus élevé selon PAR élevé

        _N_leaf_actual_content = _leaf_N_potential_offer #update feuille puit
        
        _N_leaf_actual_content = _N_leaf_actual_content -_N_leaf_remobilized #update feuille senescent

    else :
        N_stress = True
        _N_remobilized = _N_pool_offer

        leaf_fall = True#update property senescent = leaf_fall
        _N_leaf_fall = _N_leaf_actual_content

        _N_leaf_actual_content = _N_leaf_actual_content + (w * _N_remobilized)# pondération des feuilles puits selon PAR le plus élevé
       
    return senescence_stage, N_remobilisation, _N_senescence_threshold, _C_senescence_threshold, _N_plant_hunger, _N_remobilized, leaf_fall, N_stress, _N_leaf_actual_content, _N_leaf_fall


#main function
def turnover_budget(
    g, vid: int,
    C_potential_leaf_offer: float,
    C_potential_leaf_demand: float,
    N_potential_leaf_offer: float,
    w): 
    """
    Calculation of the actual growth by stress of the leaf
    
    Parameters
    ----------
    g : MTG
        above MTG
    vid : int
        Vertices id
    C_potential_leaf_offer: float
        Actual available carbon for the plant development
    C_potential_leaf_demand: float
        Actual necessary carbon for the plant development
    N_potential_leaf_offer : float
        N potential leaf offer for a given leaf
    w : float
        Weight ???

    Returns
    -------
    senescence_initiation : bool
        senescence initiation indicator
    senescence_threshold = bool
        senescence threshold indicator
    leaf_fall : bool
        viable or non viable state
    stress : bool
        stress (nutrient deficit) state
    _C_actual_uptake : float
        Final imported carbon quantity
    growth : bool
        growth state
    empty : bool
        pool state
    source : bool
        organe state
    storage : float
        C storage value

    Algorithm
    ---------
    Let g be the shoot.
    1. Query all the leaves that are gowing
    2. Query all the others
    3. Is the plant is C stressed
        * get C available in all organs (C_pool_budget)
        * distribute depending on w
    4. Else
        *  continuer

    """  
    # TODO : for all leaves
    growth, sink, source, storage, dC = leaf_statut(
        C_potential_leaf_offer,
        C_potential_leaf_demand)

    if not growth:
        # Missing : g
        leaf_offer, empty, dA = C_pool_budget(
            C_potential_leaf_offer,
            N_potential_leaf_offer,
            w)

        if not empty:

            stress = C_importation_budget(
                C_potential_leaf_offer,
                C_potential_leaf_demand,
                N_potential_leaf_offer,
                dC, w)
                

            if not stress:
                _C_actual_uptake = _C_potential_leaf_demand

            else:
                _senescence_initiation = True
                _N_senescence_threshold, _C_senescence_threshold, _leaf_fall, _N_leaf_fall = leaf_senescence()

            return _senescence_initiation, _N_senescence_threshold, _C_senescence_threshold, _leaf_fall, _N_leaf_fall  

        else : # empty == True
            _senescence_initiation = True
            _N_senescence_threshold, _C_senescence_threshold, leaf_fall, _N_leaf_fall = leaf_senescence()

        return stress, _C_actual_uptake, _senescence_initiation, _N_senescence_threshold, _C_senescence_threshold, leaf_fall, _N_leaf_fall
    
    else : # growth == True
        source, storage = leaf_statut(
            _C_potential_leaf_offer,
            _C_potential_leaf_demand)
        
        _C_actual_uptake = _C_potential_leaf_demand

    return growth, empty,_N_senescence_threshold, _C_senescence_threshold, source, storage, _C_actual_uptake

            
                
                
















   




