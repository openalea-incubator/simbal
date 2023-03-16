# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 10:44:53 2022

@author: FARRA
"""

from math import (exp, pi)

#Photosynthesis
def N_potential_shoot_offer(g,
    N_potential_plant_offer: float, 
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
        N shoot potential offer in gN
    
    Hypothesis
    ----------
    H1 Roots take precedence on the nitrogen uptake
    """
    max_scale = g.max_scale()

    leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) if g.class_name(lid)=='Leaf']
    alive = [l for l in leaves if not l.leaf_fall and l.senescence]
    senescence = [l for l in alive if l.senescence]

    N_shoot = N_potential_plant_offer - N_potential_root_demand
    PAR_total = sum(l.PAR for l in leaves)

    for l in alive:
        l.weight = l.PAR / PAR_total
        l.N_available = l.N * l.weight
    
    for l in senescence:
        l.N_available = l.N
        
    N_pool = sum(l.N_available for l in alive) 
    N_potential_shoot_offer = N_shoot + N_pool
    
    return PAR_total,N_pool, N_potential_shoot_offer


def leaf_area(g):
    """
    Calculation of the mtg leaf surface 

    Parameters
    ----------
    g : MTG
        the above part
    
    Returns
    -------
    g : MTG
        the above part

    State Variables (on MTG)
    ------------------------
    length : float
    radius : float
    area : float
        leaf surface acording to length and radius
    """
    max_scale = g.max_scale()

    length = g.property('length')
    radius = g.property('radius')

   

    leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) if g.class_name(lid)=='Leaf']
    alive = [l for l in leaves if not l.leaf_fall]

    for l in alive :
        l.area = 2 * pi * l.length  * l.radius
    
    return g
    

def photosynthetic_leaf_area(g):
    """
    Calculation of the photosynthetic leaf area (m2) 

    Parameters
    ----------
    g : MTG
        the above part
    
    Returns
    -------
    g : MTG
        the above part

    State Variables (on MTG)
    ------------------------
    area : float
        leaf surface acording to length and radius
    necrotic_area : float
        non photosynthetic leaf area
    photosynthetic_area : float
        green leaf area
    
    Hypothesis
    ----------
    H2 leaf necrotic surface doesn't have photosynthetic activities anymore
    """
    max_scale = g.max_scale()

    leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) if g.class_name(lid)=='Leaf']
    alive = [l for l in leaves if not l.leaf_fall]

    for l in alive:
        l.photosynthetic_area = l.area - l.necrotic_area

    return g


def potential_leaf_PAR(g,
    k: int,
    e_amax: float):
    """
    Intercepted light by the leaf

    Parameters
    ----------
    g : MTG
        the above part
    
    Returns
    -------
    g : MTG
        the above part

    State Variables (on MTG)
    ------------------------
    photosynthetic_area : float
        green leaf area
    PAR : float
        leaf intercepted light

    Hypothesis
    ----------
    H3 The Monteith law is valid at the leaf scale (with MTG)
    
    References
    ----------
    Bertheloot et al. 2008
    """
    max_scale = g.max_scale()

    leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) if g.class_name(lid)=='Leaf']
    alive = [l for l in leaves if not l.leaf_fall]

    for l in alive:
        l.PAR = k * e_amax * (1 - exp(k * l.photosynthetic_area))
  
    return  g


def leaf_photosynthesis(g,
    RUE: float):
    """ Transform radiation intercepted into biomass by the leaf

   Parameters
    ----------
    g : MTG
        the above part
    RUE : float
        ?????   
    
    Returns
    -------
    g : MTG
        the above part

    State Variables (on MTG)
    ------------------------
    PAR : float
        leaf intercepted light
    C_producted : float
        leaf producted biomass by photosynthesis 
    
    Hypothesis
    ----------
    H3 The Monteith law is valid at the leaf scale (with MTG)
    """
    max_scale = g.max_scale()

    leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) if g.class_name(lid)=='Leaf']
    alive = [l for l in leaves if not l.leaf_fall]

    for l in alive:
        l.C_producted = RUE * l.PAR

    return g


def leaf_N_potential_offer(g,
    PAR_total,
    N_potential_shoot_offer):
    """
    Calculation of leaves N potential demand (in gN) in function of the N by light unit necessary

    Parameters
    ----------
    g : MTG
        the above part
    N_potential_plant_offer : float
        N shoot potential offer in gN

    Returns
    -------
    g : MTG
        the above part
    PAR_total : float
        total leaves PAR
    plant_N_by_lightunit : float
        nitrogen light use efficiency ?

    State Variables (on MTG)
    ------------------------
    PAR : float
        leaf intercepted light
    N_potential_offer : float
        nitrogen potential leaf offer
    
    Hypothesis
    ----------
    H4 Nitrogen plant distribution accordind to PAR (light)

    H5 Nutrient contents and PAR are homogeneous in the leaf area

    References
    ----------
    Takashima (2004)
    Hikosaka (2005)
    Yasumura (2007)
    """
    max_scale = g.max_scale()

    leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) if g.class_name(lid)=='Leaf']
    alive = [l for l in leaves if not l.leaf_fall]

    plant_N_by_lightunit = N_potential_shoot_offer / PAR_total

    for l in alive:
        l.N_potential_offer = plant_N_by_lightunit * l.PAR
   

    return g, PAR_total, plant_N_by_lightunit


def potential_leaf_surface(g):
    """
    Calculation of potential formed leaf area  (cm2) in function of potential N leaf offer (gN) and the coefficient of the relationship of leaf area - leaf N quantity (cm2.g-1)

    Parameters
    ----------
    g : MTG
        the above part

    Returns
    -------
    g : MTG
        the above part

    State Variables (on MTG)
    ------------------------
    N : float
        nitrogen leaf content
    area : float
        leaf surface
    area_by_Nunit : float
        surface nitrogen use efficiency ?
    N_potential_offer : float
        nitrogen potential leaf offer
    potential_area : float
        potential leaf surface
    
    Hypothesis
    ----------
    H5 Nutrient contents and PAR are homogeneous in the leaf area
    """
    max_scale = g.max_scale()

    leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) if g.class_name(lid)=='Leaf']
    alive = [l for l in leaves if not l.leaf_fall]

    for l in alive:
        l.area_by_Nunit = l.photosynthetic_area / l.N
        l.potential_area = l.area_by_Nunit * l.N_potential_offer
        
    return g

def potential_leaf_light_interception(g, time_step,
    k: int,
    e_amax: float):
    """
    Intercepted light by the leaf

   Parameters
    ----------
    g : MTG
        the above part
    k: int
        extinction coefficient
    Light : int
        global incident light

    Returns
    -------
    g : MTG
        the above part
    
    State Variables (on MTG)
    ------------------------
    potential_area : float
        potential leaf surface
    PAR_potential : float
        potential leaf intercepted light
    
    Hypothesis
    ----------
    H3 The Monteith law is valid at the leaf scale (with MTG) 
    
    References
    ----------
    Bertheloot et al. 2008
    """
    max_scale = g.max_scale()

    leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) if g.class_name(lid)=='Leaf']
    alive = [l for l in leaves if not l.leaf_fall]

    for l in alive:
        l.PAR_potential = k * e_amax * (1 - exp(k * l.potential_area))
  
    return  g


def leaf_potential_photosynthesis(g,
    RUE: float):
    """ Transform radiation intercepted into biomass by the leaf

    Parameters
    ----------
    g : MTG
        the above part
    RUE : float
        radiation use efficiency ?

    Returns
    -------
    g : MTG
        the above part
    
    State Variables (on MTG)
    ------------------------
    PAR_potential : float
        potential leaf intercepted light
    C_potential_producted : float
        potential leaf producted biomass
    
    Hypothesis
    ----------
    H3 The Monteith law is valid at the leaf scale (with MTG)
    """
    max_scale = g.max_scale()

    leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) if g.class_name(lid)=='Leaf']
    alive = [l for l in leaves if not l.leaf_fall]

    for l in alive:
        l.C_potential_producted = RUE * l.PAR_potential

    return g


def FSPM_Allirand_law(g):
    """
    Allirand law ; relationship between the light intercepted by the leaf and the nitrogen content

    Parameters
    ----------
    g : MTG
        the above part

    Returns
    -------
    g : MTG
        the above part

    State Variables (on MTG)
    ------------------------
    PAR_potential : float
        potential leaf intercepted light
    N_potential_offer : float
        nitrogen potential leaf offer
    FSPM_Allirand : float
        relationship between the light intercepted by the leaf and the nitrogen content
    
    Hypothesis
    ----------
    H6 The Allirand law is valid at the leaf scale (with MTG)

    References
    ----------
    Allirand (2009)
    """
    max_scale = g.max_scale()

    leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) if g.class_name(lid)=='Leaf']
    alive = [l for l in leaves if not l.leaf_fall]

    for l in alive:
        l.FSPM_Allirand = l.PAR_potential * l.N_potential

    return g


def C_potential_leaf_demand(g):
    """
    Calculation of C potential leaf demand (in gC) in function of the potential photosynthetic area and leaf N density
    
    Parameters
    ----------
    g : MTG
        the above part

    Returns
    -------
    g : MTG
        the above part
    
    State Variables (on MTG)
    ------------------------
    area : float
        leaf surface
    potential_are : flooat
        potential leaf surface
    C_by_areaunit : float
        carbon surface use efficiency ?
    C_potential_demand : float
        leaf potential carbon demand
    
    Hypothesis
    ----------
    H5 Nutrient contents and PAR are homogeneous in the leaf area
    """
    max_scale = g.max_scale()

    leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) if g.class_name(lid)=='Leaf']
    alive = [l for l in leaves if not l.leaf_fall]

    for l in alive:
        l.C_by_areaunit = l.C / l.area
        l.C_potential_demand = l.C_by_areaunit * l.potential_area
        
    return g


def C_potential_leaf_offer(g):
    """
    Calculation of C available to the leaf (in gC), in function of C uptake (in gC) 
    
    Parameters
    ----------
    g : MTG
        the above part

    Returns
    -------
    g : MTG
        the above part
    
    State Variables (on MTG)
    ------------------------
    C_producted : float
        leaf producted biomass
    C_potential_offer : float
        leaf potential carbon offer
    """
    max_scale = g.max_scale()

    leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) if g.class_name(lid)=='Leaf']
    alive = [l for l in leaves if not l.leaf_fall]

    for l in alive:
        l.C_potential_offer = l.C + l.C_producted - l.C_respiration         

    return g



##Turnover Budget
# internal functions
def leaf_status(g):
    """
    Evaluation of the leaf state to compute his value (weigth) in the pool 
    
    Parameters
    ----------
    g: MTG
        the above part

    Returns
    -------
    g : MTG
        the above part

    State Variables (on MTG)
    ------------------------
    PAR : float
        leaf light intercepted
    C_potential_leaf_offer : float
        leaf potential C content
    C_potential_leaf_demand : float
        leaf potential C demand
    dC : float
        C deficit in the leaf
    growth : bool
        growth state
    sink : bool
        organe state
    source : bool
        organe state
    storage : bool
        C storage value
    
    Hypothesis
    ----------
    H7 Source/sink status according to leaves c stress

    H8 A leaf cannot be a source and a sink organ at the same time step

    References
    ----------
    Khudairi (1970)
    Wingler (2004)
    """
    max_scale = g.max_scale()

    leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) if g.class_name(lid)=='Leaf']

    for l in leaves :
        l.dC = l.C_potential_leaf_offer - l.C_potential_leaf_demand

        if l.dC < 0:
            l.storage = False
            l.growth = False
            l.sink = True
            l.source = False
    
        else:
            l.storage = True
            l.growth = True
            l.sink = False
            l.source = True

    return g

 

def C_pool_budget(g):
    """
    C distribution between leaves acording to their status and value (weight)
    
    Parameters
    ----------
    g: MTG
        the above part

    Returns
    -------
    g : MTG
        the above part
    pool_C : float
        C labile available in the plant
    PAR_total :
        total light intercepted by the plant
    C_plant_hunger :
        C plant stress state

    State Variables (on MTG)
    ------------------------
    sink : bool
        organe state
    leaf_fall : bool
        viability state of the leaf
    C_storage : float
        ????
    C_producted : float
        leaf producted biomass
    PAR : float
        leaf light intercepted
    weight : float
        value of the leaf participation to the pool accordint to his PAR
    C : float
        leaf carbon content
    C_exported : float
        C exported from the source leaf
    C_imported : float
        C imported to the sink leaf
    C_deficit : float
        leaf C stress
    growth : bool
        growth state
    C_potential_demand : float
        carbon potential leaf demand
    
    Hypothesis
    ----------
    H9 leaf participation to the C pool according to their PAR (higher PAR = greater weight)
    """
    max_scale = g.max_scale()

    leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) if g.class_name(lid)=='Leaf']
    alive = [l for l in leaves if not l.leaf_fall]
    sink = [l for l in alive if l.sink]

    pool_C = sum(l.C_storage + l.C_producted for l in alive)
    PAR_total = sum(l.PAR for l in alive)

    for l in alive :
        l.weight = l.PAR / PAR_total
        l.C_exported = pool_C * l.weight
        gC = l.C - l.C_exported

        if pool_C < 0 :
            C_plant_hunger = True
                 
            for l in sink :
                l.weight = l.PAR / PAR_total
                l.C_imported = pool_C * l.weight
                    
                l.dC = leaf_status()
                l.C_deficit = l.dC + l.C_imported

        else :
            C_plant_hunger = False
                
            for l in sink :
                l.growth = True
                gC = l.C_potential_demand
                    
    
    return g, pool_C, PAR_total, C_plant_hunger



def leaf_senescence(g,
    N_potential_shoot_offer) :
    """
    Leaf senescence mecanism and leaf falling
    
    Parameters
    ----------
    g: MTG
        the above part
    N_potential_plant_offer : float
        plant potential nitrogen content

    Returns
    -------
    g : MTG
        the above part
    plant_area : float
        plant total area 
    PAR_total :
        total light intercepted by the plant
    plant_N : float
        plant nitrogen content
    pool_N : float
        senescence leaves nitrogen pool
    plant_N_by_lightunit : float
        potential nitrogen quantity by light unit
    N_plant_demand : float
        plant nitrogen demand
    N_plant_deficit : float
        plant nitrogen deficit
    N_plant_hunger :
        N plant stress state
    N_plant_remobilized : float
        plant remobilized nitrogen

    State Variables (on MTG)
    ------------------------
    growth : bool
        growth state
    leaf_fall : bool
        viability state of the leaf
    senescence : bool
        leaf development status
    photosynthetic_area : float
        photosynthetic leaf surface
    N_potential_offer : float
        potential leaf nitrogen content
    N_deficit : float
        leaf nitrogen deficit
    weight : float
        value of the leaf participation to the pool accordint to his PAR
    N_exported : float
        leaf exported nitrogen from the source/senescent leaves to the N pool
    N : float
        leaf nitrogen content
    C : float
        leaf carbon content
    N_senescnence_threshold : float
        leaf senescence initiation nitrogen threshold according to the leaf state (MTG)
    C_senescnence_threshold : float
        leaf senescence initiation carbon threshold according to the leaf state (MTG)
    PAR : float
        leaf light intercepted
    
    Hypothesis
    ----------
    H10 Senescence initiation of leaf choose according to their PAR (lowest PAR = sacrificial leaf) // Samurai law

    H11 Senescent leaf participation to the N pool according to their PAR (higher PAR = greater weight)

    References
    ----------
    Skulachev (2002)
    Thomas (2013)
    """
    max_scale = g.max_scale()
    
    leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) if g.class_name(lid)=='Leaf']
    alive = [l for l in leaves if not l.leaf_fall]
    source = [l for l in alive if l.growth]
    senescence = [l for l in alive if l.senescence]

    plant_area = sum(l.photosynthetic_area for l in alive)

    for l in alive :
        l.N = (N_potential_shoot_offer / plant_area)*l.photosynthetic_area

        plant_N = sum(l.N for l in alive)
        pool_N = sum(l.N for l in source if l.senescence)
        PAR_total = leaf_N_potential_offer(g, N_potential_shoot_offer)

        l.N_deficit = l.N - l.N_potential_offer

        N_plant_demand = sum(l.N_deficit for l in alive)
        N_plant_deficit = pool_N - N_plant_demand

        if N_plant_deficit > 0 :
            N_plant_hunger = False
            N_plant_remobilized = N_plant_demand

            for l in senescence :
                l.weight = l.PAR / PAR_total
                l.N_exported = pool_N * l.weight
                l.N = l.N - l.N_exported
                l.C = l.C
                   
            for l in alive:
                l.N = l.N_potential_offer
                l.C = l.C

        else :
            N_plant_hunger = True

            alive.sort(key = lambda l:l.PAR)

            for l in alive[0]:
                l.senescence = True
                l.N_senescence_threshold = l.N
                l.C_senescence_threshold = l.C

                for l in senescence :
                    l.weight = l.PAR / PAR_total
                    l.N_exported = pool_N * l.weight
                    l.N = l.N_senescence_threshold - l.N_exported
                    l.C = l.C_senescence_threshold
                    gN = l.N
                    gC = l.C

                N_plant_remobilized = sum(l.N_exported for l in senescence)

                senescence.sort(key = lambda l:l.PAR)

                for l in senescence[0] : 
                    l.leaf_fall = True
                    l.N_leaf_fall = l.N - l.N_exported
                    l.C_leaf_fall = l.C
                    gN = l.N_leaf_fall
                    gC = l.C_leaf_fall

                survivors = [l for l in alive if not l.leaf_fall]

                for l in survivors:
                    l.weight = l.PAR / PAR_total
                    l.N_imported = pool_N * l.weight 
                    l.N_deficit = l.N + l.N_imported
                

        return g, plant_area, plant_N, pool_N, N_plant_demand, N_plant_deficit, N_plant_hunger, N_plant_remobilized
        
        
def leaf_necrotic_area(g):
    """
    Leaf senescence necrotic surface after N remobilization (yellowing area)
    
    Parameters
    ----------
    g: MTG
        the above part

    Returns
    -------
    g : MTG
        the above part

    State Variables (on MTG)
    ------------------------
    
    leaf_fall : bool
        viability state of the leaf
    senescence : bool
        leaf development status
    area_by_Nunit : float
        surface nitrogen use efficiency ?
    N_exported : float
        leaf exported nitrogen from the source/senescent leaves to the N pool
    necrotic_area : float
        leaf necrotic surface
    """
   
    max_scale = g.max_scale()
    
    leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) if g.class_name(lid)=='Leaf']
    alive = [l for l in leaves if not l.leaf_fall]
    senescence = [l for l in alive if l.senescence]

    for l in senescence:
        l.necrotic_area = l.N_exported * l.area_by_Nunit

    return g



#main function
def turnover_plant_budget(g): 
    """
    Leaf senescence mecanism and leaf falling
    
    Parameters
    ----------
    g: MTG
        the above part

    Returns
    -------
    g : MTG
        the above part
    PAR_total :
        total light intercepted by the plant
    plant_N : float
        plant nitrogen content
    pool_N : float
        senescence leaves nitrogen pool
    plant_N_by_lightunit : float
        potential nitrogen quantity by light unit
    N_plant_demand : float
        plant nitrogen demand
    N_plant_deficit : float
        plant nitrogen deficit
    N_plant_hunger :
        N plant stress state
    N_plant_remobilized : float
        plant remobilized nitrogen
    senescence_initiation : bool
        scenescence initiation indicator
    leaf_fall : bool
        leaf falling indicator

    State Variables (on MTG)
    ------------------------
    N : float
        leaf nitrogen content
    C : float
        leaf carbon content
    leaf_fall : bool
        viability state of the leaf
    dC : float
        C deficit in the leaf
    growth : bool
        growth state
    sink : bool
        organe state
    source : bool
        organe state
    storage : bool
        C storage value
    
    Hypothesis
    ----------
    H12 Carbon stress definit nitrogen remobilization necessities (photosynthetic system turnover to optimise the biomass production according to the environmental factors)

    References
    ----------
    Thornley (1998)
    Kull (2002)
    """
    max_scale = g.max_scale()
    
    leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) if g.class_name(lid)=='Leaf']
    alive = [l for l in leaves if not l.leaf_fall]   
    
    for l in alive :
        l.growth, l.sink, l.source, l.storage, l.dC = leaf_status(g)

        if not l.growth :
            g, C_plant_hunger = C_pool_budget(g)

            if C_plant_hunger :
                g, plant_area, plant_N, pool_N, plant_N_by_lightunit, N_plant_demand, N_plant_deficit, N_plant_hunger, N_plant_remobilized = leaf_senescence(g)
                g = leaf_necrotic_area(g)
                senescence_initiation = True

                if N_plant_hunger:
                    leaf_fall = True

                else :
                    leaf_fall = False
                
            else :
                senescence_initiation = False

    return g, plant_area, plant_N, pool_N, plant_N_by_lightunit, N_plant_demand, N_plant_deficit, N_plant_hunger, N_plant_remobilized, senescence_initiation, leaf_fall

            
            
                
    

##############################################################################################################################

#MODULE SENESCENCE

def leaf_senescence_module(g):
    """
    Leaf senescence mecanism and leaf falling
    
    Parameters
    ----------
    g: MTG
        the above part

    Returns
    -------
    g : MTG
        the above part
    PAR_total :
        total light intercepted by the plant
    plant_N : float
        plant nitrogen content
    pool_N : float
        senescence leaves nitrogen pool
    plant_N_by_lightunit : float
        potential nitrogen quantity by light unit
    N_plant_demand : float
        plant nitrogen demand
    N_plant_deficit : float
        plant nitrogen deficit
    N_plant_hunger :
        N plant stress state
    N_plant_remobilized : float
        plant remobilized nitrogen
    senescence_initiation : bool
        scenescence initiation indicator
    leaf_fall : bool
        leaf falling indicator

    State Variables (on MTG)
    ------------------------
    N : float
        leaf nitrogen content
    C : float
        leaf carbon content
    PAR : float
        leaf light intercepted
    C_potential : float
        potential leaf carbon demand
    dC : float
        C deficit in the leaf
    growth : bool
        growth state
    sink : bool
        organe state
    source : bool
        organe state
    storage : bool
        C storage value
    leaf_fall : bool
        viability state of the leaf
    weight : float
        value of the leaf participation to the pool accordint to his PAR
    C_exported : float
        C exported from the source leaf
    C_imported : float
        C imported to the sink leaf
    C_deficit : float
        leaf C stress
    senescence : bool
        leaf development status
    area : float
        leaf area (surface)
    N_potential : float
        potential leaf nitrogen content
    N_deficit : float
        leaf nitrogen deficit
    N_exported : float
        leaf exported nitrogen from the source/senescent leaves to the N pool
    N_senescnence_threshold : float
        leaf senescence initiation nitrogen threshold according to the leaf state (MTG)
    C_senescnence_threshold : float
        leaf senescence initiation carbon threshold according to the leaf state (MTG)
    N_leaf_fall : float
        fall leaf nitrogen content
    C_leaf_fall : float
        fall leaf carbon content
    N_imported : float
        leaf imported nitrogen from the N pool (senescent leaves)
    
    Hypothesis
    ----------
    H4 Nitrogen plant distribution accordind to PAR (light)

    H5 Nutrient contents and PAR are homogeneous in the leaf area

    H7 Source/sink status according to leaves c stress

    H8 A leaf cannot be a source and a sink organ at the same time step

    H10 Senescence initiation of leaf choose according to their PAR (lowest PAR = sacrificial leaf) // Samurai law

    H11 Senescent leaf participation to the N pool according to their PAR (higher PAR = greater weight)

    H12 Carbon stress definit nitrogen remobilization necessities (photosynthetic system turnover to optimise the biomass production according to the environmental factors)

    References
    ----------
    Takashima (2004)
    Hikosaka (2005)
    Yasumura (2007)

    Khudairi (1970)
    Wingler (2004)

    Skulachev (2002)
    Thomas (2013)

    Thornley (1998)
    Kull (2002)
    """
    max_scale = g.max_scale()
    
    leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) if g.class_name(lid)=='Leaf']
    alive = [l for l in leaves if not l.leaf_fall]
    source = [l for l in alive if l.growth]

    plant_N = sum(l.N for l in alive)
    pool_N = sum(l.N for l in source if l.senescence)
    PAR_total = sum(l.PAR for l in alive)
    plant_N_by_lightunit = plant_N / PAR_total

    for l in alive :
        l.N_potential_offer = plant_N_by_lightunit * l.PAR
        l.N_deficit = l.N - l.N_potential_offer

        N_plant_demand = sum(l.N_deficit for l in alive)
        N_plant_deficit = pool_N - N_plant_demand

        if N_plant_deficit > 0 :
            N_plant_hunger = False
            N_plant_remobilized = N_plant_demand

            for l in alive:
                l.N = l.N_potential_offer
                l.C = l.C

        else :
            N_plant_hunger = True

            alive.sort(key = lambda l:l.PAR)

            for l in alive[0]:
                g.add_property('senescence')
                senescence = True
                l.N_senescence_threshold = l.N
                l.C_senescence_threshold = l.C

            senescence = [l for l in alive if l.senescence]

            for l in senescence :
                l.weight = l.PAR / PAR_total
                l.N_exported = pool_N * l.weight
                l.N_leaf_fall = l.N - l.N_exported
                l.C_leaf_fall = l.C
                l.N = l.N_leaf_fall
                l.C = l.C_leaf_fall

            N_plant_remobilized = sum(l.N_exported for l in senescence)

            leaf_fall = [l for l in alive if l.senescence]

            for l in leaf_fall :
                g.add_property('leaf_fall')
                leaf_fall = True

            survivors = [l for l in alive if not l.leaf_fall]

            for l in survivors:
                l.weight = l.PAR / PAR_total
                l.N_imported = pool_N * l.weight 
                l.N_deficit = l.N + l.N_imported
                l.N = l.N_deficit
                    
            alive = survivors   

        return g, plant_N, pool_N, PAR_total, plant_N_by_lightunit, pool_N, N_plant_demand, N_plant_deficit, N_plant_hunger, N_plant_remobilized








   




