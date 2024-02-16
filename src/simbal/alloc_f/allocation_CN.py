""" License Cecill-C

This module compute the allocation of C and N on a shoot and root model.
allocation_CN.py is part of simbal model

:TODO: Modify the documentation
A module that 
    (1) Check available N: store, acquired, remobilized
    (2) Distribute N among organs according to rules and hypothesis
        (H1) lux --> N --> dS
        (H2) dvt --> dS --> N
    (3) Check available C : store, acquired, remobilized
    (4) Distribute C among organs according to rules

#TODO : define the roadmap
* Separate parameters from functions (C/N) from script.
"""
### Parameters:
Ndensity = 1
Cdensity = 1 
SNU = 1 # speed of N uptake by roots in gN.cm-2.s-1
N_remob_rate = 1 # N remobilization rate in gN.gN-1
coef_N_leaf_surface = 0.046 # leaf surface = f(leaf N quantity) coefficient in cm2.g-1
coef_N_stem_surface = 0.005 # stem surface = f(stem N quantity) coefficient in cm2.g-1
alpha = 0.002 # increase of leaf area per unit of time in cm2.day-1
e_amax = 0.95 # maximum of radiation interception
k = -0.88 # extinction coefficient in Beer Lambert law
RUE = 1 # Radiation use efficiency in gC.MJ-1.cm-2
C_remob_rate = 1 # C remobilization rate in gC.gC-1

###Costs
"""
_C_respiration
_C_rhizodeposition
_N_metabolic
"""

#import numpy #for loading the entire package
from math import (pi, exp) #for loading only pi function from package numpy


### Functions:

##Root architecture
#Calculation of potential variation of volume for roots
def delta_volume(below):

    # Traverse all root segments
    g = below
    length = g.property('length')
    dlength = g.property('dlength')
    radius = g.property('radius')
    dradius = g.property('dradius')

    def seg_delta_volume(vid, below=below):
        dV = pi*(2 * radius.get(vid,0) * length.get(vid,0) * dradius.get(vid,0) + dlength.get(vid,0) * (radius.get(vid,0))**2)
        return dV
    
    max_scale = g.max_scale()
    segments = [vid for vid in g.vertices_iter(scale=max_scale) if g.label(vid)=='Segment']
    
    _delta_volume = sum(seg_delta_volume(segment) for segment in segments)
    return _delta_volume

def root_surface(below):
    g = below
    max_scale = g.max_scale()
    length = g.property('length')
    radius = g.property('radius')
    def seg_surface(vid):
        dS = 2 * pi * length.get(vid,0)  * radius.get(vid,0)
        return dS
    
    segments = [vid for vid in g.vertices_iter(scale=max_scale) if g.label(vid)=='Segment']
    _root_surface = sum(seg_surface(vid) for vid in segments)
    return _root_surface

##Shoot architecture
#Calculation of leaf surface
def leaf_surface(above):
    g = above
    max_scale = g.max_scale()
    length = g.property('length')
    radius = g.property('radius')
    def leaf_surface(vid):
        dS = 2 * pi * length.get(vid,0)  * radius.get(vid,0)
        return dS
    
    leaves = [vid for vid in g.vertices_iter(scale=max_scale) if g.label(vid)=='Leaf']
    _leaf_surface = sum(leaf_surface(vid) for vid in leaves)
    return _leaf_surface

#Calculation of stem surface
def stem_surface(above):
    g = above
    max_scale = g.max_scale()
    length = g.property('length')
    radius = g.property('radius')
    def stem_surface(vid):
        dS = 2 * pi * length.get(vid,0)  * radius.get(vid,0)
        return dS
    
    stems = [vid for vid in g.vertices_iter(scale=max_scale) if g.label(vid)=='Stem']
    _stem_surface = sum(stem_surface(vid) for vid in stems)
    return _stem_surface

def leaf_stem_surface_ratio(_leaf_surface,_stem_surface):
    """
    Calculation of the surface ratio between leaf and stem (in cm².cm-²)
    """
    tau = _leaf_surface / (_leaf_surface + _stem_surface)
    return tau


## N uptake by root
def potential_N_root_uptake(_root_surface, Nsoil, time_step, SNU):
    """
    Calculation of N potential uptake by roots
    :param _root_surface: root object
    :param Nsoil: N quantity in soil in gN
    :param time_step: time step of the simulation in day
    :param SNU:
    :return: potential N uptake by roots
    """
    _potential_N_root_uptake = SNU * Nsoil * time_step * _root_surface
    print("N root uptake is", _potential_N_root_uptake, "gN")
    return _potential_N_root_uptake


## Photosynthesis by shoot
# Intercepted radiation by shoot
def intercepted_PAR(projected_leaf_area, k, e_amax):
    """
    Calculation of PAR intercepted by the plant with Beer Lambert law
    :param projected_leaf_area: projected leaf area at the beginning of the time step in cm2
    :param k: extinction coefficient
    :param e_amax: maximum intercepted radiation
    :return: intercepted PAR in MJ.cm-2
    Intercepted PAR is printed
    """
    _intercepted_PAR = e_amax * (1 - exp(k * projected_leaf_area)) #exp is defined in math
    return _intercepted_PAR

# Transform radiation intercepted into biomass
def producted_biomass(_intercepted_PAR, RUE):
    """

    :param PAR_intercepted: PAR intercepted by the plant
    :param RUE: Radiation use efficiency
    :return: 
    """
    _producted_biomass = RUE * _intercepted_PAR
    return _producted_biomass


## N remobilized from storage by the plant 
def potential_N_remobilized(N_remob_rate, Nstock):
    """
    Calculation of N remobilized from stock (in g)
    :param Nstock: N stock in the plant in gN
    :param N_remob_rate:
    :return: N remobilized by the plant in gN
    """   
    _N_plant_remobilized = N_remob_rate * Nstock
    print("N remobilized is", _N_plant_remobilized, "gN")
    return _N_plant_remobilized


## C remobilized from storage by the plant 
def potential_C_remobilized(C_remob_rate, Cstock):
    """
    Calculation of C remobilized from storage (in g)
    :param Cstock: N storage in the plant in gC
    :param C_remob_rate:
    :return: C remobilized by the plant in gC
    """   
    _C_plant_remobilized = C_remob_rate * Cstock
    print("C remobilized is", _C_plant_remobilized, "gC")
    return _C_plant_remobilized


## CN allocation

# Calculation of potential C and N demand for the root
def N_potential_root_demand(_delta_volume, Ndensity):
    """
    Calculation of N potential root demand (in gN) in function of the root volume and root N density
    :param delta_volume: potential increase of root volume in cm3
    :param Ndensity: N root density in gN.cm-3
    :return: N root potential demand in gN
    """
    _N_potential_root_demand = _delta_volume * Ndensity
    print("N potential root demand is", _N_potential_root_demand, "gN")
    return _N_potential_root_demand

def C_potential_root_demand(_delta_volume, Cdensity):
    """
    Calculation of C potential root demand (in gC) in function of the root segment volume and root C density
    :param delta_volume: potential increase of root volume in cm3
    :param Cdensity: C root density in gC.cm-3
    :return: C root potential demand in gC
    """
    _C_potential_root_demand = _delta_volume * Cdensity
    print("C potential root demand is", _C_potential_root_demand, "gC")
    return _C_potential_root_demand


# Calculation of N offer for the plant
def N_potential_offer(_potential_N_root_uptake, _N_plant_remobilized):
    """
    Calculation of N available to the plant (in gN), in function of N uptake by roots (in g) and N remobilized
    from stock (in g) and cost of respiration
    :return: _N_potential_offer to the plant in gN
    """
    _N_potential_offer = _potential_N_root_uptake + _N_plant_remobilized
    print("N available to the plant is", _N_potential_offer, "gN")
    return _N_potential_offer


# Calculation of potential N demand for shoot, leaf and stem
def N_potential_shoot_demand(_N_potential_offer, _N_potential_root_demand):
    """
    Calculation of the potential shoot demand (gN) in function of N available to the plant (gN)
    and potential root N demand (gN)
    :param _N_potential_offer: N offer to the plant in gN
    :param N_potential_root_demand: potential root N demand in gN
    :return: N shoot potential demand in gN
    """
    _N_potential_shoot_demand = _N_potential_offer - _N_potential_root_demand
    print("N potential shoot demand is", _N_potential_shoot_demand, "gN")
    return _N_potential_shoot_demand

def N_potential_leaves_stem_demand(tau, _N_potential_shoot_demand):
    """
    Calculation of leaves and stems N potential demand (in gN) in function of the ratio of leaf/stem surface (in cm².cm-²)
    and potential shoot N demand (in gN)
    :param tau: leaf stem surface ratio in cm².cm-²
    :param N_shoot_potential_demand: potential shoot N demand in gN
    :return: leaf and stem potential N demand in gN
    """
    _N_potential_leaves_demand = tau * _N_potential_shoot_demand
    print("N potential leaves demand is", _N_potential_leaves_demand, "gN")
    _N_potential_stems_demand = (1-tau) * _N_potential_shoot_demand
    print("N potential stems demand is", _N_potential_stems_demand, "gN")
    return _N_potential_leaves_demand, _N_potential_stems_demand

def N_potential_leaf_demand(_N_potential_leaves_demand):
    """
    Calculation of leaves and stems N potential demand (in gN) in function of the ratio of leaf/stem surface (in cm².cm-²)
    and potential shoot N demand (in gN)
    :param tau: leaf stem surface ratio in cm².cm-²
    :param N_shoot_potential_demand: potential shoot N demand in gN
    :return: leaf and stem potential N demand in gN
    """
    g = above
    # get light interception as a dictionary for all organs
    light_intercepted = g.property('light_intercepted')

    # select leaves only 
    max_scale = g.max_scale() # I select the organ scale 
    leaves = [vid for vid in g.vertices_iter(scale=max_scale) if g.label(vid)=='Leaf']
    """leaves = []
    for vid in g:
        if g.label(vid)=='Leaf':
            leaves.append(vid)
    """ # same code then previous instruction

    total_light = sum(light_intercepted[lid] for lid in leaves)
    N_by_lightunit = _N_potential_leaves_demand / total_light
    _N_potential_leaf_demands = dict((lid, N_by_lightunit*light_intercepted[lid]) for lid in leaves)
    
    g.add_property('N_potential_leaf_demands')
    Nps = g.property('N_potential_leaf_demands')
    Nps.update(_N_potential_leaf_demands)

    #  add to the MTG the N_potential_leaf_demand

    #_N_potential_leaf_demand = light_intercepted_leaf * _N_potential_leaves_demand
    _N_potential_leaf_demand = _N_potential_leaf_demands
    return _N_potential_leaf_demand

### N_potential_stems_demand distribitued equally between every segments = N_potential_stem_demand
def N_potential_stem_demand(_N_potential_stems_demand):
    g = above
    # select internodes only
    max_scale = g.max_scale() # I select the organ scale 
    internodes = [vid for vid in g.vertices_iter(scale=max_scale) if g.label(vid)=='Internode']
    n = len(internodes)
    _N_potential_stem_demands = dict((sid, _N_potential_stems_demand/n) for sid in internodes)
    return _N_potential_stem_demands


# Potential C demand of leaf, stem and shoot calculated from N available to each compartment
def potential_leaf_area(_N_potential_leaf_demand, coef_N_leaf_surface):
    """
    Calculation of potential leaf area (cm2) in function of potential N leaf demand (gN)
    and the coefficient of the relationship of leaf area - leaf N quantity (cm2.g-1)
    :param N_leaf_potential_demand: leaf potential N demand in gN
    :param coef_N_leaf_area: coefficient between leaf area and leaf N quantity in cm2.g-1
    :return: potential leaf area in cm2
    Potential leaf area is printed
    """
    _potential_leaf_surface = {}
    for lid in _N_potential_leaf_demand:
        _potential_leaf_surface = coef_N_leaf_surface * _N_potential_leaf_demand[lid] 
    print("potential leaf surface is", _potential_leaf_surface, "cm2")
    return _potential_leaf_surface

def potential_stem_surface(_N_potential_stem_demand, coef_N_stem_surface):
    """
    Calculation of potential stem area (cm2) in function of potential N stem demand (gN)
    and the coefficient of the relationship of stem area - stem N quantity (cm2.g-1)
    :param N_stem_potential_demand: stem potential N demand in gN
    :param coef_N_stem_area: coefficient between stem area and leaf N quantity in cm2.g-1
    :return: potential stem area in cm2
    Potential stem area is printed
    """
    _potential_stem_surface = coef_N_stem_surface * _N_potential_stem_demand
    print("potential stem surface is", _potential_stem_surface, "cm2")
    return _potential_stem_surface

### should be projected surface####
def C_potential_shoot_demand(_potential_leaf_surface, _potential_stem_surface):
    """
    Calculation of C potential shoot demand (in gC) in function of the potential photosynthetic area and shoot N density
    :param potential_photosynthetic_area: potential photosynthetic area in cm2
    :param Cdensity: C shoot density in gC.cm-2
    :return: C root potential demand in gC
    C shoot potential demand is printed
    """
    _C_potential_shoot_demand = _potential_leaf_surface + _potential_stem_surface
    print("C potential shoot demand is", _C_potential_shoot_demand, "gC")
    return _C_potential_shoot_demand

# Calculation of C actual offer for the plant

def C_actual_offer(_producted_biomass, _C_plant_remobilized, _C_respiration):
    """
    Calculation of C available to the plant (in gC), in function of C uptake (in g) and C remobilized
    from stock (in g)
    :param PAR_intercepted: PAR intercepted by the plant
    :param RUE: Radiation use efficiency
    """
    _C_actual_offer = _producted_biomass + _C_plant_remobilized - _C_respiration
    print("C available to the plant is", _C_actual_offer, "gC")
    return _C_actual_offer

# Actual C uptake resulting in confrontation between C actual offer and potential demand
def C_actual_uptake(_C_actual_offer, _C_potential_shoot_demand):
    """
    Calculation of C available to the plant (in gC), in function of C uptake (in g) and C remobilized
    from stock (in g)
    :param PAR_intercepted: PAR intercepted by the plant
    :param RUE: Radiation use efficiency
    """
    _C_actual_uptake = _C_actual_offer - _C_potential_shoot_demand
    return _C_actual_uptake

#if C actual uptake < 0, then C actual uptake = _C_actual_offer
#



# simulation
#Calculate the C and N potential demand of root (gC and gN)
#C_pot_demand, N_pot_demand=CN_potential_root_demand(delta_volume=2)
#print("The potential demand in C is", C_pot_demand, "gC.")

#list_root=[]
#list_root.append(root_element1)
#list_root_segments.append(root_element2)
#C_potential_root_total_demand(list_elements=list_root_segments)


def root_potential_growth(below, dl, dr):
    """ External function
    Can be define by the modeller for each plant? or a varity or all the plants?
    """

    # 1. Query all the segments
    g = below
    max_scale = g.max_scale()
    segments = [vid for vid in g.vertices_iter(scale=max_scale) if g.label(vid) =='Segment']
    growing_segments = [vid for vid in segments if g.is_leaf(vid)]
    
    # Init properties
    length = g.properties().setdefault('length',{}) # state
    dlength = g.properties().setdefault('dl',{})    # rate
    
    radius = g.properties().setdefault('radius',{}) # state
    dradius = g.properties().setdefault('dr',{})    # rate

    # Select all the segments that are connected to a bud
    
    for sid in growing_segments:
        length[sid] = length.get(sid, 0.)
        dlength[sid] = dl
        

    for sid in segments:
        radius[sid] = radius.get(sid, 0.)
        dradius[sid] = dr
        
        #volume[sid] = 2 * pi * radius[sid] * length[sid]
    

def shoot_potential_growth(above, dl, dS):
    """ The same function but for above ground. 
    """
    
    # 1. Query all the segments
    g = above
    max_scale = g.max_scale()

    shoots = [vid for vid in g.vertices_iter(scale=max_scale)]
    internodes = [vid for vid in shoots if g.label(vid) == 'Internode']
    leaves = [vid for vid in shoots if g.label(vid) == 'Leaf']
    growing_internodes = [g.parent(vid) for vid in shoots if (g.label(vid) == 'Bud') and (g.edge_type(vid)=='<')]

    # Init properties
    length = g.properties().setdefault('length',{})
    dlength = g.properties().setdefault('dlength',{})

    area = g.properties().setdefault('area',{}) # or surface
    darea = g.properties().setdefault('darea',{}) # rate
    radius = g.properties().setdefault('radius',{})
 
    for vid in internodes:
        radius[vid] = 1.

    for vid in growing_internodes:
        length[vid] = length.get(vid, 0.)
        dlength[vid] = dl

    for vid in leaves:
        area[vid] = area.get(vid, 0.)
        darea[vid] = dS


def main():
    """ Run the above"""

    # initializing

    starting_time=0
    end_time=10
    time_step=2

    # Create the plant
    rapeseed = Plant(name='Rapeseed')
    above = rapeseed.above
    below = rapeseed.below

    for t in range(starting_time,end_time,time_step):
        print("Time is", t)
        # Elements are growing:
        
        root_potential_growth(dl=time_step/10, dr=time_step/50)
        shoot_potential_growth(dl=time_step/10, dS=alpha * time_step)



        # Calculate potential C and N demand for total root
        dV = delta_volume(below)
        _N_potential_root_demand = N_potential_root_demand(dV, Ndensity=1)
        _C_potential_root_demand = C_potential_root_demand(dV, Cdensity=1)
        
        # Calculate N available to shoot and root
        _N_available = N_available(below, Nsoil=20, Nstock=5, time_step=time_step, SNU=SNU)


        # Calculate N potential demand of shoot, leaf and stem
        _N_shoot_potential_demand = N_shoot_potential_demand(_N_available=_N_available, _N_potential_root_demand=_N_potential_root_demand)
        tau = leaf_stem_surface_ratio(element1=leaf_element, element2=stem_element)
        _N_leaf_potential_demand = N_leaf_potential_demand(tau=tau, _N_shoot_potential_demand=_N_shoot_potential_demand)
        _N_stem_potential_demand = N_stem_potential_demand(tau=tau, _N_shoot_potential_demand=_N_shoot_potential_demand)
        # Calculate potiential leaf and stem area
        _potential_leaf_area = potential_leaf_area(_N_leaf_potential_demand=_N_leaf_potential_demand, coef_N_leaf_area=coef_N_leaf_area)
        _potential_stem_area = potential_stem_area(_N_stem_potential_demand=_N_stem_potential_demand, coef_N_stem_area=coef_N_stem_area)
        # Calculate potential total photosynthetic area
        potential_photosynthetic_area = _potential_leaf_area + _potential_stem_area
        # Calculate C potential shoot demand
        _C_potential_shoot_demand = C_potential_shoot_demand(potential_photosynthetic_area=potential_photosynthetic_area, Cdensity=1)
        # Calculate C available to the plant
        _intercepted_PAR = intercepted_PAR(initial_leaf_area=leaf_area0, k=k, e_amax=e_amax)
        _C_available = C_available(_intercepted_PAR=_intercepted_PAR, RUE=RUE, Cstock=3, time_step=time_step)






#N_seg_available=N_available(segment_length=1,segment_radius=1, Nsoil=10, N_remobilized=2, time_step=1)
#print(N_seg_available)

if __name__ == '__main__':
    main()