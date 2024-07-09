#TODO
#liste des paramètres
# define mtg for leaves
#indicate units
#simulation loop
#les priorités d’allocation des ressources dans chaque organe : dépend du PARi qqsoit l'organe
#le calcul des demandes en N des feuilles et des tiges
#tout passer en flux (rate) et pas quantité
#couplage avec aptos ? 
#N stem based on length and not area like leaves for growth ? i.e. Archisimple stem
#pour le végétatif: EN court sont tous pareils, ils sont empilés avec l'émergence de la feuille avec diamètre fixe = N potential stem est la même pour chaque EN court, 
# vitesse d'élongation maximale entre l'émergence de 2 feuilles consécutives


###Calculating photosynthesis and maintenance
# Intercepted radiation by Caribu for polygons
intercepted_PAR_per_polygon : number of photons (micromol.s-1.m-2)
#calculated by Caribu (mean intercepted PAR per organ), this value integrate the area of the polygon
# check the unit particularly .s-1

#Calculate intercepted radiation by each organ
def intercepted_PAR_per_organ(intercepted_PAR_per_polygon, polygons_of_each_organ):
    """
    :return: 
    """
    #read Excel file of Caribu outputs
    #sum for each polygon per leaf (i.e. mean intercepted PAR per leaf)
    #check the output of Caribu (micromol.s-1.m-2, link to organ area)

# Transform radiation intercepted into C
def produced_C_rate(_intercepted_PAR, SCA, organ_area):
    """

    :param PAR_intercepted in micromol.s-1.m-2: PAR intercepted by organs
    :param SCA : gC.micromol-1 Specific carbon assimilation (similar to RUE but not net photosynthesis and at the organ scale), do not vary with N (cf experiments of F.Brun)
    :param organ_area in m2: area of photosynthetic organ
    :return: produced_C_rate in gC.s-1 
    """
    _produced_C_rate = SCA * _intercepted_PAR * organ_area
    return _produced_C_rate 

# C cost for maintenance (i.e. respiration) in opposition to C for growth
def maintenance_rate(_produced_C_rate, k, organ):
    """
    k in gC/gC 
    :return: _C_maintenance_cost_rate  in gC.s-1.organ
    """
    _C_maintenance_cost_rate = k * _produced_C_rate
    return _C_maintenance_cost_rate 
# Calculation of C actual offer for the plant

def C_shoot_actual_offer_rate(_produced_C_rate, remobilized):
    """
    Calculation of C available to the plant (in gC), in function of C uptake (in g) and C remobilized
    from stock (in g)
    """
    #sum of C produced rate by each organ 
    #calculation of C net for plant produced + remobilized
    print("C available to the plant is", _C_actual_offer, "gC")
    return _C_shoot_actual_offer_rate

#either C allocation (1) proportionnal to PARi for each organ or (2) per phytomer then distributed between organs depending on area

def C_shoot_organ_allocation_rate(__C_shoot_actual_offer_rate):
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
    C_by_lightunit = _C_shoot_actual_offer_rate / total_light
    _C_shoot_actual_offer_rate = dict((lid, C_by_lightunit*light_intercepted[lid]) for lid in leaves)
    
    g.add_property('C_shoot_actual_offer_rate')
    g.property('C_shoot_actual_offer_rate') = _C_shoot_actual_offer_rate.values()

    #  add to the MTG the N_potential_leaf_demand

    _C_shoot_organ = light_intercepted_leaf * _N_potential_leaves_demand
    _C_shoot_organ = _N_potential_leaf_demands
    return _C_shoot_organ

#
def C_N_shoot_balance():
    #for each organ, priority of maintenance for shoot, option 2 : retirer apres l'allocation 
    produced_C_rate(_intercepted_PAR, SCA, organ_area)
    maintenance_rate(_produced_C_rate, k, organ) 
    C_contribution_rate = produced_C_rate - maintenance_rate
    C_shoot_organ_allocation_rate(__C_shoot_actual_offer_rate)
    C_consumption_by_leaf(_C_actual_offer, _C_leaf_demand)
    # C mobile has been updated
    leaf_actual_structural_growth(_C_actual_offer, _C_leaf_demand, external_C, Cdensity_leaf, minimal_thickness,yield_growth)
    

 



###Allocating N to each organ
#N absorbed by the root system # from Root Simple, available for entire root system or root segment
#def N_absorbed_by_root(absorbed_N_per_root_segment, number_of_root_segments):
    return(_N_absorbed_by_root)

#N remobilised from other organs calculated by Root Simple
d#ef N_remobilization():
    return _N_remobilized 

# Calculation of potential C and N demand for the root calculated by Root Simple
#def N_consumed_by_root():
    return _N_consumed_by_root

#def C_potential_root_demand():#
    return _C_potential_root_demand

# Calculation of N offer for the plant from outputs of root Simple
def N_available_for_shoot():
    return _N_available_for_shoot # add N_remobilized

# Calculation of potential N demand for shoot
#def N_potential_shoot_demand(_N_potential_offer, _N_potential_root_demand):
    """
    Calculation of the potential shoot demand (gN) in function of N available to the plant (gN)
    and potential root N demand (gN)
    :param _N_potential_offer: N offer to the plant in gN
    :param N_potential_root_demand: potential root N demand in gN
    :return: N shoot potential demand in gN
    """
    _N_potential_shoot_demand = _N_potential_offer - _N_consumption_by_root
    print("N potential shoot demand is", _N_potential_shoot_demand, "gN")
    return _N_potential_shoot_demand
    ef N_potential_leaves_stem_demand(tau, _N_potential_shoot_demand):
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

#leaves_N_demand ? stem N Demand ? 

def N_potential_leaf_allocation(_N_potential_leaves_demand):
    """
    Calculation of leaves and stems N potential demand (in gN) in function of the ratio of leaf/stem surface (in cm².cm-²)
    and potential shoot N demand (in gN)
    :intercepted light in micromol.s-1
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
    #Nps = g.property('N_potential_leaf_demands')
    #Nps.update(_N_potential_leaf_demands)
    g.property('N_potential_leaf_demands') = _N_potential_leaf_demands.values()

    #  add to the MTG the N_potential_leaf_demand

    _N_potential_leaf_allocation = light_intercepted_leaf * _N_potential_leaves_demand
    _N_potential_leaf_allocation = _N_potential_leaf_demands
    return _N_potential_leaf_allocation

### N_potential_stems_demand distribitued equally between every segments = N_potential_stem_demand
def N_potential_stem_demand(_N_potential_stems_demand):
    g = above
    # select internodes only
    max_scale = g.max_scale() # I select the organ scale 
    internodes = [vid for vid in g.vertices_iter(scale=max_scale) if g.label(vid)=='Internode']
    n = len(internodes)
    _N_potential_stem_demands = dict((sid, _N_potential_stems_demand/n) for sid in internodes)
    return _N_potential_stem_demands


###Potential C demand of leaf, stem and shoot calculated from N available to each compartment
def potential_leaf_area(_N_leaf_offer, coef_N_leaf_surface):
    """
    Calculation of potential leaf area (cm2) in function of potential N leaf demand (gN)
    and the coefficient of the relationship of leaf area - leaf N quantity (cm2.g-1)
    :param N_leaf_potential_demand: leaf potential N demand in gN
    :param coef_N_leaf_area: coefficient between leaf area and leaf N quantity in cm2.g-1
    :return: potential leaf area in cm2
    Potential leaf area is printed
    """
    _potential_leaf_surface = {}
    _potential_leaf_surface = coef_N_leaf_surface * _N_leaf_offer
    print("potential leaf surface is", _potential_leaf_surface, "cm2")
    return _potential_leaf_surface

#leaf C demand
def leaf_potential_growth(_potential_leaf_surface, minimal_thickness, LA_max):
    delta_leaf_volume = min(_potential_leaf_surface, LA_max) * minimal_thickness #LA_max vitesse max de croissance en surface
    return(delta_leaf_volume)

def C_leaf_demand(delta_leaf_volume, Cdensity_leaf, yield_growth):
    _C_leaf_demand = delta_leaf_volume * Cdensity_leaf / yield_growth
    return(_C_leaf_demand)

def stem_potential_growth(growth_rate, stem_length, stem_diameter):
    _potential_length = stem_length + growth_rate * time_step
    delta_stem_volume = (_potential_length - stem_length) * stem_diameter **2 /4 * 3.1416
    return(delta_stem_volume)

def C_stem_demand(delta_stem_volume, Cdensity_stem, yield_growth):
    """
    yield_growth: 0.8 proportion of total C used for biomass production, the rest is respired
    """
    _C_stem_demand = delta_stem_volume * Cdensity_stem / yield_growth
    return(_C_stem_demand)


# Actual C uptake resulting in confrontation between C actual offer and potential demand
def C_consumption_by_leaf(_C_actual_offer, _C_leaf_demand):
    """
    Calculation of C available to the plant (in gC), in function of C uptake (in g) and C remobilized
    from stock (in g)
    :param PAR_intercepted: PAR intercepted by the plant
    :param RUE: Radiation use efficiency
    :_C_mobile : C remaining after structural leaf growth
    """
    _C_consumption_by_leaf = min(_C_actual_offer,_C_leaf_demand)
    _C_mobile = _C_actual_offer - _C_consumption_by_leaf
    return _C_consumption_by_leaf

def leaf_actual_structural_growth(_C_actual_offer, _C_leaf_demand, external_C, Cdensity_leaf, minimal_thickness,yield_growth):
    if _C_actual_offer < _C_leaf_demand:
        actual_leaf_C_consumption = min(_C_actual_offer + external_C, _C_leaf_demand)
    else:
        actual_leaf_C_consumption = _C_leaf_demand
    actual_C_structural_mass_increment = actual_leaf_C_consumption * yield_growth 
    actual_leaf_area_increment = actual_C_structural_mass_increment / (Cdensity_leaf * minimal_thickness)
    actual_leaf_area = actual_leaf_area + actual_leaf_area_increment
    actual_C_structural_mass = actual_C_structural_mass + actual_C_structural_mass_increment 
    return()

def actual_N_leaf_consumption(actual_leaf_area_increment, coef_N_leaf_surface, _N_leaf_offer):
    _actual_N_leaf_consumption = actual_leaf_area_increment / coef_N_leaf_surface
    N_mobile =_N_leaf_offer - _actual_N_leaf_consumption
    return()

def N_remobilization(time_step):
    remobilization_rate = k * N_mobile
    N_remobilized_leaf = remobilization_rate * time_step
    N_mobile = N_mobile - N_remobilized_leaf 
    #organ_N_quantity = N_mobile

def C_shoot_mobile(_C_mobile):
    _C_shoot_mobile = max(_C_mobile-C_root_demand,0) #C_root_demand (structural + metabolique + rhizodeposition + respiration ?) is calculated from previous time step to call root model then shoot model then end of time step
    #allocation C shoot mobile to each shoot organ depending on PARi (cc function above for N or C allocation)
    organ_C_quantity = actual_C_structural_mass + _C_organ_mobile 
    organ_biomass = organ_C_quantity / coefBiomass #coefBiomass depends on organ type (env 0.4)
    #pas de simulation d'épaisseur, uniquement biomasse, QC et QN

#que fait-on de N _remobilized ?
