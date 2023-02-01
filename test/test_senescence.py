from caracol import mtg_senescence
from caracol import senescence as snc

#Scenarios
def init_scenario1_HLHN(
    Light: float = 450,
    Nitrate: float = 8,
    N: float = 1.0, 
    C: float = 100,
    Surface: float = 5,
    Length: float = 10):
    """ High Light - High Nitrate (witness)

    Properties
    ----------
    Light : float (umol.m-2.s-1)
        PPFD Photosynthetic Photon Flux Density = The amount of micromoles of photons striking a square meter per second
    Nitrate : float (mM NO3-)
        Fertilisation input in soil solution
    N : float (gN.m-2)
        Nitrogen leaf content
    C : float (mmolC.m-2)
        Carbohydrates leaf content (glucose + sucrose + starch)
    

    Hypothesis
    -------
    H1  C and N distribution drive by light
    H2  Senescence initiation by Cand N threshold
   
    """



    above = mtg_senescence.init_above_mtg()
    below = mtg_senescence.init_below_mtg()
    
    leaves = [vid for vid in above.vertices(scale=4) if above.label(vid) =='Leaf']
    segs = [vid for vid in below.vertices(scale=3) if below.label(vid) =='Root']

    for vid in leaves: 
        nid = above.node(vid)
        nid.Light = Light
        nid.N = N
        nid.C = C
        nid.Surface = Surface


    for vid in segs: 
        nid = below.node(vid)
        nid.Nitrate = Nitrate
        nid.Length = Length




    #Nitrogen Budget
    _potential_N_root_uptake = snc.potential_N_root_uptake(Nitrate, _root_surface, time_step, SNU )
    _N_potential_offer = snc.N_potential_offer(_potential_N_root_uptake)
    _N_potential_root_demand = snc.N_potential_root_demand(_delta_volume, Ndensity)
    _N_potential_plant_offer = snc.N_potential_plant_offer(_N_potential_offer, _N_potential_root_demand)

    #Light Budget
    _leaf_actual_surface = snc.leaf_actual_surface(above, length, radius)
    _photosynthetic_leaf_surface = snc.photosynthetic_leaf_surface(leaf_actual_surface, leaf_necrotic_surface)
    _potential_leaf_light_interception = snc.potential_leaf_light_interception(k, Light, _photosynthetic_leaf_surface)
    _leaf_producted_biomass = snc.leaf_photosynthesis(_leaf_intercepted_PAR, RUE)
    _leaf_N_potential_offer = snc.leaf_N_potential_offer(_N_potential_plant_offer, Light, _potential_leaf_light_interception, above)
    _potential_leaf_surface = snc.potential_leaf_surface(_leaf_N_potential_offer, coef_N_leaf_surface)
    _FSPM_Allirand_cst = snc.FSPM_Allirand(_potential_leaf_light_interception, _leaf_N_potential_offer)

    #Carbon Budget
    _C_potential_root_demand = snc.C_potential_root_demand(_delta_volume, Cdensity)
    _C_potential_leaf_demand = snc.C_potential_leaf_demand(_C_potential_leaf_surface, Cdensity)
    _C_potential_leaf_offer = snc.C_potential_leaf_offer(_leaf_producted_biomass, _C_respiration)

    #Turnover Budget
    _C_actual_uptake = snc.C_actual_uptake(_C_potential_leaf_offer, _C_potential_leaf_demand)
    _C_actual_uptake = snc.C_actual_uptake_import(_C_actual_uptake, _C_source_offer)
    _C_source_offer, _N_source_offer = snc.leaf_senescence(_C_potential_leaf_offer, _N_potential_leaf_offer, _C_sink_demand, _N_sink_demand, w)

























##########################################################################################################################################################
def init_scenario2_HLLN(g,
    Light: float = 450,
    Nitrate: float = 0.8, 
    N: float = 0.8, 
    C: float = 250):
    """ High Light - Low Nitrate 

    Properties
    ----------
    Light : float (umol.m-2.s-1)
        PPFD Photosynthetic Photon Flux Density = The amount of micromoles of photons striking a square meter per second
    Nitrate : float (mM NO3-)
        Fertilisation input in soil solution
    N : float (gN.m-2)
        Nitrogen leaf content
    C : float (mmolC.m-2)
        Carbohydrates leaf content (glucose + sucrose + starch)
    

    Hypothesis
    -------
    H1  More repressed growth than witness
    H2  Shorter increase period of plant N content than HN-
    
    """
    
    leaves = [vid for vid in g.vertices(scale=4) if g.label(vid) =='Leaf']
    segs = [vid for vid in g.vertices(scale=3) if g.label(vid) =='Root']

    for vid in leaves: 
        nid = g.node(vid)
        nid.Light = Light
        nid.N = N
        nid.C = C

    for vid in segs: 
        nid = g.node(vid)
        nid.Nitrate = Nitrate


def init_scenario3_LLHN(g,
    Light: float = 150,
    Nitrate: float = 8,
    N: float = 1.1, 
    C: float = 30):
    """ Low Light - High Nitrate

    Properties
    ----------
    Light : float (umol.m-2.s-1)
        PPFD Photosynthetic Photon Flux Density = The amount of micromoles of photons striking a square meter per second
    Nitrate : float (mM NO3-)
        Fertilisation input in soil solution
    N : float (gN.m-2)
        Nitrogen leaf content
    C : float (mmolC.m-2)
        Carbohydrates leaf content (glucose + sucrose + starch)
    

    Hypothesis
    -------
    H1  More repressed growth than witness
    H2  Slower increase of plan N content than HL-
    H3  Slower N remobilisation than HL- (remobilisation = decrease of leaf N content)
    
    """
    
    leaves = [vid for vid in g.vertices(scale=4) if g.label(vid) =='Leaf']
    segs = [vid for vid in g.vertices(scale=3) if g.label(vid) =='Root']

    for vid in leaves: 
        nid = g.node(vid)
        nid.Light = Light
        nid.N = N
        nid.C = C

    for vid in segs: 
        nid = g.node(vid)
        nid.Nitrate = Nitrate


def init_scenario4_LLLN(g,
    Light: float = 150,
    Nitrate: float = 0.8, 
    N: float = 0.7, 
    C: float = 10):
    """ Low Light - Low Nitrate 

    Properties
    ----------
    Light : float (umol.m-2.s-1)
        PPFD Photosynthetic Photon Flux Density = The amount of micromoles of photons striking a square meter per second
    Nitrate : float (mM NO3-)
        Fertilisation input in soil solution
    N : float (gN.m-2)
        Nitrogen leaf content
    C : float (mmolC.m-2)
        Carbohydrates leaf content (glucose + sucrose + starch)
    

    Hypothesis
    -------
    H1  More repressed growth than witness
    H2  Tardive first leave senesced
    H3  Slower increase of plant N content than HL-
    H4  Faster N remobilisation than HN- (remobilisation = decrease of leaf N content)
    
    """
    
    leaves = [vid for vid in g.vertices(scale=4) if g.label(vid) =='Leaf']
    segs = [vid for vid in g.vertices(scale=3) if g.label(vid) =='Root']

    for vid in leaves: 
        nid = g.node(vid)
        nid.Light = Light
        nid.N = N
        nid.C = C

    for vid in segs: 
        nid = g.node(vid)
        nid.Nitrate = Nitrate