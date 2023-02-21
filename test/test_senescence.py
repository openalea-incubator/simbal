from caracol import mtg_senescence
from caracol import senescence as snc

scenario = mtg_senescence.PlantScenario
g = above = scenario.init_above_mtg()

max_scale = g.max_scale()

#Test de Ono (1996)
def init_scenario1_HLHN(g):
    """ High Light - High Nitrate (positive witness)

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
    """
    PAR1, PAR2 = 200, 450
    s1, s2 = 5., 2.
    C1, C2 = 10, 10
    N1, N2 = 10, 10

    leaves = [lid for lid in g.vertices(scale=max_scale) if g.class_name(lid)=='Leaf']
    assert len(leaves)==2

    l1, l2 = leaves
    l1, l2 = g.node(l1), g.node(l2)

    l1.PAR = PAR1
    l2.PAR = PAR2

    l1.surface = s1    
    l2.surface = s2

    l1.C = C1    
    l2.C = C2

    l1.N = N1    
    l2.N = N2

    l1.growth = False
    l2.growth = True

    return g


def init_scenario2_HLLN(g):
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
    """
    PAR1, PAR2 = 200, 450
    s1, s2 = 5., 2.
    C1, C2 = 10, 10
    N1, N2 = 10, 10

    leaves = [lid for lid in g.vertices(scale=max_scale) if g.class_name(lid)=='Leaf']
    assert len(leaves)==2

    l1, l2 = leaves
    l1, l2 = g.node(l1), g.node(l2)

    l1.PAR = PAR1
    l2.PAR = PAR2

    l1.surface = s1    
    l2.surface = s2

    l1.C = C1    
    l2.C = C2

    l1.N = N1    
    l2.N = N2

    l1.growth = False
    l2.growth = True

    return g


def init_scenario3_LLHN(g):
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
    """
    PAR1, PAR2 = 200, 450
    s1, s2 = 5., 2.
    C1, C2 = 10, 10
    N1, N2 = 10, 10

    leaves = [lid for lid in g.vertices(scale=max_scale) if g.class_name(lid)=='Leaf']
    assert len(leaves)==2

    l1, l2 = leaves
    l1, l2 = g.node(l1), g.node(l2)

    l1.PAR = PAR1
    l2.PAR = PAR2

    l1.surface = s1    
    l2.surface = s2

    l1.C = C1    
    l2.C = C2

    l1.N = N1    
    l2.N = N2

    l1.growth = False
    l2.growth = True

    return g


def init_scenario4_LLLN(g):
    """ Low Light - Low Nitrate (negative witness)

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
    """

    PAR1, PAR2 = 200, 450
    s1, s2 = 5., 2.
    C1, C2 = 10, 10
    N1, N2 = 10, 10

    leaves = [lid for lid in g.vertices(scale=max_scale) if g.class_name(lid)=='Leaf']
    assert len(leaves)==2

    l1, l2 = leaves
    l1, l2 = g.node(l1), g.node(l2)

    l1.PAR = PAR1
    l2.PAR = PAR2

    l1.surface = s1    
    l2.surface = s2

    l1.C = C1    
    l2.C = C2

    l1.N = N1    
    l2.N = N2

    l1.growth = False
    l2.growth = True

    return g






##################################################################################################################################################################################


time_step = 10
for i in range(time_step):
    g = init_scenario1_HLHN(g, time_step)
    g, plant_N, pool_N, PAR_total, plant_N_by_lightunit, pool_N, N_plant_demand, N_plant_deficit, N_plant_hunger, N_plant_remobilized = snc.leaf_senescence_module(g, time_step)










































