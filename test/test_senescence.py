#Scenarios
def init_scenario1_HLHN(g,
    Light: float = 450,
    Nitrate: float = 8,
    N: float = 1.0, 
    C: float = 100):
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
    H1  Fast increase of dry mass
    H2  Optimal growth
    H3  fast increase of plant N content
    H4  Faster N remobilisation than LL- (remobilisation = decrease of leaf N content)
    
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