#MTG
from openalea.mtg import MTG, fat_mtg
from .plant_mtg import Plant

class PlantScenario(Plant):
    """ Define a plant as an above and below ground MTG.

    Some usefull functions can be defined as proxies.
    
    """

    @staticmethod
    def init_above_mtg():
        """ Simple function to create a simple aerial MTG
        Return a seedling.
        This seedling is composed of 4 scales:
            - Plant: scale 1
            - Axis: scale 2
            - Phytomer: scale 3
            - Organ: scale 4
                - Internode
                - Leaf
                - Bud
                - Petiole
        
        Note: 
            - < : successor
            + : ramification
        """

        g = MTG()

        # Scale Plant
        plant1 = g.add_component(g.root, label='Plant')

        # Scale Axis
        axis1 = g.add_component(plant1, label='Axis')

        # Scale Phytomer
        phyto1 = g.add_component(axis1, label='Phytomer')

        # Scale Organ
        inter1 = g.add_component(phyto1, label = 'Internode')
        leaf1 = g.add_child(inter1, edge_type='+', label = 'Leaf')

        #Add a new leaf
        phyto2 = g.add_child(phyto1, edge_type='<', label='Phytomer')
        inter2, phyto2 = g.add_child_and_complex(parent=inter1,complex=phyto2, edge_type='<', label = 'Internode')
        leaf2 = g.add_child(inter2, edge_type='+', label = 'Leaf')

        fat_mtg(g)

        """
        assert g.nb_vertices(scale=1) == 1
        assert g.nb_vertices(scale=2) == 5
        assert g.complex(v4) == plant1
        """
        return g


    @staticmethod
    def init_below_mtg():
        """ Simple function to create a simple belowground MTG or RootSystem Architecture
            Return a seedling.
        This seedling is composed of 3 scales:
            - RootSystem Architecture RSA : scale 1
            - Root: scale 2   
            -Segment: scale 3
                
        Note: 
            - < : successor
            + : ramification
        """
        g = MTG()
        
        # Add a plant (scale is 1)
        rsa = g.add_component(g.root, label='RSA')
        root = g.add_component(rsa, label='Root')
        
        # Last scale
        seg = g.add_component(root, label='Segment')
    
        fat_mtg(g)
        return g



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
    
    leaves = [vid for vid in g.vertices() if g.label(vid) =='Leaf']
    roots = [vid for vid in g.vertices() if g.label(vid) =='Root']

    for vid in leaves: 
        g.property('Light', 'N', 'C')

    for vid in roots: 
        g.property('Nitrate')


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
    
    leaves = [vid for vid in g.vertices() if g.label(vid) =='Leaf']
    roots = [vid for vid in g.vertices() if g.label(vid) =='Root']

    for vid in leaves: 
        g.property('Light', 'N', 'C')

    for vid in roots: 
        g.property('Nitrate')


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
    
    leaves = [vid for vid in g.vertices() if g.label(vid) =='Leaf']
    roots = [vid for vid in g.vertices() if g.label(vid) =='Root']

    for vid in leaves: 
        g.property('Light', 'N', 'C')

    for vid in roots: 
        g.property('Nitrate')


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
    
    leaves = [vid for vid in g.vertices() if g.label(vid) =='Leaf']
    roots = [vid for vid in g.vertices() if g.label(vid) =='Root']

    for vid in leaves: 
        g.property('Light', 'N', 'C')

    for vid in roots: 
        g.property('Nitrate')