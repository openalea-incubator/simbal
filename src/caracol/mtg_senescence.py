#MTG
from openalea.mtg import MTG

class Plant:
    """ Define a plant as an above and below ground MTG.

    Some usefull functions can be defined as proxies.
    
    """
    def __init__(self, name='', above=None, below=None):
        self.name = name
        self.above = self.init_above_mtg() if above is None else above
        self.below = self.init_below_mtg() if below is None else below



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
    phyto1 = g.add_component(a1, label='Phytomer')

    # Scale Organ
    inter1 = g.add_component(phyto1, label = 'Internode')
    leaf1 = g.add_child(inter1, edge_type='+', label = 'Leaf')

    #Add a new leaf
    phyto2 = g.add_child(phyto1, edge_type='<', label='Phytomer')
    inter2 = g.add_component(phyto2, label = 'Internode')
    leaf2 = g.add_child_and_complex(parent=inter1,complex=phyto1, edge_type='+', label = 'Leaf')

    fat_mtg(g)
    return g

assert g.nb_vertices(scale=1) == 1
assert g.nb_vertices(scale=2) == 5
assert g.complex(v4) == plant1


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

