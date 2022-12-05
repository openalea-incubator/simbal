#MTG
from openalea.mtg import MTG

def aerial_plant():
    """ Return a seedling.
    This seedling is composed of 4 scales:
        - Plant: scale 1
        - Axis: scale 2
        - Phytomer: scale 3
        - Organ: scale 4
            - S: stem
            - L: leaf
            - B: bud
            - PL: petiole
    
    Note: 
        - < : successor
        + : ramification
    """

g = MTG()

# Scale Plant
plant1 = g.add_component(g.root)

# Scale Axis
v1 = g.add_component(plant1, edge_type='<', label='A')

# Scale Phytomer
phytomer = g.add_component(v1, edge_type='<', label='P')

# Scale Organ
stem = g.add_component(v2, label = 'S')
leaf = g.add_child(stem, edge_type='+', label = 'L')

#Add a new leaf
v2 = g.add_child(v2, edge_type='<', label='P')
leaf_new = g.add_child_and_complex(parent=stem,complex=phytomer, edge_type='+', label = 'L')


assert g.nb_vertices(scale=1) == 1
assert g.nb_vertices(scale=2) == 5
assert g.complex(v4) == plant1
