""" plant.py

Define the topological structure, the vaious scales and growth and development functions.
This module is responsible of the morphogenesis program.

"""

from openalea.mtg import *

def aerial_plant():
    """ Return a seedling.
    This seedling is composed of 4 scales:
        - Plant: scale 1
        - Axis: scale 2
        - Metamer: scale 3
        - Organ: scale 4
            - S: stem
            - L: leaf
            - B: bud
    
    Note: 
        - < : successor
        + : ramification
    """
    
    g = MTG()
    vid = g.add_component(g.root, label='P')
    vid = g.add_component(vid, label='A')
    mid = g.add_component(vid, label='M')
    # organ
    stem_id = g.add_component(mid, label='S') 
    leaf = g.add_child(stem_id, label = 'L', edge_type='+')

    # Add a new metamer with a bud
    mid = g.add_child(mid, label='M', edge_type='<')
    bud = g.add_child_and_complex(parent=stem_id, complex=mid, edge_type='<', label='B')

    fat_mtg(g)

    return g

def simple_axis(n=10):
    """ Build a mtg composed of a sequence of metamers.
    """

    g = MTG()
    vid = g.add_component(g.root, label='S')
    lid = g.add_child(vid, label='L', edge_type='+')

    for i in range(n-1):
        vid = g.add_child(vid, label='S', edge_type='<')
        lid = g.add_child(vid, label='L', edge_type='+')

 
    return g

def complex_plant(n=10):
    g = MTG()
    vid = g.add_component(g.root, label='S')
    lid = g.add_child(vid, label='L', edge_type='+')

    for i in range(n-1):
        vid = g.add_child(vid, label='S', edge_type='<')
        lid = g.add_child(vid, label='L', edge_type='+')

    vids = [int(n/4), int(n/2)]

    for vid in vids:
        position = 1+vid*2
        vid = g.add_child(position, label='S', edge_type='+')
        lid = g.add_child(vid, label='L', edge_type='+')

        for i in range(2):
            vid = g.add_child(vid, label='S', edge_type='<')
            lid = g.add_child(vid, label='L', edge_type='+')

    return g

    
def all_scales(g):
    """ insert different scales from a simple MTG at one scale

    g1= all_scales(simple_axis())
    g2= all_scales(complex_plant())
    """   

    if g.max_scale() != 1:
        print('ERROR')
        return g
    
    labels =  set(g.property('label').values())
    if (len(labels) != 2) or (not labels.issuperset(set(('L','S')))):
        print('ERROR', labels)
        return g
    
    # Add metamers: define partitions
    def partition_metamer(v):
        return g.label(v)=='S'

    def partition_axis(v):
        return ((g.parent(v) is None) or (g.edge_type(v)=='+'))

    partition_plant = lambda v: g.parent(v) is None

    # insert scales
    g = g.insert_scale(inf_scale=1, partition=partition_metamer, default_label='M')
    g = g.insert_scale(inf_scale=1, partition=partition_axis, default_label='A')
    g = g.insert_scale(inf_scale=1, partition=partition_plant, default_label='P')
    return g


def init_N_randomly(g, N=10.):
    max_scale = g.max_scale()
    stems = [vid for vid in g.vertices(scale=max_scale) if g.label(vid) =='S']
    leaves = [vid for vid in g.vertices(scale=max_scale) if g.label(vid) =='L']

    ns = len(stems)
    nl = len(leaves)

    Natomic = N/(2*nl + ns)
    for vid in stems:
        nid = g.node(vid)
        nid.N = Natomic

    for vid in leaves:
        nid = g.node(vid)
        nid.N = 2*Natomic

    print('N total', sum(g.property('N').values()))
    print('N init', N)

    return g

def phytomer():
    """ Return a phytomer.
    Note: 
        - < : successor
        + : ramification
    """
    g = MTG()
    stem_id = g.add_component(g.root, label='Internode') 
    pet_id = g.add_child(stem_id, label = 'Petiole', edge_type='+')
    leaf_id = g.add_child(pet_id, label = 'Leaf', edge_type='<')
    bud = g.add_child(parent=stem_id, edge_type='<', label='Apex')
    fat_mtg(g)
    return g


###############################################################################

class Plant:
    """ Define a plant as an above and below ground MTG.

    Some usefull functions can be defined as proxies.
    
    """
    def __init__(self, above=None, below=None):
        
        self.above = self.init_above_mtg() if above is None else above
        self.below = self.init_below_mtg() if below is None else below


    @staticmethod
    def init_above_mtg():
        """ Simple function to create a simple aerial MTG"""
        g = MTG()
        
        # Add a plant (scale is 1)
        vid = g.add_component(g.root, label='Plant')
        vid = g.add_component(vid, label='Axis')
        vid = g.add_component(vid, label='Phytomer')
        
        # Last scale
        vid = g.add_component(vid, label='Internode')
        leaf_id = g.add_child(vid, edge_type='+', label='Leaf')
        bud_id = g.add_component(vid, edge-type='<', label='Bud')

        fat_mtg(g)
        return g

    @staticmethod
    def init_below_mtg():
        """ Simple function to create a simple belowground MTG or RootSystem Architecture"""
        g = MTG()
        
        # Add a plant (scale is 1)
        vid = g.add_component(g.root, label='RSA')
        vid = g.add_component(vid, label='Root')
        
        # Last scale
        vid = g.add_component(vid, label='Segment')

        fat_mtg(g)
        return g

    