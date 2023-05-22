from openalea.mtg import *
from openalea.mtg.turtle import *
from openalea.plantgl import all as pgl
from openalea.plantgl.all import *
from openalea.plantgl.all import Vector3, Color3, Scaled, Viewer

leafshape = NurbsPatch([[(-0.0241355, -0.0549754, -0.00804807, 1), (-0.00877065, -0.00933309, 0.00100194, 1), (0.0221632, 0.0405477, -0.00898557, 1), (-0.00223729, -0.00247516, 0.000385988, 1)], 
                        [(0.314084, -0.472921, 0.392759, 1), (0.000957513, 0.000400347, 0.33376, 1), (-0.00092555, 0.00517653, 0.332921, 1), (0.408856, 0.483731, 0.413612, 1)], 
                        [(0.329337, -0.501345, 0.730665, 1), (0.000326354, 0.0164005, 0.665257, 1), (-0.000118909, 0.019707, 0.66817, 1), (0.414001, 0.558038, 0.74373, 1)], 
                        [(0.00223079, 0.0223592, 0.999352, 1), (0.00112035, 0.028812, 1.0005, 1), (-0.000767699, 0.0327184, 0.999658, 1), (-0.0011293, 0.0361733, 1.00279, 1)]])

#feuille fine 
#leafshape = NurbsPatch([[(-0.0241355, -0.0549754, -0.00804807, 1), (0.0221632, 0.0405477, -0.00898557, 1), (-0.00223729, -0.00247516, 0.000385988, 1)], 
#                        
#                       [(0.329337, -0.501345, 0.730665, 1), (-0.000118909, 0.019707, 0.66817, 1), (0.414001, 0.558038, 0.74373, 1)], 
#                       [(0.00223079, 0.0223592, 0.999352, 1),  (-0.000767699, 0.0327184, 0.999658, 1), (-0.0011293, 0.0361733, 1.00279, 1)]])


colors = dict()
colors['Internode'] = (100, 100, 80)
colors['Petiole'] = (100, 100, 80)
colors['Leaf'] = (9, 82, 40)



def rapeseed_visitor(g, v, turtle):
    
    nb_phytomer = g.nb_vertices(scale=3)


    nid = g.node(v)
    label = g.label(v)
    phyllo = 137.5


    petiole_length =nb_phytomer/(v)

    turtle.rollL(phyllo)
    
    if g.edge_type(v) == '+':
        short = False
        if label == 'Leaf':
            short = nid.parent().short

            rank = g.Rank(g.parent(v))

            #s = (nb_phytomer-(0.3*v))
            s = (nb_phytomer-rank)+1
            s = max(1, s/4.)
            scale = Vector3(1, s, s)

            
            
            #scale = Vector3(1, s,2*s)

        if short:
            turtle.down(80)
        else:
            turtle.down(30)

        turtle.setColor(2)

        turtle.F(petiole_length)

        turtle.setColor(3)
        my_leaf = Scaled(scale, leafshape)
        turtle.customGeometry(my_leaf)


    if label == 'Internode':
        
        turtle.setColor(1)
        turtle.F(nid.length)

    turtle.setId(v)



def scene3d(g):
    t = pgl.PglTurtle()
    # set colors 
    color = Color3(*colors['Internode'])
    t.setColorAt(1, color)
    color = Color3(*colors['Petiole'])
    t.setColorAt(2, color)
    color = Color3(*colors['Leaf'])
    t.setColorAt(3, color)



    max_scale = g.max_scale()

    vid = next(g.component_roots_at_scale_iter(g.root, scale=max_scale))

    scene = turtle.traverse_with_turtle(
        g, vid, 
        visitor=rapeseed_visitor, turtle=t, gc=False)

    return scene



class Plant:
    """ Define a plant as an above and below ground MTG.

    Some usefull functions can be defined as proxies.
    
    """
    def __init__(self, name='', above=None, below=None):
        self.name = name
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
        bud_id = g.add_component(vid, edge_type='<', label='Bud')

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


# Create an MTG 
def vegetative_rapeseed(
    n_short: int = 10,
    n_long: int= 0,
    short_length: float= 0.04,
    long_length: float= 2.5,
    ):
    """ Build a plant with n_short internodes and n_long ones.
    """

    above = g = MTG()
    vid = g.add_component(g.root, label='Plant')
    vid = g.add_component(vid, label='Axis')

    pid = g.add_component(vid, label='Phytomer')
    vid = g.add_component(pid, label='Internode')
    g.node(vid).length = short_length
    g.node(vid).short = True
    leaf_id = g.add_child(vid, edge_type='+', label='Leaf')
    
    for i in range(n_short-1):
        pid = g.add_child(pid, edge_type='<', label='Phytomer')
        vid, pid = g.add_child_and_complex(vid, complex=pid, edge_type='<', label='Internode')
        g.node(vid).length = short_length
        g.node(vid).short = True

        leaf_id = g.add_child(vid, edge_type='+', label='Leaf')

    for i in range(n_long):
        pid = g.add_child(pid, edge_type='<', label='Phytomer')
        vid, pid = g.add_child_and_complex(vid, complex=pid, edge_type='<', label='Internode')
        g.node(vid).length = long_length
        g.node(vid).short = False
        leaf_id = g.add_child(vid, edge_type='+', label='Leaf')

    return Plant(name='rapeseed', above = above)

g=vegetative_rapeseed().above
#print(display(g))
#print(g.properties)

#Viewer.display(scene3d(g))