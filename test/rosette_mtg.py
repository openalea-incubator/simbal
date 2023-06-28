from openalea.mtg import *
from openalea.mtg.turtle import *
from openalea.plantgl import all as pgl
from openalea.plantgl.all import *
from openalea.plantgl.all import Vector3, Color3, Scaled, Viewer
from numpy import *

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


def vegetative_rapeseed(
    n_nodes: int = 10,
    node_length: float= 0.04,
    ):
    """ Build a plant with n_nodes internodes
    """

    above = g = MTG()
    vid = g.add_component(g.root, label='Plant')
    vid = g.add_component(vid, label='Axis')

    pid = g.add_component(vid, label='Phytomer')
    vid = g.add_component(pid, label='Internode')
    g.node(vid).NodeLength = node_length
    #g.node(vid).short = True
    leaf_id = g.add_child(vid, edge_type='+', label='Leaf')
    
    for i in range(n_nodes-1):
        pid = g.add_child(pid, edge_type='<', label='Phytomer')
        vid, pid = g.add_child_and_complex(vid, complex=pid, edge_type='<', label='Internode')
        g.node(vid).NodeLength = n_nodes
        #g.node(vid).short = True

        leaf_id = g.add_child(vid, edge_type='+', label='Leaf')

    return Plant(name='rapeseed', above = above)

def vegetative_rapeseed_DJ(
    DJ: int = 950,
    node_length: float= 0.02,
    phyllochrone = 0.015

    ):
    """ Build a plant with n_nodes internodes
    """
    n_nodes=round(DJ*phyllochrone)
    above = g = MTG()
    vid = g.add_component(g.root, label='Plant')
    vid = g.add_component(vid, label='Axis')

    pid = g.add_component(vid, label='Phytomer')
    vid = g.add_component(pid, label='Internode')
    g.node(vid).NodeLength = node_length
    #g.node(vid).short = True
    leaf_id = g.add_child(vid, edge_type='+', label='Leaf')
    
    for i in range(n_nodes-1):
        pid = g.add_child(pid, edge_type='<', label='Phytomer')
        vid, pid = g.add_child_and_complex(vid, complex=pid, edge_type='<', label='Internode')
        g.node(vid).NodeLength = node_length
        #g.node(vid).short = True

        leaf_id = g.add_child(vid, edge_type='+', label='Leaf')

    return Plant(name='rapeseed', above = above)



def update_MTG(g, outputs):
    """
    Update the MTG 
    :param dict outputs: {'param1': { vid1: , vid2, ...}, 'param2': { vid1: , vid2, ...}}
    """
     # add the missing property
    for param in outputs.keys():
        if param not in g.properties():
            g.add_property(param)

        # update the MTG
        g.property(param).update(outputs[param])   

    

def rapeseed_visitor(
        g, 
        v, 
        turtle,
        phyllo : int = 137.5,
        angle_insert : int = 70,
        coeff_leaf: int=20,
        coeff_width: int =0.831562,
        coeff_petiole: int =0.51017 
        ):
    
    #nb_phytomer = g.nb_vertices(scale=3)


    nid = g.node(v)
    label = g.label(v)
    surfaces=g.property('LeafSurface')
    #diam = g.property('diam') #a definir 
    turtle.rollL(phyllo)
    

    if g.edge_type(v) == '+':
        #short = False
        #if v in diam:
            #turtle.setWidth(diam[v]/2.)
        if label == 'Leaf':
            
            """
            #short = nid.parent().short

            #rank = g.Rank(g.parent(v))

            #s = (nb_phytomer-(0.3*v))
           # s = (nb_phytomer-rank)+1
            #s = max(1, s/4.)
            """
            bell_s=surfaces.get(v)*coeff_leaf
            petiole_length =bell_s*coeff_petiole 
            leafscale = Vector3(1, coeff_width*bell_s, bell_s)


        
        turtle.down(angle_insert)
        
        turtle.setColor(2)

        turtle.F(petiole_length)

        turtle.setColor(3)
        my_leaf = Scaled(leafscale, leafshape)
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
        visitor=rapeseed_visitor_v2, turtle=t, gc=False)

    return scene

def rapeseed_visitor_v2(
    g, 
    v, 
    turtle,
    ):
    
    #nb_phytomer = g.nb_vertices(scale=3)


    nid = g.node(v)
    label = g.label(v)
    turtle.rollL(137.5)
    

    if g.edge_type(v) == '+':
        
        turtle.setWidth(nid.PetioleDiam)
        if label == 'Leaf':
            
            leafscale = Vector3(1, nid.LeafWidth, nid.LeafLength)


        
        turtle.down(nid.InsertionAngle)
        
        turtle.setColor(2)

        turtle.F(nid.PetioleLength)

        turtle.setColor(3)
        my_leaf = Scaled(leafscale, leafshape)
        turtle.customGeometry(my_leaf)


    if label == 'Internode':
        turtle.setWidth(nid.StemDiam)
        turtle.setColor(1)
        turtle.F(nid.NodeLength)

    turtle.setId(v)




colors = dict()
colors['Internode'] = (100, 100, 80)
colors['Petiole'] = (100, 100, 80)
colors['Leaf'] = (9, 82, 40)


leafshape = NurbsPatch([[(-0.0241355, -0.0549754, -0.00804807, 1), (-0.00877065, -0.00933309, 0.00100194, 1), (0.0221632, 0.0405477, -0.00898557, 1), (-0.00223729, -0.00247516, 0.000385988, 1)], 
                        [(0.314084, -0.472921, 0.392759, 1), (0.000957513, 0.000400347, 0.33376, 1), (-0.00092555, 0.00517653, 0.332921, 1), (0.408856, 0.483731, 0.413612, 1)], 
                        [(0.329337, -0.501345, 0.730665, 1), (0.000326354, 0.0164005, 0.665257, 1), (-0.000118909, 0.019707, 0.66817, 1), (0.414001, 0.558038, 0.74373, 1)], 
                        [(0.00223079, 0.0223592, 0.999352, 1), (0.00112035, 0.028812, 1.0005, 1), (-0.000767699, 0.0327184, 0.999658, 1), (-0.0011293, 0.0361733, 1.00279, 1)]])



import numpy

def bell_shaped_dist(total_area=1, nb_phy=15, rmax=.7, skewness=5):
    """ returns leaf area of individual leaves along bell shaped model """

    r = numpy.linspace(1./nb_phy, 1, nb_phy)
    k = skewness
    relative_surface = numpy.exp(-k / rmax * ( 2 * (r - rmax)**2 + (r - rmax)**3))
    leaf_area = relative_surface / relative_surface.sum() * total_area
    return leaf_area.tolist()


def phenotype(g,
    leaf_max: int=0.5,
    skew: int=5,
    phyllot : int = 137.5,
    ins_angle : int = 70,
    coeff_leaf: int=20, #surface to length
    coeff_width: int =0.831562, #length to width
    coeff_petiole: int =0.31017, #length to petiole length 0.51017
    stem_d: int=0.1,
    coeff_petiole_d : int=0.5
    ):
    sortie=dict()
    leaves=[k for k, v in g.properties()['label'].items()  if v=='Leaf']
    internodes=[k for k, v in g.properties()['label'].items()  if v=='Internode']
    
    surface=bell_shaped_dist(1,len(leaves),leaf_max,skew)
    sortie['LeafSurface']=dict(zip(leaves,surface))

    sortie['LeafLength']=dict(zip(leaves,[element * coeff_leaf for element in surface]))
    sortie['LeafWidth']=dict(zip(leaves,[element * coeff_leaf * coeff_width for element in surface]))
    sortie['PetioleLength']=dict(zip(leaves,[element * coeff_leaf * coeff_petiole for element in surface]))
    
    sortie['InsertionAngle']=dict(zip(leaves,[ins_angle]*len(leaves)))
    sortie['Phyllotaxy']=dict(zip(leaves,[phyllot]*len(leaves)))

    sortie['StemDiam']=dict(zip(internodes,[stem_d]*len(internodes)))
    sortie['PetioleDiam']=dict(zip(leaves,[stem_d*coeff_petiole_d]*len(leaves)))

    update_MTG(g,sortie)



#1 create an MTG for a plant at a certain age (DD), for a given phyllochron (n_leaves/DD)
g=vegetative_rapeseed_DJ(950).above 
#2 attribute architectural properties to the MTG 

#3 generate and display a 3D representation of the architecture
#%gui qt
phenotype(g)
Viewer.display(scene3d(g))

""" OBSOLETE leaves=[k for k, v in g.properties()['label'].items()  if v=='Leaf']
ref=bell_shaped_dist(1,len(leaves),.7,5)
sortie=dict()
sortie['LeafSurface']=dict(zip(leaves,ref))
update_MTG(g,sortie)
"""

################ CARIBU

import openalea.plantgl.all as pgl
from alinea.caribu.CaribuScene import CaribuScene
from alinea.caribu.data_samples import data_path


#can = str(data_path('f331s1_100plantes.can'))
can = scene3d(g)
sky = str(data_path('Turtle16soc.light'))
#opts = list(map(str, [data_path('par.opt'), data_path('nir.opt')]))
pattern = str(data_path('filter.8'))
    
# complete set of files
cs = CaribuScene(scene=can, light=sky,pattern=pattern)
raw,agg=cs.run(simplify=True)

#runtime ~30s 

#scene,values = cs.plot(raw['par']['Ei'],display=False)

#cs.plot(agg['par']['Ei'],display=True)

cs.plot(raw['Ei'],display=True)
"""


import numpy as np
v99 = np.plenercentile(values, 99)
nvalues=np.array(values)
nvalues[nvalues>v99]=v99
values = nvalues.tolist()

from oawidgets.plantgl import *
PlantGL(scene, group_by_color=False, property=values) 

from openalea.plantgl.all import Viewer

#%gui qt
scene, values = cs.plot(raw['par']['Ei'], display=True)
"""
