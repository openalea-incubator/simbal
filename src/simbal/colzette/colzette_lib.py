"""
Colzette Library
Clean python file containing the necessary functions for 'fullcycle.py' to operate
"""

from openalea.mtg import *
from openalea.mtg.turtle import *
from openalea.plantgl import all as pgl
from openalea.plantgl.all import *
from openalea.plantgl.all import Vector3, Color3, Scaled, Viewer
from numpy import *
import numpy as np

import openalea.plantgl.all as pgl

#color dict of the 3D elements
colors = dict()
colors['Internode'] = (100, 100, 80)
colors['Petiole'] = (100, 100, 80)
colors['Leaf'] = (9, 82, 40)

# 3D shape of the leaves

#Nurbs mesh that looks like a standard leaf 
#leafshape = NurbsPatch([[(-0.0241355, -0.0549754, -0.00804807, 1), (-0.00877065, -0.00933309, 0.00100194, 1), (0.0221632, 0.0405477, -0.00898557, 1), (-0.00223729, -0.00247516, 0.000385988, 1)], 
#                        [(0.314084, -0.472921, 0.392759, 1), (0.000957513, 0.000400347, 0.33376, 1), (-0.00092555, 0.00517653, 0.332921, 1), (0.408856, 0.483731, 0.413612, 1)], 
#                        [(0.329337, -0.501345, 0.730665, 1), (0.000326354, 0.0164005, 0.665257, 1), (-0.000118909, 0.019707, 0.66817, 1), (0.414001, 0.558038, 0.74373, 1)], 
#                        [(0.00223079, 0.0223592, 0.999352, 1), (0.00112035, 0.028812, 1.0005, 1), (-0.000767699, 0.0327184, 0.999658, 1), (-0.0011293, 0.0361733, 1.00279, 1)]])


#"circle" leaf, made of 8 triangles.
def unileaflet(length=1., width=1.):
    """Generates a unileaflet shape

    :param length: length of the shape, defaults to 1.
    :type length: float, optional
    :param width: width of the shape, defaults to 1.
    :type width: float, optional
    :return: the shape
    :rtype: pgl.Group
    """  
    disc = pgl.Translated((-0.5,0,0), pgl.Disc())
    disc = pgl.Scaled((length, width,1), disc)
    disc = pgl.AxisRotated(axis=(0,1,0), angle=radians(90.), geometry=disc)

    d3 = pgl.AxisRotated(axis=(1,0,0), angle=0., geometry=disc)

    shape = pgl.Group([d3])
    return shape

leafshape=unileaflet()


#MTG generation and edition

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
    """ Build an mtg of a single plant. The number of internodes is based on its age (DD) and phyllochrone 
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


##Mtg builder with mutliple plants
# functions from Adel used in vegetative_rapeseed_field()
def adel_regular(nb_plants, nb_rank, dx, dy, nx=None):
    
    if nx is None:
        nx = int(nb_plants/nb_rank)
    ny = nb_rank
    domain = ((0,0),(nx*dx, ny*dy))
    return [(i*dx+dx/2., j*dy+dy/2., 0.) for j in range(ny) for i in range(nx)], domain

def agronomicplot(nplants=None, length=None, width=None, sowing_density=None, plant_density=None, inter_row=None, noise=0, convunit=100, center_scene=True):
    """ Returns the number of plants, the positions, the domain (scene units), the domain area (square meter) and the conversion coefficient for meter to scene unit (1/convunit) of a micro-plot specified with agronomical variables
    nplants is the number of plants of the micro-plot (mandatory if length or width is None)
    length (m) is plot dimension along row direction (mandatory if nplants is None)
    width (m) is plot dimension perpendicular to row direction (mandatory if nplants is None)
    sowing density is the density of seeds sawn (mandatory)
    plant_density is the density of plants that are present (after loss due to bad emergence, early death...) (mandatory)
    inter_row (m) is for the  distance between rows (mandatory)
    noise (%), indicates the precision of the sowing for the inter plant spacing
    convunit is the conversion factor from meter to scene unit
    center_scene allows to center the position around origin. If False, the scene is in the x+,y+ sector, the origin being at the lower left corner of the domain
    
    Rows are parallel to x-axis
    Length and Width are adjusted to produce a canopy centered in its domain and compliant with infinitisation
    """
    from operator import itemgetter
    import random
    
    
    inter_plant = 1. / inter_row / sowing_density
    
    if nplants is None: # length and width are NOT None
        nrow = max(1, int(float(width) / inter_row))
        plant_per_row = max(1, int(float(length) / inter_plant))
        nplants = nrow * plant_per_row
        
    else: # nplants is NOT None ; length and width are None 
        import math
        nrow = max(1, int(math.sqrt(nplants)))
        
    dx = inter_plant * convunit
    dy = inter_row * convunit
    positions, domain = adel_regular(nplants, nrow, dx, dy)
    n_emerged = int(nplants * plant_density / sowing_density)
    positions = random.sample(positions, n_emerged)
    # sorting by ranks
    positions = sorted(positions, key= itemgetter(1,0))
    domain_area = abs(domain[1][0] - domain[0][0]) / convunit * abs(domain[1][1] - domain[0][1]) / convunit
    if center_scene:
        xc = float(domain[1][0] + domain[0][0]) / 2
        yc = float(domain[1][1] + domain[0][1]) / 2
        positions = [(x - xc, y - yc, z) for x, y, z in positions]
        domain = ((domain[0][0] - xc, domain[0][1] - yc), (domain[1][0] - xc, domain[1][1] - yc))
    
    return n_emerged, positions, domain, domain_area, 1. / convunit

def vegetative_rapeseed_field(
    DJ: int = 950,
    node_length: float= 0.02,
    phyllochrone = 0.015,
    coord = [(0,0,0)]
    
    ):
    """ Build a field with nrows(coord) plants that have n_nodes internodes

    """
    n_nodes=round(DJ*phyllochrone)
    #print('plantes à',n_nodes,'entrenoeuds')
    above = g = MTG()


    for pl in coord:
       
        #print(pl)
        vid = g.add_component(g.root, label='Plant',edge_type='/',position=pl) #vid = vertex_id
        vid = g.add_component(vid, label='Axis')

        pid = g.add_component(vid, label='Phytomer') #pid= phytomer_id
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

def bell_shaped_dist(total_area=1, nb_phy=15, rmax=.7, skewness=5):
    """ OLD, use generate_normalized_bell_curve_points() instead
    returns the leaf area of individual leaves along a bell shaped model, as a proportion of the total leaf area """

    r = np.linspace(1./nb_phy, 1, nb_phy)
    k = skewness
    relative_surface = np.exp(-k / rmax * ( 2 * (r - rmax)**2 + (r - rmax)**3))
    leaf_area = relative_surface / relative_surface.sum() * total_area
    return leaf_area.tolist()



def generate_normalized_bell_curve_points(mean, std_dev, amplitude, N):
    """ returns the leaf area of individual leaves along a bell shaped model, as a proportion of the total leaf area """
    x_values = np.linspace(mean - 3 * std_dev, mean + 3 * std_dev, N)
    y_values = amplitude * np.exp(-0.5 * ((x_values - mean) / std_dev)**2)
    
    # Normalize y values to ensure the sum is 1
    normalized_y_values = y_values / np.sum(y_values)
    
    points = list(normalized_y_values)
    return points



def phenotype(g,
    leaf_max: float=0.6587057,
    skew: float=5,  #
    sd: float= 0.246378, 
    amp: float=0.1188,
    phyllot : float = 137.5,
    ins_angle : float = 60,
    leaf_angle : float = 80, #à mesurer expérimentallement 
    total_surface: float=300, 
    coeff_width: float =0.831562, 
    coeff_petiole: float =0.51017, #limb length to petiole length , 0.31017 length/whole leag length
    stem_d: float=0.1,
    coeff_petiole_d : float=0.5
    ): 
    """
    applies phenotypical attributes to an MTG based on parameters. 
        leaf_max : [float between 0 and 1] Expectation of the bell curve => Position of the biggest leaf along the stem.
        sd: [float] Standard deviation of the bell curve => repartition of the leaf surface
        amp: [float] amplitude of the bell curve 
        phyllot: [float, degrees] phyllotaxy angle 
        ins_angle: [float, degrees] insertion angle of the petiole
        leaf_angle: [float, degrees] angle of the flat of the leaf
        total_surface: [int, cm2] surface to be shared following the bell proportion
        coeff_width: [float] Leaf Width / Leaf Length ratio (limb only)
        coeff_petiole: [float] Petiole Length / Leaf Length 
        stem_d: [float, cm] stem diameter
        coeff_petiole_d : [float] petiole diameter / stem diameter ration 
    """

    sortie=dict()
    plants=[k for k, v in g.properties()['label'].items()  if v=='Plant']
    leaves=[k for k, v in g.properties()['label'].items()  if v=='Leaf']
    internodes=[k for k, v in g.properties()['label'].items()  if v=='Internode']
    
    #surface=bell_shaped_dist(1,int(len(leaves)/len(plants)),leaf_max,skew) Old functio
    surface=generate_normalized_bell_curve_points(mean=leaf_max,std_dev=sd,amplitude=amp,N=int(len(leaves)/len(plants)))
    sortie['SurfaceRepartition']=dict(zip(leaves,surface*len(leaves)))
    sortie['LeafSurface']=dict(zip(leaves,[element * total_surface for element in surface]*len(leaves)))


    sortie['LeafLength']= dict(zip(leaves,
             [2*(element/(coeff_width*pi))**(0.5) for element in sortie['LeafSurface'].values()]*len(leaves)))

    sortie['LeafWidth']=dict(zip(leaves,
                     [element *coeff_width for element in sortie['LeafLength'].values()]*len(leaves)))
    
    sortie['SurfTheo']=dict(zip(leaves,[pi*(element/2)*(coeff_width*element/2) for element in sortie['LeafLength'].values()]*len(leaves)))

    sortie['PetioleLength']=dict(zip(leaves,[element * coeff_petiole for element in sortie['LeafLength'].values()]*len(leaves)))
    
    sortie['InsertionAngle']=dict(zip(leaves,[ins_angle]*len(leaves)))
    sortie['LeafAngle']=dict(zip(leaves,[leaf_angle]*len(leaves)))
    sortie['Phyllotaxy']=dict(zip(leaves,[phyllot]*len(leaves)))

    sortie['StemDiam']=dict(zip(internodes,[stem_d]*len(internodes)))
    sortie['PetioleDiam']=dict(zip(leaves,[stem_d*coeff_petiole_d]*len(leaves)))

    update_MTG(g,sortie)




###3D scene generation    

#Visitors
def rapeseed_visitor(
    g, 
    v, 
    turtle,
    ):
    """
    OpenAlea visitor
    This function is called by the scene3D() function when creating the 3D scene
    
    g is an MTG tree. In contains all the geometrical information needed by the visitor 
    v is the ID of the vertex being generated (int)
    turtle is the turtle used to create the 3D space
    """

    nid = g.node(v) #retrieves the current vertex based on its ID
    label = g.label(v) #retrieves the current vertex's label
    turtle.rollL(137.5) #rotates the turtle laterally by 137.5°
    

    if g.edge_type(v) == '+': #if '+', an element is added to the scene
        
        turtle.setWidth(nid.PetioleDiam)
        if label == 'Leaf':
            
            leafscale = Vector3(1, nid.LeafWidth, nid.LeafLength) #(bending,width,length) of a leaf


        
        turtle.down(nid.InsertionAngle) #turtle rotates vertically to match the insertion angle
        
        turtle.setColor(2)

        turtle.F(nid.PetioleLength) #turtle moves forward and draws a petiole

        turtle.setColor(3)
        my_leaf = Scaled(leafscale, leafshape) 
        turtle.customGeometry(my_leaf) #a leaf is generated at the end of the petiole


    if label == 'Internode': 
        turtle.setWidth(nid.StemDiam) 
        turtle.setColor(1)
        turtle.F(nid.NodeLength) #turtle moves forward and draws an internode

    turtle.setId(v) 


def ColzetteVisitor(
    g, 
    v, 
    turtle,
    ):
    """
    Visitor that can handle an MTG file with multiple plants 
    This function is called by the scene3D() function when creating the 3D scene
    
    g is an MTG tree. In contains all the geometrical information needed by the visitor 
    v is the ID of the vertex being generated (int)
    turtle is the turtle used to create the 3D space
    """
    
   
    #retrieve the node and its label 
    nid = g.node(v)
    label = g.label(v)

    if nid.parent() is None:#this is a new plant base
        p = nid.complex_at_scale(scale=1)
        #print('plant no',p)
        if 'position' in p.properties():
            #print (p.label, 'moving to ', p.position)
            turtle.move(list(map(float,p.position)))
        else:
            print('standing still')
            turtle.move(0,0,0)

        if 'azimuth' in p.properties():
            turtle.rollR(p.azimuth)
    
    turtle.rollL(137.5)

    if g.edge_type(v) == '+':
        
        turtle.setWidth(nid.PetioleDiam)
        if label == 'Leaf':
            
            leafscale = Vector3(1, nid.LeafWidth, nid.LeafLength)


        
        turtle.down(nid.InsertionAngle)
        
        turtle.setColor(2)

        turtle.F(nid.PetioleLength)

        turtle.setColor(3)
        
        turtle.down(nid.LeafAngle-nid.InsertionAngle)

        my_leaf = Scaled(leafscale, leafshape)
        turtle.customGeometry(my_leaf)


    if label == 'Internode':
        turtle.setWidth(nid.StemDiam)
        turtle.setColor(1)
        turtle.F(nid.NodeLength)

    turtle.setId(v)

#Scene creation

def scene3d(g,select_visitor=ColzetteVisitor):
    """
    calls a PlantGL turtle to generate a 3D scene
    g is an MTG file
    """
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

    scene = turtle.TurtleFrame(
        g, 
        visitor=select_visitor, turtle=t, gc=False, all_roots=True)
    
    return scene

### Display Results

def show_result(res):
    print ("Index of the triangles : ", res['index'])
    print ("The internal barcode : ", res['label'])
    print ("The individual areas of triangles : ", res['area'])
    print ("Surfacic density of energy : \n")
    print ("\t- absorbed by the triangles : ", res['Eabs'])
    print ("\t- incoming on the triangles :", res['Ei'])
    print ("\t- incoming on the inferior face of the triangle", res['Ei_inf'])
    print ("\t- incoming on the superior face of the triangle", res['Ei_sup'])
