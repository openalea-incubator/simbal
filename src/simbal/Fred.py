from openalea.mtg import *
from openalea.mtg.turtle import *
from openalea.mtg import turtle as turt
from openalea.mtg.plantframe import color
from openalea.plantgl import all as pgl
from openalea.plantgl.all import Vector3, Color3, Scaled

from numpy import *
import numpy as np
from pandas import *
import pandas as pd
#import os ?

# GENERAL FUNCTIONS ########################################################################################################################
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# This function enables to record the properties of each element of the MTG:
def recording_MTG_properties(g, file_name='g_properties.csv'):
    """This function records the properties of each node of the MTG "g" inside a csv file.
    
    Parameters
    ---------
    g: MTG 
        MTG where properties are recorded
    file_name: csv fil name
        the name of the csv file where properties of each node will be recorded
    
    Returns
    ------
    [no return]

    Authors
    -------
    Frederic REES

    """
    # We define the list of all properties of the MTG:
    initial_list_of_properties = list(g.properties().keys())
    # We verify that the list does not contain any "None" property:
    list_of_properties = []
    for property in initial_list_of_properties:
        if property is None:
            print("WARNING: a None property was detected in the MTG!")
            continue
        else:
            list_of_properties.append(property)
    # We sort the list by alphabetical order
    list_of_properties.sort()

    # We create an empty list of node indices:
    node_index = []
    # We create an empty list that will contain the properties of each node:
    g_properties = []

    # We cover all the vertices in the MTG:
    for vid in g.vertices_iter(scale=1):
        # Initializing an empty list of properties for the current node:
        node_properties = []
        # Adding the index at the beginning of the list:
        node_properties.append(vid)
        # n represents the vertex:
        n = g.node(vid)
        # For each possible property:
        for property in list_of_properties:
            # We add the value of this property to the list:
            node_properties.append(getattr(n, property, "NA"))
        # Finally, we add the new node's properties list as a new item in g_properties:
        g_properties.append(node_properties)
    # We create a list containing the headers of the dataframe:
    column_names = ['node_index']
    column_names.extend(list_of_properties)
    # We create the final dataframe:
    data_frame = pd.DataFrame(g_properties, columns=column_names)
    # We record the dataframe as a csv file:
    data_frame.to_csv(file_name, na_rep='NA', index=False, header=True)

    return

def create_MTG_from_csv_file(_init_, csv_filename):
    """This function use datas from the csv file to create each node of the MTG "g".
    
    Parameters
    ---------
    g: MTG 
        MTG where properties are recorded
    csv_file_name: csv fil name
        the name of the csv file where properties of each node will be recorded
    
    Returns
    ------
    [no return]

    AUTHORS
    -------
    Frederic REES

    """
    # We first read the csv file where all the properties of each vertex has been previously recorded:
    dataframe = pd.read_csv(csv_filename, sep=',', header=0)
    # We initialize an empty MTG:
    g = MTG()

    # We define the list of vertices and properties to be added to the MTG:
    list_of_vid = dataframe['node_index']
    list_of_properties = list(dataframe.columns.values)
    list_of_properties.remove('node_index')

    # We initialize the first vid, to which a new element will be added:
    vid = g.root

    for vid in list_of_vid:
        n = g.node(vid-1)
        n = n.add_child()
        for property in list_of_properties:
            g.property(property)[vid] = dataframe[property].loc[dataframe['node_index'] == vid]

    return g



# COLZETTE ########################################################################################################################
#color dict of the 3D elements
colors = dict()
colors['Internode'] = (100, 100, 80)
colors['Petiole'] = (100, 100, 80)
colors['Leaf'] = (9, 82, 40)
#####################################

def unileaflet(length=1., width=1.):
    """Generates a unileaflet shape (circle leaf) made with 8 triangles.
    
    Parameters
    ---------
    length: float
        length of the shape defaults to 1.
    width: float
        width of the shape, defaults to 1.
    
    Returns
    ------
    shape : pgl.Group

    AUTHORS
    -------
    Armand GRUMEL and Laurène PERTHAME
        Colzette model

    """
    disc = pgl.Translated((-0.5,0,0), pgl.Disc())
    disc = pgl.Scaled((length, width,1), disc)
    disc = pgl.AxisRotated(axis=(0,1,0), angle=radians(90.), geometry=disc)
    d3 = pgl.AxisRotated(axis=(1,0,0), angle=0., geometry=disc)

    shape = pgl.Group([d3])

    return shape

######################
leafshape = unileaflet()
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# This function enables to set a MTG with properties:
class Plant:
    """ Define a plant as an above and below ground MTG.

    Some usefull functions can be defined as proxies.

    Authors
    -------
    Arman GRUMEL and Laurène PERTHAME
        Colzette model

    """
    def __init__(self, name='Rapessed', above=None, below=None):
        self.name = name
        self.above = self.init_above_mtg() if above is None else above
        self.below = self.init_below_mtg() if below is None else below

    @staticmethod
    def init_above_mtg():
        """ Simple function to create a simple aerial MTG

        Authors
        -------
        Arman GRUMEL and Laurène PERTHAME
            Colzette model

        """
        g = MTG()
        
        # Add a plant (scale is 1)
        vid = g.add_component(g.root, label='Plant')
        vid = g.add_component(vid, label='Axis')
        vid = g.add_component(vid, label='Phytomer')
        
        # Last scale
        vid = g.add_component(vid, label='Internode')
        lid = g.add_child(vid, edge_type='+', label='Leaf')
        bid = g.add_component(vid, edge_type='<', label='Bud')

        fat_mtg(g)

        return g

    @staticmethod
    def init_below_mtg():
        """ Simple function to create a simple belowground MTG or RootSystem Architecture

        Authors
        -------
        Arman GRUMEL and Laurène PERTHAME
            Colzette model

        """
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

def agronomicplot(nplants=None, length=None, width=None, sowing_density=None, 
                  plant_density=None, inter_row=None, noise=0, convunit=100, center_scene=True):
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

def metabolism(g,
    fall : bool,
    photosynthesis = 0.1,
    nitrogen_content = 0.1,
    carbon_content= 0.1,
    potential_photosynthesis= 0.1,
    potential_area= 0.1,
    nitrogen_import= 0.1,
    nitrogen_export= 0.1,
    weight_import= 0.1,
    weight_export= 0.1,
    carbon_import= 0.1,
    carbon_export= 0.1,
    carbon_spend= 0.1,
    carbon_nitrogen_ratio= 0.1,
    nitrogen_storage= 0.1,
    N_fall_threshold_empirical= 0.1,
    PAR_fall_threshold_empirical= 0.1,
    N_fall_threshold_mecanism= 0.1): 
    """
    """
    sortie=dict()
    plants=[k for k, v in g.properties()['label'].items()  if v=='Plant']
    leaves=[k for k, v in g.properties()['label'].items()  if v=='Leaf']
    internodes=[k for k, v in g.properties()['label'].items()  if v=='Internode']

    #surface=bell_shaped_dist(1,int(len(leaves)/len(plants)),leaf_max,skew) Old functio
    
    sortie['NitrogenContent']=dict(zip(leaves,nitrogen_content))
    sortie['CarbonContent']=dict(zip(leaves,carbon_content))

    sortie['Fall']= dict(zip(leaves,fall))

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

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#DISPLAYING AND VIEW THE MTG:

g = vegetative_rapeseed().above
phenotype(g)
#print(display(g))
#print(g.properties)
#pgl.Viewer.display(scene3d(g))

# DEFINING METABOLIC FUNCTIONS ##############################################################################################################
class Photosynthesis:
    """Photosynthesis methods.

    Authors
    -------
    Solen FARRA

    """
    print("Checking the carbon produced of leaves according to incident light...")

    def __init__(self, g, photosynthesis = None) :

        self.photosynthesis = self.leaf_photosynthesis(g) if photosynthesis is None else photosynthesis

        return
    
    # This function calculates the carbon produced by a leaf:
    @classmethod
    def leaf_photosynthesis(self, g) :
        """Calculation of a potential (optimum) quantity of carbon produced of each organ.
    
        Parameters
        ---------
        g: MTG 
            MTG where properties are recorded
        
        alpha: MTG proporties
            quantical yield (slope coefficient) from the light response-curve of photosynthesis

        PAR: MTG proporties
            Incident light on the organ (W/m2)
        
        Respiration: MTG proporties
            Carbon quantity spend for respiration activity (g)
        
        Returns
        ------
        g: MTG 
            MTG where properties are recorded

        Authors
        -------
        Solen FARRA

        References
        ----------
        Drouet (1998)

        """

        max_scale = g.max_scale()

        leaves=[l for l, v in g.properties()['label'].items()  if v=='Leaf']
        alive = [l for l in leaves if not l.fall]

        for l in alive :
            l.photosynthesis = (l.alpha * l.PAR) - l.respiration

        return g
class Growth:
    """Growth methods. 

    Authors
    -------
    Solen FARRA

    """
    print("Checking the potentiality, demand and cost in nutrient of growth processus...")

    def __init__(self, g, potential_area = None, potential_photosynthesis = None, area = None) :
        
        self.g = g
        self.potential_area = self.leaf_potential_area(g) if potential_area is None else potential_area
        self.potential_photosynthesis = self.potential_photosynthesis(g) if potential_photosynthesis is None else potential_photosynthesis
        self.area = self.leaf_area(g) if area is None else area

        return 
    
    # This function calculates the potential growth of the leaf surface according to RUE:
    @classmethod
    def leaf_potential_area(self, g) :
        """Estimation of the potential leaf area based on a radiation use efficiency of each organ.
    
        Parameters
        ---------
        g: MTG 
            MTG where properties are recorded
    
        RUE: MTG proporties
            Radiation Use Efficiency coefficient

        PAR : MTG proporties
            Incident light on the organ (W/m2)
        
        Returns
        ------
        g: MTG 
            MTG where properties are recorded

        Authors
        -------
        Solen FARRA

        """
        g = self.g

        max_scale = g.max_scale()

        leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) 
                  if g.class_name(lid) =='Leaf']
        alive = [l for l in leaves if not l.leaf_fall]

        for l in alive :
            l.potential_area = l.RUE * l.PAR

        return g
    
    # This function calculates the potential net photosynthesis according to the potential leaf area:
    @classmethod
    def leaf_potential_photosynthesis(self, g,
                                      epsilon_a: float, epsilon_b: float ):
        """Calculation of a potential (optimum) quantity of carbon produced of each organ.
    
        Parameters
        ---------
        g: MTG 
            MTG where properties are recorded
        
        alpha: MTG proporties
            quantical yield (slope coefficient) from the light response-curve of photosynthesis

        PAR: MTG proporties
            Incident light on the organ (W/m2)
        
        epsilon_b: float
            Biological efficiency coefficient (g/MJ)

        epsilon_a: float
            Absorption efficiency coefficient

        
        Returns
        ------
        g: MTG 
            MTG where properties are recorded

        Authors
        -------
        Solen FARRA

        References
        ----------
        Monteith (1981)
            Penman-Monteith equation

        """
        g = self.g

        max_scale = g.max_scale()

        leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) 
                  if g.class_name(lid) =='Leaf']
        alive = [l for l in leaves if not l.leaf_fall]

        for l in alive :
            l.potential_photosynthesis = (l.PAR*epsilon_a)*epsilon_b

        return g
    
    @classmethod
    def leaf_area(self, g) :
        """leaf area extend according to final carbon content.
        
        Parameters
        ----------
        g: MTG
            the above part of the plant

        carbon_demand : MTG proporties
            Nitrogen content of each organs (g)

        potential_area : MTG proporties
            Potential area of each organs according to a radiation use efficiency (m2)

        carbon_growth : MTG proporties
            Carbon quantity left to growth after nitrogen exchange among leaves (g)

        carbon_content : MTG proporties 
            Carbon content of each organs (g)

        Returns
        -------
        g : MTG
            the above part of the plant

        """
        g = self.g

        max_scale = g.max_scale()

        leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) 
                if g.class_name(lid) =='Leaf']
        alive = [l for l in leaves if not self.leaf_fall]
        growing = [l for l in alive if self.growth]

        for l in growing :
            
            if self.carbon_demand < 0 :
                l.area = self.potential_area
            
            else :
                l.area = (self.carbon_growth*self.area)/self.carbon_content

        return g
    
class Competition:
    """Nutrient competition status methods. 

    Authors
    -------
    Solen FARRA

    """
    print("Checking whether leaves are in source or sink status in nitrogen or not...")

    def __init__(self, g, 
                 carbon_deficit = None, source_N = None, sink_N = None) :
    
        self.g = g
        self.carbon_deficit = self.leaf_carbon_deficit() if carbon_deficit is None else carbon_deficit
        self.source_N = self.leaf_nitrogen_competition_status() if source_N is None else source_N
        self.sink_N = self.leaf_nitrogen_competition_status() if sink_N is None else sink_N
    
        return 
    
    @classmethod
    def leaf_carbon_deficit(self, g):
        """Calculation of the carbon deficit according to the actual photosynthesis and the carbon stress for each organ.
    
        Parameters
        ---------
        g: MTG 
            MTG where properties are recorded
    
        carbon demand: MTG properties
            Actual quantity of carbon produced by photosynthesis (g)
        
        potential_area: MTG properties
            Carbon quantity spend for nitrogen redistribution among organs (g)

        area: MTG properties
            Actual leaf area 
        
        Returns
        ------
        g: MTG 
            MTG where properties are recorded

        Authors
        -------
        Solen FARRA

        """
        g = self.g

        leaves = [g.node(lid) for lid in g.vertices(scale=1) 
                  if g.class_name(lid) =='Leaf']
        alive = [l for l in leaves if not l.leaf_fall]

        for l in alive :
            l.carbon_deficit = l.photosynthesis - l.carbon_demand

        return g
    
    @classmethod
    def leaf_nitrogen_competition_status(self, g):
        """Evaluation of the nitrogen source-sink relationship between organs.
    
        Parameters
        ---------
        g: MTG 
            MTG where properties are recorded
    
        carbon demand: MTG properties
            Actual quantity of carbon produced by photosynthesis (g)
        
        potential_area: MTG properties
            Carbon quantity spend for nitrogen redistribution among organs (g)

        area: MTG properties
            Actual leaf area 
        
        Returns
        ------
        g: MTG 
            MTG where properties are recorded

        Authors
        -------
        Solen FARRA

        """
        print("Cheking the nitrogen source-sink relationship between organs") 

        g = self.g

        leaves = [g.node(lid) for lid in g.vertices(scale=1) 
                  if g.class_name(lid) =='Leaf']
        alive = [l for l in leaves if not l.leaf_fall]

        for l in alive :
            if l.carbon_deficit < 0 :
                l.source_N = True
                l.sink_N = False
            else :
                l.source_N = False
                l.sink_N = True
            
        return g

class Nitrogen:
    """Nitrogen exchange methods. 

    Authors
    -------
    Solen FARRA

    """
    print("Checking whether leaves are in senescence state or not...")
    print("Checking nitrogen quantities re-allocate between leaves...")
    print("Cheking leaf yellowing area evolution...")

    def __init__(self, g, nitrogen_demand = None, weight_export = None, 
               weight_import = None, nitrogen_import = None, nitrogen_export = None,
               nitrogen_storage = None, RUE = None, carbon_expend = None, nitrogen_content = None,
               carbon_nitrogen_ratio = None) :
    
        self.g = g

        self.nitrogen_demand = self.leaf_nitrogen_demand() if nitrogen_demand is None else nitrogen_demand
        self.weight_export = self.N_exchange_within_shoots() if weight_export is None else weight_export
        self.weight_import = self.N_exchange_within_shoots() if weight_import is None else weight_import
        self.nitrogen_import = self.N_exchange_within_shoots() if nitrogen_import is None else nitrogen_import
        self.nitrogen_export = self.N_exchange_within_shoots() if nitrogen_export is None else nitrogen_export
        self.nitrogen_storage = self.N_exchange_within_shoots() if nitrogen_storage is None else nitrogen_storage
        self.RUE = self.leaf_RUE() if RUE is None else RUE
        self.carbon_expend = self.leaf_senescence_carbon_cost() if carbon_expend is None else carbon_expend
        self.nitrogen_content = self.leaf_nitrogen_content if nitrogen_content is None else nitrogen_content
        self.carbon_nitrogen_ratio = self.leaf_carbon_nitrogen_ratio() if carbon_nitrogen_ratio is None else carbon_nitrogen_ratio
    
        return 
    

    # This function calculates the leaf nitrogen demand according to the potential leaf area:
    @classmethod
    def leaf_nitrogen_demand(self, g):
        """Calculation of the nitrogen demand of each organ.
    
        Parameters
        ---------
        g: MTG 
            MTG where properties are recorded
        
        potential_area: MTG proporties
            Potential organ area according to the RUE (m2)

        N_content : MTG proporties
            Actual nitrogen content of the oran (g)
        
        Returns
        ------
        g: MTG 
            MTG where properties are recorded

        Authors
        -------
        Solen FARRA

        """
        g = self.g

        leaves = [g.node(lid) for lid in g.vertices(scale=1) 
                  if g.class_name(lid) =='Leaf']
        alive = [l for l in leaves if not l.leaf_fall]

        for l in alive :
            l.nitrogen_demand = l.potential_area - l.nitrogen_content
        
        return g
    

    # This function calculates the nitrogen remobilised quantity according to ponderation,
    # the redistribution of N among leaves according to ponderation,
    # and the final leaf nitrogen content:
    @classmethod
    def N_exchange_within_shoots(self, g,
                                 PAR: float = None) :
        """Leaf senescence mecanism 
        
        Parameters
        ----------
        g: MTG
            the above part of the plant

        PAR : float
            Global incident light on each organs (W/m2)
        
        PAR : MTG proporties
            Incident light on the organ (W/m2)
            
        N_content : MTG proporties
            Nitrogen content of each organs (g)

        nitrogen_demand : MTG proporties
            Nitrogen content of each organs (g)

        Returns
        -------
        g : MTG
            the above part of the plant

        References
        ----------
        Lemaire et al. (2008)

        Examples
        --------
        >>> leaf_senescence()

        """
        g = self.g

        leaves = [g.node(lid) for lid in g.vertices(scale=1) 
                if g.class_name(lid) =='Leaf']
        alive = [l for l in leaves if not l.leaf_fall]
        senescence = [l for l in alive if l.source_N]

        for l in senescence :
            l.weight_export = l.PAR/PAR
            l.nitrogen_exported = l.nitrogen_content * l.weight_export
            

            pool_N = sum(l.N_exported)
        
            for l in alive :
                l.weight_import = l.PAR/PAR

                if l.nitrogen_demand < 0 :
                    l.nitrogen_import = pool_N * l.weight_import

                else :
                    l.nitrogen_storage = l.nitrogen_demand
                    l.nitrogen_content = l.nitrogen_content

        return g

    # This function calculates the final radiation use efficiency leaf coefficient:
    @classmethod
    def leaf_RUE(self, g):
        """Calculation of the radiation use efficiency for each organs.
    
        Parameters
        ---------
        g: MTG 
            MTG where properties are recorded
    
        PAR : MTG proporties
            Incident light on the organ (W/m2)
            
        N_content : MTG proporties
            Nitrogen content of each organs (g) 
        
        Returns
        ------
        g: MTG 
            MTG where properties are recorded

        Authors
        -------
        Solen FARRA

        References
        ----------
        Luquet et al. (2006)
            EcoMeristem model

        """
        g = self.g

        leaves = [g.node(lid) for lid in g.vertices(scale=1) 
                if g.class_name(lid) =='Leaf']
        alive = [l for l in leaves if not l.leaf_fall]
        
        for l in alive : 
            l.RUE = l.nitrogen_content/l.PAR

        return g
    
    # This function calculates the carbon cost of the redistribution of N among leaves:
    @classmethod
    def leaf_senescence_carbon_cost(self, g,
                                     kr : float = 0.08):
        """Calculation of the carbon cost of nitrogen exchange among organs.
    
        Parameters
        ---------
        g: MTG 
            MTG where properties are recorded
            
        N_import : MTG proporties
            Nitrogen quantity import of each organs (g) 
        
        N_content : MTG proporties
            Actual nitrogen content of each organs (g)
        
        kr : float
            Carbon cost coefficient of nitrogen resorption for each organs (kg/C/m2)
        
        Returns
        ------
        g: MTG 
            MTG where properties are recorded

        Authors
        -------
        Solen FARRA

        References
        ----------
        Fisher et al. (2010)
        Brzostek et al. (2014)
            FUN 2.0 model

        """
        g = self.g

        leaves = [g.node(lid) for lid in g.vertices(scale=1) 
                if g.class_name(lid) =='Leaf']
        alive = [l for l in leaves if not l.leaf_fall]

        for l in alive : 
            l.carbon_expend = (kr/l.nitrogen_content)*l.nitrogen_import
            
        return g
    
    @classmethod
    def leaf_nitrogen_content(self, g):
        """Calculation of the final nitrogen content after exchange among organs.
    
        Parameters
        ---------
        g: MTG 
            MTG where properties are recorded
            
        N_import : MTG proporties
            Nitrogen quantity import of each organs (g) 
        
        N_export : MTG proporties
            Nitrogen quantity export of each organs (g) 
        
        Returns
        ------
        g: MTG 
            MTG where properties are recorded

        Authors
        -------
        Solen FARRA

        References
        ----------
        Fisher et al. (2010)
        Brzostek et al. (2014)
            FUN 2.0 model

        """
        g = self.g

        leaves = [g.node(lid) for lid in g.vertices(scale=1) 
                if g.class_name(lid) =='Leaf']
        alive = [l for l in leaves if not l.leaf_fall]
        senescence = [l for l in alive if l.source_N]

        for l in alive : 
            l.nitrogen_content = l.N_content - l.N_import

        for l in senescence :
            l.nitrogen_content = l.N_content - l.N_export

    @classmethod
    def leaf_carbon_nitrogen_ratio(self, g):
        """Evaluation of the C:N ratio of each organ for growth regulation and shoot:root balance.
    
        Parameters
        ---------
        g: MTG 
            MTG where properties are recorded
            
        N_content : MTG proporties
           Final nitrogen content of each organs (g)

        carbon_growth : MTG properties
            Carbon quantity left to growth after nitrogen exchange among organs.
        
        Returns
        ------
        g: MTG 
            MTG where properties are recorded

        Authors
        -------
        Solen FARRA

        References
        ----------
        Fisher et al. (2010)
        Brzostek et al. (2014)
            FUN 2.0 model

        """
        g = self.g

        leaves = [g.node(lid) for lid in g.vertices(scale=1) 
                if g.class_name(lid) =='Leaf']
        alive = [l for l in leaves if not l.leaf_fall]

        for l in alive : 
            l.carbon_nitrogen_ratio = l.carbon_growth/l.N_import

        return g
    
class Carbon:
    """Carbon balance methods. 

    Authors
    -------
    Solen FARRA

    """
    def __init__(self, g, carbon_growth = None, carbon_demand = None, carbon_export = None, 
               carbon_import = None, growth = None) :

        self.g = g
        self.carbon_growth = self.leaf_carbon_growth() if carbon_growth is None else carbon_growth
        self.carbon_demand = self.leaf_carbon_demand() if carbon_demand is None else carbon_demand
        self.carbon_export = self.C_exchange_within_shoots() if carbon_export is None else carbon_export
        self.carbon_import = self.C_exchange_within_shoots() if carbon_import is None else carbon_import
        self.growth = self.C_exchange_within_shoots() if growth is None else growth
            
        return 
    
    # This function calculates the carbon finally available for growth after nutrient redistribution:
    @classmethod
    def leaf_carbon_growth(self, g) :
        """Calculation of a potential (optimum) quantity of carbon produced of each organ.
    
        Parameters
        ---------
        g: MTG 
            MTG where properties are recorded
    
        photosynthesis: MTG properties
            Actual quantity of carbon produced by photosynthesis (g)
        
        carbon_spend: MTG properties
            Carbon quantity spend for nitrogen redistribution among organs (g)
        
        Returns
        ------
        g: MTG 
            MTG where properties are recorded

        Authors
        -------
        Solen FARRA

        """
        g = self.g

        leaves = [g.node(lid) for lid in g.vertices(scale=1) 
                  if g.class_name(lid) =='Leaf']
        alive = [l for l in leaves if not l.leaf_fall]

        for l in alive :
            l.carbon_growth = l.photosynthesis - l.carbon_spend
        
        return g
    
    # This function calculates the carbon demand according to the potential carbon produced associated with the potential leaf area:
    @classmethod
    def leaf_carbon_demand(self, g) :
        """Calculation of a potential (optimum) quantity of carbon produced of each organ.
    
        Parameters
        ---------
        g: MTG 
            MTG where properties are recorded
        
        potential_photosynthesis: MTG properties
            Potential (optimum) quantity of carbon produced by photosynthesis according to a potential leaf area (g)
        
        carbon_growth: MTG properties
            Carbon quantity left for growth of each organs (g)
        
        Returns
        ------
        g: MTG 
            MTG where properties are recorded

        Authors
        -------
        Solen FARRA

        """
        g = self.g

        leaves = [g.node(lid) for lid in g.vertices(scale=1) 
                  if g.class_name(lid) =='Leaf']
        alive = [l for l in leaves if not l.leaf_fall]
    
        for l in alive :
            l.carbon_demand = l.potential_photosynthesis - l.carbon_growth

        return g
    
    @classmethod
    def C_exchange_within_shoots(self, g) :
        """Carbon storage mecanism and carbon stress for growth state for each organs.
        
        Parameters
        ----------
        g: MTG
            the above part of the plant

        C_storage: MTG properties
            Carbon quantity storage after carbon balance
        
        weight_export: MTG properties
            Ponderation coefficient of each organs for nutrient pool export

        weight_import: MTG properties
            Ponderation coefficient of each organs for nutrient pool import

        carbon_growth : MTG proporties
            Carbon quantity left to growth after nitrogen exchange among leaves (g)

        Returns
        -------
        g : MTG
            the above part of the plant

        References
        ----------
        Lemaire et al. (2008)

        """
        g = self.g

        leaves = [g.node(lid) for lid in g.vertices(scale=1) 
                if g.class_name(lid) =='Leaf']
        alive = [l for l in leaves if not l.leaf_fall]
        sink = [l for l in alive if l.sink_N]

        for l in alive :
            l.carbon_export = l.carbon_storage*l.weight_export
            
            pool_C = sum(l.carbon_export)
        
        for l in sink :

            if l.carbon_growth < 0 :
                l.carbon_import = pool_C*l.weight_import
                l.carbon_growth = l.carbon_growth + l.carbon_import

                if l.carbon_growth > 0 :
                    l.growth = True

                else :
                    l.growth = False

            else : l.growth = True    

        return g

    
class Fall:
    """Leaves falling methods. 

    Authors
    -------
    Solen FARRA

    """
    print("Checking whether leaves should fall or not...")
    print("Checking nitrogen and carbon quantity thresholds of leaves falling...")

    def __init__(self, g, fall = None, nitrogen_fall_threshold = None, nitrogen_fall_threshold_empirical = False,
               PAR_fall_threshold_empirical = False, carbon_fall_threshold = None) :

        self.g = g
        self.fall = self.leaf_falling() if fall is None else fall
        self.nitrogen_fall_threshold_empirical = self.leaf_falling_empirical_threshold() if nitrogen_fall_threshold_empirical is None else nitrogen_fall_threshold
        self.nitrogen_fall_threshold = self.leaf_falling_mecanism_threshold() if nitrogen_fall_threshold is None else nitrogen_fall_threshold_empirical
        self.PAR_fall_threshold_empirical = self.leaf_falling_empirical_threshold() if PAR_fall_threshold_empirical is None else carbon_fall_threshold
        self.carbon_fall_threshold = self.leaf_falling_mecanism_threshold if carbon_fall_threshold is None else PAR_fall_threshold_empirical
            
        return 

    # This function considers whether a leaf should fall or not:
    @classmethod
    def leaf_falling(self, g):
        """Leaf fall mecanism.
    
        Parameters
        ---------
        g: MTG 
            MTG where properties are recorded
    
        carbon_growth : MTG proporties
            Carbon left for growth after nitrogen exchange accoring to photosynthesis of each organs (g) 
        
        Returns
        ------
        g: MTG 
            MTG where properties are recorded

        Authors
        -------
        Solen FARRA

        References
        ----------
        Field (1983)

        """
        g = self.g

        leaves = [g.node(lid) for lid in g.vertices(scale=1) 
                if g.class_name(lid) =='Leaf']
        alive = [l for l in leaves if not l.fall]
        senescence = [l for l in alive if l.source_N]
        
        for l in senescence :
            if l.carbon_growth < 0:
                l.fall = True

            else:
                l.fall = False

        return g
    
    # This function considers the final carbon and nitrogen mecanistical threshold for leaf fall:
    @classmethod
    def leaf_falling_mecanism_threshold(self, g):
        """Leaf fall nutrient threshold estimation.
    
        Parameters
        ---------
        g: MTG 
            MTG where properties are recorded
    
        N_content : MTG proporties
            Actual nitrogen content of each organs (g)
        
        C_content : MTG propoerties
            Actual carbon content of each organs (g)

        Returns
        ------
        g: MTG 
            MTG where properties are recorded

        Authors
        -------
        Solen FARRA

        References
        ----------
        Field (1983)

        """
        g = self.g

        leaves = [g.node(lid) for lid in g.vertices(scale=1) 
                if g.class_name(lid) =='Leaf']
        dead = [l for l in leaves if l.fall]

        for l in dead :
            l.C_fall_threshold = l.C_content
            l.N_fall_threshold = l.N_content

        return g
    
    # This function considers the final nitrogen and light empirical threshold for leaf fall:
    @classmethod
    def leaf_falling_empirical_threshold(self, g):
        """Leaf fall empirical (light) threshold estimation.
    
        Parameters
        ---------
        g: MTG 
            MTG where properties are recorded
    
        N_content : MTG proporties
            Actual nitrogen content of each organs (g)

        PAR : MTG proporties
            Incident light on the organ (W/m2)
        
        Returns
        ------
        g: MTG 
            MTG where properties are recorded

        Authors
        -------
        Solen FARRA

        References
        ----------
        Thomas (2013)

        """
        g = self.g

        leaves = [g.node(lid) for lid in g.vertices(scale=1) 
            if g.class_name(lid) =='Leaf']
        dead = [l for l in leaves if l.fall]

        for l in dead :
            l.PAR_fall_threshold = l.PAR
            l.N_fall_threshold = l.N_content

        return g
    
# MAIN SIMULATION ##############################################################################################################################
class plot_MTG:
    
    def __init__(self, csv_filename = False, input_dataframe = None) :
        
        self.csv_filename = self.create_MTG_from_csv_file() if csv_filename is None else input_dataframe
        self.input_dataframe = self.create_MTG_from_scratch() if input_dataframe is None else csv_filename

        return 

    # OPTION 1: We create a MTG from its known properties (NOT COMPLETELY FINISHED!):
    @staticmethod
    def create_MTG_from_csv_file():
        
        csv_filename='MTG_00003.csv'

        # We first read the csv file where all the properties of each vertex has been previously recorded:
        dataframe = pd.read_csv(csv_filename, sep=',', header=0)
    
        # We initialize an empty MTG:
        g = MTG()

        # We define the list of vertices and properties to be added to the MTG:
        list_of_vid = dataframe['node_index']
        list_of_properties = list(dataframe.columns.values)
        list_of_properties.remove('node_index')

        # We cover each new vertex to be added to the MTG:
        for vid in list_of_vid:
            # We get the previous element on which to add the new element:
            n = g.node(vid-1)
            # We now define the new element as the child of the previous element:
            n = n.add_child()
            # For this new element, we cover all the properties defined in the csv file:
            for property in list_of_properties:
                # For the specific property, we get the single value corresponding to the current vertex ID:
                selected_series = dataframe[property].loc[dataframe['node_index'] == vid].iloc[0]
                # NB: we need to add 'iloc[0]' as the previous expression otherwise returns a series, not a single value.
        
        return g

    # OPTION 2: We create before a generic MTG:
    @staticmethod
    def create_MTG_from_scratch():
    # We read the environemental/initial conditions:
        #input_dataframe = pd.read_excel(r'C:\Users\FARRA\AppData\Local\GitHubDesktop\app-3.2.3\src\caracol\src\caracol\APTOS_input_1.xlsx', header=0, sheet_name="scenario")

        #We plot the MTG:
        g = vegetative_rapeseed().above
        g = phenotype(g)

        # We initialize the content of N in each leaf/internode:
        #for vid in g.vertices_iter(scale=1):
            #n = g.node(vid)
            # n.N_content = input_dataframe['Initial_N_content'].loc[vid-1]
            #n.N_content =input_dataframe['Initial_N_content'].loc[input_dataframe['element_ID'] == vid]

        # We initialize the PAR in each leaf/internode:
        #for vid in g.vertices_iter(scale=1):
            #n = g.node(vid)
            #n.PAR = input_dataframe['PAR'].loc[input_dataframe['element_ID'] == vid]

        # We initialize a new property of the MTG g that corresponds to the status "Fallen / Not fallen":
        #for vid in g.vertices_iter(scale=1):
            #n = g.node(vid)
            #n.fallen = False

        # We record the properties in its initial state:
        recording_MTG_properties(g, file_name = 'MTG_0000.csv')

        return g


class Simulation:
    # We define an initial and final step for covering a certain period:
    def __init__(self, g, initial_step = 1, final_step = 3) :
        
        self.g = g
        self.initial_step = initial_step
        self.final_step = final_step 

    @classmethod
    def main(self, g, initial_step = 1, final_step = 3):

        self.g = g

        # For each time step covering the simulation period:
        for step in range(initial_step, final_step+1):
            print("")
            print("For step", step,":")

            # We successively call the different class and methods:
            Main_photosynthesis = Photosynthesis(g, photosynthesis = None)
            
            Main_growth = Growth(g, potential_photosynthesis = None,
                                 potential_area = None, 
                                 area = None)
            
            Main_competition = Competition(g, carbon_deficit = None, 
                                           source_N = None, 
                                           sink_N = None)
            
            Main_nitrogen = Nitrogen(g, nitrogen_demand = None, 
                                     weigth_import = None, 
                                     weigth_export = None,
                                     nitrogen_import =None, 
                                     nitrogen_export = None,
                                     nitrogen_storage = None,
                                     carbon_expend = None, 
                                     nitrogen_content = None,
                                     RUE = None, 
                                     carbon_nitrogen_ratio = None)
            
            Main_carbon = Carbon(g, carbon_growth = None, 
                                 carbon_demand = None,
                                 carbon_import = None, 
                                 carbon_export = None)
            
            Main_fall = Fall(g, fall = None,  
                             nitrogen_fall_threshold = None, 
                             nitrogen_fall_threshold_empirical = False,
                             PAR_fall_threshold_empirical = False, 
                             carbon_fall_threshold = None)

            g = Main_photosynthesis.leaf_photosynthesis(g)
            
            g = Main_growth.leaf_potential_photosynthesis(g)
            g = Main_growth.leaf_potential_area(g)
            
            g = Main_competition.leaf_carbon_deficit(g)
            g = Main_competition.leaf_nitrogen_competition_status(g)
            
            g = Main_nitrogen.leaf_nitrogen_demand(g)
            g = Main_nitrogen.N_exchange_within_shoots(g)
            g = Main_nitrogen.leaf_senescence_carbon_cost(g)
            g = Main_nitrogen.leaf_nitrogen_content(g)
            g = Main_nitrogen.leaf_RUE(g)
            g = Main_nitrogen.leaf_carbon_nitrogen_ratio(g)

            g = Main_carbon.leaf_carbon_growth(g)
            g = Main_carbon.leaf_carbon_demand(g)
            g = Main_carbon.C_exchange_within_shoots(g)

            g = Main_growth.leaf_area(g)

            g = Main_fall.leaf_falling(g)
            g = Main_fall.leaf_falling_empirical_threshold(g)
            g = Main_fall.leaf_falling_mecanism_threshold(g)

            # We record the main properties of the plant:
            properties_filename = 'MTG_%.5d.csv' % step
            recording_MTG_properties(g, properties_filename)

            return g

        print("")
        print("Simulation is done!")

scenario_0 = Simulation(g, initial_step = 1.0, final_step = 3.0)
scenario_0.main(g)

pgl.Viewer.display(scene3d(g))
