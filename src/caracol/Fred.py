
from openalea.mtg import *
from openalea.mtg.turtle import *
from openalea.mtg import turtle as turt
from openalea.mtg.plantframe import color
from openalea.plantgl import all as pgl
from openalea.plantgl.all import Vector3, Color3, Scaled, Viewer

from numpy import *
import numpy as np
from pandas import *
import pandas as pd
import os

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

    AUTHORS
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

def create_MTG_from_csv_file(csv_filename):
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
    Armand GRUMEL (colzette module)

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
    
    """
    def __init__(self, name='Rapessed', above=None, below=None):
        self.name = name
        self.above = self.init_above_mtg() if above is None else above
        self.below = self.init_below_mtg() if below is None else below


    @staticmethod
    def init_above_mtg():
        """ Simple function to create a simple aerial MTG
        
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
    n_nodes: int = 6,
    length_node: float= 2.0):
    """ Build a rapeseed plant with n_nodes internodes and other properties

    """
    above = g = MTG()

    vid = g.add_component(g.root, label='Plant')
    vid = g.add_component(vid, label='Axis')
    pid = g.add_component(vid, label='Phytomer')
    vid = g.add_component(pid, label='Internode')
    
    g.node(vid).NodeLength = length_node 
    #g.node(vid).short = True

    lid = g.add_child(vid, 
                          edge_type='+', 
                          label='Leaf')
    
    for i in range(n_nodes-1):
        pid = g.add_child(pid, 
                          edge_type='<', 
                          label='Phytomer')
        vid, pid = g.add_child_and_complex(vid, complex=pid, edge_type='<', label='Internode')
        
        g.node(vid).NodeLength = n_nodes
        #g.node(vid).short = True

        lid = g.add_child(vid, edge_type='+', label='Leaf')

    return Plant(name='rapeseed', above = above)

def update_MTG(g, outputs):
    """Update the MTG 

    Parameters
    ----------
    outputs: dict
        {'param1': { vid1: , vid2, ...}, 'param2': { vid1: , vid2, ...}}

    """
     # add the missing property
    for param in outputs.keys():
        if param not in g.properties():
            g.add_property(param)

        # update the MTG
        g.property(param).update(outputs[param])  


def generate_normalized_bell_curve_points(mean, std_dev, amplitude, N):
    """ returns the leaf area of individual leaves along a bell shaped model, as a proportion of the total leaf area 
    
    """
    x_values = np.linspace(mean - 3 * std_dev, mean + 3 * std_dev, N)
    y_values = amplitude * np.exp(-0.5 * ((x_values - mean) / std_dev)**2)
    
    # Normalize y values to ensure the sum is 1
    normalized_y_values = y_values / np.sum(y_values)
    
    points = list(normalized_y_values)
    return points

def phenotype(g,
    leaf_max: float=0.6587057,
    skew: float=5,
    sd: float= 0.246378, 
    amp: float=0.1188,
    phyllot : float = 137.5,
    ins_angle : float = 60,
    leaf_angle : float = 80, 
    total_surface: float=300, 
    coeff_width: float =0.831562, 
    coeff_petiole: float =0.51017, #limb length to petiole length , 0.31017 length/whole leag length
    stem_d: float=0.1,
    coeff_petiole_d : float=0.5): 
    """applies phenotypical attributes to an MTG based on parameters. 
        
    Parameters
    ----------
    leaf_max : float [between 0 and 1]
        Expectation of the bell curve => Position of the biggest leaf along the stem.
    sd: float
        Standard deviation of the bell curve => repartition of the leaf surface
    amp: float 
        amplitude of the bell curve 
    phyllot: float [degrees] 
        phyllotaxy angle 
    ins_angle: float[degrees] 
        insertion angle of the petiole
    leaf_angle: float [degrees] 
        angle of the flat of the leaf
    total_surface: int [cm2] 
        surface to be shared following the bell proportion
    coeff_width: float
        Leaf Width / Leaf Length ratio (limb only)
    coeff_petiole: float
        Petiole Length / Leaf Length 
    stem_d: float [cm] 
        stem diameter
    coeff_petiole_d : float 
        petiole diameter / stem diameter ration 

    """
    outputs=dict()

    plants=[k for k, v in g.properties()['label'].items()  if v=='Plant']
    leaves=[k for k, v in g.properties()['label'].items()  if v=='Leaf']
    internodes=[k for k, v in g.properties()['label'].items()  if v=='Internode']

    surface=generate_normalized_bell_curve_points(mean=leaf_max,std_dev=sd,amplitude=amp,N=int(len(leaves)/len(plants)))
   
    outputs['SurfaceRepartition']= dict(zip(leaves,
                                            surface*len(leaves)))
    outputs['LeafSurface']=dict(zip(leaves,
                                    [element * total_surface for element in surface]*len(leaves)))

    outputs['LeafLength']= dict(zip(leaves,
             [2*(element/(coeff_width*pi))**(0.5) for element in outputs['LeafSurface'].values()]*len(leaves)))
    
    outputs['LeafWidth']=dict(zip(leaves,
                     [element *coeff_width for element in outputs['LeafLength'].values()]*len(leaves)))
    
    outputs['SurfTheo']=dict(zip(leaves,
                                 [pi*(element/2)*(coeff_width*element/2) for element in outputs['LeafLength'].values()]*len(leaves)))
    
    outputs['PetioleLength']=dict(zip(leaves,
                                      [element * coeff_petiole for element in outputs['LeafLength'].values()]*len(leaves)))
    
    outputs['InsertionAngle']=dict(zip(leaves,
                                       [ins_angle]*len(leaves)))
    
    outputs['LeafAngle']=dict(zip(leaves,
                                  [leaf_angle]*len(leaves)))
    
    outputs['Phyllotaxy']=dict(zip(leaves,
                                   [phyllot]*len(leaves)))
    
    outputs['StemDiam']=dict(zip(internodes,
                                 [stem_d]*len(internodes)))
    
    outputs['PetioleDiam']=dict(zip(leaves,
                                    [stem_d*coeff_petiole_d]*len(leaves)))
    
    update_MTG(g, outputs) 

    
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# This function enables to plot an above MTG with PlantGL:

def rapeseed_visitor(g, v, turtle):
    """OpenAlea visitor
    This function is called by the scene3D() function when creating the 3D scene
    
    g is an MTG tree. In contains all the geometrical information needed by the visitor 
    v is the ID of the vertex being generated (int)
    turtle is the turtle used to create the 3D space

    """
    nid = g.node(v) #retrieves the current vertex based on its ID
    label = g.label(v) #retrieves the current vertex's label
    turtle.rollL(137.5) #rotates the turtle laterally by 137.5° = leaf insertion angle
    

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


def scene3d(g, select_visitor=rapeseed_visitor):
    """calls a PlantGL turtle to generate a 3D scene g is an MTG file
    
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
        visitor=select_visitor, turtle=t, gc=False, all_roots=False)
    
    return scene

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#DISPLAYING AND VIEW THE MTG:
g=vegetative_rapeseed().above
phenotype(g)
#print(display(g))
#print(g.properties)
pgl.Viewer.display(scene3d(g))


# DEFINING METABOLIC FUNCTIONS ##############################################################################################################
# This function considers whether a leaf should start senescing or not:

class Growth:

    def _init_(self, g) :
        """Leaf scale turnover module initial parameters
          
        """
        for vid in g.vertices(scale=1):
            l = g.node(vid)

        self = l 
            
        return self
    
    def leaf_potential_area(_init_, g) :
        return
    
    def leaf_nitrogen_demand(_init_, g) :
        return
    
    def leaf_potential_photosynthesis(_init_, g) :
        return
    
    def leaf_actual_growth(_init_, g) :
        return
    
    def leaf_carbon_demand(_init_, g) :
        return
    
    def leaf_actual_area(_init_, g) :
        return
    

class Competition:
    print("Checking whether leaves are in source or sink status in carbon or not...")
    print("Checking whether leaves are in source or sink status in nitrogen or not...")

    def _init_(self, g) :
        """Leaf scale turnover module initial parameters
          
        """
        for vid in g.vertices(scale=1):
            l = g.node(vid)

        self = l 
            
        return self
    
    def leaf_carbon_competition_status(_init_, g) :
        return
    
    def leaf_nitrogen_competition_status(_init_, g) :
        return

class Senescence:
    print("Checking whether leaves are in senescence state or not...")
    print("Checking nitrogen quantities re-allocation between leaves...")
    print("Cheking leaf yellowing area evolution...")

    def _init_(self, g) :
        """Leaf scale turnover module initial parameters
          
        """
        for vid in g.vertices(scale=1):
            l = g.node(vid)

        self = l 
            
        return self
    
    # This function calculates the nitrogen remobilised quantity according to ponderation:
    def leaf_senescence(_init_, g,
                                N_potential_shoot_offer) :
            """Leaf senescence mecanism 
        
        Parameters
        ----------
        g: MTG
            the above part of the plant
        N_potential_plant_offer : float
            plant potential nitrogen content

        Returns
        -------
        g : MTG
            the above part of the plant
        plant_area : float
            plant total area 
        PAR_total :
            total light intercepted by the plant
        plant_N : float
            plant nitrogen content
        pool_N : float
            senescence leaves nitrogen pool
        plant_N_by_lightunit : float
            potential nitrogen quantity by light unit
        N_plant_demand : float
            plant nitrogen demand
        N_plant_deficit : float
            plant nitrogen deficit
        N_plant_hunger :
            N plant stress state
        N_plant_remobilized : float
            plant remobilized nitrogen

        State Variables (on MTG)
        ------------------------
        growth : bool
            growth state
        leaf_fall : bool
            viability state of the leaf
        senescence : bool
            leaf development status
        photosynthetic_area : float
            photosynthetic leaf surface
        N_potential_offer : float
            potential leaf nitrogen content
        N_deficit : float
            leaf nitrogen deficit
        weight : float
            value of the leaf participation to the pool accordint to his PAR
        N_exported : float
            leaf exported nitrogen from the source/senescent leaves to the N pool
        N : float
            leaf nitrogen content
        C : float
            leaf carbon content
        N_senescnence_threshold : float
            leaf senescence initiation nitrogen threshold according to the leaf 
            state (MTG)
        C_senescnence_threshold : float
            leaf senescence initiation carbon threshold according to the leaf 
            state (MTG)
        PAR : float
            leaf light intercepted
        
        Hypothesis
        ----------
        > Senescence initiation of leaf choose according to their PAR 
            (lowest PAR = sacrificial leaf) // Samurai law

        > Senescent leaf participation to the available N pool according to their PAR 
            (higher PAR = greater weight)

        References
        ----------
        Skulachev (2002)
        Lemaire et al. (2008)
        Thomas (2013)

        Examples
        --------
        >>> leaf_senescence()

        """
            max_scale = g.max_scale()
            self = l
        
            leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) 
                    if g.class_name(lid) =='Leaf']
            alive = [l for l in leaves if not self.leaf_fall]
            source = [l for l in alive if self.growth]
            senescence = [l for l in alive if self.senescence]

            plant_area = sum(self.leaf_area for l in alive)

            for l in alive :

                plant_N = sum(self.N for l in alive)
                pool_N = sum(self.N for l in source if self.senescence)
                PAR_total = N_potential_shoot_offer
                l.N_demand = self.N_potential_offer
                l.N_deficit = self.N - self.N_demand
                N_plant_demand = sum(l.N_deficit for l in alive)
                N_plant_deficit = pool_N - N_plant_demand

                if N_plant_deficit > 0 :
                    N_plant_hunger = False
                    N_plant_remobilized = N_plant_demand

                    for l in senescence :
                        l.weight = self.PAR / PAR_total
                        l.N_exported = pool_N * l.weight
                        l.N = self.N - l.N_exported
                        l.C = self.C
                    
                    for l in alive :
                        l.N = self.N_potential_offer
                        l.C = self.C

                else :
                    N_plant_hunger = True

                    alive.sort(key=lambda l : l.PAR)

                    for l in alive[0]:
                        l.senescence = True
                        l.N_senescence_threshold = self.N
                        l.C_senescence_threshold = self.C

                        for l in senescence :
                            l.weight = self.PAR / PAR_total
                            l.N_exported = pool_N * l.weight
                            l.N = l.N_senescence_threshold - l.N_exported
                            l.C = l.C_senescence_threshold

                        N_plant_remobilized = sum(l.N_exported for l in senescence)

                        senescence.sort(key=lambda l : l.PAR)

                        for l in senescence[0] : 
                            l.leaf_fall = True
                            l.N_leaf_fall = l.N - l.N_exported
                            l.C_leaf_fall = l.C
                            l.N = l.N_leaf_fall
                            l.C = l.C_leaf_fall

                        survivors = [l for l in alive if not l.leaf_fall]

                        for l in survivors :
                            l.weight = l.PAR / PAR_total
                            l.N_imported = pool_N * l.weight 
                            l.N_deficit = self.N + l.N_imported
                    

                return g, plant_area, plant_N, pool_N, N_plant_demand, 
            N_plant_deficit, N_plant_hunger, N_plant_remobilized


    # This function calculates the redistribution of N among leaves according to ponderation:
    def N_exchange_within_shoots(_init_, g):
        return
    
    # This function calculates the final leaf nitrogen content:
    def leaf_nitrogen_content(_init_, g):
        return
    
    # This function calculates the final radiation use efficiency leaf coefficient:
    def leaf_RUE(_init_, g):
        return
    
    # This function calculates the carbon cost of the redistribution of N among leaves:
    def leaf_senescence_carbon_cost(_init_, g):
        return

class Fall:
    print("Checking whether leaves should fall or not...")
    print("Checking nitrogen and carbon quantity thresholds of leaves falling...")

    def _init_(self, g) :

        for vid in g.vertices(scale=1):
            l = g.node(vid)

        self = l 
            
        return self

    # This function considers whether a leaf should fall or not:
    def leaf_falling(g):
        # for vid in g.vertices_iter(scale=1):
        #     print("The vid is", vid)
        #     n = g.node(vid)
        #      if n.N_content < N_treshold:
        #          n.fallen = True

        return
    
    # This function considers the final carbon and nitrogen mecanistical threshold for leaf fall:
    def leaf_falling_mecanism_threshold(g):
        return
    
    # This function considers the final nitrogen and light empirical threshold for leaf fall:
    def leaf_falling_empirical_threshold(g):
        return


# MAIN SIMULATION ##############################################################################################################################
# We read the environemental/initial conditions:
input_dataframe = pd.read_excel(r'C:\Users\FARRA\AppData\Local\GitHubDesktop\app-3.2.3\src\caracol\src\caracol\APTOS_input_1.xlsx', header=0, sheet_name="scenario")

# OPTION 1: We create a MTG from its known properties (NOT COMPLETELY FINISHED!):
#g = create_MTG_from_csv_file()

# OPTION 2: We create an initial MTG corresponding to the shoot (stem and leaves):
g = vegetative_rapeseed().above

# We initialize the content of N in each leaf/internode:
for vid in g.vertices_iter(scale=1):
    n = g.node(vid)
    # n.N_content = input_dataframe['Initial_N_content'].loc[vid-1]
    n.N_content =input_dataframe['Initial_N_content'].loc[input_dataframe['element_ID'] == vid]

# We initialize the PAR in each leaf/internode:
for vid in g.vertices_iter(scale=1):
    n = g.node(vid)
    n.PAR = input_dataframe['PAR'].loc[input_dataframe['element_ID'] == vid]

# We initialize a new property of the MTG g that corresponds to the status "Fallen / Not fallen":
for vid in g.vertices_iter(scale=1):
    n = g.node(vid)
    n.fallen = False

# We record the properties in its initial state:
recording_MTG_properties(g, 'MTG_0000.csv')

# We define an initial and final step for covering a certain period:
initial_step = 1
final_step = 3

print("Simulation has started...")

# For each time step covering the simulation period:
for step in range(initial_step, final_step+1):

    print("")
    print("For step", step,":")

    # We successively call the different modules:
    starting_senescence(g)
    N_exchange_within_shoots(g)
    leaves_falling(g, N_treshold = 0.02)

    # We record the main properties of the plant:
    properties_filename = 'MTG_%.5d.csv' % step
    recording_MTG_properties(g, properties_filename)

print("")
print("Simulation is done!")