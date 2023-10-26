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
class Files:
    """Files management methods

    Authors
    -------
    Solen FARRA

    """
    def _init_(self, g) :

        for vid in g.vertices(scale=g.max_scale):
            l = g.node(vid)

        self = l 
            
        return self

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

class Rapeseed:
    """Define a MTG with rapeseed properties.

    Authors
    -------
    Solen FARRA

    """
    def _init_(g) :
        
        g = MTG()
            
        return g
    
    def vegetative_rapeseed(_init_, g, n_nodes: int = 6, length_node: float= 2.0):
        """ Build a rapeseed plant with n_nodes internodes and other properties

        Authors
        -------
        Arman GRUMEL and Laurène PERTHAME
            Colzette model

        """
        above = g 

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

    def update_MTG(_init_, g, outputs: dict):
        """Update the MTG 

        Parameters
        ----------
        outputs: dict
            {'param1': { vid1: , vid2, ...}, 'param2': { vid1: , vid2, ...}}
        
        Authors
        -------
        Arman GRUMEL and Laurène PERTHAME
            Colzette model

        """
        # add the missing property
        for param in outputs.keys():
            if param not in g.properties():
                g.add_property(param)

            # update the MTG
            g.property(param).update(outputs[param])  

        return
    
    def generate_normalized_bell_curve_points(_init_, mean, std_dev, amplitude, N):
        """ returns the leaf area of individual leaves along a bell shaped model, as a proportion of the total leaf area 

        Authors
        -------
        Arman GRUMEL and Laurène PERTHAME
            Colzette model

        """
        x_values = np.linspace(mean - 3 * std_dev, mean + 3 * std_dev, N)
        y_values = amplitude * np.exp(-0.5 * ((x_values - mean) / std_dev)**2)
        
        # Normalize y values to ensure the sum is 1
        normalized_y_values = y_values / np.sum(y_values)
        
        points = list(normalized_y_values)

        return points

    def phenotype(_init_, self, g,
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

        Authors
        -------
        Arman GRUMEL and Laurène PERTHAME
            Colzette model
        """
        outputs=dict()

        plants=[k for k, v in g.properties()['label'].items()  if v=='Plant']
        leaves=[k for k, v in g.properties()['label'].items()  if v=='Leaf']
        internodes=[k for k, v in g.properties()['label'].items()  if v=='Internode']

        surface = self.generate_normalized_bell_curve_points(mean=leaf_max,std_dev=sd,amplitude=amp,N=int(len(leaves)/len(plants)))
    
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
        
        self.update_MTG(g, outputs) 

    
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# This function enables to plot an above MTG with PlantGL:

def rapeseed_visitor(g, v, turtle):
    """OpenAlea visitor
    This function is called by the scene3D() function when creating the 3D scene
    
    g is an MTG tree. In contains all the geometrical information needed by the visitor 
    v is the ID of the vertex being generated (int)
    turtle is the turtle used to create the 3D space

    Authors
    -------
    Arman GRUMEL and Laurène PERTHAME
        Colzette model

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

    Authors
    -------
    Arman GRUMEL and Laurène PERTHAME
        Colzette model

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
g = Rapeseed.vegetative_rapeseed().above
Rapeseed.phenotype(g)
#print(display(g))
#print(g.properties)
pgl.Viewer.display(scene3d(g))


# DEFINING METABOLIC FUNCTIONS ##############################################################################################################
class Photosynthesis:
    """Photosynthesis methods.

    Authors
    -------
    Solen FARRA

    """
    print("Checking the carbon produced of leaves according to incident light...")

    def _init_(self, g) :
        
        for vid in g.vertices(scale=1):
            l = g.node(vid)

        self = l 
            
        return self
    
    # This function calculates the carbon produced by a leaf:
    def leaf_photosynthesis(_init_, g) :
        return
    
class Growth:
    """Growth methods. 

    Authors
    -------
    Solen FARRA

    """
    print("Checking the potentiality, demand and cost in nutrient of growth processus...")

    def _init_(self, g) :

        for vid in g.vertices(scale=1):
            l = g.node(vid)

        self = l 
            
        return self
    
    # This function calculates the potential growth of the leaf surface according to RUE:
    def leaf_potential_area(_init_, g, self) :
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
        leaves = [g.node(lid) for lid in g.vertices(scale=1) 
                  if g.class_name(lid) =='Leaf']
        alive = [l for l in leaves if not l.leaf_fall]

        for l in alive :
            l.potential_area = self.RUE * self.PAR

        return g
    
    # This function calculates the leaf nitrogen demand according to the potential leaf area:
    def leaf_nitrogen_demand(_init_, g, self) :
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
        leaves = [g.node(lid) for lid in g.vertices(scale=1) 
                  if g.class_name(lid) =='Leaf']
        alive = [l for l in leaves if not l.leaf_fall]

        for l in alive :
            l.N_demand = self.potential_area - self.N_content
        
        return g
    
    # This function calculates the potential net photosynthesis according to the potential leaf area:
    def leaf_potential_photosynthesis(_init_, g, self) :
        """Calculation of a potential (optimum) quantity of carbon produced of each organ.
    
        Parameters
        ---------
        g: MTG 
            MTG where properties are recorded
        
        photosynthesis_slope: MTG proporties
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
        leaves = [g.node(lid) for lid in g.vertices(scale=1) 
                  if g.class_name(lid) =='Leaf']
        alive = [l for l in leaves if not l.leaf_fall]

        for l in alive :
            l.potential_photosynthesis = (self.photosynthesis_slope * self.PAR) - self.respiration

        return g
    
    # This function calculates the carbon finally available for growth after nutrient redistribution:
    def leaf_carbon_growth(_init_, g, self) :
        """Calculation of a potential (optimum) quantity of carbon produced of each organ.
    
        Parameters
        ---------
        g: MTG 
            MTG where properties are recorded
    
        photosynthesis: MTG proporties
            Actual quantity of carbon produced by photosynthesis (g)
        
        carbon_spend: MTG proporties
            Carbon quantity spend for nitrogen redistribution among organs (g)
        
        Returns
        ------
        g: MTG 
            MTG where properties are recorded

        Authors
        -------
        Solen FARRA

        """
        leaves = [g.node(lid) for lid in g.vertices(scale=1) 
                  if g.class_name(lid) =='Leaf']
        alive = [l for l in leaves if not l.leaf_fall]

        for l in alive :
            l.carbon_growth = self.photosynthesis - self.carbon_spend
        
        return g
    
    # This function calculates the carbon demand according to the potential carbon produced associated with the potential leaf area:
    def leaf_carbon_demand(_init_, g, self) :
        """Calculation of a potential (optimum) quantity of carbon produced of each organ.
    
        Parameters
        ---------
        g: MTG 
            MTG where properties are recorded
        
        potential_photosynthesis: MTG proporties
            Potential (optimum) quantity of carbon produced by photosynthesis according to a potential leaf area (g)
        
        carbon_growth: MTG proporties
            Carbon quantity left for growth of each organs (g)
        
        Returns
        ------
        g: MTG 
            MTG where properties are recorded

        Authors
        -------
        Solen FARRA

        """
        leaves = [g.node(lid) for lid in g.vertices(scale=1) 
                  if g.class_name(lid) =='Leaf']
        alive = [l for l in leaves if not l.leaf_fall]
    
        for l in alive :
            l.carbon_demand = self.potential_photosynthesis * self.carbon_growth

        return g
    
    # This function calculates the final leaf area according to the carbon demand:
    def leaf_actual_area(_init_, g, self) :
        """Calculation of a potential (optimum) quantity of carbon produced of each organ.
    
        Parameters
        ---------
        g: MTG 
            MTG where properties are recorded
    
        carbon_demand: MTG properties
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
        leaves = [g.node(lid) for lid in g.vertices(scale=1) 
                  if g.class_name(lid) =='Leaf']
        alive = [l for l in leaves if not l.leaf_fall]

        for l in alive :
            if self.carbon_demand < 0 :
                l.area = self.potential_area
            else :
                l.area = self.area
                print("The plant is in carbon stress, the growth is inhibate.") 
                print("It needs to produced more carbon to get an optimal growth and leaf area")
            
        return
    
class Competition:
    """Nutrient competition status methods. 

    Authors
    -------
    Solen FARRA

    """
    print("Checking whether leaves are in source or sink status in nitrogen or not...")

    def _init_(self, g) :
        for vid in g.vertices(scale=1):
            l = g.node(vid)

        self = l 
            
        return self
    
    def leaf_carbon_deficit(_init_, g, self) :
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
        leaves = [g.node(lid) for lid in g.vertices(scale=1) 
                  if g.class_name(lid) =='Leaf']
        alive = [l for l in leaves if not l.leaf_fall]

        for l in alive :
            l.carbon_deficit = self.photosynthesis - self.carbon_demand

        return g
    
    def leaf_nitrogen_competition_status(_init_, g, self) :
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

        leaves = [g.node(lid) for lid in g.vertices(scale=1) 
                  if g.class_name(lid) =='Leaf']
        alive = [l for l in leaves if not l.leaf_fall]

        for l in alive :
            if self.carbon_deficit < 0 :
                l.source_N = True
                l.sink_N = False
            else :
                l.source_N = False
                l.sink_N = True
            
        return g

class Senescence:
    """Nitrogen exchange methods. 

    Authors
    -------
    Solen FARRA

    """
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
    
    # This function calculates the nitrogen remobilised quantity according to ponderation,
    # the redistribution of N among leaves according to ponderation,
    # and the final leaf nitrogen content:
    def N_exchange_within_shoots(_init_, g, self, PAR) :
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
        leaves = [g.node(lid) for lid in g.vertices(scale=1) 
                if g.class_name(lid) =='Leaf']
        alive = [l for l in leaves if not self.leaf_fall]
        senescence = [l for l in alive if self.source_N]

        for l in senescence :
            l.weight_export = self.PAR/PAR
            l.N_exported = self.N_content * l.weight_export
            l.N_content = self.N_content - l.N_exported

            pool_N = sum(l.N_exported)
        
            for l in alive :
                l.weight_import = self.PAR/PAR

                if self.nitrogen_demand < 0 :
                    l.nitrogen_import = pool_N * l.weight_import

                else :
                    l.N_storage = self.nitrogen_demand
                    l.N_content = self.N_content

        return g

    # This function calculates the final radiation use efficiency leaf coefficient:
    def leaf_RUE(_init_, g, self):
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
        leaves = [g.node(lid) for lid in g.vertices(scale=1) 
                if g.class_name(lid) =='Leaf']
        alive = [l for l in leaves if not self.leaf_fall]
        senescence = [l for l in alive if self.source_N]
        
        for l in alive : 
            l.RUE = self.N_content/self.PAR

        return g
    
    # This function calculates the carbon cost of the redistribution of N among leaves:
    def leaf_senescence_carbon_cost(_init_, g, self, cost_N_exchange):

        leaves = [g.node(lid) for lid in g.vertices(scale=1) 
                if g.class_name(lid) =='Leaf']
        alive = [l for l in leaves if not self.leaf_fall]
        senescence = [l for l in alive if self.source_N]
        
        for l in alive : 
            l.carbon_spend = self.N_import*cost_N_exchange

        return

class Fall:
    """Leaves falling methods. 

    Authors
    -------
    Solen FARRA

    """
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
g = Rapeseed.vegetative_rapeseed().above

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
Files.recording_MTG_properties(g, 'MTG_0000.csv')

# We define an initial and final step for covering a certain period:
initial_step = 1
final_step = 3

print("Simulation has started...")

# For each time step covering the simulation period:
for step in range(initial_step, final_step+1):

    print("")
    print("For step", step,":")

    # We successively call the different methods:
    Growth.leaf_actual_area(g)
    Growth.leaf_carbon_demand(g)
    Growth.leaf_carbon_growth(g)
    Growth.leaf_nitrogen_demand(g)
    Growth.leaf_potential_area(g)
    Growth.leaf_potential_photosynthesis(g)
    Senescence.N_exchange_within_shoots(g)
    Senescence.leaf_RUE(g)
    Senescence.leaf_senescence_carbon_cost(g, cost_N_exchange = None)
    Fall.leaf_falling(g, N_treshold = 0.02)

    # We record the main properties of the plant:
    properties_filename = 'MTG_%.5d.csv' % step
    Files.recording_MTG_properties(g, properties_filename)

print("")
print("Simulation is done!")