import openalea.mtg
from openalea.mtg import *
from openalea.mtg import turtle as turt
from openalea.mtg.traversal import pre_order, post_order
from openalea.mtg.plantframe import color
import openalea.plantgl.all as pgl

from alinea.caribu.CaribuScene import CaribuScene
from alinea.caribu.data_samples import data_path

import pyvista as pv

import numpy as np
import pandas as pd
from math import sqrt, pi, floor, trunc, exp, isnan, cos, acos
import os

import simbal_parameters as param

########################################################################################################################
########################################################################################################################

# DEFINING GENERAL FUNCTIONS:
#############################

#-----------------------------------------------------------------------------------------------------------------------
# General function for easily adding a new child element onto a parent:
def ADDING_A_CHILD_ORGAN(parent_element,
                         edge_type = '<',
                         label = 'Internode',
                         properties_identical_to_parent = False,
                         nil_properties = False):
    """
    This function enables to add a child element on the parent element, with properties equal to those of the parent,
    or to 0.
    :param parent_element: the element of the MTG onto which we add a new child element
    :param edge_type: the type of connection between the child and the parent ('<' or '+')
    :param label: the label of the child
    :param properties_identical_to_parent: if True, most properties of the child will be identical to those of the parent
    :param nil_properties: if True, most properties of the child will be 0
    :return: the new child element
    """

    # We add a new element on the parent:
    child = parent_element.add_child(edge_type = edge_type,
                                     label = label,
                                     radius = parent_element.radius)

    # We define the list of properties of the child that need to be updated, based on the properties of the parent element:
    list_of_attributes_to_set = list(parent_element.properties().keys())
    # We remove unwanted properties that have already been defined:
    list_of_attributes_to_set.remove("label")
    list_of_attributes_to_set.remove("radius")
    if "edge_type" in list_of_attributes_to_set:
        list_of_attributes_to_set.remove("edge_type")

    # We cover each attribute of the parent to be updated:
    for attribute in list_of_attributes_to_set:
        # If the properties of the child are expected to be identical to those of the parent:
        if properties_identical_to_parent:
            # Then we assign the values of the parent:
            setattr(child, attribute, parent_element.properties()[attribute])
        # Otherwise, the properties of the child are set to 0:
        elif nil_properties:
            setattr(child, attribute, 0)
        else:
            setattr(child, attribute, np.NaN)

    return child
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# General function for easily adding a new phytomer onto a parent element:
def ADDING_A_PHYTOMER(parent_element,
                      struct_mass_density=1e5, minimal_leaf_thickness=1e-3,
                      internode_radius=0.,
                      internode_length=0.,
                      leaf_radius=0.,
                      leaf_length=0.,
                      bud_radius=0.,
                      bud_length=0.,
                      properties_identical_to_parent = True,
                      nil_properties = False
                      ):

    # 1) We first create a new internode:
    # ------------------------------------
    internode = ADDING_A_CHILD_ORGAN(parent_element,
                                     edge_type='<',
                                     label='Internode',
                                     properties_identical_to_parent=properties_identical_to_parent,
                                     nil_properties=nil_properties)
    # We specified some properties that may be distinct from the parent element:
    internode.radius = internode_radius
    internode.length = internode_length
    internode.angle_down = 0.
    internode.angle_roll = param.phyllotactic_angle
    # We calculate the surface, volume and structural mass of the new organ:
    internode.surface = internode.radius * 2. * pi * internode.length
    internode.volume = internode.radius ** 2. * pi * internode.length
    internode.struct_mass = internode.volume * struct_mass_density
    internode.photosynthetic_area = internode.surface

    # 2) We then create a new leaf that is actually a "branch" on the previous internode:
    # ------------------------------------------------------------------------------------
    leaf = ADDING_A_CHILD_ORGAN(parent_element=internode,
                                edge_type='+',
                                label='Leaf',
                                properties_identical_to_parent=properties_identical_to_parent,
                                nil_properties=nil_properties)
    # We specified some properties that may be distinct from the parent element:
    leaf.radius = leaf_radius
    leaf.length = leaf_length
    leaf.angle_down = param.leaf_insertion_angle
    leaf.angle_roll = 0.
    # We calculate the surface, assuming the leaf is an ellipse characterized by its length and width:
    leaf.surface = leaf.length / 2. * leaf.radius * pi
    # We calculate the volume, taking into account a minimal thickness:
    leaf.volume = leaf.surface * minimal_leaf_thickness
    leaf.struct_mass = leaf.volume * struct_mass_density
    leaf.photosynthetic_area = leaf.surface

    # 3) We finally create an axilliary bud that is actually a "branch" on the previous internode:
    # ---------------------------------------------------------------------------------------------
    bud = ADDING_A_CHILD_ORGAN(parent_element=internode,
                               edge_type='+',
                               label='Bud',
                               properties_identical_to_parent=False,
                               nil_properties=True)
    # We specified some properties that may be distinct from the parent element:
    bud.radius = bud_radius
    bud.length = bud_length
    bud.angle_down = param.bud_insertion_angle
    bud.angle_roll = 0.
    # We assign the surfaces and volumes:
    bud.surface = 0.
    bud.volume = 0.
    bud.struct_mass = 0.
    bud.photosynthetic_area = 0.

    # Eventually, we redefine the starting coordinates of the internode, leaf and bud as the last point of the parent:
    for element in [internode, leaf, bud]:
        element.x1 = parent_element.x2
        element.x1 = parent_element.x2
        element.y1 = parent_element.y2
        element.z1 = parent_element.z2

    return [internode, leaf, bud]
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function for starting the formation of new phytomers on the shoot MTG:
def new_phytomers_initiation(g, new_initiation = True, initial_internode_length = 0.0005, initial_leaf_length = 0.003):
    """
    This function defines how new phytomers are formed over time.
    :param g: the MTG corresponding to the shoot system
    :return:
    """

    # We initialize the current list of apices in the shoot system that may produce new phytomers:
    list_of_apices = [g.node(vid) for vid in g.vertices_iter(scale=1) if g.node(vid).label == "Apex"]
    # We initialize a condition defining whether new phytomers should develop from a given apex:
    apparition_condition=new_initiation # TODO: Define how and when new phytomers are formed!

    # As long as the list of apices to treat is not empty:
    while list_of_apices != []:
        # We cover each apex present in the list:
        for apex in list_of_apices:
            # We remove the current apex from the list, so that it is not considered next time:
            list_of_apices.remove(apex)
            # If the conditions for developing a new phytomer are OK:
            if apparition_condition is True:
                # Then we add a new phytomer at the top of the current apex:
                [internode, leaf, bud] = ADDING_A_PHYTOMER(apex,
                                                           internode_radius=apex.radius,
                                                           internode_length=initial_internode_length,
                                                           leaf_radius=initial_leaf_length * param.leaf_width_to_length_ratio / 2.,
                                                           leaf_length=initial_leaf_length,
                                                           bud_radius=0.,
                                                           bud_length=0.,
                                                           properties_identical_to_parent=False,
                                                           nil_properties=True
                                                           )
                # print("We have added a new phytomer on the apex", apex.index(), "!")
                # We remove the previous apex:
                g.remove_vertex(apex.index(), reparent_child=True)
                # And we add a new apex at the top of the new internode, of length 0:
                new_apex = ADDING_A_CHILD_ORGAN(parent_element=internode,
                                                edge_type="<",
                                                label="Apex",
                                                properties_identical_to_parent = True,
                                                nil_properties=False)
                # print("We have moved the previous apex", apex.index(), "on the node", internode.index(),"; its index is now", new_apex.index())
                new_apex.length = 0.
                # TODO: WATCH OUT BELOW, IN CASE YOU WANT TO CONTINUE ADDING PHYTOMERS
                # # Finally, we update the list of apices by adding the new apex that has been formed,
                # # in case it could itself form a new phytomer:
                # list_of_apices.append(new_apex)

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# This function enables to generate an initial MTG with a certain number of leaves and corresponding internodes:
def create_shoot_MTG(number_of_phytomers=3,
                     initial_leaf_length=0.003,
                     initial_internode_radius=0.001,
                     initial_internode_height = 0.0005,
                     struct_mass_density=1e5, minimal_leaf_thickness=1e-3):

    # We initialize an empty MTG:
    g = MTG()

    # We define the base of the shoot system, ADDING ALL THE POSSIBLE PROPERTIES THAT A SHOOT ORGAN WILL LATER CONTAIN:
    id_base = g.add_component(g.root,
                              label='Base_of_the_shoot_system',
                              # Geometry:
                              radius=initial_internode_radius,
                              length=0.,
                              angle_down=0.,
                              angle_roll=0.,
                              x1=0.,
                              x2=0.,
                              y1=0.,
                              y2=0.,
                              z1=0.,
                              z2=0.,
                              # C & N labile pools:
                              C_labile=0.,
                              C_labile_deficit=0.,
                              N_labile=0.,
                              N_labile_deficit=0.,
                              # Rates:
                              photosynthesized_C_rate = 0.,
                              C_input_rate = 0.,
                              shoot_maintenance_respiration_rate = 0.,
                              actual_C_consumption_rate_by_shoot_growth = 0.,
                              N_input_rate = 0.,
                              actual_N_consumption_rate_by_shoot_growth = 0.,
                              shoot_N_mobilization_rate = 0.,
                              net_shoot_N_immobilization_rate = 0.,
                              # Other attributes:
                              senescence_status = False,
                              PAR_intercepted=0.,
                              Caribu_Eabs = 0.,
                              Caribu_Ei = 0.,
                              Caribu_area = 0.
                              )

    # We access this first element, and initialize the internode on which the next elements will be attached:
    new_internode = g.node(id_base)
    # We calculate the surface, volume and structural mass of this first element:
    new_internode.surface = new_internode.radius * 2. * pi * new_internode.length
    new_internode.volume = new_internode.radius ** 2. * pi * new_internode.length
    new_internode.struct_mass = new_internode.volume * struct_mass_density
    new_internode.photosynthetic_area = new_internode.surface

    # For each new phytomer:
    for i in range(0, number_of_phytomers):
        # We first create an internode, and then a connected leaf and a bud:
        [new_internode, new_leaf, new_bud] = ADDING_A_PHYTOMER(parent_element=new_internode,
                                                               struct_mass_density=struct_mass_density,
                                                               minimal_leaf_thickness=minimal_leaf_thickness,
                                                               internode_radius=new_internode.radius * 0.9,
                                                               internode_length=initial_internode_height,
                                                               leaf_radius=initial_leaf_length * param.leaf_width_to_length_ratio,
                                                               leaf_length=initial_leaf_length,
                                                               bud_radius=0.,
                                                               bud_length=0.)

    # Finally, on the last internode, we add a new apex that may continue to form new phytomers over the simulation:
    new_apex = ADDING_A_CHILD_ORGAN(parent_element=new_internode,
                                    edge_type="<",
                                    label="Apex",
                                    properties_identical_to_parent=True)
    new_apex.length = 0.

    return g
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function for loading a MTG from CSV file containing its properties:
def create_MTG_from_csv_file(csv_filename='MTG_00003.csv'):

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
#-----------------------------------------------------------------------------------------------------------------------
# Function for moving the turtle and create specific shapes for stems and leaves:

# Function for moving the turtle and create specific shapes for stems and leaves:
def get_shoot_visitor(leaf_shape="oval", upscaling_factor=1):
    """
    This function describes the movement of the 'turtle' along the shoot MTG for creating a scene for Caribu or PlantGL.
    :leaf_shape: if "oval", the leaves will be displayed as an ellipse based on their length and width; else, the leaves
    are displayed with a typical leaf curved shape, based on their length only.
    :return: shoot_visitor
    """

    def oval_shape(length=1., width=1.):
        """
        This simple function generates a PlantGL shape of an ellipse, defined by a length and a width.
        """
        # We create first a round disc:
        disc = pgl.Disc()
        # We specify the number of triangles to build this shape (8 by default):
        disc.slices = 12
        # And we now 'deform' the disc so that it becomes an ellipse characterized by a length and a width:
        shape = pgl.Scaled((length, width, 1), disc)
        return shape

    def unileaflet(length=1., width=1.):
        """Generates a unileaflet shape [taken from Colzette].

        :param length: length of the shape, defaults to 1.
        :type length: float, optional
        :param width: width of the shape, defaults to 1.
        :type width: float, optional
        :return: the shape
        :rtype: pgl.Group
        """
        disc = pgl.Translated((-0.5, 0, 0), pgl.Disc())
        disc = pgl.Scaled((length, width, 1), disc)
        disc = pgl.AxisRotated(axis=(0, 1, 0), angle=pi/2., geometry=disc)

        d3 = pgl.AxisRotated(axis=(1, 0, 0), angle=0., geometry=disc)

        shape = pgl.Group([d3])
        return shape

    def shoot_visitor(g, v, turtle):

        # We access the node corresponding to the vertex ID "v" of the MTG "g":
        n = g.node(v)

        # We look at the geometrical properties already defined within the root element:
        radius = n.radius * upscaling_factor
        length = n.length * upscaling_factor
        angle_down = n.angle_down
        angle_roll = n.angle_roll

        # We get the x,y,z coordinates from the beginning of the root segment, before the turtle moves:
        position1 = turtle.getPosition()
        n.x1 = position1[0] * upscaling_factor
        n.y1 = position1[1] * upscaling_factor
        n.z1 = position1[2] * upscaling_factor
        # The turtle is moved:
        turtle.setId(v)
        # The direction of the turtle is changed:
        turtle.down(angle_down)
        turtle.rollL(angle_roll)

        # 1) CASE OF A LEAF:
        #-------------------
        # In the current element corresponds to a leaf:
        if n.label == 'Leaf':
            # A. OPTION WITH OVAL SHAPE:
            if leaf_shape == "oval":
                # First, we need to move the turtle (without showing it) so that it is located in the middle of the oval
                # to be displayed:
                turtle.f(length/2.)
                # We also need to rotate the turtle so that the oval is correctly oriented:
                turtle.down(90)
                # Finally, we print the oval shape according to the length and radius of the leaf:
                turtle.customGeometry(oval_shape(length = length, width = 2 * radius))
            # B. OPTION WITH 'TYPICAL' LEAF SHAPE:
            elif leaf_shape == "typical_leaf_shape":
                scale_factor = 1. * length
                leafshape = pgl.NurbsPatch([[(-0.0241355*scale_factor, -0.0549754*scale_factor, -0.00804807*scale_factor, 1),
                                         (-0.00877065*scale_factor, -0.00933309*scale_factor, 0.00100194*scale_factor, 1),
                                         (0.0221632*scale_factor, 0.0405477*scale_factor, -0.00898557*scale_factor, 1),
                                         (-0.00223729*scale_factor, -0.00247516*scale_factor, 0.000385988*scale_factor, 1)],
                                        [(0.314084*scale_factor, -0.472921*scale_factor, 0.392759*scale_factor, 1),
                                         (0.000957513*scale_factor, 0.000400347*scale_factor, 0.33376*scale_factor, 1),
                                         (-0.00092555*scale_factor, 0.00517653*scale_factor, 0.332921*scale_factor, 1),
                                         (0.408856*scale_factor, 0.483731*scale_factor, 0.413612*scale_factor, 1)],
                                        [(0.329337*scale_factor, -0.501345*scale_factor, 0.730665*scale_factor, 1),
                                         (0.000326354*scale_factor, 0.0164005*scale_factor, 0.665257*scale_factor, 1),
                                         (-0.000118909*scale_factor, 0.019707*scale_factor, 0.66817*scale_factor, 1),
                                         (0.414001*scale_factor, 0.558038*scale_factor, 0.74373*scale_factor, 1)],
                                        [(0.00223079*scale_factor, 0.0223592*scale_factor, 0.999352*scale_factor, 1),
                                         (0.00112035*scale_factor, 0.028812*scale_factor, 1.0005*scale_factor, 1),
                                         (-0.000767699*scale_factor, 0.0327184*scale_factor, 0.999658*scale_factor, 1),
                                         (-0.0011293*scale_factor, 0.0361733*scale_factor, 1.00279*scale_factor, 1)]])
                turtle.customGeometry(leafshape)
                # We move the turtle by the length of the leaf, without modifying its shape:
                turtle.f(length)
            elif leaf_shape == "unileaflet":
                turtle.customGeometry(unileaflet(length=length, width=2 * radius))
            # C. OPTION WITH CLASSICAL CYLINDERS:
            else:
                # We define the radius of the cylinder to be displayed:
                turtle.setWidth(radius)
                # We move the turtle by the length of the root segment:
                turtle.F(length)

        # 2) CASE OF ANOTHER SHOOT ORGAN:
        # -------------------------------
        # Otherwise, if the current element is not a leaf, it is defined as a cylinder:
        else:
            # We define the radius of the cylinder to be displayed:
            turtle.setWidth(radius)
            # We move the turtle by the length of the root segment:
            turtle.F(length)

        # We get the x,y,z coordinates from the end of the root segment, after the turtle has moved:
        position2 = turtle.getPosition()
        n.x2 = position2[0] * upscaling_factor
        n.y2 = position2[1] * upscaling_factor
        n.z2 = position2[2] * upscaling_factor

    return shoot_visitor

def create_a_scene_with_plants_popuplation(g,
                                           number_of_plants_per_square_meter = 9,
                                           inter_row_distance = 0.45,
                                           scene_length=1.,
                                           scene_width=1.,
                                           true_color=True,
                                           surrounding_plants_transparency=0.,
                                           leaf_shape="oval",
                                           upscaling_factor=1):
    """
    This function generates a scene containing a population of plants, all cloned from the same plant MTG "g",
    and distributed within rows with a random orientation. The plant in the middle of the scene is the last set of shapes
    given by the scene. The surrounding plants can be displayed as (semi-) transparent or not.

    Parameters
    ----------
    g: the MTG corresponding to the shoot system
    number_of_plants_per_square_meter: plant density (plants per m2)
    inter_row_distance: distance between two rows (m)
    scene_length: length of the rows (m)
    scene_width: width of the scene (m)
    true_color: if so, ...
    surrounding_plants_transparency: level of transparency of plants adjacent to the central one in the scene (0 = non transparency, 1 = full trnasparency)
    leaf_shape: if "oval", the leaves will be displayed as an ellipse based on their length and width; else, the leaves
    are displayed with a typical leaf curved shape, based on their length only.
    upscaling_factor: a factor by which all dimensions should be multiplied (e.g. useful to convert dimensions from meters to centimeters)

    Returns
    -------
    final_scene: the final scene with the population of plants

    """

    # 1.COMPUTING PLANT POSITIONS:
    ##############################

    # We calculate the number of rows and plants to be show in the scene:
    number_of_rows = floor(scene_width/float(inter_row_distance)) + 1
    number_of_plants = int(number_of_plants_per_square_meter * (scene_length*scene_width))
    # We make sure that there is a least one row containing plants:
    if number_of_plants < number_of_rows:
        number_of_rows = number_of_plants

    # We calculate the number of plants per row:
    number_of_plants_per_row = floor(number_of_plants / float(number_of_rows))
    # We calculate the distance between two plants within a row (in meters):
    inter_plant_distance = float(scene_length) / float(number_of_plants_per_row)
    # We finally calculate the length and width of the lateral borders, so that the population of plants is centered
    # within the square scene:
    lateral_border_width = (scene_width - (number_of_rows - 1) * inter_row_distance) / 2.
    lateral_border_length = (scene_length - (number_of_plants_per_row - 1) * inter_plant_distance) / 2.

    # We initialize the list that will contain the x,y,z coordinates of each plant in the scene:
    list_of_positions = []
    # We move to the position of the first plant:
    x = lateral_border_length * upscaling_factor
    y = lateral_border_width * upscaling_factor
    z = 0. * upscaling_factor
    # We cover each plant of the scene:
    for plant_index in range(0, number_of_plants):
        # We add the new position (x,y,z) to the list:
        list_of_positions.append((x,y,z))
        # And we move from one plant to another on the same row:
        x += inter_plant_distance * upscaling_factor
        # If we have reached the maximal length of a given row:
        if x > scene_length * upscaling_factor:
            # Then we switch to the next row:
            y +=  inter_row_distance * upscaling_factor
            # And we start at the beginning of the new row:
            x = lateral_border_length * upscaling_factor
    # We define the index of the position of the plant that stays in the middle of the scene :
    index_middle_plant = floor(number_of_plants/ 2.)
    # And we move the position of the middle plant to the last position of the list, so that it is easily recognizable:
    list_of_positions.append(list_of_positions[index_middle_plant])
    list_of_positions.remove(list_of_positions[index_middle_plant])

    # CREATING THE FINAL SCENE:
    ###########################

    # We initialize a turtle in PlantGL:
    turtle = turt.PglTurtle()
    # We initialize the final scene:
    final_scene = pgl.Scene()

    # We create a "soil surface" corresponding to the limits of the scene:
    points = [(0, 0, 0),
              (scene_length * upscaling_factor, 0, 0),
              (scene_length * upscaling_factor, scene_width * upscaling_factor, 0),
              (0, scene_width * upscaling_factor, 0)]
    normals = [(0, 0, 1) for i in range(4)]
    indices = [(0, 1, 2, 3)]
    square = pgl.QuadSet(points, indices, normals, indices)
    # We define the same color for each of the 4 points:
    square.colorList = [pgl.Color4(100,100,0,200) for i in range(4)]
    final_scene += square

    # If the option "true_color" has been required:
    if true_color:
        # We try to access the dictionnary containing the colors of each vertex:
        colors = g.property('color')
        # However, if the dictionnary is empty, we reset the option to False:
        if colors == {}:
            true_color = False

    # We cover each plant in the population:
    for plant_index in range(0,len(list_of_positions)):

        # We first define the seed of random, depending on the index of the plant:
        np.random.seed(int(param.random_choice * plant_index))
        # We then modify the angle of the plant, so that leaves are not oriented the same way within the population:
        angle_roll = abs(np.random.normal(180, 180))

        # We reset the turtle:
        turtle.reset()
        # And we reposition the turtle from there:
        turtle.move(list_of_positions[plant_index])
        turtle.rollL(angle_roll)
        # We create a scene for one plant by moving the turtle along the MTG:
        scene = turt.TurtleFrame(g, visitor=get_shoot_visitor(leaf_shape=leaf_shape, upscaling_factor=upscaling_factor), turtle=turtle, gc=False)

        # We get a dictionnary of all shapes in the scene:
        shapes = dict((sh.id, sh) for sh in scene)

        # We cover each shape of the current MTG:
        for shape in scene:
            # If the true color has to be activated:
            if true_color:
                # If the plant corresponds to the central plant (which we have moved to the last position in the list):
                if plant_index == len(list_of_positions) - 1:
                    # Then we add the color of the corresponding vertex, with no transparency:
                    shape.appearance = pgl.Material(colors[shape.id], transparency=0.0)
                else:
                    # Then we add the color of the corresponding vertex, with semi-transparency:
                    shape.appearance = pgl.Material(colors[shape.id], transparency=surrounding_plants_transparency)
            else:
                # If the plant corresponds to the central plant (which we have moved to the last position in the list):
                if plant_index == len(list_of_positions) - 1:
                    # Then we add the color of the corresponding vertex, with no transparency:
                    shape.appearance = pgl.Material([0, 250, 0], transparency=0.0)
                else:
                    # For all other plants, we set nearly transparent plants with greenish colors:
                    shape.appearance = pgl.Material([50, 200, 50], transparency=surrounding_plants_transparency)
            # Finally, we add the shape to the final scene:
            final_scene += shape

    return final_scene

# Function for creating a color map for a MTG:
def my_colormap(g, property_name, cmap='jet', vmin=None, vmax=None, lognorm=True):
    """
    This function computes a property 'color' on a MTG based on a given MTG's property.
    :param g: the investigated MTG
    :param property_name: the name of the property of the MTG that will be displayed
    :param cmap: the type of color map
    :param vmin: the min value to be displayed
    :param vmax: the max value to be displayed
    :param lognorm: a Boolean describing whether the scale is logarithmic or not
    :return: the MTG with the corresponding color
    """

    # We make sure that the user did not accidently switch between vmin and vmax:
    if vmin >= vmax:
        raise Exception("Sorry, the vmin and vmax values of the color scale of the graph are wrong!")
    if lognorm and (vmin <= 0 or vmax <= 0):
        raise Exception(
            "Sorry, it is not possible to represent negative values in a log scale - check vmin and vmax!")

    prop = g.property(property_name)
    keys = prop.keys()
    values = list(prop.values())
    # m, M = int(values.min()), int(values.max())

    _cmap = color.get_cmap(cmap)
    norm = color.Normalize(vmin, vmax) if not lognorm else color.LogNorm(vmin, vmax)
    values = norm(values)

    colors = (_cmap(values)[:, 0:3]) * 255
    colors = np.array(colors, dtype=int).tolist()

    g.properties()['color'] = dict(zip(keys, colors))

    return g

# Function for setting the center anc camera position in a scene:
def setting_PGLViewer(width=1200, height=1200,
                      x_center=0., y_center=0., z_center=0.,
                      x_cam=0., y_cam=0., z_cam=-0.2,
                      grid=True,
                      background_color=[255, 255, 255]):
    """
    This function sets the center of the graph and the relative position of the camera for displaying a scene in
    PGL.
    :param width: width of the window (pixels)
    :param height: height of the window (pixels)
    :param x_center: x-coordinate of the center at which the camera looks
    :param y_center: y-coordinate of the center at which the camera looks
    :param z_center: z-coordinate of the center at which the camera looks
    :param x_cam: x-coordinate of the camera position
    :param y_cam: y-coordinate of the camera position
    :param z_cam: z-coordinate of the camera position
    :param grid: if True, all grids will be displayed
    :param background_color: a RGB vector corresponding to the background color of the graph
    :return:
    """

    # We define the coordinates of the point cam_target that will be the center of the graph:
    cam_target = pgl.Vector3(x_center, y_center, z_center)
    # We define the coordinates of the point cam_pos that represents the position of the camera:
    cam_pos = pgl.Vector3(x_cam, y_cam, z_cam)
    # We position the camera in the scene relatively to the center of the scene:
    pgl.Viewer.camera.lookAt(cam_pos, cam_target)
    # We define the dimensions of the graph:
    pgl.Viewer.frameGL.setSize(width, height)
    # We define the background of the scene:
    pgl.Viewer.frameGL.setBgColor(background_color[0], background_color[1], background_color[2])
    # We define whether grids are displayed or not:
    pgl.Viewer.grids.set(grid, grid, grid, grid)

    return

#-----------------------------------------------------------------------------------------------------------------------
# This function enables to plot a shoot MTG with PlantGL:
def plotting_shoot_MTG(g,
                       width=1200, height=900,
                       x_center=0, y_center=-0, z_center=0.,
                       x_cam=5., y_cam=5., z_cam=-1,
                       grid=True,
                       background_color=[94, 76, 64],
                       single_color=None,
                       property_name='radius', cmap='jet', vmin=1e-4, vmax=1e-2, lognorm=True,
                       recording_plot=True,
                       plot_name='plot.png'):

    # CREATING THE PLOT:
    #-------------------

    # We first compute the colors of the graph:
    my_colormap(g, property_name=property_name, cmap=cmap, vmin=vmin, vmax=vmax, lognorm=lognorm)
    # We use the property 'color' of the MTG calculated by the function 'my_colormap':
    colors = g.property('color')

    # We initialize a turtle in PlantGL:
    turtle = turt.PglTurtle()

    # We initialize a scene:
    new_scene = pgl.Scene()

    # SETTING THE NEW SCENE WITH THE TURTLE:
    # We create the scene by moving the turtle along the current MTG:
    scene = turt.TurtleFrame(g, visitor=get_shoot_visitor(leaf_shape="oval"), turtle=turtle, gc=False)

    # We get a dictionnary of all shapes in the scene:
    shapes = dict((sh.id, sh) for sh in scene)
    # We cover each node of the MTG:
    for vid in colors:
        if vid in shapes:
            n = g.node(vid)
            if single_color is None:
                shapes[vid].appearance = pgl.Material(colors[vid], transparency=0.0)
            else:
                # Otherwise, we can print it with a unique color:
                shapes[vid].appearance = pgl.Material(single_color, transparency=0.0)

    # We add each shape of the current MTG to the final scene:
    for shape in scene:
        new_scene += shape

    # SETTING THE VIEWER OPTIONS:
    setting_PGLViewer(width=width, height=height,
                      x_center=x_center, y_center=y_center, z_center=z_center,
                      x_cam=x_cam, y_cam=y_cam, z_cam=z_cam,
                      grid=grid,
                      background_color=background_color)

    # DISPLAYING THE MTG:
    pgl.Viewer.display(new_scene)

    # RECORDING THE MTG:
    if recording_plot:
        print("Saving the PlantGL plot...")
        pgl.Viewer.saveSnapshot(plot_name)
        print("Plot saved!")

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
def plotting_shoot_MTG_population(g,
                                  number_of_plants_per_square_meter=33,
                                  inter_row_distance=0.45,
                                  scene_length=1.,
                                  scene_width=1.,
                                  leaf_shape="oval",
                                  width=1200, height=900,
                                  x_center=0, y_center=-0, z_center=0.,
                                  x_cam=5., y_cam=5., z_cam=-1,
                                  grid=True,
                                  background_color=[94, 76, 64],
                                  single_color=None,
                                  property_name='PAR_intercepted', cmap='jet', vmin=0, vmax=400, lognorm=False,
                                  recording_plot=True,
                                  plot_filepath='plot.png'
                                  ):

    my_colormap(g, property_name=property_name, cmap=cmap, vmin=vmin, vmax=vmax, lognorm=lognorm)
    complete_scene = create_a_scene_with_plants_popuplation(g,
                                                            number_of_plants_per_square_meter=number_of_plants_per_square_meter,
                                                            inter_row_distance=inter_row_distance,
                                                            scene_length=scene_length,
                                                            scene_width=scene_width,
                                                            leaf_shape=leaf_shape)
    setting_PGLViewer(width=width, height=height,
                      x_center=x_center, y_center=y_center, z_center=z_center,
                      x_cam=x_cam, y_cam=y_cam, z_cam=z_cam,
                      grid=grid,
                      background_color=background_color)
    # background_color=[150, 200, 250])
    pgl.Viewer.display(complete_scene)
    pgl.Viewer.saveSnapshot(plot_filepath)

    return


# Function that methodically creates each shape (takes a long time for large root systems):
def plotting_shoots_with_pyvista(g,
                                showing=True, recording_image=True, image_file='plot_.png',
                                factor_of_higher_resolution = 1,
                                background_color = [94, 76, 64],
                                plot_width=1000, plot_height=750,
                                camera_x=2., camera_y=0., camera_z=0.5,
                                focal_x=0., focal_y=0., focal_z=0.,
                                closing_window=False):
    """
    This functions aims to plot a shoot system with Pyvista, creating step by step every shape.
    :param g: the root MTG to be displayed
    :param showing: if True, a Pyvista window will be opened
    :param recording_image: if True, an image file of the plot will be recorded
    :param image_file: name of the image file
    :param closing_window: if True, the Pyvista window will be closed after having shown the plot
    :param background_color: a list containing three integers corresponding to Red, Green and Blue values (0-255)
    :param plot_width: width of the plot (in pixels)
    :param plot_height: height of the plot (in pixels)
    :param camera_x: x-coordinate of the camera position
    :param camera_y: y-coordinate of the camera position
    :param camera_z: z-coordinate of the camera position
    :param focal_x: x-coordinate of the focal point (from which plot is centered and allowed to rotate from)
    :param focal_y: y-coordinate of the focal point (from which plot is centered and allowed to rotate from)
    :param focal_z: z-coordinate of the focal point (from which plot is centered and allowed to rotate from)
    :return:
    """

    # We initialize a plot:
    p = pv.Plotter()
    # p = pvqt.BackgroundPlotter() # If we want to update the graph when opened
    if recording_image and not showing:
        p.off_screen=True

    # We initialize a list of blocks and empty lists:
    blocks = pv.MultiBlock()
    colors=[]
    opacity=[]

    # CREATING THE SHAPES:
    ######################

    # The creation of the cylinders with no directions raises a Warning in numpy, but we can ignore it.
    np.seterr(invalid='ignore')

    print("    > Creating the shapes of the plot...")
    # We cover all the vertices in the MTG g:
    for vid in g.vertices_iter(scale=1):

        # We define the current root element as n:
        n = g.node(vid)
        # We skip this element if it does not have a positive length:
        if n.length <= 0.:
            continue

        # We calculate the geometric features required by pyvista to display a cylinder at the right spot:
        x1 = n.x1
        x2 = n.x2
        y1 = n.y1
        y2 = n.y2
        z1 = n.z1
        z2 = n.z2
        center_x = (x1 + x2) / 2.
        center_y = (y1 + y2) / 2.
        center_z = (z1 + z2) / 2.
        direction_x = x2 - x1
        direction_y = y2 - y1
        direction_z = z2 - z1
        # normal_x = direction_x * np.cos(90) - direction_y * np.sin(90)
        # normal_y = direction_x * np.sin(90) + direction_y * np.cos(90)
        # normal_z = direction_z
        height = ((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2) ** 0.5
        radius = n.radius
        # We add the color and opacity to a list:
        color = n.color
        colors.append(color)
        opacity.append(1) # Here, all the elements are supposed to be fully opaques

        if n.label == "Internode":
            # We create a new cylinder with the proper dimensions and position:
            shape = pv.Cylinder(center=(center_x,center_y,center_z), direction=(direction_x, direction_y, direction_z),
                                radius=radius, height=height)
        elif n.label == "Leaf":
            shape = pv.Disc(center=(center_x, center_y, center_z), inner=0., outer=n.length/2., c_res=12,
                            normal=(-direction_x, -direction_y, direction_z))
            # shape = pv.Triangle([[x1,y1,z1], [x2-0.1,y2-0.1,z2], [x2+0.1,y2+0.1,z2]])
            # shape = shape.rotate_vector(vector=(x2, -y2, 0), angle=180, point=(center_x, center_y, center_z))
            # shape = shape.rotate_vector(vector=(direction_x, direction_y, 0), angle=180, point=(center_x, center_y, center_z))
        # And we add this shape to the list of blocks:
        blocks.append(shape)

    # We reset the warning option in Numpy:
    np.seterr(invalid=None)

    # CREATING THE SCENE:
    # Adding the multi-blocks to a plot:
    actor, mapper = p.add_composite(blocks)
    print("    > Adding the colors and transparencies of the plot...")
    # We cover each block and set its color and opacity:
    for i in range(1,blocks.n_blocks+1):
        mapper.block_attr[i].color = colors[i-1]
        mapper.block_attr[i].opacity = opacity[i-1]

    # We set the color of the background:
    p.background_color = background_color

    # ADJUSTING THE CAMERA:
    # We can define the target position and the camera position:
    p.camera.focal_point = (focal_x, focal_y, focal_z)
    p.camera.position = (camera_x, camera_y, camera_z)
    # Or we can define the relative distance and elevation of the camera compared to the middle of the plot:
    # p.camera.distance = 0.2
    # p.camera_elevation = -0.2
    # p.reset_camera()
    # p.view_yz()
    # pv.Axes.origin=(0,0,0.0)
    p.camera_set = True

    # SETTING AXES:
    # p.add_axes()
    # p.add_axes_at_origin()

    # DIVERSE SETTINGS:
    # We set the background color of the plot:
    p.background_color = background_color
    # We set the final size of the window:
    p.window_size = [plot_width, plot_height]
    # We can improve the resolution of the image:
    p.image_scale = factor_of_higher_resolution

    # p.line_smoothing = True
    # p.add_bounding_box(line_width=3, color='black')
    # p.enable_terrain_style(mouse_wheel_zooms=True, shift_pans=True)

    print("> The plot has been created!")

    # SHOWING/RECORDING THE PLOT:
    if recording_image:
        # p.screenshot(filename='test_MTG.png')
        p.show(screenshot=image_file, auto_close=closing_window)
        # # p.close()
    elif showing:
        p.show(auto_close=closing_window)

    if closing_window:
        pv.close_all()

    return p
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# This function enables to record the properties of each element of the MTG:
def recording_MTG_properties(g, file_name='g_properties.csv'):
    """
    This function records the properties of each node of the MTG "g" inside a csv file.
    :param g: the MTG where properties are recorded
    :param file_name: the name of the csv file where properties of each node will be recorded
    :return: [no return]
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

    return data_frame
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
def computing_light_interception(g,
                                 incoming_light_density = 400,
                                 isolated_plant=False,
                                 zenithal_light=False,
                                 infinite=False,
                                 number_of_plants_per_square_meter=33,
                                 inter_row_distance=0.45,
                                 scene_length=0.3,
                                 scene_width=1.,
                                 leaf_shape="oval"
                                 ):
    """

    Parameters
    ----------
    g: the MTG corresponding to the shoot system
    incoming_light_density: the zenital light density input (umol of photons per m2 per second)

    Returns
    -------

    """

    # can = scene3d(g)
    # zenith = [(100, (0, 0, -1))]
    # (xmin,ymin,xmax,ymax)=(0,0,1,1)
    # cs = CaribuScene(scene=can, light=zenith, pattern=(0,0,1,1))
    # raw, agg = cs.run(infinite=True, simplify=True, direct=False)

    # print("Computing light interception...")

    # We initialize a turtle in PlantGL:
    turtle = turt.PglTurtle()
    # We create the scene with the plant(s).
    # If we want to study the light interception by the isolated plant:
    if isolated_plant:
        # Then we create a scene with only one plant with the classical PlantGL turtle:
        my_scene = turt.TurtleFrame(g, visitor=get_shoot_visitor(leaf_shape=leaf_shape, upscaling_factor=100),
                                    turtle=turtle, gc=False, all_roots=True)
    else:
        # Otherwise, we simulate a population of plants, from which the central plant will be the one on which Caribu
        # will compute the light interception:
        my_scene = create_a_scene_with_plants_popuplation(g,
                                                          number_of_plants_per_square_meter=number_of_plants_per_square_meter,
                                                          inter_row_distance=inter_row_distance,
                                                          scene_length=scene_length,
                                                          scene_width=scene_width,
                                                          leaf_shape=leaf_shape,
                                                          upscaling_factor=100)

    # # EXAMPLES FROM ALINEA:
    # can = str(data_path('f331s1_100plantes.can'))
    # opts = map(str, [data_path('par.opt'), data_path('nir.opt')])
    # pattern = str(data_path('filter.8'))
    # pattern_file = pd.read_csv(pattern, delim_whitespace=True, header=None)

    # If only the zenithal light is considered
    if zenithal_light:
       # We set only one vertical light direction:
       light_direction = (0, 0, -1)
       LIGHT = [(incoming_light_density,light_direction)]
    else:
        # We open the standard light input file as a dataframe:
        sky = str(data_path('Turtle16soc.light'))
        reference_light_file = pd.read_csv(sky, delim_whitespace=True, header=None)
        # The file is made of one line by collimated source, i.e. by direction of the sky hemisphere.
        # Each line contains the source-induced irradiance (i.e. the irradiance of a horizontal surface located just above
        # the canopy), followed by the 3D vector of the direction of light. Ex: (650,0,0,-1) for the sun located at zenith.
        # We transform the dataframe into a new "light" file by multiplying the first column corresponding to the reference
        # radiation by the incoming light density:
        reference_light_file.iloc[:,0] = reference_light_file.iloc[:,0] * incoming_light_density
        # We initialize the final list corresponding to the light inputs:
        LIGHT = []
        # For each light source:
        for row_index in range(0,len(reference_light_file)):
            # We get the coordinates of the light vecto in the last 3 rows:
            coordinates = reference_light_file.loc[row_index, 1:4].values.flatten().tolist()
            # We create a tuple containing the irradiance, then the tuple of coordinates:
            line = (reference_light_file.loc[row_index, 0], tuple(coordinates))
            # And we add this new line to the light file:
            LIGHT.append(line)

    # for shape in my_scene:
    #     print(shape.geometry)
    #     if shape.geometry
    #     for point in shape.geometry.pointList:
    #         new_point = point * 100
    #         print(new_point)

    # In case we want to replicate the scene to infinity:
    if infinite:
        # We create a pattern based on the size of the scene:
        pattern = tuple([0*unit_factor,0*unit_factor,scene_length*unit_factor,scene_width*unit_factor])
    else:
        pattern = None

    # Finally, we create the Caribu scene and run the light interception:
    cs = CaribuScene(scene=my_scene, light=LIGHT, scene_unit='cm', pattern=pattern) #,opt=opts)
    raw, agg = cs.run(infinite=infinite, simplify=True, direct=True)
    # The resulting dictionnary "agg" contains three dictionnaries for each "active" vertex if the MTG:
    #   area (float): the individual areas (m2)
    #   Eabs (float): the surfacic density of energy absorbed (m-2)
    #   Ei (float): the surfacic density of energy incoming (m-2)
    # We update the properties of the MTG with these dictionnaries:
    g.properties()["Caribu_area"].update(agg['area'])
    g.properties()["Caribu_Eabs"].update(agg['Eabs'])
    g.properties()["Caribu_Ei"].update(agg['Ei'])

    # We update the property "PAR_intercepted" with the results coming from Caribu's run:
    g.properties()["PAR_intercepted"].update(agg['Eabs'].items())
    # TODO: Check this once for all: is it the total light per m2 of soil surface? per m2 of leaf surface?
    # Note that the PAR_intercepted is expressed here in umol of photons per m2 of organ (?) per second.

    # print("Ending light interception!")

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function for calculating the actual photosynthesis rate:
def computing_photosynthesis_rate(g, SCA=1e-6):
    """
    This function calculates the actual photosythesis rate (in gC per second) for each organ of the plant g, based on its
    intercepted PAR (micromol of photons per m2 per second), its photosynthetic surface (m2) and a parameter SCA, i.e.
    its specific carbon assimilation (gC per micromol of photons).
    :param g: the shoot MTG containing each organ as node
    :param SCA: specific carbon assimilation (gC per micromol of photon) - this is similar to RUE at the organ scale
                but not for net photosynthesis; it does not vary with the N content of the organ
    :return: [no return - the MTG g has been updated with the property 'photosynthesized_C_rate' (gC per second)]
    """

    # We cover each organ n of the shoot:
    for vid in g.vertices_iter(scale=1):
        n = g.node(vid)
        n.photosynthetic_area = n.surface
        # We calculate the photosynthesis rate for this organ:
        n.photosynthesized_C_rate = SCA * n.PAR_intercepted * n.photosynthetic_area

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function for computing the C cost for maintenance respiration (as opposed to growth respiration):
def computing_maintenance_rate(g, Vmax_maintenance=1., Km_maintenance=1.):
    """
    This function calculates the respiration rate (in gC per second) for the maintenance of each organ of the plant g.
    :param g: the shoot MTG containing each organ as node
    :param Vmax_maintenance: maximal rate of respiration (gC per gram of dry structural mass per second)
    :param Km_maintenance: Michaelis_Menten constant for maintenance (gC per gram of dry structural mass),
                           i.e. the concentration of labile C for which the respiration rate is half of the maximal rate
    :return: [no return - the MTG g has been with the updated property 'maintenance_respiration_rate' (gC per second)]
    """
    # We cover each organ n of the shoot:
    for vid in g.vertices_iter(scale=1):
        n = g.node(vid)

        # If the element does not have any C left, it can't respire!
        if n.C_labile <= 0.:
            n.shoot_maintenance_respiration_rate = 0.
            continue

        # The maintenance rate is based on a maximal rate of respiration per gram of structural mass, modulated by the
        # amount of labile C in the organ:
        n.shoot_maintenance_respiration_rate = Vmax_maintenance * n.C_labile / (Km_maintenance + n.C_labile) * n.struct_mass

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function for redistributing the C among shoot organs:
# either C allocation (1) proportionnal to PARi for each organ or (2) per phytomer then distributed between organs depending on area
def computing_shoot_C_redistribution(g, shoot_available_C = None, intercepted_light = None, time_step_in_seconds=1.):

    # 1) We calculate the amount of shoot C available and the total intercepted light:
    #---------------------------------------------------------------------------------
    if shoot_available_C is None:
        # We initialize the available C in the shoot system to 0:
        available_C = 0.
    else:
        # Otherwise, the available C is directly the provided shoot available C:
        available_C = shoot_available_C
    if intercepted_light is None:
        # We initialize the intercepted light by the shoot system to 0:
        intercepted_light_by_the_plant = 0.
    else:
        # Otherwise, the intercepted light is directly the light that has been provided:
        intercepted_light_by_the_plant = intercepted_light

    # If either shoot_available_C or intercepted_light have not been provided:
    if shoot_available_C is None or intercepted_light is None:
        # Then we cover each organ n of the shoot:
        for vid in g.vertices_iter(scale=1):
            n = g.node(vid)
            # We add the initial amount of C_labile and the net gain in C of the current organ:
            available_C += n.C_labile + (n.photosynthesized_C_rate - n.shoot_maintenance_respiration_rate) * time_step_in_seconds
            # We add the light intercepted by the organ:
            intercepted_light_by_the_plant += n.PAR_intercepted

    # 2) We redistribute the C according to specific rules:
    #------------------------------------------------------

    # We cover each organ n of the shoot:
    for vid in g.vertices_iter(scale=1):
        n = g.node(vid)
        # The corresponding rate of C redistribution for this organ is calculated:
        # (NOTE: if positive, it corresponds to a C gain; if negative, it corresponds to a C source)
        n.shoot_C_redistribution_rate = available_C * n.PAR_intercepted / intercepted_light_by_the_plant / time_step_in_seconds

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function for distributing the N coming from the roots into each shoot organ:
def computing_shoot_N_distribution(g, shoot_available_N=None, time_step_in_seconds=60.*60.*24.):
    """
    This function distributes the whole amount of labile N available to the whole plant to each organ which has a
    positive surface, proportionnaly to the mean surfacic PAR intercepted by the organ compared to the sum of all
    surfacic PARs of shoot elements.

    :param g: the MTG corresponding to the shoot system
    :param time_step_in_seconds: the time step of the simulation (in seconds)
    return: [no return - the MTG g has been updated with the property 'N_labile']

    """

    # 1) Computing the intercepted light (and possibly the available N to distribute):
    #---------------------------------------------------------------------------------

    # We initialize the total amount of light intercepted by the whole plant:
    intercepted_light_by_the_plant = 0.
    sum_of_mean_surfacic_PAR = 0.

    # If the amount of N available in the shoot system has not been provided yet:
    if shoot_available_N is None:
        # We initialize the available N in the shoot system to 0:
        available_N = 0.
    else:
        # Otherwise, the available N is directly the provided shoot available N:
        available_N = shoot_available_N

    # Then we cover each organ n of the shoot:
    for vid in g.vertices_iter(scale=1):
        n = g.node(vid)
        if n.surface <= 0.:
            continue
        # We compute not only the total intercepted light, but also the sum of the mean surfacic PARs:
        intercepted_light_by_the_plant += n.PAR_intercepted
        n.mean_surfacic_PAR = n.PAR_intercepted / n.surface
        sum_of_mean_surfacic_PAR +=  n.mean_surfacic_PAR
        # If the amount of shoot available N has not been provided, we add the amount of available N of the current element:
        if shoot_available_N is None:
            available_N += n.N_labile - n.N_labile_deficit

    # print("Before N redistribution among the shoot organs, the total amount of labile N was", available_N)

    # 2) Computing the rate of N distribution in each shoot organ
    # -----------------------------------------------------------

    # We cover each organ of the shoot:
    for vid in g.vertices_iter(scale=1):
        # We access the current organ"n":
        n = g.node(vid)
        if n.surface <= 0.:
            continue
        # The new target for C_labile in this organ is calculated according to the relative intercepted light:
        # N_labile_target = plant_available_N * n.PAR_intercepted / intercepted_light_by_the_plant
        fraction = n.mean_surfacic_PAR / sum_of_mean_surfacic_PAR
        # print("   xxx For element", vid, "the fraction of light was", fraction)
        N_labile_target = available_N * fraction
        n.N_input_rate = N_labile_target / time_step_in_seconds
        # print("                    Its N input rate was", n.N_input_rate)

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
def computing_shoot_N_remobilization(g, k_N_mobilization = 1/(24.*60.*60.)):
    """
    This function calculates the rate of N mobilization (gN per second) in each shoot organ, based on first-order kinetics.
    :param g: the MTG corresponding to the shoot system
    :param k_N_mobilization: the time constant of the N mobilization (second-1)
    :return: [no return - the MTG g has been updated with the property 'N_mobilization_rate']
    """

    # We cover each organ of the shoot:
    for vid in g.vertices_iter(scale=1):
        # We access the current organ"n":
        n = g.node(vid)

        # We exclude elements where a deficit of labile N has been previously recorded:
        if n.N_labile_deficit >= 0.:
            n.shoot_N_mobilization_rate = 0.
            continue

        # We compute the rate N mobilization:
        n.shoot_N_mobilization_rate = k_N_mobilization * n.N_labile

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
def computing_shoot_N_immobilization(g):
    """
    This function calculates the rate of N immobilization (gN per second) in each shoot organ.
    :param g: the MTG corresponding to the shoot system
    :return: [no return - the MTG g has been updated with the property 'net_shoot_N_immobilization_rate']
    """

    # We cover each organ of the shoot:
    for vid in g.vertices_iter(scale=1):
        # We access the current organ"n":
        n = g.node(vid)
        # We compute the rate N mobilization:
        n.net_shoot_N_immobilization_rate = 0.

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
def computing_shoot_C_storage(g, time_step_in_seconds=60.*60.*24.):
    """
    This function calculates the rate of labile C reallocation to each leaf.
    :param g: the MTG corresponding to the shoot system
    :param time_step_in_seconds: the time step of the simulation (seconds)
    :return: [no return - the MTG g has been updated with the property 'C_storage_rate' (gC per second)]
    """

    # 1) Computing the total intercepted light (if not calculated elsewhere)
    #-----------------------------------------------------------------------
    # OPTION 1: if it has not been calculated yet
    # We initialize the total amount of light intercepted by the whole plant:
    intercepted_light_by_the_plant = 0.
    # We cover each organ n of the shoot:
    for vid in g.vertices_iter(scale=1):
        n = g.node(vid)
        intercepted_light_by_the_plant += n.PAR_intercepted

    # # OPTION 2: if it has already been calculated and stored in the base of the shoot system:
    # intercepted_light_by_the_plant = g.node(1).total_PAR_intercepted_by_the_plant

    # 2) Computing the rate of C distribution in each shoot organ
    # -----------------------------------------------------------

    # We access the value of the total amount of available C, stored in the base of the shoot system:
    plant_available_C = g.node(1).shoot_available_C

    # We cover each organ of the shoot:
    for vid in g.vertices_iter(scale=1):
        # We access the current organ"n":
        n = g.node(vid)
        # The new C storage rate in this organ is calculated according to the relative intercepted light:
        n.shoot_C_storage_rate = plant_available_C * n.PAR_intercepted / intercepted_light_by_the_plant / time_step_in_seconds

    # We reset the value of the plant available C:
    g.node(1).shoot_available_C = 0.

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

def potential_shoot_growth(g,
                           SFQN=1,
                           LA_max=10,
                           LA_max_relative_increase_per_day = 0.5,
                           leaf_minimal_thickness = 1e-3,
                           elong_internode=5e-7,
                           time_step_in_seconds = 60.*60.*24.):
    """
    This function computes the potential increase in volume of each organ in the shoot system.
    :param g: the shoot system
    :param SFQN: the slope of the relationship of Leaf area = f(N_labile) (m2 per gN)
    :param LA_max: the maximal area for a leaf (m2)
    :param leaf_minimal_thickness: the minimal thickness of a leaf (m)
    :param elong_internode: a potential rate of elongation (m per second)
    :param time_step_in_seconds: the time step of the simulation (seconds)
    :return:
    """

    # We cover each organ n of the shoot:
    for vid in g.vertices_iter(scale=1):
        # We define the current organ as "n":
        n = g.node(vid)
        # For a leaf:
        if n.label == "Leaf":
            # We calculate the potential surface based on the available N in the leaf:
            potential_leaf_surface = SFQN * n.N_labile
            # We calculate the potential surface increase (m2), knowing that the surface can't become higher than LA_max:
            potential_surface_increase = min(potential_leaf_surface - n.surface, LA_max - n.surface)
            # We calculate the maximal surface increase (m2), knowing that the surface increase can't become higher than the maximal increase:
            maximal_surface_increase = n.surface * LA_max_relative_increase_per_day * (time_step_in_seconds / (60*60.*24.))
            # The potential volume increment is based on a minimal thickness:
            n.potential_volume_increment = min(potential_surface_increase, maximal_surface_increase) * leaf_minimal_thickness
            # TODO: what happens if the leaf has already a thickness higher than the minimal one?
        # For an internode:
        elif n.label == "Internode":
            # We calculate a potential elongation rate:
            potential_elongation_rate = elong_internode * n.radius * 2.
            potential_elongation_rate = elong_internode * n.radius * 2.
            # TODO: the elongation rate of the internode should react to a signal or plant stage!
            # TODO: an internode should also increase in diameter over time!
            # And we compute the potential volume increment:
            n.potential_volume_increment = potential_elongation_rate * (n.radius ** 2 * pi) * time_step_in_seconds
        else:
            n.potential_volume_increment = 0.
        # TODO: define a case for the transformation of axillary buds into active apices!

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
def computing_C_redistribution_for_shoot_growth(g,
                                                yield_growth=0.8,
                                                structural_mass_density=1e5,
                                                structural_C_content=0.44,
                                                time_step_in_seconds = 60.*60.*24.
                                                ):

    # We initialize the stock of C that can be redistributed to sink leaves and the one that will be eventually used:
    total_available_C = 0.
    C_used_by_sink_leaves = 0.
    # We initialize the total intercepted light:
    # We add the light intercepted by the organ:
    intercepted_light = 0.
    # We initialize a list that will contain the indices of the leaves that will need additional C:
    list_of_sink_leaves = []

    # 1) DEFINING THE AMOUNT OF AVAILABLE C:
    ########################################

    # We cover each organ n of the shoot:
    for vid in g.vertices_iter(scale=1):
        # We define the current organ as "n":
        n = g.node(vid)

        # # We exclude any organ that is not a leaf:
        # if n.label != "Leaf":
        #     continue
        # TODO: should other organs be excluded in the C redistribution for growth?

        # We add the light intercepted by the organ:
        intercepted_light += n.PAR_intercepted

        # We calculate the C that may be consumed by the potential growth, and the C that may remain in a "source" leaf:
        C_consumption_by_growth = (n.potential_volume_increment * structural_mass_density
                                                 * structural_C_content / yield_growth)
        remaining_C_labile = n.C_labile - C_consumption_by_growth

        # If there is some C left after the potential growth:
        if remaining_C_labile > 0.:
            # The leaf is a source, and it increases the total amount of available C:
            n.C_sink_status = False
            total_available_C += remaining_C_labile
        else:
            # The leaf is a sink, and may receive additional C for its growth:
            n.C_sink_status = True
            list_of_sink_leaves.append(vid)

    # 2) RE-ALLOCATING C TO SINK LEAVES:
    ####################################

    # We cover each organ n of the shoot:
    for vid in list_of_sink_leaves:
        # We define the current organ as "n":
        n = g.node(vid)
        # The new amount of C_labile in this organ is calculated according to the relative intercepted light:
        additional_C_for_sink_leaf = total_available_C * n.PAR_intercepted / intercepted_light
        # TODO: Here, only the fraction of C corresponding to (light_intercepted_by_sink_leaves) / (total_light_intercepted_by all_leaves) is used...
        # TODO: Watch out - leaves that have not emerged don't get any C!!!

        # n.C_labile += additional_C_for_sink_leaf
        n.C_input_rate = additional_C_for_sink_leaf / time_step_in_seconds
        C_used_by_sink_leaves += additional_C_for_sink_leaf

    # Eventually, we store the remaining available C in the first node:
    base = g.node(1)
    base.shoot_available_C = total_available_C - C_used_by_sink_leaves # The remaining amount of C here may be allocated to roots!

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
def actual_shoot_growth(g,
                        yield_growth = 0.8,
                        structural_mass_density = 1e5,
                        structural_C_content = 0.44,
                        structural_N_content = 0.02,
                        leaf_minimal_thickness = 1e-3,
                        time_step_in_seconds = 60.*60.*24.):
    """
    This function computes the actual growth of the shoot system, depending on the C and N available in each organ.
    :param g: the shoot system
    :param yield_growth: the fraction of the total C used for growth that will be incorporated to the new mass (the rest is respired) (gC per gC)
    :param struct_mass_density: the density of the structural mass (gram of structural maqss per m3)
    :param structural_C_content: the C content of the structural mass (gC per gram of structural mass)
    :param structural_N_content: the N content of the structural mass (gN per gram of structural mass)
    :param leaf_minimal_thickness: the minimal thickness of a leaf (m)
    :param time_step_in_seconds: the time step (seconds)
    :return:
    """

    # We cover each organ n of the shoot:
    for vid in g.vertices_iter(scale=1):

        # We reinitialize temporary variables:
        final_volume_increment = 0.
        length_increment = 0.

        # We define the current organ as "n":
        n = g.node(vid)

        # If the element does not have any potential growth, the we move to the next element:
        if n.potential_volume_increment <=0.:
            continue
        # If the element does not have any C or N left, it cannot grow:
        if n.C_labile <=0. or n.N_labile <=0.:
            continue

        # The maximal volume increment is defined according to the amount of labile C or N available for growth:
        max_volume_increment_by_C = n.C_labile / (structural_mass_density * structural_C_content / yield_growth)
        max_volume_increment_by_N = n.N_labile / (structural_mass_density * structural_N_content)

        # The final volume is the minimal volume between the potential volume calculated from potential growth
        # and the maximal volumes that can be created from the available C and N:
        final_volume_increment = min(n.potential_volume_increment, max_volume_increment_by_C, max_volume_increment_by_N)

        # The new volume of the organ and the C and N used for creating the corresponding structural mass is calculated:
        n.volume += final_volume_increment
        n.actual_C_consumption_rate_by_shoot_growth = (final_volume_increment * structural_mass_density
                                                       * structural_C_content / yield_growth / time_step_in_seconds)
        n.actual_N_consumption_rate_by_shoot_growth = (final_volume_increment * structural_mass_density
                                                       * structural_N_content / time_step_in_seconds)
        # The C respired because of growth is also calculated:
        n.shoot_growth_respiration_rate = n.actual_C_consumption_rate_by_shoot_growth * (1 - yield_growth)

        # The final dimensions of the organ are recorded:
        if n.label == "Internode":
            length_increment = final_volume_increment / (n.radius ** 2 * pi)
            n.length += length_increment
            n.surface = n.length * n.radius * 2 * pi
        elif n.label == "Leaf":
            surface_increment = final_volume_increment / leaf_minimal_thickness
            n.surface += surface_increment
            # We assume that the leaf is an ellipse characterized by its length and its width, which is calculated as
            # the product of the length with an allometric coefficient. The surface of the ellipse is then defined as:
            # surface = (length / 2.) * (param.leaf_width_to_length_ratio * length) / 2. * pi
            # i.e. surface = length ^ 2 * param.leaf_width_to_length_ratio / 4. * pi
            # Hence:
            length_increment = sqrt(n.surface * 4 / (param.leaf_width_to_length_ratio * pi)) - n.length
            # The length is incremented with the corresponding result:
            n.length += length_increment
            # The radius, which is half of the width, is also calculated based on the allometric coefficient
            n.radius += length_increment * param.leaf_width_to_length_ratio / 2.

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
def updating_coordinates(g):

    # We define "root" as the starting point of the loop below:
    root_gen = g.component_roots_at_scale_iter(g.root, scale=1)
    root = next(root_gen)

    # We travel in the MTG from the base to the extremity:
    for vid in pre_order(g, root):

        # We define the current element as n:
        n = g.node(vid)
        # We make sure that the current element has a positive length:
        if n.length<=0.:
            continue

        # We get the current coordinates of the vector corresponding to the element:
        x1 = n.x1
        x2 = n.x2
        y1 = n.y1
        y2 = n.y2
        z1 = n.z1
        z2 = n.z2

        # The direction of the vector from point 1 to point 2 will be unchanged:
        direction_x = x2 - x1
        direction_y = y2 - y1
        direction_z = z2 - z1

        # We need to change the coordinates of the new point 2 based on the updated length.
        # The distance between point 1 and the original point 2 can be calculated by:
        d = ((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)**0.5

        try:
            x = 1/d
        except:
            print("ERROR: the direction of the element", n.index(), "could not be calculated!")
        # The new coordinates of the second point of the current element are then:
        n.x2 = x1 + n.length / d * direction_x
        n.y2 = y1 + n.length / d * direction_y
        n.z2 = z1 + n.length / d * direction_z

        # We try then to update the x1,y1,z1 coordinates of the next element, if any:
        son_id = g.Successor(vid)
        if son_id:
            son = g.node(son_id)
            son.x1 = n.x2
            son.y1 = n.y2
            son.z1 = n.z2

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
def computing_senescence(g):

    # We cover each organ of the shoot:
    for vid in g.vertices_iter(scale=1):

        # We define the current element as n:
        n = g.node(vid)
        # We make sure that the current element has a positive length:
        if n.length<=0.:
            continue

        if n.C_labile_deficit >= 0.:
            n.senescence_status = True

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
def shoot_organ_CN_balance(shoot_MTG, time_step_in_seconds = 60.*60.*24., resetting_rates=False):

    g = shoot_MTG

    # We cover each organ of the shoot:
    for vid in g.vertices_iter(scale=1):

        # We define the current element as n:
        n = g.node(vid)
        # We make sure that the current element has a positive length:
        if n.length<=0.:
            continue

        # COMPUTING THE NEW RESULTING NET INCREASE RATES IN C & N:
        # We compute the net C increase in the labile C pool associated to the organ:
        net_C_increase_rate = (n.photosynthesized_C_rate
                               + n.C_input_rate # TODO: should C_storage_rate be included in the balance of the structural mass or here?
                               - n.shoot_maintenance_respiration_rate
                               - n.actual_C_consumption_rate_by_shoot_growth)
        # We compute the net N increase in the labile N pool associated to the organ:
        net_N_increase_rate = (n.N_input_rate
                               - n.actual_N_consumption_rate_by_shoot_growth
                               - n.shoot_N_mobilization_rate)

        # COMPUTING THE NEW RESULTING NET GAIN IN LABILE C & N:
        n.C_labile += net_C_increase_rate * time_step_in_seconds - n.C_labile_deficit
        if n.C_labile < 0.:
            print("    !!! Warning: for the shoot element", n.index(), "a deficit in labile C was obtained! We therefore set this pool to 0 gC and will record its deficit for the next time!")
            n.C_labile_deficit = - n.C_labile
            n.C_labile = 0.

        n.N_labile += net_N_increase_rate * time_step_in_seconds - n.N_labile_deficit
        if n.N_labile < 0.:
            print("    !!! Warning: for the shoot element", n.index(), "a deficit in labile N was obtained! We therefore set this pool to 0 gN and will record its deficit for the next time!")
            n.N_labile_deficit = - n.N_labile
            n.N_labile = 0.

        # # TODO: Check whether the C storage increases the biomass without changing the leaf thickness and without C:N stoechiometry!
        # n.struct_mass += n.shoot_C_storage_rate / param.struct_mass_C_content * time_step_in_seconds

    if resetting_rates:
        n.photosynthesized_C_rate = 0.
        n.C_input_rate= 0.
        n.shoot_maintenance_respiration_rate = 0.
        n.actual_C_consumption_rate_by_shoot_growth = 0.
        N_input_rate = 0.
        n.actual_N_consumption_rate_by_shoot_growth = 0.
        n.shoot_N_mobilization_rate = 0.

    return

def summing_in_shoot_system(shoot_system, shoot_MTG, resetting_rates=False):

    # We rename the shoot system as "s" and the shoot MTG as "g":
    s = shoot_system
    g = shoot_MTG

    # TODO: WATCH OUT !!!!!!!!
    s.shoot_labile_C_pool = 0.
    s.shoot_labile_N_pool = 0.

    # print("Before summing, the shoot labile N pool is", s.shoot_labile_N_pool)

    # We cover each organ of the shoot:
    for vid in g.vertices_iter(scale=1):

        # We define the current organ as "n":
        n = g.node(vid)

        # We add its "additive" properties to the general properties of the shoot system:
        s.shoot_labile_C_pool += n.C_labile
        s.shoot_labile_N_pool += n.N_labile
        # print("   > For the element", vid, "we just added the following amounts of C & N:", n.C_labile, n.N_labile)

        s.photosynthesized_C_rate += n.photosynthesized_C_rate
        s.shoot_maintenance_respiration_rate += n.shoot_maintenance_respiration_rate
        s.actual_C_consumption_rate_by_shoot_growth += n.actual_C_consumption_rate_by_shoot_growth
        s.net_shoot_N_immobilization_rate += n.net_shoot_N_immobilization_rate
        s.actual_N_consumption_rate_by_shoot_growth += n.actual_N_consumption_rate_by_shoot_growth

        if resetting_rates:
            n.photosynthesized_C_rate = 0.
            n.shoot_maintenance_respiration_rate = 0.
            n.actual_C_consumption_rate_by_shoot_growth = 0.
            n.net_shoot_N_immobilization_rate = 0.
            n.actual_N_consumption_rate_by_shoot_growth = 0.

    # print("After summing, the shoot labile N pool is", s.shoot_labile_N_pool)

    return

def shoot_CN_balance(shoot_system, time_step_in_seconds = 60.*60):

    # We rename the shoot system as "s":
    s = shoot_system

    # We update the pool of labile C, taking into account the possible deficit in C previously recorded,
    # and all processes affecting the C pool:
    s.shoot_labile_C_pool += (-s.shoot_labile_C_deficit
                             + (- s.net_C_transfer_from_shoot_to_root_rate
                                + s.photosynthesized_C_rate
                                - s.shoot_maintenance_respiration_rate
                                - s.actual_C_consumption_rate_by_shoot_growth) * time_step_in_seconds)
    # We reset the shoot labile C deficit to 0:
    s.shoot_labile_C_deficit = 0.

    # print("Before N balance, the labile N pool in shoots was", s.shoot_labile_N_pool)
    # print("The deficit in N is", s.shoot_labile_N_deficit, "mol of N and here are the fluxes involved:", s.net_N_transfer_from_root_to_shoot_rate, s.net_shoot_N_immobilization_rate,
    #       s.actual_N_consumption_rate_by_shoot_growth)
    # We update the pool of labile N, taking into account the possible deficit in N previously recorded,
    # and all processes affecting the N pool: n.N_input_rate - n.actual_N_consumption_rate_by_growth - n.N_mobilization_rate
    s.shoot_labile_N_pool += (- s.shoot_labile_N_deficit +
                             (s.net_N_transfer_from_root_to_shoot_rate
                              - s.net_shoot_N_immobilization_rate
                              - s.actual_N_consumption_rate_by_shoot_growth) * time_step_in_seconds)
    # print("After N balance, the labile N pool in shoots is", s.shoot_labile_N_pool)
    # We reset the shoot labile N deficit to 0:
    s.shoot_labile_N_deficit = 0.

    # If the pools of C & N are negative, we reset them to 0 and adjust the deficits that will be taken into account
    # at the next step:
    if s.shoot_labile_C_pool < 0.:
        s.shoot_labile_C_deficit = - s.shoot_labile_C_pool
        s.shoot_labile_C_pool = 0.
    if s.shoot_labile_N_pool < 0.:
        s.shoot_labile_N_deficit = - s.shoot_labile_N_pool
        s.shoot_labile_N_pool = 0.

    return

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
def calling_the_different_shoot_modules(shoot_g,
                                        time_step_in_days=1,
                                        phytomers_initiation = True,
                                        number_of_plants_per_square_meter=33,
                                        inter_row_distance=0.45,
                                        scene_length=1.,
                                        scene_width=1.,
                                        C_allocation_option="Full_priority_to_leaves"
                                        ):

    # DIFFERENT POSSIBILITES OF C ALLOCATION CAN BE TESTED:
    # -"Full_priority_to_leaves": Labile C is used by all leaves before C allocation to roots
    # -"No_priority_to_leaves": We get a pool of labile C available after photosynthesis and before growth,
    # and we allocate a fraction to roots, then to all leaves according to PAR intercepted
    # -"C_allocation_to_roots_before_C_storage": We first achieve the actual growth of source leaves and of the sink
    # leaves before C reallocation, then we allocate C to roots, and the remaining is used for reallocating C to sink leaves.

    g = shoot_g

    # WE SUCCESSIVELY CALL THE DIFFERENT MODULES:

    # 1) Inputs of C & N in each organ, and maintenance:
    #---------------------------------------------------
    # We calculate light interception based on the architecture:
    computing_light_interception(g,
                                 number_of_plants_per_square_meter=number_of_plants_per_square_meter,
                                 inter_row_distance=inter_row_distance,
                                 scene_length=scene_length,
                                 scene_width=scene_width
                                 )

    # We calculate the C production based on PAR:
    computing_photosynthesis_rate(g,
                                  SCA = param.SCA)

    # We calculate the C respiration cost associated to maintenance:
    computing_maintenance_rate(g,
                               Vmax_maintenance=param.Vmax_maintenance,
                               Km_maintenance = param.Km_maintenance)

    # We allocate the N available from the base of the shoot up to each leaf, based on PAR:
    computing_shoot_N_distribution(g,
                                   time_step_in_seconds=time_step_in_days*60.*60.*24.)

    # This calculates the total amount of available C, and distributes it in each organ according to PAR:
    # computing_shoot_C_redistribution(g)

    # 2) Growth of organs:
    #---------------------
    # We calculate the possible increase in volume of each organ:
    potential_shoot_growth(g,
                           SFQN=param.SFQN, #TODO: is SFQN constant over the whole range of (0, LA_max)? Is it two distinct curves?
                           LA_max=param.LA_max,
                           leaf_minimal_thickness=param.leaf_minimal_thickness,
                           elong_internode=param.elong_internode,
                           time_step_in_seconds=time_step_in_days*60.*60.*24.
                           )

    # We calculate the fraction of C that can be redistributed by source leaves, and distribute a fraction of it in sink leaves:
    computing_C_redistribution_for_shoot_growth(g,
                                                yield_growth=param.yield_growth,
                                                structural_mass_density=param.struct_mass_density,
                                                structural_C_content=param.structural_C_content,
                                                time_step_in_seconds=time_step_in_days*60.*60.*24.
                                                )

    # We calculate the actual increase in size and the corresponding C & N costs:
    actual_shoot_growth(g,
                        yield_growth=param.yield_growth,
                        structural_mass_density=param.struct_mass_density,
                        structural_C_content=param.structural_C_content,
                        structural_N_content=param.structural_N_content,
                        leaf_minimal_thickness=param.leaf_minimal_thickness,
                        time_step_in_seconds=time_step_in_days*60.*60.*24.
                        )

    # We initiate the structure of new phytomers (if any):
    new_phytomers_initiation(g, new_initiation=phytomers_initiation)

    # We update all dimensions of the organs:
    updating_coordinates(g)

    # 3) Mobilization/immobilization:
    #--------------------------------

    # This function reallocates all the remaining C available to the different organs according to their PAR,
    # and stores it as additional structural mass:
    computing_shoot_C_storage(g,time_step_in_seconds=time_step_in_days*60.*60.*24.)

    # This function mobilizes N from the organs and stores it at the base of the shoot system:
    computing_shoot_N_remobilization(g, k_N_mobilization=param.k_N_mobilization)

    # This function immobilizes N in each shoot organ:
    computing_shoot_N_immobilization(g)

    # 4) C&N mass balance:
    #---------------------
    # This function performs the actual C & N balance and updates the labile pools of C & N in each organ:
    shoot_organ_CN_balance(g)

    # 5) Senescence and leaf fall:
    #-----------------------------
    computing_senescence(g)

    return

def loading_parameters_from_file(parameters_filepath = "parameters_AllocF.xlsx",
                                 outputs_path = "outputs"):

    # We read the scenario to be run:
    scenarios_df = pd.read_excel(parameters_filepath, header=0, sheet_name="scenarios_as_columns")
    # We remove the columns containing unnecessary details about parameters:
    useless_columns = ["Explanation", "Type / Unit", "Reference value", "Justification of the reference value"]
    scenarios_df.drop(useless_columns, axis=1, inplace=True)
    # Before transposing the dataframe, we rename the first column as 'Scenario':
    scenarios_df.rename(columns={'Parameter': 'Scenario'}, inplace=True)
    # Because the way of dealing with data type is not well suited with the current Excel file, we transpose the
    # dataframe, record it as a CSV file in the outputs folder, and reload it with the proper type for each
    # parameter:
    transposed_df = scenarios_df.T
    transposed_df.to_csv(os.path.join(outputs_path, 'scenarios_list.csv'), na_rep='NA', header=False)
    new_scenarios_df = pd.read_csv(os.path.join(outputs_path, 'scenarios_list.csv'), index_col='Scenario', header=0)
    # We convert the dataframe into a dictionnary:
    scenario = new_scenarios_df.to_dict(orient='records')[0]
    # And we finally update the parameters from the original python file with the parameters taken from this dictionnary:
    param.__dict__.update(scenario)

    return

def ramification(g):

    # On part sur une structure potentielle maximale paramétrée.
    # Le virage floral est atteint au bout d'un certain temps thermique (entre mi-décembre et mi-janvier).

    # Au virage floral, le méristème principal ne fera plus de feuilles et que une fleur (ou plusieurs).
    # Le nombre total de phytomères sur la tige principal sera:
    # n_phytomers_formés (y compris non émergés) = k * (n_feuilles_rosette_émergées au virage floral)
    # Modèle "sapin" pour définir le nombre d'entre-noeuds sur les ramifs primaires (+1 en descendant du haut de la tige vers le bas)
    # Pour les secondaires : plutôt en haut...

    # Le profil de longueur d'entre-noeuds sur la tige principale est plutôt ventru (max au milieu de la tige, pas en bas).
    # Elongation forcée + limitation trophique + dominance apicale (depuis le haut)

    # La date d'élongation est différente du virage floral (vers février-mars). Temps thermique constant ?

    return

########################################################################################################################
########################################################################################################################

# MAIN SIMULATION:
##################

if __name__ == '__main__':

    # my_path = os.getcwd()
    my_path = "C:/Users/frees/SIMBAL/"

    # 1) INITIAL CONDITIONS
    #----------------------
    # OPTION 1: We create a MTG from its known properties (NOT COMPLETELY FINISHED!):
    # # We read the input file:
    # input_dataframe = pd.read_excel('APTOS_input_1.xlsx', header=0, sheet_name="scenario")
    # g = create_MTG_from_csv_file()

    # OPTION 2: We create an initial MTG corresponding to the shoot (stem and leaves):
    g = create_shoot_MTG(number_of_phytomers=5)

    # We delete the previous files containing plots or MTG properties:
    plot_directory = os.path.join(my_path, "simulation", "outputs","plots")
    # If this directory doesn't exist:
    if not os.path.exists(plot_directory):
        # Then we create it:
        os.mkdir(plot_directory)
    else:
        # Otherwise, we delete all the images that are already present inside:
        for root, dirs, files in os.walk(plot_directory):
            for file in files:
                os.remove(os.path.join(root, file))
    MTG_directory = os.path.join(my_path, "simulation", "outputs", "MTG_properties")
    # If this directory doesn't exist:
    if not os.path.exists(MTG_directory):
        # Then we create it:
        os.mkdir(MTG_directory)
    else:
        # Otherwise, we delete all the images that are already present inside:
        for root, dirs, files in os.walk(MTG_directory):
            for file in files:
                os.remove(os.path.join(root, file))

    # UPDATING PARAMETERS FROM FILE:
    loading_parameters_from_file(parameters_filepath = os.path.join(my_path, "simulation", "inputs", "simbal_parameters.xlsx"),
                                     outputs_path = os.path.join(my_path, "simulation", "outputs"))

    # # We display the MTG in its initial state with PlantGL:
    # empty_scene = pgl.Scene()
    # pgl.Viewer.display(empty_scene)
    # plotting_shoot_MTG(g, width=1200, height=900,
    #                x_center=0, y_center=0, z_center=2.,
    #                x_cam=5., y_cam=5., z_cam=5.,
    #                grid=False,
    #                background_color=[0, 0, 0],
    #                single_color=None,
    #                property_name='PAR_intercepted', cmap='jet', vmin=0.1, vmax=100, lognorm=False,
    #                recording_plot=True,
    #                plot_name=os.path.join(plot_directory,'plot_0000.png'))

    # # Plotting the MTG with Pyvista:
    # plotting_shoots_with_pyvista(g)

    # We record the properties in its initial state:
    recording_MTG_properties(g, os.path.join(MTG_directory, 'MTG_00000.csv'))

    # We define an initial and final step for covering a certain period:
    initial_time = 0.
    simulation_period_in_days = 50.
    time_step_in_days = 1.
    n_steps = trunc(simulation_period_in_days / time_step_in_days)
    elapsed_time_in_days_since_last_leaf_emergence = 0.

    plant_density_in_square_meter = 33
    inter_row_distance = 0.45
    row_length = 1.
    scene_width = 1.

    # 2) OVER TIME:
    #--------------
    # For each time step covering the scenarios period:
    for step in range(1, n_steps+1):

        # We calculate the current time in hours:
        initial_time_in_days = initial_time + (step-1) * time_step_in_days
        final_time_in_days = initial_time + step * time_step_in_days
        print("")
        print("From", initial_time_in_days,"to", final_time_in_days,"days:")
        print("-----------------------")

        # CALLING THE EXCHANGE WITH ROOTS:
        # TODO: change this forced inputs with the actual connection with RootSimple
        g.node(1).shoot_available_N = 10.

        # CALLING SUCCESSIVELY THE DIFFERENT SHOOT MODULES:

        # We first examine whether new phytomers will be produced at this time step.
        # We reset the authorization to add a new phytomer:
        phytomer_formation = False
        # We increment the time elapsed since the last phytomer formation with the current time step:
        elapsed_time_in_days_since_last_leaf_emergence += time_step_in_days
        # If the elapsed time since last phytomer formation is higher than the prescribed phyllochrone:
        if elapsed_time_in_days_since_last_leaf_emergence >= param.phyllochron_in_days:
            # We authorize the formation of a new phytomer:
            phytomer_formation = True
            # And we reset the elapsed time:
            elapsed_time_in_days_since_last_leaf_emergence = elapsed_time_in_days_since_last_leaf_emergence - param.phyllochron_in_days

        # We call the different modules:
        calling_the_different_shoot_modules(shoot_g = g,
                                            phytomers_initiation=phytomer_formation,
                                            number_of_plants_per_square_meter=plant_density_in_square_meter,
                                            inter_row_distance=inter_row_distance,
                                            scene_length=row_length,
                                            scene_width=scene_width)


        # PLOTTING ONE SINGLE PLANT OVER TIME:
        # We display the MTG in its new state:
        # Option 1A - With PlantGL:
        #--------------------------
        graph = plotting_shoot_MTG(g, width=1200, height=900,
                               x_center=0, y_center=0, z_center=0.,
                               x_cam=3, y_cam=3, z_cam=1.,
                               grid=False,
                               background_color=[0, 0, 0],
                               # single_color=[100,150,50],
                               # property_name='C_labile', cmap='jet', vmin=1e-5, vmax=1e1, lognorm=True,
                               property_name='PAR_intercepted', cmap='jet', vmin=0, vmax=600, lognorm=False,
                               # property_name='Caribu_Eabs', cmap='jet', vmin=0., vmax=1., lognorm=False,
                               # property_name='z1', cmap='jet', vmin=0., vmax=0.5, lognorm=False,
                               recording_plot=True,
                               plot_name = os.path.join(plot_directory, 'plot_%.5d.png' % step))

        # # Option 1B - With Pyvista:
        # #--------------------------
        # plotting_shoots_with_pyvista(g)

        # # PLOTTING A POPULATION OF PLANTS:
        # # Option 2A - With PlantGL:
        # #--------------------------
        # plotting_shoot_MTG_population(g,
        #                               number_of_plants_per_square_meter=plant_density_in_square_meter,
        #                               inter_row_distance=inter_row_distance,
        #                               scene_length=row_length,
        #                               scene_width=scene_width,
        #                               leaf_shape="oval",
        #                               # leaf_shape="typical_leaf_shape",
        #                               width=800, height=800,
        #                               x_center=0., y_center=-0., z_center=0.,
        #                               x_cam=3., y_cam=3., z_cam=1.,
        #                               grid=False,
        #                               background_color=[0, 0, 0],
        #                               single_color=None,
        #                               property_name='PAR_intercepted', cmap='jet', vmin=0, vmax=400, lognorm=False,
        #                               recording_plot=True,
        #                               plot_filepath=os.path.join(plot_directory, 'plot_%.5d.png' % step))

        # RECORDING THE PROPERTIES:
        # We can record the main properties of the plant:
        dataframe = recording_MTG_properties(g, file_name = os.path.join(MTG_directory, 'MTG_%.5d.csv' % step))

        print("   > The shoot system now contains", g.nb_vertices(),"elements!")

    print("")
    print("Simulation is done!")

