#  -*- coding: utf-8 -*-

"""
    The script 'AllocF' contains all functions for running the development anf functionning of the aerial parts of the plant.

    :copyright: see AUTHORS.
    :license: see LICENSE for details.
"""

import openalea.mtg
from openalea.mtg import *
from openalea.mtg import turtle as turt
from openalea.mtg.traversal import pre_order, post_order
from openalea.mtg.plantframe import color
import openalea.plantgl.all as pgl

from alinea.caribu.CaribuScene import CaribuScene
from alinea.caribu.data_samples import data_path

import simbal_parameters as param
from tools import (loading_parameters_from_file, temperature_modification, recording_MTG_properties,
                   circle_coordinates, updating_coordinates, create_a_scene_with_plants_popuplation,
                   plotting_shoot_MTG, plotting_shoots_with_pyvista)

import os
import numpy as np
import pandas as pd
from math import sqrt, pi, floor, trunc, exp, isnan, cos, sin, acos

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
                      leaf_width=0.,
                      leaf_length=0.,
                      bud_length=0.,
                      properties_identical_to_parent = True,
                      nil_properties = False
                      ):
    """
    This function adds one phytomer (internode + leaf + bud) on the same axis as the parent element. Note that the leaf
     is defined as the direct child of the internod, and the bud is the direct child of the leaf (not the internode).

    Parameters
    ----------
    parent_element: the parent element above which the phytomer is added
    struct_mass_density: the density of the structural mass (gDW per cubic meter)
    minimal_leaf_thickness: the tickness of the leaf to be created (meters)
    internode_radius: the radius of the internode to be created (meters)
    internode_length: the length of the internode to be created (meters)
    leaf_width: the width of the leaf to be created (meters)
    leaf_length: the length of the leaf to be created (meters)
    bud_length: the length of the bud to be created (meters)
    properties_identical_to_parent: if True, the non-geometrical properties of the phytomers will be derived from those of the parent element
    nil_properties: if True, the non-geometrical properties of the phytomers will be set to 0

    Returns
    -------
    A list containing the three objects of the new phytomer, i.e. the internode, the leaf and the bud
    """


    # 1) We first create a new internode:
    # ------------------------------------
    internode = ADDING_A_CHILD_ORGAN(parent_element,
                                     edge_type='<',
                                     label='Internode',
                                     properties_identical_to_parent=properties_identical_to_parent,
                                     nil_properties=nil_properties)
    # We specified some properties that may be distinct from the parent element:
    internode.radius = internode_radius
    internode.potential_radius = internode_radius
    internode.length = internode_length
    internode.potential_length = internode_length
    internode.angle_down = 0.
    internode.angle_roll = param.phyllotactic_angle
    # We calculate the surface, volume and structural mass of the new organ:
    internode.surface = internode.radius * 2. * pi * internode.length
    internode.volume = internode.radius ** 2. * pi * internode.length
    internode.struct_mass = internode.volume * struct_mass_density
    internode.photosynthetic_area = internode.surface
    # We specify the ramification order :
    internode.ramification_order = parent_element.ramification_order
    # And we reset the time since formation:
    internode.actual_time_since_formation = 0.
    internode.thermal_time_since_formation = 0.

    # 2) We then create a new leaf that is actually a "branch" on the previous internode:
    # ------------------------------------------------------------------------------------
    leaf = ADDING_A_CHILD_ORGAN(parent_element=internode,
                                edge_type='+',
                                label='Leaf',
                                properties_identical_to_parent=properties_identical_to_parent,
                                nil_properties=nil_properties)
    # We specified some properties that may be distinct from the parent element:
    leaf.radius = leaf_width
    leaf.length = leaf_length
    leaf.angle_down = param.leaf_insertion_angle
    leaf.angle_roll = 0.
    # We calculate the surface, assuming the leaf is an ellipse characterized by its length and width:
    leaf.surface = leaf.length / 2. * leaf.radius * pi
    # We calculate the volume, taking into account a minimal thickness:
    leaf.volume = leaf.surface * minimal_leaf_thickness
    leaf.struct_mass = leaf.volume * struct_mass_density
    leaf.photosynthetic_area = leaf.surface
    # We specify the ramification order :
    leaf.ramification_order = parent_element.ramification_order
    # And we reset the time since formation:
    leaf.actual_time_since_formation = 0.
    leaf.thermal_time_since_formation = 0.

    # 3) We finally create an axilliary bud that is actually a "branch" on the previous internode:
    # ---------------------------------------------------------------------------------------------
    bud = ADDING_A_CHILD_ORGAN(parent_element=leaf, # WATCH OUT: here, we set the bud as the child of the leaf, not the internode!
                               edge_type='+',
                               label='Bud',
                               properties_identical_to_parent=False,
                               nil_properties=True)
    # We specified some properties that may be distinct from the parent element:
    bud.radius = param.bud_radius_ratio * parent_element.radius
    bud.potential_radius = bud.radius
    bud.length = bud_length
    bud.potential_length = bud_length
    # The angle_down of the bud is defined according to the angle down of the turtle
    # (we go back on the axis direction, then we apply the insertion angle of the bud):
    bud.angle_down = - param.leaf_insertion_angle + param.bud_insertion_angle
    bud.angle_roll = 0.
    # We assign the surfaces and volumes:
    bud.surface = 0.
    bud.volume = 0.
    bud.struct_mass = 0.
    bud.photosynthetic_area = 0.
    # We specify the ramification order :
    bud.ramification_order = parent_element.ramification_order + 1
    # We set the emergence possibility to False by default:
    bud.emergence_possibility = False
    # And we reset the time since formation:
    bud.actual_time_since_formation = 0.
    bud.thermal_time_since_formation = 0.

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
def new_phytomers_initiation(shoot_system,
                             shoot_MTG, new_initiation = True,
                             initial_internode_length = 0.002,
                             initial_leaf_length = 0.002):
    """
    This function creates new phytomers (internode + leaf + bud) wherever neceassry in the shoot MTG.
    Parameters
    ----------
    shoot_system: the object corresponding to the whole shoot system (not the MTG)
    shoot_MTG: the shoot MTG containing all shoot elements
    new_initiation: a condition for allowing the creation of new phytomers (True by default)
    initial_internode_length: the initial length of the internode of the phytomers to be created (meters)
    initial_leaf_length: the initial length of the leaf of the phytomers to be created (meters)

    Returns
    -------
    [Nothing is returned, the shoot system has been updated with new elements corresponding to the new phytomers]
    """

    # We initialize the current list of apices in the shoot system that may produce new phytomers:
    list_of_apices = [shoot_MTG.node(vid) for vid in shoot_MTG.vertices_iter(scale=1)
                      if shoot_MTG.node(vid).label == "Apex"]

    # We initialize a condition defining whether new phytomers should develop from a given apex:
    apparition_condition = new_initiation

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
                                                           leaf_width=initial_leaf_length * param.leaf_width_to_length_ratio / 2.,
                                                           leaf_length=initial_leaf_length,
                                                           bud_length=0.,
                                                           properties_identical_to_parent=False,
                                                           nil_properties=True
                                                           )
                # We remove the previous apex:
                shoot_MTG.remove_vertex(apex.index(), reparent_child=True)
                # And we add a new apex at the top of the new internode, of length 0:
                new_apex = ADDING_A_CHILD_ORGAN(parent_element=internode,
                                                edge_type="<",
                                                label="Apex",
                                                properties_identical_to_parent = True,
                                                nil_properties=False)
                # print("We have moved the previous apex", apex.index(), "on the node", internode.index(),"; its index is now", new_apex.index())
                new_apex.length = 0.
                # TODO: WATCH OUT BELOW, IN CASE YOU WANT TO CONTINUE ADDING PHYTOMERS ON THE TOP OF EACH OTHERS
                # # Finally, we update the list of apices by adding the new apex that has been formed,
                # # in case it could itself form a new phytomer:
                # list_of_apices.append(new_apex)

                # REGISTERING THE COST OF THE FORMATION FOR THE NEW PHYTOMER:
                # We (re)initialize the structural mass produced:
                new_structural_mass = 0.
                # We add the structural mass of each new element:
                for element in [internode, leaf, bud]:
                    new_structural_mass += element.struct_mass
                # We calculate the corresponding C (gC) and N (gN) cost for the plant:
                C_cost = new_structural_mass * param.structural_C_content / param.yield_growth
                N_cost = new_structural_mass * param.structural_N_content
                # We include these costs in the labile C and N pools in the shoot system:
                shoot_system.shoot_labile_C_pool -= C_cost
                shoot_system.shoot_labile_N_pool -= N_cost

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# This function enables to generate an initial MTG with a certain number of leaves and corresponding internodes:
def create_shoot_MTG(number_of_phytomers=3,
                     initial_leaf_length=0.003,
                     initial_internode_radius=0.001,
                     initial_internode_height = 0.0005,
                     struct_mass_density=1e5, minimal_leaf_thickness=1e-3):
    """
    This function is used to create the shoot part of a plant with a succession of different phytomers on the same branch.

    Parameters
    ----------
    number_of_phytomers: the number of successive phytomers to be created
    initial_leaf_length: the length of the leaves to be created (meters)
    initial_internode_radius: the radius of the internodes to be created (meters)
    initial_internode_height: the height of the internodes to be created (meters)
    struct_mass_density: the density of the structural mass (gDW per cubic meter)
    minimal_leaf_thickness: the thickness of the leaves to be created (meters)

    Returns
    -------
    The new shoot MTG with the appropriate number of phytomers.
    """

    # We initialize an empty MTG:
    g = MTG()

    # We define the base of the shoot system, ADDING ALL THE POSSIBLE PROPERTIES THAT A SHOOT ORGAN WILL LATER CONTAIN:
    id_base = g.add_component(g.root,
                              label='Base_of_the_root_system',
                              # Geometry:
                              radius=initial_internode_radius,
                              potential_radius=initial_internode_radius,
                              length=0.,
                              potential_length=0.,
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

    # We define the ramification order of this first internode as 0:
    new_internode.ramification_order = 0

    # We set the time since formation:
    new_internode.actual_time_since_formation = 0.
    new_internode.thermal_time_since_formation = 0.

    # For each new phytomer to be created:
    for i in range(0, number_of_phytomers):
        # We first create an internode, and then a connected leaf and a bud:
        [new_internode, new_leaf, new_bud] = ADDING_A_PHYTOMER(parent_element=new_internode,
                                                               struct_mass_density=struct_mass_density,
                                                               minimal_leaf_thickness=minimal_leaf_thickness,
                                                               internode_radius=new_internode.radius * 0.9,
                                                               internode_length=initial_internode_height,
                                                               leaf_width=initial_leaf_length * param.leaf_width_to_length_ratio,
                                                               leaf_length=initial_leaf_length,
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
#-----------------------------------------------------------------------------------------------------------------------
def computing_light_interception(shoot_MTG,
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
    This function computes the light intercepted by a plant within a specific population design defined by a scene.

    Parameters
    ----------
    shoot_MTG: the MTG structure of a plant
    incoming_light_density: the density of incoming light above the plant (in micromoles of photons per square meter of land per second)
    isolated_plant: if True, the plant will be considered isolated from any other plants; otherwise, the light will be computed for a population of plants
    zenithal_light: if True, light will be considered coming as only a vertical vector; otherwise, it will come from all directions
    infinite: if True, the scene will be mutliplied on each side to precisely define the intercepted light
    number_of_plants_per_square_meter: the density of plants in the scene (plants per square meter of land)
    inter_row_distance: the distance between two rows (meters)
    scene_length: the length of the scene (meters)
    scene_width: the width of the scene (meters)
    leaf_shape: option to design the shape of the leaves to be displayed ('oval', 'typical_leaf_shape', 'unileaflet', 'drop_shape')

    Returns
    -------
    [Nothing is returned; the shoot MTG has been updated with the properties 'PAR_intercepted' and 'Caribu_area', 'Caribu_Ei', 'Caribu_Eabs']
    """

    # We rename the shoot MTG as 'g':
    g = shoot_MTG

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
            # We get the coordinates of the light vector in the last 3 rows:
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

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function for calculating the actual photosynthesis rate:
def computing_photosynthesis_rate(shoot_MTG, air_temperature_in_degree_Celsius=10.):
    """
    This function calculates the actual photosythesis rate (in gC per second) for each organ of the shoot MTG, based on
    its intercepted PAR (micromol of photons per m2 per second), its photosynthetic surface (m2)
    and a temperature-dependent parameter SCA, i.e. its specific carbon assimilation (gC per micromol of photons).

    Parameters
    ----------
    shoot_MTG: the shoot MTG containing each organ as node
    air_temperature_in_degree_Celsius: the average temperature of the air experienced by the shoots

    Returns
    -------
    [Nothing is returned; the MTG has been updated with the property 'photosynthesized_C_rate' (gC per second)]
    """

    # We correct the specific carbon assimilation, based on air temperature:
    corrected_SCA = temperature_modification(temperature_in_Celsius=air_temperature_in_degree_Celsius,
                                             process_at_T_ref=param.SCA,
                                             T_ref=param.SCA_T_ref,
                                             A=param.SCA_A,
                                             B=param.SCA_B,
                                             C=param.SCA_C)

    # We initialize the total intercepted PAR (in micromoles of photons per second per plant):
    total_PAR_intercepted = 0.

    # We cover each organ n of the shoot:
    for vid in shoot_MTG.vertices_iter(scale=1):
        n = shoot_MTG.node(vid)
        n.photosynthetic_area = n.surface
        # We calculate the photosynthesis rate for this organ:
        n.photosynthesized_C_rate = corrected_SCA * n.PAR_intercepted * n.photosynthetic_area
        # TODO: have the SCA vary with the level of "green" within the leaf? with the level of N? with the senescence status?
        # We increment the total PAR intercepted by the plant:
        total_PAR_intercepted += n.PAR_intercepted * n.photosynthetic_area

    return total_PAR_intercepted # in micromoles of photons per second per plant
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function for computing the C cost for maintenance respiration (as opposed to growth respiration):
def computing_maintenance_rate(shoot_MTG, air_temperature_in_degree_Celsius=10.):
    """
    This function calculates the respiration rate (in gC per second) for the maintenance of each organ of the shoot MTG.

    Parameters
    ----------
    shoot_MTG: the shoot MTG containing each organ as node
    air_temperature_in_degree_Celsius: the average temperature of the air experienced by the shoots

    Returns
    -------
    [Nothing is returned; the MTG has been updated with the property 'maintenance_respiration_rate' (gC per second)]
    """

    # We compute a modifying factor based on temperature:
    temperature_modifier = temperature_modification(temperature_in_Celsius=air_temperature_in_degree_Celsius,
                                                    process_at_T_ref=1,
                                                    T_ref=param.resp_maintenance_max_T_ref,
                                                    A=param.resp_maintenance_max_A,
                                                    B=param.resp_maintenance_max_B,
                                                    C=param.resp_maintenance_max_C)

    # We convert the maximal rate (originally expressed in mol of C per gDW per second, now in gC per gDW per second),
    # and we adjust it with temperature:
    corrected_Vmax = param.resp_maintenance_max * 12.01 * temperature_modifier

    # We also adjust the unit of the Km inherited from RhizoDep (originally expressed in moles of equivalent-hexose):
    Km_maintenance = param.Km_maintenance * 6 * 12.01  # now in gC per gram of structural mass per second

    # We cover each organ n of the shoot:
    for vid in shoot_MTG.vertices_iter(scale=1):
        n = shoot_MTG.node(vid)

        # If the element does not have any C left, it cannot respire!
        if n.C_labile <= 0.:
            n.shoot_maintenance_respiration_rate = 0.
            continue

        # The maintenance rate is based on a maximal rate of respiration per gram of structural mass, modulated by the
        # amount of labile C in the organ:
        n.shoot_maintenance_respiration_rate = corrected_Vmax * n.C_labile / (Km_maintenance + n.C_labile) * n.struct_mass

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function for redistributing the C among shoot organs, with a C allocation. Different possibilities can be considered:
# - proportionnal to PARi for each organ compared to the total PAR intercepted by the plant
# - proportionnal to PARi/surface for each organ compared to the sum of PARi/surface of all organs
# - proportionnal to PARi/surface of the leaf of the phytomer for each organ compared to the sum of PARi/surface of all leaves
def computing_shoot_C_redistribution(shoot_MTG,
                                     shoot_available_C = None,
                                     intercepted_light = None,
                                     time_step_in_seconds=60*60,
                                     option="classic_PAR"):

    """

    Parameters
    ----------
    shoot_MTG
    shoot_available_C
    intercepted_light
    time_step_in_seconds
    option: "classic_PAR" (micromol of photons per second),
            "surfacic_PAR" (micromol of photons per m2 of organ per second)

    Returns
    -------

    """

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

    # We initialize the sum of surfacic PAR, if needed:
    if option == "surfacic_PAR":
        sum_of_surfacic_PAR = 0.

    # If either shoot_available_C or intercepted_light have not been provided:
    if shoot_available_C is None or intercepted_light is None:
        # Then we cover each organ n of the shoot:
        for vid in shoot_MTG.vertices_iter(scale=1):
            n = shoot_MTG.node(vid)
            # We add the initial amount of C_labile and the net gain in C of the current organ:
            # TODO: WATCH OUT - Here, we already add the C from photosynthetis and withdraw the C for maintenance to calculate the remaining C!
            available_C += n.C_labile + (n.photosynthesized_C_rate - n.shoot_maintenance_respiration_rate) * time_step_in_seconds
            # We add the light (or light density) intercepted by the organ:
            if option == "classic_PAR":
                intercepted_light_by_the_plant += n.PAR_intercepted
            elif option == "surfacic_PAR" and n.surface > 0.:
                sum_of_surfacic_PAR += n.PAR_intercepted / n.surface

    # 2) We redistribute the C according to specific rules:
    #------------------------------------------------------

    # We cover each organ n of the shoot:
    for vid in shoot_MTG.vertices_iter(scale=1):
        n = shoot_MTG.node(vid)
        # The corresponding rate of C redistribution for this organ is calculated:
        # (NOTE: if positive, it corresponds to a C gain; if negative, it corresponds to a C source)
        if option == "classic_PAR":
            n.shoot_C_redistribution_rate = (available_C * n.PAR_intercepted / intercepted_light_by_the_plant
                                             / time_step_in_seconds)
        elif option == "surfacic_PAR" and n.surface > 0. and sum_of_surfacic_PAR > 0.:
            n.shoot_C_redistribution_rate = (available_C * n.PAR_intercepted / n.surface / sum_of_surfacic_PAR
                                             / time_step_in_seconds)
        else:
            n.shoot_C_redistribution_rate = 0.

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function for distributing the N coming from the roots into each shoot organ:
def computing_shoot_N_distribution(shoot_MTG, shoot_available_N=None, time_step_in_seconds=60.*60.*24.):
    """
    This function distributes the whole amount of labile N available to the whole plant to each organ which has a
    positive surface, proportionnaly to the mean surfacic PAR intercepted by the organ compared to the sum of all
    surfacic PARs of shoot elements.

    :param shoot_MTG: the MTG corresponding to the shoot system
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
    for vid in shoot_MTG.vertices_iter(scale=1):
        n = shoot_MTG.node(vid)
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
    for vid in shoot_MTG.vertices_iter(scale=1):
        # We access the current organ"n":
        n = shoot_MTG.node(vid)
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
def computing_shoot_N_remobilization(shoot_MTG, k_N_mobilization = 1/(24.*60.*60.)):
    """
    This function calculates the rate of N mobilization (gN per second) in each shoot organ, based on first-order kinetics.
    :param g: the MTG corresponding to the shoot system
    :param k_N_mobilization: the time constant of the N mobilization (second-1)
    :return: [no return - the MTG g has been updated with the property 'N_mobilization_rate']
    """

    # We cover each organ of the shoot:
    for vid in shoot_MTG.vertices_iter(scale=1):
        # We access the current organ"n":
        n = shoot_MTG.node(vid)

        # We exclude elements where a deficit of labile N has been previously recorded:
        if n.N_labile_deficit >= 0.:
            n.shoot_N_mobilization_rate = 0.
            continue

        # We compute the rate N mobilization:
        n.shoot_N_mobilization_rate = k_N_mobilization * n.N_labile

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
def computing_shoot_N_immobilization(shoot_MTG):
    """
    This function calculates the rate of N immobilization (gN per second) in each shoot organ.
    :param shoot_MTG: the MTG corresponding to the shoot system
    :return: [no return - the MTG g has been updated with the property 'net_shoot_N_immobilization_rate']
    """

    # We cover each organ of the shoot:
    for vid in shoot_MTG.vertices_iter(scale=1):
        # We access the current organ"n":
        n = shoot_MTG.node(vid)
        # We compute the rate N mobilization:
        n.net_shoot_N_immobilization_rate = 0.

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
def computing_shoot_C_storage(shoot_system, shoot_MTG, time_step_in_seconds=60.*60.*24.):
    """
    This function calculates the rate of labile C reallocation to each leaf.
    :param shoot_system: the object corresponding to the whole shoot system (not the MTG!)
    :param shoot_MTG: the MTG corresponding to the shoot system
    :param time_step_in_seconds: the time step of the simulation (seconds)
    :return: [no return - the MTG g has been updated with the property 'C_storage_rate' (gC per second)]
    #TODO: Add the increase in thickness corresponding to the additional C storage.
    """

    # We rename the shoot MTG for simplification:
    g=shoot_MTG

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
    # intercepted_light_by_the_plant = plant.total_PAR_intercepted_by_the_plant

    # 2) Computing the rate of C distribution in each shoot organ
    # -----------------------------------------------------------

    # We access the value of the total amount of available C, stored in the base of the shoot system:
    plant_available_C = shoot_system.shoot_labile_C_pool

    # We cover each organ of the shoot:
    for vid in g.vertices_iter(scale=1):
        # We access the current organ"n":
        n = g.node(vid)
        # The new C storage rate in this organ is calculated according to the relative intercepted light:
        n.shoot_C_storage_rate = plant_available_C * n.PAR_intercepted / intercepted_light_by_the_plant / time_step_in_seconds

    # We reset the value of the plant available C:
    shoot_system.shoot_labile_C_pool = 0.

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
def potential_shoot_growth(shoot_MTG,
                           SFQN=1,
                           LA_max=10,
                           leaf_minimal_thickness = 1e-3,
                           elong_internode=5e-7,
                           plant_age=0.,
                           time_step_in_seconds = 60.*60.*24.,
                           air_temperature_in_degree_Celsius = 10.):
    """
    This function computes the potential increase in volume of each organ in the shoot system.
    :param shoot_MTG: the MTG corresponding to the shoot system
    :param SFQN: the slope of the relationship of Leaf area = f(N_labile) (m2 per gN)
    :param LA_max: the maximal area for a leaf (m2)
    #TODO: should the leaf stop growing after a predetermined duration, instead of after reaching the LA_max?
    :param leaf_minimal_thickness: the minimal thickness of a leaf (meter)
    :param elong_internode: a potential rate of elongation (meter per second)
    :param time_step_in_seconds: the time step of the simulation (seconds)
    :param air_temperature_in_degree_Celsius: the mean temperature experienced by the shoots (degree Celsius)
    :return:
    """

    # We rename the MTG for simplification:
    g=shoot_MTG

    # We consider a temperature-modifier of time, in order to calculate the equivalent of thermal time:
    temperature_adjustement = temperature_modification(temperature_in_Celsius=air_temperature_in_degree_Celsius,
                                                                 process_at_T_ref=1,
                                                                 T_ref=param.shoot_growth_T_ref,
                                                                 A=param.shoot_growth_A,
                                                                 B=param.shoot_growth_B,
                                                                 C=param.shoot_growth_C)
    # We calculate the equivalent time step in thermal time:
    thermal_time_step_in_seconds = time_step_in_seconds * temperature_adjustement
    thermal_time_step_in_days = thermal_time_step_in_seconds/(60.*60.*24.)

    # To cover the MTG in a specific order, we need to define its "root":
    root_gen = g.component_roots_at_scale_iter(g.root, scale=1)
    root = next(root_gen)
    # We cover each organ n of the shoot, from the tip of axes up to the base:
    for vid in post_order(g, root):
        # We define the current organ as "n":
        n = g.node(vid)

        # CASE 1: the shoot element is a leaf
        #------------------------------------
        if n.label == "Leaf":
            # We calculate the potential surface based on the available N in the leaf:
            potential_leaf_surface = SFQN * n.N_labile
            # We calculate the potential surface increase (m2), knowing that the surface can't become higher than LA_max:
            potential_surface_increase = min(potential_leaf_surface - n.surface, LA_max - n.surface)
            # We calculate the maximal surface increase (m2), knowing that the surface increase cannot become higher
            # than the maximal increase rate allowed:
            maximal_surface_increase = n.surface * param.LA_max_relative_increase_per_day * thermal_time_step_in_days
            # The potential volume increment is based on a minimal thickness:
            n.potential_volume_increment = min(potential_surface_increase, maximal_surface_increase) * leaf_minimal_thickness
            # TODO: what happens if the leaf has already a thickness higher than the minimal one?

        # CASE 2: the shoot element is an internode
        #------------------------------------------
        elif n.label == "Internode":

            # 1) We calculate a potential elongation rate.

            # We check whether the lateral bud has emerged or is about to emerge.
            try:
                # We initialize a Boolean that will determine whether elongation is fast or not:
                fast_elongation = False
                # We try to get the ID of the first element of the ramification, i.e. the lateral axis of the leaf:
                leaf_vid = g.Sons(vid, EdgeType='+')[0]
                lateral_vid = g.Sons(leaf_vid, EdgeType='+')[0]
                bud = g.node(lateral_vid)
                # We then examine whether the bud has emerged or is about to emerge:
                if bud.emergence_possibility is True:
                    fast_elongation = True
            except:
                fast_elongation = False

            # If the plant has reached a minimal age above which it can fully elongate:
            if plant_age > param.minimal_age_for_internode_elongation and fast_elongation:
                # Then we set the potential elongation rate to the full potential:
                potential_elongation_rate = elong_internode * temperature_adjustement
            else:
                # Otherwise, we only allow only a very slow elongation rate:
                potential_elongation_rate = elong_internode * temperature_adjustement / 100.
            # Eventually, we record the potential new length:
            n.potential_length = n.length + potential_elongation_rate * time_step_in_seconds

            # 2) We calculate the potential radial increase.
            # We first need to consider all children of the current internode in order to verifiy whether the radius
            # of the internode should increase according to the pipe model rule.
            # We initialize the sections of the children:
            son_section=0.
            sum_of_lateral_sections=0.
            # We cover each child of the internode (either the neighbouring internode, or a lateral leaf):
            for child in n.children():
                # If the child belongs to the same axis, it is an internode:
                if child.properties()['edge_type'] == '<':
                    # Then we record the POTENTIAL section of this child:
                    son_section = child.potential_radius ** 2 * pi
                # Otherwise if the child belongs to a lateral axis, it is a leaf:
                elif child.properties()['edge_type'] == '+':
                    # The we go the the child of this leaf, which corresponds to an unemerged bud or the internode of a ramification:
                    grand_child = child.children()[0]
                    # If its length is positive, it is a lateral ramification:
                    if grand_child.length  > 0.:
                        # We add the POTENTIAL section of this lateral ramification to a sum of lateral sections:
                        sum_of_lateral_sections += grand_child.potential_radius ** 2 * pi
            # Eventually, we calculate the theoretical new radius of the initial internode according to its children:
            theoretical_radius = sqrt(son_section / pi + param.SGC_shoot * sum_of_lateral_sections / pi)
            # If the net difference is below 0.1% of the initial radius:
            if (theoretical_radius - n.radius) <= 0.001 * n.radius:
                # Then we reset the theoretical radius value to the initial value to avoid computational problems:
                theoretical_radius = n.radius
            # We also estimate the maximal thickening rate (meter per second) according to air temperature:
            thickening_rate = param.max_internode_thickening * n.radius * temperature_adjustement
            # The maximal possible new radius according to this regulation is therefore:
            new_radius_max = (1 + thickening_rate * time_step_in_seconds) * n.radius
            # If the potential new radius is higher than the maximal new radius:
            if theoretical_radius > new_radius_max:
                # Then potential thickening is limited up to the maximal new radius:
                n.potential_radius = new_radius_max
            # Otherwise, the potential radius to achieve is equal to the theoretical one:
            else:
                n.potential_radius = theoretical_radius

            #  Finally, we compute the potential volume increment:
            potential_volume = n.potential_length * (n.potential_radius ** 2 * pi)
            n.potential_volume_increment = potential_volume - n.length * (n.radius ** 2 * pi)

        # CASE 3: the shoot element is anything else
        #-------------------------------------------
        else:
            # Then we set the potential volume increment to 0:
            n.potential_volume_increment = 0.

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
def computing_C_redistribution_for_shoot_growth(shoot_system,
                                                shoot_MTG,
                                                yield_growth=0.8,
                                                structural_mass_density=1e5,
                                                structural_C_content=0.44,
                                                time_step_in_seconds = 60.*60.*24.,
                                                option="classic_PAR"
                                                ):

    # We rename the shoot MTG for simplification:
    g=shoot_MTG

    # We initialize the stock of C that can be redistributed to sink leaves and the one that will be eventually used:
    total_available_C = 0.
    C_used_by_sink_organs = 0.
    # We initialize the total intercepted light:
    if option == "classic_PAR":
        intercepted_light = 0.
    elif option == "surfacic_PAR":
        sum_of_surfacic_PAR = 0.
    # We initialize a list that will contain the indices of the leaves that will need additional C:
    list_of_sink_organs = []

    # 1) DEFINING THE AMOUNT OF AVAILABLE C:
    ########################################

    # We cover each organ n of the shoot:
    for vid in g.vertices_iter(scale=1):
        # We define the current organ as "n":
        n = g.node(vid)

        # We add the light intercepted by the organ:
        if option == "classic_PAR":
            intercepted_light += n.PAR_intercepted
        elif option == "surfacic_PAR" and n.surface > 0.:
            sum_of_surfacic_PAR += n.PAR_intercepted / n.surface

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
            list_of_sink_organs.append(vid)

    # 2) RE-ALLOCATING C TO SINK ORGANS:
    ###################################

    # We cover each organ n of the shoot:
    for vid in list_of_sink_organs:
        # We define the current organ as "n":
        n = g.node(vid)
        # The new amount of C_labile in this organ is calculated, after its initialization:
        additional_C_for_sink_organ = 0.
        if option == "classic_PAR":
            additional_C_for_sink_organ = total_available_C * n.PAR_intercepted / intercepted_light
        elif option == "surfacic_PAR" and n.surface > 0. and sum_of_surfacic_PAR > 0.:
            additional_C_for_sink_organ = total_available_C * n.PAR_intercepted / n.surface / sum_of_surfacic_PAR
        # TODO: Note that in either case, only the fraction of C corresponding to (light_intercepted_by_sink_organs)
        #  / (total_light_intercepted_by all_organs) is used...
        # We adjust the properties of the organ:
        # n.C_labile += additional_C_for_sink_organ
        n.C_input_rate = additional_C_for_sink_organ / time_step_in_seconds
        C_used_by_sink_organs += additional_C_for_sink_organ

    # Eventually, we store the remaining available C:
    shoot_system.shoot_labile_C_pool = total_available_C - C_used_by_sink_organs # The remaining amount of C here may be allocated to roots!

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
def actual_shoot_growth(shoot_MTG,
                        yield_growth = 0.8,
                        structural_mass_density = 1e5,
                        structural_C_content = 0.44,
                        structural_N_content = 0.02,
                        leaf_minimal_thickness = 1e-3,
                        time_step_in_seconds = 60.*60.*24.,
                        air_temperature_in_degree_Celsius = 10.):
    """
    This function computes the actual growth of the shoot system, depending on the C and N available in each organ,
    as well as the ageing of each shoot element.
    :param shoot_MTG: the shoot system
    :param yield_growth: the fraction of the total C used for growth that will be incorporated to the new mass (the rest is respired) (gC per gC)
    :param struct_mass_density: the density of the structural mass (gram of structural maqss per m3)
    :param structural_C_content: the C content of the structural mass (gC per gram of structural mass)
    :param structural_N_content: the N content of the structural mass (gN per gram of structural mass)
    :param leaf_minimal_thickness: the minimal thickness of a leaf (m)
    :param time_step_in_seconds: the time step (seconds)
    :param air_temperature_in_degree_Celsius: the mean temperature experienced by the shoots (degree Celsius)
    :return:
    """

    # We rename the shoot MTG for simplification:
    g=shoot_MTG

    # We consider a temperature-modifier of time, in order to calculate the equivalent of thermal time:
    temperature_adjustment = temperature_modification(temperature_in_Celsius=air_temperature_in_degree_Celsius,
                                                                 process_at_T_ref=1,
                                                                 T_ref=param.shoot_growth_T_ref,
                                                                 A=param.shoot_growth_A,
                                                                 B=param.shoot_growth_B,
                                                                 C=param.shoot_growth_C)
    # We calculate the equivalent time step in thermal time:
    thermal_time_step_in_seconds = time_step_in_seconds * temperature_adjustment

    # We cover each organ n of the shoot:
    for vid in g.vertices_iter(scale=1):

        # We define the current organ as "n":
        n = g.node(vid)

        # First, we increment the time since the shoot element has been formed:
        n.actual_time_since_formation +=  time_step_in_seconds
        n.thermal_time_since_formation += thermal_time_step_in_seconds

        # If the element does not have any potential growth, the we move to the next element:
        if n.potential_volume_increment <=0.:
            continue
        # If the element does not have any C or N left, it cannot grow:
        if n.C_labile <=0. or n.N_labile <=0.:
            continue

        # COMPUTING THE MAXIMAL VOLUME INCREMENT:
        #########################################

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

        # COMPUTING THE FINAL DIMENSIONS OF THE SHOOT ELEMENT
        #####################################################

        # CASE 1: The shoot element is an internode
        #------------------------------------------
        if n.label == "Internode":
            # We compute the maximal length increment and radius increment based on the final volume increment:
            possible_length_increment = final_volume_increment / (n.radius ** 2 * pi)
            possible_radius_increment = sqrt(n.radius ** 2 + (final_volume_increment / (n.length * pi))) - n.radius
            # If there is no expected radial increase:
            if n.potential_radius == n.radius:
                # Then we simply elongated the internode according to the final volume increment:
                n.length += possible_length_increment
            # If there is no expected elongation:
            elif n.potential_length == n.length:
                # Then we simply increase the radius according to the final volume increment:
                n.radius += possible_radius_increment
            # Otherwise we can increase both the length and the radius of the internode:
            else:
                # If the potential length previously computed cannot be created in totality:
                if n.potential_length > n.length + possible_length_increment:
                    # Then we give priority to elongation and elongate the internode according to the final volume increment:
                    n.length += possible_length_increment
                    # The radius is however not increased.
                # Otherwise, we can perform the elongation up to the full potential:
                else:
                    # First, we calculate the volume increment corresponding to this elongation:
                    volume_increment_by_elongation = (n.radius ** 2 * pi) * (n.potential_length - n.length)
                    # We increase the actual length of the internode up to its potential length:
                    n.length =  n.potential_length
                    # And we increase the radius according to the remaining "volume_increment" that is left:
                    remaining_volume_increment = final_volume_increment - volume_increment_by_elongation
                    n.radius = sqrt(n.radius ** 2 + (remaining_volume_increment / (n.length * pi)))

            # Eventually, we recompute the external surface of the internode:
            n.surface = n.length * n.radius * 2 * pi

        # CASE 2: The shoot element is a leaf
        #------------------------------------
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
def ageing(shoot_MTG, time_step_in_seconds=24 * 60. * 60., air_temperature_in_degree_Celsius=20):
    """
    This functions increases the age of each organ in the shoot MTG. The age is expressed in seconds, either as actual
    time or as time equivalent to the reference growth temperature.

    Parameters
    ----------
    shoot_MTG: the shoot MTG containing each organ as node
    time_step_in_seconds: the actual time step of the simulation (in seconds)
    air_temperature_in_degree_Celsius: the average temperature of the air experienced by the shoots (in degree Celsius)

    Returns
    -------
    [Nothing is returned; the MTG has been updated with the properties 'actual_time_since_formation' and 'thermal_time_since_formation' (seconds)]
    """

    # We consider a temperature-modifier of time, in order to calculate the equivalent of thermal time:
    temperature_adjustment = temperature_modification(
        temperature_in_Celsius=air_temperature_in_degree_Celsius,
        process_at_T_ref=1,
        T_ref=param.shoot_growth_T_ref,
        A=param.shoot_growth_A,
        B=param.shoot_growth_B,
        C=param.shoot_growth_C)

    for vid in g.vertices_iter(scale=1):
        element = g.node(vid)
        element.actual_time_since_formation += time_step_in_seconds
        element.thermal_time_since_formation += time_step_in_seconds * temperature_adjustment

    return
# -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------
def computing_ramification(shoot_MTG, initial_internode_length=0.0005, initial_leaf_length=0.001):
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

    # ***

    # La floraison est indéterminée. A partir du virage floral, on a élongation tige principale, puis apparition
    # inflorescence à la fin des axes + élongation des ramifications.
    # Idée de l’archi potentielle : x feuilles au virage floral conditionnent un nombre total de feuilles n
    # (sachant que les entrenoeuds courts ne donnent pas de ramifications).. A partir du premier entrenoeud qui ramifie,
    # tous les entrenoeuds au-dessus vont former des ramifications. Chaque axe latéral va avoir un nombre max d’entrenoeud :
    # 1 au sommet, puis 2 sur l’axe latéral en dessous, puis 3, etc. Pour les secondaires, en théorie on aurait la même chose.
    # Comment est régulé le numéro de la ramification la plus basse ? Les axes ramifiés démarrent à peu près tous en même temps,
    # donc on les fait partir en même temps, à partir du moment où une part. Le RC/RS diminue en partant du haut dans le couvert.
    # Il y aurait donc un seuil critique qui déclenche la ramif à partir de cet entrenoeud  ça détermine donc
    # la position de l’entrenoeud de départ.
    # Quand déclencher ? Avec un temps thermique donné de départ pour la plante, ou pour l’entrenoeud de départ.
    # Ou photopériode ? Ou température critique ? Faut-il tenir compte aussi de l’ensoleillement global ?
    # Pourquoi les ramif du haut s’allongeraient plus lentement ? Dominance apicale à coder en dur par rapport à la distance
    # à l’inflorescence de l’axe principal ?

    # We rename the shoot_MTG for simplification:
    g = shoot_MTG

    # DEFINING THE CONDITIONS OF TURNING DORMANT BUDS INTO ACTIVE BUDS:
    for vid in g.vertices_iter(scale=1):
        if g.node(vid).label != "Bud":
            continue
        else:
            bud = g.node(vid)
            # We get the leaf below the bud (as the bud has been defined as the child of the leaf):
            leaf_index = g.Father(vid, EdgeType='+')
            leaf = g.node(leaf_index)
            # We calculate the density of light intercepted by this leaf (micromoles of photons per second per square meter of leaf):
            light_interception = leaf.PAR_intercepted / leaf.surface
            # We define the conditions of bud emergence according to its thermal age and the PAR intercepted by the leaf nearby:
            if (light_interception >= param.bud_light_threshold
                    and bud.thermal_time_since_formation >= param.bud_minimal_dormancy_period
                    and bud.ramification_order <= 1):
                bud.emergence_possibility = True

    # TRANSFORMING THE ACTIVE BUDS INTO A PHYTOMER AND PHYTOMER:
    # We select the vertices corresponding to the buds that are allowed to emerge:
    list_of_apices = [g.node(vid) for vid in g.vertices_iter(scale=1) if g.node(vid).label == "Bud"
                      and g.node(vid).emergence_possibility == True]
    # As long as the list of emerging buds to treat is not empty:
    while list_of_apices != []:
        # We cover each bud present in the list:
        for apex in list_of_apices:
            print("For the bud", apex.index(),"the initial angle down is", apex.angle_down)
            # We add a new phytomer at the top of the current apex:
            [internode, leaf, bud] = ADDING_A_PHYTOMER(apex,
                                                       internode_radius=apex.radius,
                                                       internode_length=initial_internode_length,
                                                       leaf_width=initial_leaf_length * param.leaf_width_to_length_ratio / 2.,
                                                       leaf_length=initial_leaf_length,
                                                       bud_length=0.,
                                                       properties_identical_to_parent=False,
                                                       nil_properties=True
                                                       )
            # We correct the angle down of the new internode, so that its direction corresponds to the one from the initial bud:
            internode.angle_down = apex.angle_down
            print("For the new internode created from the bud", apex.index(), "the angle down is", internode.angle_down)
            # By default, we set the bud as dormant:
            bud.emergence_possibility = False
            # We remove the previous apex, as a new one will now be formed on the top of the new phytomer:
            g.remove_vertex(apex.index(), reparent_child=True)
            # And we add a new apex at the top of the new internode, of length 0:
            new_apex = ADDING_A_CHILD_ORGAN(parent_element=internode,
                                            edge_type="<",
                                            label="Apex",
                                            properties_identical_to_parent=True,
                                            nil_properties=False)
            new_apex.length = 0.
            new_apex.angle_down = 0.
            # We remove the former apex from the list, so that it is not considered next time:
            list_of_apices.remove(apex)
            # TODO: WATCH OUT BELOW, IN CASE YOU WANT TO CONTINUE ADDING PHYTOMERS
            # # Finally, we update the list of apices by adding the new apex that has been formed,
            # # in case it could itself form a new phytomer:
            # list_of_apices.append(new_apex)

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
def computing_senescence(shoot_MTG,
                         time_step_in_seconds=1/(60.*60.*24.),
                         air_temperature_in_degree_Celsius=10.):

    # We consider a temperature-modifier of time, in order to calculate the equivalent of thermal time:
    temperature_adjustement = temperature_modification(temperature_in_Celsius=air_temperature_in_degree_Celsius,
                                                                 process_at_T_ref=1,
                                                                 T_ref=param.shoot_growth_T_ref,
                                                                 A=param.shoot_growth_A,
                                                                 B=param.shoot_growth_B,
                                                                 C=param.shoot_growth_C)
    # We calculate the equivalent time step in thermal time:
    thermal_time_step_in_seconds = time_step_in_seconds * temperature_adjustement
    thermal_time_step_in_days = thermal_time_step_in_seconds/(60.*60.*24.)

    # We cover each organ of the shoot:
    for vid in shoot_MTG.vertices_iter(scale=1):

        # We define the current element as n:
        n = shoot_MTG.node(vid)
        # We make sure that the current element has a positive length and is not an axilliary bud:
        if n.length<=0. or n.label=="Bud":
            continue
        # We verify whether the current element shoud start senescing:
        if n.senescence_status == False:
            # If there is no labile C to sustain the element and if it has reached a maximal "youth" age:
            if (n.C_labile_deficit > 0.
                    and (n.thermal_time_since_formation + thermal_time_step_in_seconds) >= param.thermal_youth_duration):
                # We declare that the element should enter senescence:
                n.senescence_status = True

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
def computing_leaf_fall(shoot_MTG,
                         time_step_in_seconds=1/(60.*60.*24.),
                         air_temperature_in_degree_Celsius=10.):

    # We rename the shoot MTG for simplicity:
    g=shoot_MTG

    # We get the list of all living leaves in the shoot MTG:
    list_of_leaves = [g.node(vid) for vid in g.vertices_iter(scale=1) if g.node(vid).label == "Leaf"]

    # We cover each leaf:
    for leaf in list_of_leaves:
        # If the leaf has to fall:
        if (leaf.senescence_status == True
                and (leaf.C_labile / leaf.struct_mass) <= param.C_labile_content_treshold_for_leaf_fall
                and (leaf.N_labile / leaf.struct_mass) <= param.N_labile_content_treshold_for_leaf_fall):
            leaf.label="Fallen_leaf"

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
def shoot_organ_CN_balance(shoot_MTG, time_step_in_seconds = 60.*60.*24., resetting_rates=False):

    # We rename the MTG for simplification:
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
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
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
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
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

def calling_the_different_shoot_modules(shoot_MTG,
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

    # We rename the MTG for simplification:
    g = shoot_MTG

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
    computing_photosynthesis_rate(g)

    # We calculate the C respiration cost associated to maintenance:
    computing_maintenance_rate(g)

    # We allocate the N available from the base of the shoot up to each leaf, based on PAR:
    computing_shoot_N_distribution(g,
                                   time_step_in_seconds=time_step_in_days*60.*60.*24.)

    # TODO: Is this step necessary or not?
    # This calculates the total amount of available C, and distributes it in each organ according to PAR:
    computing_shoot_C_redistribution(g)

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
    computing_C_redistribution_for_shoot_growth(shoot_system=plant,
                                                shoot_MTG=g,
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
    new_phytomers_initiation(shoot_system = plant,
                             shoot_MTG = g,
                             new_initiation=phytomers_initiation,
                             initial_internode_length=param.initial_internode_length,
                             initial_leaf_length=param.initial_leaf_length)

    # We update all dimensions of the organs:
    updating_coordinates(g)

    # 3) Mobilization/immobilization:
    #--------------------------------

    # This function reallocates all the remaining C available to the different organs according to their PAR,
    # and stores it as additional structural mass:
    computing_shoot_C_storage(shoot_system=plant, shoot_MTG=g, time_step_in_seconds=time_step_in_days*60.*60.*24.)

    # This function mobilizes N from the organs and stores it at the base of the shoot system:
    computing_shoot_N_remobilization(g, k_N_mobilization=param.k_N_mobilization)

    # This function immobilizes N in each shoot organ:
    computing_shoot_N_immobilization(g)

    # 4) C&N mass balance:
    #---------------------
    # This function performs the actual C & N balance and updates the labile pools of C & N in each organ:
    shoot_organ_CN_balance(g)

    # 5) Ageing of the plant:
    # -----------------------
    ageing(g, time_step_in_seconds=time_step_in_days*24*60.*60., air_temperature_in_degree_Celsius=10.)

    # 6) Ramification:
    #-----------------
    computing_ramification(g,
                           initial_internode_length=param.initial_internode_length,
                           initial_leaf_length=param.initial_leaf_length)

    # 7) Senescence and leaf fall:
    #-----------------------------
    computing_senescence(g)
    computing_leaf_fall(g)

    return
#-----------------------------------------------------------------------------------------------------------------------
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
    g = create_shoot_MTG(number_of_phytomers=3,
                         initial_internode_height=0.002,
                         initial_internode_radius=0.0008,
                         initial_leaf_length=0.02)

    # By default, we consider the object "shoot_system" representing the whole shoots as the first node of the MTG:
    plant = g.node(1)
    # We initialize the labile content of C and N:
    plant.shoot_labile_C_pool = 1e-2
    plant.shoot_labile_N_pool = 1e-3

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

    # PREPARING THE COORDINATES OF CAMERA ROTATION:
    # We calculate the x,y,z coordinates of the different positions of the camera:
    [x_coordinates, y_coordinates, z_coordinates] = circle_coordinates(x_center=param.x_center,
                                                                       y_center=param.y_center,
                                                                       z_center=param.z_camera,
                                                                       radius=sqrt(
                                                                           (param.x_camera ** 2 + param.y_camera ** 2)),
                                                                       n_points=param.camera_rotation_n_points)
    # We initialize the index in the list of the coordinates:
    index_camera = 0

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
        # We artificially replenish the pool of labile N in the shoots:
        plant.shoot_labile_N_pool = 0.1

        # CALLING SUCCESSIVELY THE DIFFERENT SHOOT MODULES:

        # We first examine whether new phytomers will be produced at this time step.
        # We reset the authorization to add a new phytomer:
        phytomer_formation = False
        # We increment the time elapsed since the last phytomer formation with the current time step:
        elapsed_time_in_days_since_last_leaf_emergence += time_step_in_days
        # TODO: WATCH OUT: this does not take into account the actual temperature!
        # If the elapsed time since last phytomer formation is higher than the prescribed phyllochrone:
        if elapsed_time_in_days_since_last_leaf_emergence >= param.phyllochron / (60.*60.*24.):
            # We authorize the formation of a new phytomer:
            phytomer_formation = True
            # And we reset the elapsed time:
            elapsed_time_in_days_since_last_leaf_emergence = (elapsed_time_in_days_since_last_leaf_emergence
                                                              - (param.phyllochron / (60.*60.*24.)))

        # We call the different modules:
        calling_the_different_shoot_modules(shoot_MTG = g,
                                            phytomers_initiation=phytomer_formation,
                                            number_of_plants_per_square_meter=plant_density_in_square_meter,
                                            inter_row_distance=inter_row_distance,
                                            scene_length=row_length,
                                            scene_width=scene_width)

        # PLOTTING ONE SINGLE PLANT OVER TIME:
        # We display the MTG in its new state:
        # Option 1A - With PlantGL:
        #--------------------------
        x_cam = x_coordinates[index_camera]
        y_cam = y_coordinates[index_camera]
        z_cam = z_coordinates[index_camera]
        # We define the index of the coordinates to read at the next step:
        index_camera = index_camera + 1
        # If this index is higher than the number of coordinates in each vector:
        if index_camera >= param.camera_rotation_n_points:
            # Then we reset the index to 0:
            index_camera = 0
        graph = plotting_shoot_MTG(g,
                                   leaf_shape="oval",
                                   # leaf_shape="typical_leaf_shape",
                                   # leaf_shape="drop_shape",
                                   width=1200, height=900,
                                   x_center=0, y_center=0, z_center=0.03,
                                   # x_cam=0.18, y_cam=0.18, z_cam=0.03,
                                   x_cam=x_cam, y_cam=y_cam, z_cam=z_cam,
                                   grid=False,
                                   # background_color=[0, 0, 0],
                                   single_color=[100,150,50],
                                   # property_name='C_labile', cmap='jet', vmin=1e-5, vmax=1e1, lognorm=True,
                                   property_name='PAR_intercepted', cmap='jet', vmin=0, vmax=600, lognorm=False,
                                   # property_name='Caribu_Eabs', cmap='jet', vmin=0., vmax=1., lognorm=False,
                                   # property_name='z1', cmap='jet', vmin=0., vmax=0.5, lognorm=False,
                                   recording_plot=True,
                                   plot_name = os.path.join(plot_directory, 'plot_%.5d.png' % step))

        # Option 1B - With Pyvista:
        #--------------------------
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

