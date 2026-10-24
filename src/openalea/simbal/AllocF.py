#  -*- coding: utf-8 -*-

"""
    The script 'AllocF' contains all functions for running the development and functionning of the aerial parts of the plant.

    :copyright: see AUTHORS.
    :license: see LICENSE for details.
"""

import openalea.mtg
from openalea.mtg import *
from openalea.mtg import turtle as turt
from openalea.mtg.traversal import pre_order, post_order
from openalea.mtg.plantframe import color
import openalea.plantgl.all as pgl

from openalea.caribu.CaribuScene import CaribuScene
from openalea.caribu.data_samples import data_path

from openalea.simbal import simbal_parameters as param
from openalea.simbal.tools import (loading_parameters_from_file, temperature_modification, recording_MTG_properties,
                                   circle_coordinates, updating_coordinates, leaf_insertion_angle_evolution,
                                   plotting_shoot_MTG, plotting_shoots_with_pyvista)

import os
import warnings
import numpy as np
import pandas as pd
from math import sqrt, pi, floor, trunc, exp, isnan, cos, sin, acos
from decimal import Decimal

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
            try:
                # If the attribute of the parent element corresponds to a list:
                if isinstance(getattr(parent_element, attribute), list):
                    # Then we set the attribute of the child to an empty list:
                    setattr(child, attribute, [])
                # If the attribute of the parent element corresponds to a dictionnary:
                elif isinstance(getattr(parent_element, attribute), dict):
                    # Then we set the attribute of the child to an empty dictionnary:
                    setattr(child, attribute, {})
                # If the attribute of the parent element corresponds to a dictionnary:
                elif isinstance(getattr(parent_element, attribute), str):
                    # Then we set the attribute of the child to "[Undetermined]":
                    setattr(child, attribute, "[Undetermined]")
                # Otherwise, we assume that this should be a numerical value and set it to 0:
                else:
                    setattr(child, attribute, 0)
            # If the above has not worked, we set the attribute to None in the child:
            except:
                setattr(child, attribute, None)
        # In the worst case scenario, everything is set to NA:
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
    internode.node_rank = parent_element.node_rank + 1
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
    # We reset the time since formation:
    internode.actual_time_since_formation = 0.
    internode.thermal_time_since_formation = 0.
    # We initialize the share of light intensity perceived by this organ
    # compared to the sum of light intensities of all other organs:
    internode.light_intensity_share = 0.
    internode.surfacic_PAR_averaged_over_time =  0.
    # We initialize a list of previously recorded light intensity shares
    # (used to calculate the average share over a certain period):
    internode.previous_light_intensity_shares = []
    internode.previous_PAR_intercepted_values = []

    # 2) We then create a new leaf that is actually a "branch" on the previous internode:
    # ------------------------------------------------------------------------------------
    leaf = ADDING_A_CHILD_ORGAN(parent_element=internode,
                                edge_type='+',
                                label='Leaf',
                                properties_identical_to_parent=properties_identical_to_parent,
                                nil_properties=nil_properties)
    # We specified some properties that may be distinct from the parent element:
    leaf.node_rank = parent_element.node_rank + 1
    leaf.length = leaf_length
    leaf.radius = leaf_length * param.leaf_width_to_length_ratio / 2.
    leaf.angle_down = param.leaf_insertion_angle
    leaf.angle_roll = 0.
    leaf.type = "Normal"
    # We calculate the surface, assuming the leaf is an ellipse characterized by its length and width:
    leaf.surface = leaf.length / 2. * leaf.radius * pi
    # We define the maximal surface of the new leaf:
    leaf.surface_max = param.LA_max
    # We calculate the volume, taking into account a minimal thickness:
    leaf.volume = leaf.surface * minimal_leaf_thickness
    leaf.struct_mass = leaf.volume * struct_mass_density
    leaf.photosynthetic_area = leaf.surface
    # We specify the ramification order :
    leaf.ramification_order = parent_element.ramification_order
    # And we reset the time since formation:
    leaf.actual_time_since_formation = 0.
    leaf.thermal_time_since_formation = 0.
    # We initialize the share of light intensity perceived by this organ
    # compared to the sum of light intensities of all other organs:
    leaf.light_intensity_share = 0.
    leaf.surfacic_PAR_averaged_over_time
    # We initialize a list of previously recorded light intensity shares
    # (used to calculate the average share over a certain period):
    leaf.previous_light_intensity_shares = []
    leaf.previous_PAR_intercepted_values = []

    # 3) We finally create an axilliary bud that is actually a "branch" on the previous internode:
    # ---------------------------------------------------------------------------------------------
    bud = ADDING_A_CHILD_ORGAN(parent_element=leaf, # WATCH OUT: here, we set the bud as the child of the leaf, not the internode!
                               edge_type='+',
                               label='Bud',
                               properties_identical_to_parent=False,
                               nil_properties=True)
    # We specified some properties that may be distinct from the parent element:
    bud.node_rank = parent_element.node_rank + 1
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
    # We reset the time since formation:
    bud.actual_time_since_formation = 0.
    bud.thermal_time_since_formation = 0.
    # We initialize the share of light intensity perceived by this organ
    # compared to the sum of light intensities of all other organs:
    bud.light_intensity_share = 0.
    bud.surfacic_PAR_averaged_over_time = 0.
    # We initialize a list of previously recorded light intensity shares
    # (used to calculate the average share over a certain period):
    bud.previous_light_intensity_shares = []
    bud.previous_PAR_intercepted_values = []

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
                             initial_leaf_length = 0.002,
                             time_step_in_seconds = 60*60):
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

    # We initialize the structural mass produced by the apparition of new phytomers:
    new_structural_mass = 0.

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
                                                properties_identical_to_parent = False,
                                                nil_properties=True)
                # print("We have moved the previous apex", apex.index(), "on the node", internode.index(),"; its index is now", new_apex.index())
                new_apex.radius = internode.radius
                # TODO: WATCH OUT BELOW, IN CASE YOU WANT TO CONTINUE ADDING PHYTOMERS ON THE TOP OF EACH OTHERS
                # # Finally, we update the list of apices by adding the new apex that has been formed,
                # # in case it could itself form a new phytomer:
                # list_of_apices.append(new_apex)

                # REGISTERING THE COST OF THE FORMATION FOR THE NEW PHYTOMER:
                # For each new element of the new phytomer:
                for element in [internode, leaf, bud]:
                    # We add the structural mass of each new element:
                    new_structural_mass += element.struct_mass
                    # We register locally the costs of C&N for their construction:
                    element.actual_C_consumption_rate_by_shoot_growth\
                        = element.struct_mass * param.structural_C_content / param.yield_growth / time_step_in_seconds
                    element.actual_N_consumption_rate_by_shoot_growth\
                        = element.struct_mass * param.structural_N_content / time_step_in_seconds

    # We calculate the corresponding C (gC) and N (gN) cost for the plant, including the C respiratory cost:
    C_cost = new_structural_mass * param.structural_C_content / param.yield_growth
    N_cost = new_structural_mass * param.structural_N_content
    # # We update the structural mass of the plant:
    # shoot_system.shoot_struct_mass += new_structural_mass

    return C_cost, N_cost
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# This function enables to generate an initial MTG with a certain number of leaves and corresponding internodes:
def create_shoot_MTG(number_of_cotyledons=2,
                     number_of_phytomers=3,
                     initial_cotyledon_length=0.003,
                     initial_leaf_length=0.003,
                     initial_internode_radius = 0.001,
                     initial_hypocotyl_length = 0.01,
                     initial_internode_height = 0.0005,
                     initial_labile_C_concentration = 0.001,
                     initial_labile_N_concentration = 0.0001,
                     struct_mass_density=1e5, minimal_leaf_thickness=1e-3,
                     sowing_depth=0.):
    """
    This function is used to create the shoot part of a plant with a succession of different phytomers on the same branch.

    Parameters
    ----------
    number_of_phytomers: the number of cotyledons to be created
    number_of_phytomers: the number of successive phytomers to be created
    initial_leaf_length: the length of the leaves to be created (meters)
    leaf_surface_max: the maximal surface that the leaves will be allowed to reach (useful for cotyledons)
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
                              type="Hypocotyl",
                              node_rank=0,
                              # Geometry:
                              radius=initial_internode_radius,
                              potential_radius=initial_internode_radius,
                              length=initial_hypocotyl_length,
                              potential_length=initial_hypocotyl_length,
                              angle_down=0.,
                              angle_roll=0.,
                              x1=0.,
                              x2=0.,
                              y1=0.,
                              y2=0.,
                              z1=-sowing_depth,
                              z2=-sowing_depth,
                              # C & N labile pools:
                              labile_C=0.,
                              labile_C_deficit=0.,
                              labile_N=0.,
                              labile_N_deficit=0.,
                              # Rates:
                              photosynthesized_C_rate = 0.,
                              shoot_C_distribution_rate = 0.,
                              shoot_N_distribution_rate = 0.,
                              C_input_rate_for_sink = 0.,
                              shoot_maintenance_respiration_rate = 0.,
                              actual_C_consumption_rate_by_shoot_growth = 0.,
                              N_input_rate = 0.,
                              actual_N_consumption_rate_by_shoot_growth = 0.,
                              shoot_N_mobilization_rate = 0.,
                              net_shoot_N_immobilization_rate = 0.,
                              # Other attributes:
                              senescence_status = False,
                              PAR_intercepted=0.,
                              surfacic_PAR_intercepted=0.,
                              surfacic_PAR_averaged_over_time=0.,
                              previous_PAR_intercepted_values=[],
                              light_intensity_share=0.,
                              previous_light_intensity_shares=[],
                              Caribu_Eabs = 0.,
                              Caribu_Ei = 0.,
                              Caribu_area = 0.
                              )

    # We access this first element (hypocotyl), and consider it as an internode on which the next upper shoot elements
    # will be attached:
    new_internode = g.node(id_base)
    # We calculate the surface, volume and structural mass of this first element:
    new_internode.surface = new_internode.radius * 2. * pi * new_internode.length
    new_internode.volume = new_internode.radius ** 2. * pi * new_internode.length
    new_internode.struct_mass = new_internode.volume * struct_mass_density
    new_internode.photosynthetic_area = new_internode.surface

    new_internode.labile_C = initial_labile_C_concentration * new_internode.struct_mass
    new_internode.labile_N = initial_labile_N_concentration * new_internode.struct_mass

    # We define the ramification order of this first internode as 0:
    new_internode.ramification_order = 0

    # We set the time since formation:
    new_internode.actual_time_since_formation = 0.
    new_internode.thermal_time_since_formation = 0.

    # For each new cotyledon to be created:
    for i in range(0, number_of_cotyledons):
        # We first create an internode of zero length, and then a connected cotyledon and a fake bud:
        [new_internode, new_leaf, new_bud] = ADDING_A_PHYTOMER(parent_element=new_internode,
                                                               struct_mass_density=struct_mass_density,
                                                               minimal_leaf_thickness=minimal_leaf_thickness,
                                                               internode_radius=new_internode.radius * 0.9,
                                                               internode_length=0.,
                                                               leaf_length=initial_cotyledon_length,
                                                               bud_length=0.,
                                                               properties_identical_to_parent = True)
        # We set the thermal time since the phytomer formed according to the phyllochron:
        new_internode.node_rank = 0
        new_internode.thermal_time_since_formation = 0.
        new_leaf.thermal_time_since_formation = 0.
        new_bud.thermal_time_since_formation = 0.
        # We define the type of the internode:
        new_internode.type = "Cotyledon_supporting_segment"
        new_internode.angle_roll = 180
        # We set its labile content to 0, as it has no length or mass:
        new_internode.labile_C = 0.
        new_internode.labile_N = 0.
        # We define the leaf as a cotyledon:
        new_leaf.type = "Cotyledon"
        new_leaf.node_rank = 0
        new_leaf.surface_max = param.cotyledon_area_max
        # We set the bud diameter to 0, so that it will never be able to emerge:
        new_bud.radius = 0.
        new_bud.node_rank = 0

    # For each new phytomer to be created:
    for i in range(0, number_of_phytomers):
        # We first create an internode, and then a connected leaf and a bud:
        [new_internode, new_leaf, new_bud] = ADDING_A_PHYTOMER(parent_element=new_internode,
                                                               struct_mass_density=struct_mass_density,
                                                               minimal_leaf_thickness=minimal_leaf_thickness,
                                                               internode_radius=new_internode.radius * 0.9,
                                                               internode_length=initial_internode_height,
                                                               leaf_length=initial_leaf_length,
                                                               bud_length=0.,
                                                               properties_identical_to_parent = True)
        # We set the thermal time since the phytomer formed according to the phyllochron:
        new_internode.thermal_time_since_formation = (number_of_phytomers - i) * param.phyllochron
        new_leaf.thermal_time_since_formation = (number_of_phytomers - i) * param.phyllochron
        new_bud.thermal_time_since_formation = (number_of_phytomers - i) * param.phyllochron
        # We set the type of the elements of this phytomer as normal:
        new_internode.type= "Normal"
        new_leaf.type = "Normal"
        new_bud.type = "Normal"
        # We also need to specify the labile contents of the internode and the leaf:
        if new_internode.length > 0.:
            new_internode.labile_C = initial_labile_C_concentration
            new_internode.labile_N = initial_labile_N_concentration
        if new_leaf.length > 0.:
            new_leaf.labile_C = initial_labile_C_concentration
            new_leaf.labile_N = initial_labile_N_concentration

    # Finally, on the last internode, we add a new apex that may continue to form new phytomers over the simulation:
    new_apex = ADDING_A_CHILD_ORGAN(parent_element=new_internode,
                                    edge_type="<",
                                    label="Apex",
                                    nil_properties=True)
    new_apex.length = 0.
    new_apex.node_rank = 0

    # Eventually, we add special properties specific to the base of the shoot system that are not expected elsewhere
    # in the shoot system:
    base = g.node(0)
    # The following property will allow ramification to start:
    base.ramification_allowed = False

    return g
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function for loading a MTG from CSV file containing its properties (not working?):
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
def geometry_evolution(shoot_MTG, plant_thermal_age_in_degree_days):
    """
    This function updates the evolution of geometrical features of the shoot MTG over time.
    Parameters
    ----------
    shoot_MTG: the shoot MTG containing each organ as node
    plant_thermal_age: the thermal age of the whole plant

    Returns
    -------
    [Nothing; the shoot MTG has been updated with some new geometrical properties for specific organs]
    """

    # We cover each organ n of the shoot:
    for vid in shoot_MTG.vertices_iter(scale=1):
        n = shoot_MTG.node(vid)
        # If the organ is a leaf:
        if n.label == "Leaf":
            leaf_age = n.thermal_time_since_formation / (24. * 60. * 60.) * param.T_ref_growth
            if n.type != "Cotyledon":
                # We compute the evolution of the insertion angle of this leaf on the stem, and modify the angle accordingly:
                n.angle_down = leaf_insertion_angle_evolution(plant_thermal_age_in_degree_days=plant_thermal_age_in_degree_days,
                                                              leaf_age_in_degree_days = leaf_age)
            else:
                n.angle_down = 75
    return shoot_MTG
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
                                 leaf_shape="oval",
                                 averaging_light_intensity_share=False,
                                 period_of_light_intensity_averaging_in_days=1,
                                 time_step_in_days=1/24.
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
    averaging_light_intensity_share: if True, the property "light_intensity_share" will be computed based on this time step and the light_intensity_share computed previously
    period_of_light_intensity_averaging_in_days: the period over which the light intensity share should be averaged (days)
    time_step_in_days: the time step of the simulation (days)

    Returns
    -------
    [Nothing is returned; the shoot MTG has been updated with the properties 'PAR_intercepted' and 'Caribu_area', 'Caribu_Ei', 'Caribu_Eabs']
    """

    # We rename the shoot MTG as 'g':
    g = shoot_MTG

    # We initialize the total PAR intercepted by the shoot (umol of photons per second)
    # and the sum of surfacic PAR of each organ (umol of photons per square meter of organ per second):
    total_PAR_intercepted = 0.
    total_surface = 0.
    sum_of_surfacic_PAR = 0.

    # LIGHT INTERCEPTION BY CARIBU:
    ###############################

    # can = scene3d(g)
    # zenith = [(100, (0, 0, -1))]
    # (xmin,ymin,xmax,ymax)=(0,0,1,1)
    # cs = CaribuScene(scene=can, light=zenith, pattern=(0,0,1,1))
    # raw, agg = cs.run(infinite=True, simplify=True, direct=False)

    # We initialize a turtle in PlantGL:
    turtle = turt.PglTurtle()
    # We create the scene with the plant(s).
    # If we want to study the light interception by the isolated plant:
    if isolated_plant:
        # Then we create a scene with only one plant with the classical PlantGL turtle:
        my_scene = turt.TurtleFrame(g,
                                    visitor=get_shoot_visitor(leaf_shape=leaf_shape, upscaling_factor=100),
                                    turtle=turtle, gc=False, all_roots=True)
    else:
        # Otherwise, we simulate a population of plants, from which the central plant will be the one on which Caribu
        # will compute the light interception:

        # # We create an updated list containing one shoot MTG of each type copied and distributed over the scene:
        # scene_information, list_of_scene_shoot_MTG, list_of_scene_root_MTG \
        #     = prepare_heterogeneous_plants_population(list_of_distinct_plant_names=["plant"],
        #                                               list_of_distinct_shoot_MTG=[shoot_g],
        #                                               list_of_distinct_root_MTG=[],
        #                                               list_of_numbers_of_plants_per_square_meter=[number_of_plants_per_square_meter],
        #                                               list_of_inter_row_distances=[inter_row_distance],
        #                                               list_of_number_of_rows_to_show=[5],
        #                                               list_of_scene_lengths=[scene_length],
        #                                               list_of_scene_widths=[scene_width],
        #                                               upscaling_factor=100)
        # # We create an updated scene, with reindexed shapes compared to the original indices in the list of MTG:
        # my_scene, scene_information = create_plants_population_scene(scene_information=scene_information,
        #                                                              list_of_shoot_MTG=list_of_scene_shoot_MTG,
        #                                                              list_of_root_MTG=[],
        #                                                              leaf_shape=leaf_shape,
        #                                                              simplified_leaf_shape_for_surface_calculation="diamond_shape",
        #                                                              upscaling_factor=100)
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

        # # As a stupid warning linked to delim_whitespace used in read_csv occurs below, we hide it using warnings:
        # warnings.filterwarnings("ignore", message=".*The 'delim_whitespace' keyword in pd.read_csv is deprecated and will be removed in a future version.*")
        reference_light_file = pd.read_csv(sky,
                                           # delim_whitespace=True,
                                           delimiter=" ",
                                           header=None)

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

    # In case we want to replicate the scene to infinity:
    if infinite:
        # We create a pattern based on the size of the scene:
        pattern = tuple([0*unit_factor,0*unit_factor,scene_length*unit_factor,scene_width*unit_factor])
    else:
        pattern = None

    # Finally, we create the Caribu scene and run the light interception:
    try:
        cs = CaribuScene(scene=my_scene, light=LIGHT, scene_unit='cm', pattern=pattern) #,opt=opts)
        raw, agg = cs.run(infinite=infinite, simplify=True, direct=True)
        # The resulting dictionnaries "agg" contains three dictionnaries for each "active" vertex if the MTG:
        #   area (float): the individual areas (m2)
        #   Eabs (float): the surfacic density of energy absorbed (m-2)
        #   Ei (float): the surfacic density of energy incoming (m-2)
        # NOTE: The other dictionnary raw contains three similar dictionnary at the scale of each triangle.
        # We update the properties of the MTG with these dictionnaries (i.e. one value for each vertex):
        g.properties()["Caribu_area"].update(agg['area'])
        g.properties()["Caribu_Eabs"].update(agg['Eabs'])
        g.properties()["Caribu_Ei"].update(agg['Ei'])
    except Exception as e:
        print("!!! ERROR: Caribu could not be run! Here is the error:")
        print(e)
        print("!!! We therefore set an arbitrary intercepted light by the plant:")
        total_PAR_intercepted = 0
        sum_of_surfacic_PAR = 0

    # ASSIGNING LIGHT INTENSITY AMONG THE ORGANS AND COMPUTING SUMS:
    ################################################################

    # We cover each organ n of the shoot:
    for vid in g.vertices_iter(scale=1):
        n = g.node(vid)
        try:
            n.PAR_intercepted = n.Caribu_Eabs * n.Caribu_area
            total_PAR_intercepted += n.Caribu_Eabs * n.Caribu_area
            total_surface += n.Caribu_area
            n.surfacic_PAR_intercepted = n.Caribu_Eabs
            sum_of_surfacic_PAR += n.Caribu_Eabs
        except:
            print("WATCH OUT: the intercepted PAR of the shoot element", vid, "could not be calculated!")
            print("We therefore set it to 0 umol m-2 s-1.")
            n.PAR_intercepted = 0.
            n.surfacic_PAR_intercepted = 0.

    # COMPUTING LIGHT INTENSITY SHARE AND AVERAGE OVER TIME:
    #######################################################

    # If we need to consider previous light intensity shares:
    if averaging_light_intensity_share:
        # Then we will need the number of previous actual shares to include in the calculation of the average share:
        number_of_previous_shares = round(period_of_light_intensity_averaging_in_days / time_step_in_days)

    # We initialize a counter for checking that the sum of light intensity share of each organ will be equal to 1:
    sum_of_light_fractions = 0.

    # We cover each organ n of the shoot:
    for vid in g.vertices_iter(scale=1):
        n = g.node(vid)

        # We first compute the actual share of light intensity for this time step, based on surfacic PAR,
        # or else based on the share of surfaces within the plant:
        if sum_of_surfacic_PAR > 0.:
            new_light_intensity_share = n.surfacic_PAR_intercepted / sum_of_surfacic_PAR
        elif total_surface > 0.:
            new_light_intensity_share = n.Caribu_area / total_surface
        else:
            new_light_intensity_share = 0.

        # If we do not take into account previous light intensity share:
        if not averaging_light_intensity_share:
            # Then the light intensity share is directly the one that has just been determined
            # based on the light intercepted at this time step:
            n.light_intensity_share = new_light_intensity_share
        # Otherwise, if we need to consider previous light intensity shares:
        else:
            # If light was intercepted at this time step, then we add the share of light intensity for this organ
            # to the list of previous light intensity shares:
            if sum_of_surfacic_PAR > 0.:
                # If the previously recorded list does not contain yet the right number of shares to be averaged:
                if len(n.previous_light_intensity_shares) < number_of_previous_shares:
                    # Then we just add the current actual share to the list:
                    n.previous_light_intensity_shares.append(new_light_intensity_share)
                    n.previous_PAR_intercepted_values.append(n.surfacic_PAR_intercepted)
                else:
                    # We remove the first element of the list of previous actual light intensity shares (i.e. the oldest one):
                    n.previous_light_intensity_shares.pop(0)
                    n.previous_PAR_intercepted_values.pop(0)
                    # And we add at the end of the list the new actual share of light intensity for this time step:
                    n.previous_light_intensity_shares.append(new_light_intensity_share)
                    n.previous_PAR_intercepted_values.append(n.surfacic_PAR_intercepted)
            # Finally, we calculate the new average light intensity share based on previously computed shares:
            if len(n.previous_light_intensity_shares) > 0.:
                n.light_intensity_share = np.mean(n.previous_light_intensity_shares)
                n.surfacic_PAR_averaged_over_time = np.mean(n.previous_PAR_intercepted_values)
            else:
                n.light_intensity_share = new_light_intensity_share
                n.surfacic_PAR_averaged_over_time = n.surfacic_PAR_intercepted

        # We add this share to the sum of fractions for eventually checking that the sum is equal to 1:
        sum_of_light_fractions += n.light_intensity_share

    # CONTROL & ADJUSTMENTS:
    # If the sum of all light fractions is not equal to 1 (this may particularly occur when new elements are created):
    if abs(sum_of_light_fractions - 1) > 1e-14:
        print("!!! WATCH OUT: the sum of the light intensities perceived by each shoot organ was not equal to 1!")
        print("We therefore adjusted each light intensity share to correct this!")
        # We will adjust the light intensity share of each organ so that their sum is now equal to 1.
        # We cover each organ n of the shoot:
        for vid in g.vertices_iter(scale=1):
            n = g.node(vid)
            # We adjust the light intensity share so that the sum of all shares now equals 0:
            n.light_intensity_share = n.light_intensity_share * 1. / sum_of_light_fractions

    return total_PAR_intercepted, sum_of_surfacic_PAR
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
def computing_light_interception_on_plant_population(scene,
                                                     scene_information,
                                                     genotype_names_list=["rapeseed", "faba_bean"],
                                                     shoot_MTG_dict={},
                                                     list_of_scene_MTG=[],
                                                     incoming_light_density = 400,
                                                     isolated_plant=False,
                                                     zenithal_light=False,
                                                     infinite=False,
                                                     averaging_light_intensity_share=False,
                                                     period_of_light_intensity_averaging_in_days=1,
                                                     time_step_in_days=1/24.
                                                     ):
    """
    This function computes the light intercepted by a plant within a specific population design defined by a scene.
    """

    # We initialize the total PAR intercepted by the shoot (umol of photons per second)
    # and the sum of surfacic PAR of each organ (umol of photons per square meter of organ per second):
    total_PAR_intercepted_dict = {}
    sum_of_surfacic_PAR_dict = {}

    # We verify that the light input is not nil:
    if incoming_light_density <= 0.:
        for genotype_name in genotype_names_list:
            total_PAR_intercepted_dict[genotype_name] = 0.
            sum_of_surfacic_PAR_dict[genotype_name] = 0.
        return total_PAR_intercepted_dict, sum_of_surfacic_PAR_dict

    # LIGHT INTERCEPTION BY CARIBU:
    ###############################

    # If only the zenithal light is considered
    if zenithal_light:
        # We set only one vertical light direction:
        light_direction = (0, 0, -1)
        LIGHT = [(incoming_light_density,light_direction)]
    else:
        # We open the standard light input file as a dataframe:
        sky = str(data_path('Turtle16soc.light'))

        # # As a stupid warning linked to delim_whitespace used in read_csv occurs below, we hide it using warnings:
        # warnings.filterwarnings("ignore", message=".*The 'delim_whitespace' keyword in pd.read_csv is deprecated and will be removed in a future version.*")
        reference_light_file = pd.read_csv(sky,
                                           # delim_whitespace=True,
                                           delimiter=" ",
                                           header=None)

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

    # In case we want to replicate the scene to infinity:
    if infinite:
        # We create a pattern based on the size of the scene:
        pattern = tuple([0*unit_factor,0*unit_factor,scene_length*unit_factor,scene_width*unit_factor])
    else:
        pattern = None

    # Finally, we create the Caribu scene and run the light interception:
    try:
        cs = CaribuScene(scene=scene, light=LIGHT, scene_unit='cm', pattern=pattern) #,opt=opts)
        raw, agg = cs.run(infinite=infinite, simplify=True, direct=True)
        # The resulting dictionnaries "agg" contains three dictionnaries for each "active" vertex if the MTG:
        #   area (float): the individual areas (m2)
        #   Eabs (float): the surfacic density of energy absorbed (m-2)
        #   Ei (float): the surfacic density of energy incoming (m-2)
        # NOTE: The other dictionnary raw contains three similar dictionnary at the scale of each triangle.

        # Now we reallocate the properties of light interception to each MTG in the population of plants within the scene:
        list_of_plant_id = scene_information["plant_id"]
        # For each plant in the population:
        for plant_id in list_of_plant_id:
            # We get the corresponding MTG in the list:
            g = list_of_scene_MTG[plant_id]
            # We get the dictionnary which allocated new id to original id of each shape corresponding to the plant:
            reindexer = scene_information["new_index_to_old_index_shoot"][plant_id]
            # We get the corresponding result of the Caribu-computed intercepted light fo this plant:
            sub_Ei = {reindexer[id]: agg["Ei"][id] for id in reindexer.keys()}
            sub_Eabs = {reindexer[id]: agg["Eabs"][id] for id in reindexer.keys()}
            sub_area = {reindexer[id]: agg["area"][id] for id in reindexer.keys()}
            # We update the properties of the current MTG with this new property that has been reindexed to correspond to the original vid:
            g.properties()["Caribu_area"].update(sub_area)
            g.properties()["Caribu_Ei"].update(sub_Ei)
            g.properties()["Caribu_Eabs"].update(sub_Eabs)
            g.properties()["PAR_intercepted"].update(sub_Eabs)

    except Exception as e:
        print("!!! ERROR: Caribu could not be run! Here is the error:")
        print(e)

    # ASSIGNING LIGHT INTENSITY AMONG THE ORGANS AND COMPUTING SUMS:
    ################################################################

    # We cover the distinct shoot MTG for which we want to update the light properties:
    for genotype_name in genotype_names_list:

        # We initialize the total properties for the plant:
        total_PAR_intercepted = 0.
        total_surface = 0.
        sum_of_surfacic_PAR = 0.

        # We access the corresponding original MTG:
        g = shoot_MTG_dict[genotype_name]
        # We get the ID number of the copied MTG in the middle of the scene for which light properties should be taken:
        middle_id = scene_information["middle_plant_id_" + genotype_name][0]
        # We access the copied MTG:
        scene_g = list_of_scene_MTG[middle_id]

        # We cover each organ n of the original shoot MTG:
        for vid in g.vertices_iter(scale=1):
            n = g.node(vid)
            scene_n = scene_g.node(vid)
            try:
                n.Caribu_area = scene_n.Caribu_area
                n.Caribu_Ei = scene_n.Caribu_Ei
                n.Caribu_Eabs = scene_n.Caribu_Eabs
                n.PAR_intercepted = scene_n.Caribu_Eabs * scene_n.Caribu_area
                total_PAR_intercepted += scene_n.Caribu_Eabs * scene_n.Caribu_area
                total_surface += scene_n.Caribu_area
                n.surfacic_PAR_intercepted = scene_n.Caribu_Eabs
                sum_of_surfacic_PAR += scene_n.Caribu_Eabs
            except:
                print("WATCH OUT: the intercepted PAR of the shoot element", vid, "could not be calculated!")
                print("We therefore set it to 0 umol m-2 s-1.")
                n.PAR_intercepted = 0.
                n.surfacic_PAR_intercepted = 0.

        # COMPUTING LIGHT INTENSITY SHARE AND AVERAGE SURFACIC PAR OVER TIME:
        #####################################################################

        # If we need to consider previous light intensity shares:
        if averaging_light_intensity_share:
            # Then we will need the number of previous actual shares to include in the calculation of the average share:
            number_of_previous_shares = round(period_of_light_intensity_averaging_in_days / time_step_in_days)

        # We initialize a counter for checking that the sum of light intensity share of each organ will be equal to 1:
        sum_of_light_fractions = 0.

        # We cover each organ n of the shoot:
        for vid in g.vertices_iter(scale=1):
            n = g.node(vid)

            # We first compute the actual share of light intensity for this time step, based on surfacic PAR,
            # or else based on the share of surfaces within the plant:
            if sum_of_surfacic_PAR > 0.:
                new_light_intensity_share = n.surfacic_PAR_intercepted / sum_of_surfacic_PAR
            elif total_surface > 0.:
                new_light_intensity_share = n.Caribu_area / total_surface
            else:
                new_light_intensity_share = 0.

            # If we do not take into account previous light intensity share:
            if not averaging_light_intensity_share:
                # Then the light intensity share is directly the one that has just been determined
                # based on the light intercepted at this time step:
                n.light_intensity_share = new_light_intensity_share
                n.surfacic_PAR_averaged_over_time = n.surfacic_PAR_intercepted
            # Otherwise, if we need to consider previous light intensity shares:
            else:
                # If the previously recorded list does not contain yet the right number of shares to be averaged:
                if len(n.previous_light_intensity_shares) < number_of_previous_shares:
                    # If light was intercepted at this time step:
                    if sum_of_surfacic_PAR > 0.:
                        # Then we just add the current actual share to the list:
                        n.previous_light_intensity_shares.append(new_light_intensity_share)
                    # Whether light was present or not at this time step:
                    n.previous_PAR_intercepted_values.append(n.surfacic_PAR_intercepted)
                else:
                    # If light was intercepted at this time step:
                    if sum_of_surfacic_PAR > 0.:
                        # We remove the first element of the list of previous actual light intensity shares (i.e. the oldest one):
                        n.previous_light_intensity_shares.pop(0)
                        # And we add at the end of the list the new actual share of light intensity for this time step:
                        n.previous_light_intensity_shares.append(new_light_intensity_share)
                    # Whether light was present or not at this time step:
                    n.previous_PAR_intercepted_values.pop(0)
                    n.previous_PAR_intercepted_values.append(n.surfacic_PAR_intercepted)

                # Finally, we calculate the new average light intensity share based on previously computed shares:
                if len(n.previous_light_intensity_shares) > 0.:
                    n.light_intensity_share = np.mean(n.previous_light_intensity_shares)
                    n.surfacic_PAR_averaged_over_time = np.mean(n.previous_PAR_intercepted_values)
                else:
                    n.light_intensity_share = new_light_intensity_share
                    n.surfacic_PAR_averaged_over_time = n.surfacic_PAR_intercepted

            # We add this share to the sum of fractions for eventually checking that the sum is equal to 1:
            sum_of_light_fractions += n.light_intensity_share

        # CONTROL & ADJUSTMENTS:
        # If the sum of all light fractions is not equal to 1 (this may particularly occur when new elements are created):
        if abs(sum_of_light_fractions - 1) > 1e-14:
            print("!!! WATCH OUT: the sum of the light intensities perceived by each shoot organ was", sum_of_light_fractions, "and not 1!")
            print("We therefore adjusted each light intensity share to correct this!")
            # We will adjust the light intensity share of each organ so that their sum is now equal to 1.
            # We cover each organ n of the shoot:
            for vid in g.vertices_iter(scale=1):
                n = g.node(vid)
                # We adjust the light intensity share so that the sum of all shares now equals 0:
                n.light_intensity_share = n.light_intensity_share * 1. / sum_of_light_fractions

        # We record the total intercepted PAR and the sum of surfacic PAR for this plant:
        total_PAR_intercepted_dict[genotype_name] = total_PAR_intercepted
        sum_of_surfacic_PAR_dict[genotype_name] = sum_of_surfacic_PAR

    return total_PAR_intercepted_dict, sum_of_surfacic_PAR_dict
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function for calculating the actual photosynthesis rate:
def computing_photosynthesis_rate(shoot_system, shoot_MTG, air_temperature_in_degree_Celsius=10.):
    """
    This function calculates the actual photosythesis rate (in gC per second) for each organ of the shoot MTG, based on
    its intercepted PAR (micromol of photons per m2 per second), its photosynthetic surface (m2)
    and a temperature-dependent parameter SCA, i.e. its specific carbon assimilation (gC per micromol of photons).

    Parameters
    ----------
    shoot_system: the object containing the main properties of the shoot system (not the shoot MTG!)
    shoot_MTG: the shoot MTG containing each organ as node
    air_temperature_in_degree_Celsius: the average temperature of the air experienced by the shoots

    Returns
    -------
    total_PAR_intercepted: the PAR intercepted by the whole plant (micromoles of photons per second per plant)
    total_photosynthesis_rate: the rate of photosynthesis of the whole plant (gC per second)
    [The MTG has been updated with the property 'photosynthesized_C_rate' in each organ.]
    """

    # We correct the specific carbon assimilation, based on air temperature:
    corrected_SCA = temperature_modification(temperature_in_Celsius=air_temperature_in_degree_Celsius,
                                             process_at_T_ref=param.SCA,
                                             T_ref=param.SCA_T_ref,
                                             A=param.SCA_A,
                                             B=param.SCA_B,
                                             C=param.SCA_C)

    # We initialize the total intercepted PAR (in micromoles of photons per second per plant)
    # and the total photosynthetic rate (gC per second per plant):
    total_PAR_intercepted = 0.
    total_photosynthesis_rate = 0.

    # We correct the specific carbon assimilation according to mean shoot labile C concentration (possible inhibition):
    shoot_C_conc = (shoot_system.shoot_labile_C_pool - shoot_system.shoot_labile_C_deficit) / shoot_system.shoot_struct_mass

    #TODO: WATCH OUT BELOW - NO PHOTOSYNTHESIS INHIBITION!
    # inhibiting_labile_C_concentration = 0.1
    # corrected_SCA = corrected_SCA * (1 - shoot_C_conc / inhibiting_labile_C_concentration)
    # if corrected_SCA < 0.:
    #     corrected_SCA = 0.
    #     print("WATCH OUT: photosynthesis was inhibited because of a too high shoot labile C concentration!")

    # We cover each organ n of the shoot:
    for vid in shoot_MTG.vertices_iter(scale=1):
        n = shoot_MTG.node(vid)
        if n.surface <= 0. or n.struct_mass <= 0.:
            continue
        n.photosynthetic_area = n.surface
        # We calculate the photosynthesis rate for this organ:
        n.photosynthesized_C_rate = corrected_SCA * n.surfacic_PAR_intercepted * n.photosynthetic_area
        #TODO: have the SCA vary with the level of "green" within the leaf? with the level of N? with the senescence status?
        # We increment the total PAR intercepted by the plant:
        total_PAR_intercepted += n.surfacic_PAR_intercepted * n.photosynthetic_area
        total_photosynthesis_rate +=  n.photosynthesized_C_rate

    return total_PAR_intercepted, total_photosynthesis_rate
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function for computing the supply of resources from the seed to the new seedling after germination:
def computing_supply_from_seed_rate(shoot_system, air_temperature_in_degree_Celsius=10.):

    # We compute a modifying factor based on temperature:
    temperature_modifier = temperature_modification(temperature_in_Celsius=air_temperature_in_degree_Celsius,
                                                    process_at_T_ref=1,
                                                    T_ref=param.shoot_growth_T_ref,
                                                    A=param.shoot_growth_A,
                                                    B=param.shoot_growth_B,
                                                    C=param.shoot_growth_C)

    # We record the rates of C and N supply by the germinating seed within the shoot system object:
    shoot_system.C_supply_from_seed_rate = param.rate_of_C_supply_from_seed * temperature_modifier
    shoot_system.N_supply_from_seed_rate = param.rate_of_N_supply_from_seed * temperature_modifier

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function for computing the C cost for maintenance respiration (as opposed to growth respiration):
def computing_maintenance_rate(shoot_MTG, air_temperature_in_degree_Celsius=10.):
    """
    This function calculates the respiration rate (in gC per second) for the maintenance of each organ of the shoot MTG.
    It is based on the residual maintenance respiration calculated by Thornley and Cannell (2000) or Barillot et al. (2016).

    Parameters
    ----------
    shoot_MTG: the shoot MTG containing each organ as node
    air_temperature_in_degree_Celsius: the average temperature of the air experienced by the shoots (degree Celsius)

    Returns
    -------
    total_maintenance_rate: total maintenance respiration rate for the whole plant (gC per second)
    [The MTG has been updated with the property 'maintenance_respiration_rate' per organ (gC per second)]
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

    # We initialize the total maintenance rate of the plant (gC per second per plant):
    total_maintenance_rate = 0.

    # We cover each organ n of the shoot:
    for vid in shoot_MTG.vertices_iter(scale=1):
        n = shoot_MTG.node(vid)

        # If the element does not have any C left, it cannot respire!
        if n.labile_C <= 0.:
            n.shoot_maintenance_respiration_rate = 0.
            continue

        # The maintenance rate is based on a maximal rate of respiration per gram of structural mass, modulated by the
        # amount of labile C in the organ:
        n.shoot_maintenance_respiration_rate = corrected_Vmax * n.labile_C / (Km_maintenance + n.labile_C) * n.struct_mass

        # We increase the total maintenance rate with the contribution of the current organ:
        total_maintenance_rate +=  n.shoot_maintenance_respiration_rate

    return total_maintenance_rate
#-----------------------------------------------------------------------------------------------------------------------
# #-----------------------------------------------------------------------------------------------------------------------
# # Function for redistributing the C among shoot organs, with a C allocation. Different possibilities can be considered:
# # - proportionnal to PARi for each organ compared to the total PAR intercepted by the plant
# # - proportionnal to PARi/surface for each organ compared to the sum of PARi/surface of all organs
# # - proportionnal to PARi/surface of the leaf of the phytomer for each organ compared to the sum of PARi/surface of all leaves
# def computing_shoot_C_redistribution(shoot_MTG,
#                                      shoot_available_C = None,
#                                      intercepted_light = None,
#                                      time_step_in_seconds=60*60,
#                                      option="classic_PAR"):
#
#     """
#     This function computes the rate of labile C reallocation (positive = source, negative = sink) for each aerial organ,
#     based on its intercepted light. It does not however updates the amount of labile C in each organ (this is done in
#     another function).
#
#     Parameters
#     ----------
#     shoot_MTG: the MTG corresponding to the aerial parts of the plant
#     shoot_available_C: the amount of available C to be spread (gC per plant) (if None, it will be recalculated based
#                        on the labile C available in each organ)
#     intercepted_light: the total PAR intercepted by the plant (micromoles of photons per second per plant)
#                       (if None, it will be recalculated based on the intercepted PAR defined for each organ)
#     time_step_in_seconds: the time step of the simulation (seconds)
#     option: either "classic_PAR" (we compare micromol of photons per second per organ vs per plant)
#             or "surfacic_PAR" (we compared micromol of photons per m2 of organ per second per organ vs per plant)
#
#     Returns
#     -------
#     available_C: the total amount of C available in the plant that will be (re)allocated (gC per plant)
#     [The MTG has been updated with the property "shoot_C_distribution_rate" in each organ]
#     """
#
#     # 1) We calculate the amount of shoot C available and the total intercepted light:
#     #---------------------------------------------------------------------------------
#     if shoot_available_C is None:
#         # We initialize the available C in the shoot system to 0:
#         available_C = 0.
#     else:
#         # Otherwise, the available C is directly the provided shoot available C:
#         available_C = shoot_available_C
#
#     if intercepted_light is None:
#         # We initialize the intercepted light by the shoot system to 0:
#         intercepted_light_by_the_plant = 0.
#     else:
#         # Otherwise, the intercepted light is directly the light that has been provided:
#         intercepted_light_by_the_plant = intercepted_light
#
#     # We initialize the sum of surfacic PAR, if needed:
#     if option == "surfacic_PAR":
#         sum_of_surfacic_PAR = 0.
#
#     # If either shoot_available_C or intercepted_light have not been provided:
#     if shoot_available_C is None or intercepted_light is None:
#         # Then we cover each organ n of the shoot:
#         for vid in shoot_MTG.vertices_iter(scale=1):
#             n = shoot_MTG.node(vid)
#             # We add the initial amount of labile_C and the net gain in C of the current organ:
#             available_C += n.labile_C
#             # We add the light (or light density) intercepted by the organ:
#             if option == "classic_PAR":
#                 intercepted_light_by_the_plant += n.PAR_intercepted
#             elif option == "surfacic_PAR" and n.surface > 0.:
#                 sum_of_surfacic_PAR += n.PAR_intercepted / n.surface
#
#     # 2) We redistribute the C according to specific rules:
#     #------------------------------------------------------
#
#     # We cover each organ n of the shoot:
#     for vid in shoot_MTG.vertices_iter(scale=1):
#         n = shoot_MTG.node(vid)
#         # The corresponding rate of C redistribution for this organ is calculated:
#         # (NOTE: if positive, it corresponds to a C gain; if negative, it corresponds to a C source)
#         if option == "classic_PAR":
#             n.shoot_C_distribution_rate = (available_C * n.PAR_intercepted / intercepted_light_by_the_plant
#                                              / time_step_in_seconds)
#         elif option == "surfacic_PAR" and n.surface > 0. and sum_of_surfacic_PAR > 0.:
#             n.shoot_C_distribution_rate = (available_C * n.PAR_intercepted / n.surface / sum_of_surfacic_PAR
#                                              / time_step_in_seconds)
#         else:
#             n.shoot_C_distribution_rate = 0.
#
#     return available_C
#-----------------------------------------------------------------------------------------------------------------------
# # Function for distributing the N coming from the roots into each shoot organ:
# def computing_shoot_N_distribution(shoot_MTG, shoot_available_N=None, time_step_in_seconds=60.*60.*24.):
#     """
#         This function distributes the whole amount of labile N available to the whole plant to each organ which has a
#         positive surface, proportionnaly to the mean surfacic PAR intercepted by the organ compared to the sum of all
#         surfacic PARs of shoot elements.
#
#         Parameters
#         ----------
#         shoot_MTG: the MTG corresponding to the shoot system
#         shoot_available_N: the amount of available N to be allocated (gN per plant) (if None, it will be recalculated based
#                            on the labile N available in each organ)
#         time_step_in_seconds: the time step of the simulation (in seconds)
#
#         Returns
#         -------
#         available_N: the total amount of N that will be allocated to the shoot organs (gC per plant)
#         [The MTG has been updated with the property 'shoot_N_distribution_rate' in each organ (gN per second per organ)]
#         """
#
#     # 1) Computing the intercepted light (and possibly the available N to distribute):
#     #---------------------------------------------------------------------------------
#
#     # We initialize the total amount of light intercepted by the whole plant:
#     intercepted_light_by_the_plant = 0.
#     sum_of_mean_surfacic_PAR = 0.
#
#     # If the amount of N available in the shoot system has not been provided yet:
#     if shoot_available_N is None:
#         # We initialize the available N in the shoot system to 0:
#         available_N = 0.
#     else:
#         #-----------------------------------------------------------------------------------------------------------------------
#         # Otherwise, the available N is initialized with the provided shoot available N:
#         available_N = shoot_available_N
#
#     # Then we cover each organ n of the shoot:
#     for vid in shoot_MTG.vertices_iter(scale=1):
#         n = shoot_MTG.node(vid)
#         if n.surface <= 0.:
#             continue
#         # We compute not only the total intercepted light, but also the sum of the mean surfacic PARs:
#         intercepted_light_by_the_plant += n.PAR_intercepted
#         n.mean_surfacic_PAR = n.PAR_intercepted / n.surface
#         sum_of_mean_surfacic_PAR +=  n.mean_surfacic_PAR
#         # If the amount of shoot available N has not been provided, we add the amount of available N of the current element:
#         if shoot_available_N is None:
#             available_N += n.labile_N - n.labile_N_deficit
#
#     # print("Before N redistribution among the shoot organs, the total amount of labile N was", available_N)
#
#     # 2) Computing the rate of N distribution in each shoot organ
#     # -----------------------------------------------------------
#
#     # We cover each organ of the shoot:
#     for vid in shoot_MTG.vertices_iter(scale=1):
#         # We access the current organ"n":
#         n = shoot_MTG.node(vid)
#         if n.surface <= 0.:
#             continue
#         # The new target for labile_C in this organ is calculated according to the relative intercepted light:
#         # labile_N_target = plant_available_N * n.PAR_intercepted / intercepted_light_by_the_plant
#         fraction = n.mean_surfacic_PAR / sum_of_mean_surfacic_PAR
#         # print("   xxx For element", vid, "the fraction of light was", fraction)
#         labile_N_target = available_N * fraction
#         n.shoot_N_distribution_rate = labile_N_target / time_step_in_seconds
#
#     return available_N
# #-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# # TODO: Complete or delete this function
# def computing_shoot_N_remobilization(shoot_MTG, k_N_mobilization):
#     """
#     This function calculates the rate of N mobilization (gN per second) from the N labile pool in each shoot organ,
#     based on first-order kinetics.
#     :param g: the MTG corresponding to the shoot system
#     :param k_N_mobilization: the time constant associated to the first order decay (second-1)
#     :return: [no return - the MTG g has been updated with the property 'N_mobilization_rate']
#     """
#
#     # We cover each organ of the shoot:
#     for vid in shoot_MTG.vertices_iter(scale=1):
#         # We access the current organ"n":
#         n = shoot_MTG.node(vid)
#
#         # We exclude elements where a deficit of labile N has been previously recorded:
#         if n.labile_N_deficit >= 0.:
#             n.shoot_N_mobilization_rate = 0.
#             continue
#
#         # We compute the rate of N mobilization:
#         n.shoot_N_mobilization_rate = k_N_mobilization * n.labile_N
#
#     return
# #-----------------------------------------------------------------------------------------------------------------------
# #-----------------------------------------------------------------------------------------------------------------------
# # TODO: Complete or delete this function - only relevant if labile N is disintguished from reserve N...
# def computing_shoot_N_immobilization(shoot_MTG):
#     """
#     This function calculates the rate of N immobilization (gN per second) in each shoot organ.
#     :param shoot_MTG: the MTG corresponding to the shoot system
#     :return: [no return - the MTG g has been updated with the property 'shoot_N_immobilization_rate']
#     """
#
#     # We cover each organ of the shoot:
#     for vid in shoot_MTG.vertices_iter(scale=1):
#         # We access the current organ"n":
#         n = shoot_MTG.node(vid)
#         # We compute the rate N mobilization:
#         n.net_shoot_N_immobilization_rate = 0.
#
#     return
# #-----------------------------------------------------------------------------------------------------------------------
# #-----------------------------------------------------------------------------------------------------------------------
# def computing_shoot_C_storage(shoot_system, shoot_MTG, time_step_in_seconds=60.*60.*24.):
#     """
#     This function calculates the rate of labile C immobilization into C storage.
#     :param shoot_system: the object corresponding to the whole shoot system (not the MTG!)
#     :param shoot_MTG: the MTG corresponding to the shoot system
#     :param time_step_in_seconds: the time step of the simulation (seconds)
#     :return: [no return - the MTG g has been updated with the property 'C_storage_rate' (gC per second)
#               and the C available in the shoot system has been reset to 0]
#     #TODO: Convert the additional C storage into an increase in organ thickness?
#     """
#
#     # We rename the shoot MTG for simplification:
#     g=shoot_MTG
#
#     # 1) Computing the total intercepted light (if not calculated elsewhere)
#     #-----------------------------------------------------------------------
#     # OPTION 1: if it has not been calculated yet
#     # We initialize the total amount of light intercepted by the whole plant:
#     intercepted_light_by_the_plant = 0.
#     # We cover each organ n of the shoot:
#     for vid in g.vertices_iter(scale=1):
#         n = g.node(vid)
#         intercepted_light_by_the_plant += n.PAR_intercepted
#
#     # # OPTION 2: if it has already been calculated and stored in the base of the shoot system:
#     # intercepted_light_by_the_plant = plant.total_PAR_intercepted_by_the_plant
#
#     # 2) Computing the rate of C distribution in each shoot organ
#     # -----------------------------------------------------------
#
#     # We access the value of the total amount of available C, stored in the base of the shoot system:
#     plant_available_C = shoot_system.shoot_labile_C_pool - shoot_system.shoot_labile_C_deficit
#
#     # We cover each organ of the shoot:
#     for vid in g.vertices_iter(scale=1):
#         # We access the current organ"n":
#         n = g.node(vid)
#         # The new C storage rate in this organ is calculated according to the relative intercepted light:
#         n.shoot_C_storage_rate = plant_available_C * n.PAR_intercepted / intercepted_light_by_the_plant / time_step_in_seconds
#
#     # We reset the value of the plant available C:
#     shoot_system.shoot_labile_C_pool = 0.
#
#     return
# #-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
def potential_shoot_growth(shoot_MTG,
                           SFQN = 1,
                           leaf_minimal_thickness = 1e-3,
                           leaf_shape="diamond_shape",
                           elong_internode = 5e-7,
                           plant_age = 0.,
                           time_step_in_seconds = 60.*60.*24.,
                           air_temperature_in_degree_Celsius = 10.):
    """
    This function computes the potential increase in volume of each organ in the shoot system.
    :param shoot_MTG: the MTG corresponding to the shoot system
    :param SFQN: the slope of the relationship of Leaf area = f(labile_N) (m2 per gN)
    :param leaf_minimal_thickness: the minimal thickness of a leaf (meter)
    :param elong_internode: a potential rate of elongation (meter per second)
    :param time_step_in_seconds: the time step of the simulation (seconds)
    :param air_temperature_in_degree_Celsius: the mean temperature experienced by the shoots (degree Celsius)
    :return:
    """

    # We rename the MTG for simplification:
    g = shoot_MTG

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

            # # We calculate the potential surface based on the available N in the leaf:
            # potential_leaf_surface = SFQN * n.labile_N

            #TODO: WATCH OUT - here we remove SFQN law and make a direct relationship with light interception!
            # We define a coefficient of surface increase per day at 20°C as a function of PAR:
            SIP = 1e-4 / 400. # We assume a potential increase by 1 cm2 per day for a PAR of 400 µmol m-2 s-1
            # potential_leaf_surface = n.surface * (1 + SIP * n.surfacic_PAR_intercepted * thermal_time_step_in_days/1.)
            potential_leaf_surface = n.surface + SIP * n.surfacic_PAR_averaged_over_time * thermal_time_step_in_days/1.
            # We calculate the potential surface increase (m2), knowing that the surface cannot become higher than a max:
            potential_surface_increase = min(potential_leaf_surface - n.surface, n.surface_max - n.surface)

            # We calculate the maximal surface increase (m2), knowing that the surface increase cannot become higher
            # than the maximal rate allowed:
            if n.type != "Cotyledon":
                maximal_surface_increase = n.surface * param.LA_max_relative_increase_per_day * thermal_time_step_in_days
            else:
                maximal_surface_increase = n.surface * param.LA_max_relative_increase_per_day/3. * thermal_time_step_in_days
            # The potential volume increment is based on a minimal thickness:
            n.potential_volume_increment = min(potential_surface_increase, maximal_surface_increase) * leaf_minimal_thickness
            #TODO: what happens if the leaf has already a thickness higher than the minimal one?

            if leaf_shape == "diamond_shape":
                # We consider the leaf as a diamond characterized by its length and its width, which is calculated as
                # the product of the length with an allometric coefficient. The surface of the diamond is then defined as:
                # surface = (length / 2.) * (param.leaf_width_to_length_ratio * length) / 2.
                # i.e. surface = length ^ 2 * param.leaf_width_to_length_ratio / 4.
                # Hence:
                n.potential_length = sqrt((n.surface+maximal_surface_increase) * 4 / (param.leaf_width_to_length_ratio))
                # The radius, which is half of the width, is also calculated based on the allometric coefficient:
                n.potential_radius = n.potential_length * param.leaf_width_to_length_ratio / 2.
            else:
                # We assume that the leaf is an ellipse characterized by its length and its width, which is calculated as
                # the product of the length with an allometric coefficient. The surface of the ellipse is then defined as:
                # surface = (length / 2.) * (param.leaf_width_to_length_ratio * length) / 2. * pi
                # i.e. surface = length ^ 2 * param.leaf_width_to_length_ratio / 4. * pi
                # Hence:
                n.potential_length = sqrt((n.surface+maximal_surface_increase) * 4 / (param.leaf_width_to_length_ratio * pi))
                # The radius, which is half of the width, is also calculated based on the allometric coefficient:
                n.potential_radius = n.potential_length * param.leaf_width_to_length_ratio / 2.

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
            if plant_age > param.minimal_age_for_internode_elongation and fast_elongation is True:
                # Then we set the potential elongation rate to the full potential:
                potential_elongation_rate = elong_internode * temperature_adjustement
            else:
                # Otherwise, we only allow only a very slow elongation rate:
                potential_elongation_rate = elong_internode * temperature_adjustement / 100. #TODO: change the parameter with the actual elongation before stem elongation
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
# #-----------------------------------------------------------------------------------------------------------------------
# # Function for computing labile C reallocation within shoot organs, based on intercepted PAR:
# def computing_C_redistribution_for_shoot_growth(shoot_system,
#                                                 shoot_MTG,
#                                                 yield_growth=0.8,
#                                                 structural_mass_density=1e5,
#                                                 structural_C_content=0.44,
#                                                 time_step_in_seconds = 60.*60.*24.,
#                                                 option="classic_PAR"
#                                                 ):
#     """
#         This function evalates whether shoot organs should get or give labile C to other organs, based on their intercepted
#         PAR, and computes the rate of exchange between each organ and a general, shared pool of labile C, which is updated.
#
#         Parameters
#         ----------
#         shoot_system: the object corresponding the general shoot system (not the MTG)
#         shoot_MTG: the MTG containing all shoot organs
#         yield_growth: growth yield as defined by Thornley and Cannell (2000), which corresponds to the fraction of C consumed
#                       by growth that is actually used for biomass production (gC per gC)
#         structural_mass_density: density of the structural mass of shoot tissue (gram of dry structural mass per cubic meter)
#         structural_C_content: carbon content of the dry structural mass (gC per gram of dry structural mass)
#         time_step_in_seconds: the time step of the simulation (seconds)
#         option: either "classic_PAR" (we compare micromol of photons per second per organ vs per plant)
#                 or "surfacic_PAR" (we compared micromol of photons per m2 of organ per second per organ vs per plant)
#
#         Returns
#         -------
#         [The MTG has been updated with the property "C_input_rate_for_sink", without modifying the labile C content of each organ]
#         [The pool of labile C in the shoot system object has been updated to account for the next exchange with all shoot organs]
#         """
#
#     # We rename the shoot MTG for simplification:
#     g=shoot_MTG
#
#     # We initialize the stock of C that can be redistributed to sink leaves and the one that will be eventually used:
#     total_available_C = shoot_system.shoot_labile_C_pool
#     C_used_by_sink_organs = 0.
#     # We initialize the total intercepted light:
#     if option == "classic_PAR":
#         intercepted_light = 0.
#     elif option == "surfacic_PAR":
#         sum_of_surfacic_PAR = 0.
#         #-----------------------------------------------------------------------------------------------------------------------
#     # We initialize a list that will contain the indices of the leaves that will need additional C:
#     list_of_sink_organs = []
#
#     # 1) DEFINING THE AMOUNT OF AVAILABLE C:
#     ########################################
#
#     # We cover each organ n of the shoot:
#     for vid in g.vertices_iter(scale=1):
#         # We define the current organ as "n":
#         n = g.node(vid)
#
#         # We add the light intercepted by the organ:
#         if option == "classic_PAR":
#             intercepted_light += n.PAR_intercepted
#         elif option == "surfacic_PAR" and n.surface > 0.:
#             sum_of_surfacic_PAR += n.PAR_intercepted / n.surface
#
#         # We calculate the C that may be consumed by the potential growth, and the C that may remain in a "source" leaf:
#         C_consumption_by_growth = (n.potential_volume_increment * structural_mass_density
#                                    * structural_C_content / yield_growth)
#         remaining_labile_C = n.labile_C - C_consumption_by_growth
#
#         # If there is some C left after the potential growth:
#         if remaining_labile_C > 0.:
#             # The leaf is a source, and it increases the total amount of available C:
#             n.C_sink_status = False
#             total_available_C += remaining_labile_C
#         else:
#             # The leaf is a sink, and may receive additional C for its growth:
#             n.C_sink_status = True
#             list_of_sink_organs.append(vid)
#
#     # 2) RE-ALLOCATING C TO SINK ORGANS:
#     ###################################
#
#     # We cover each organ n of the shoot:
#     for vid in list_of_sink_organs:
#         # We define the current organ as "n":
#         n = g.node(vid)
#         # The new amount of labile_C in this organ is calculated, after its initialization:
#         additional_C_for_sink_organ = 0.
#         if option == "classic_PAR":
#             additional_C_for_sink_organ = total_available_C * n.PAR_intercepted / intercepted_light
#         elif option == "surfacic_PAR" and n.surface > 0. and sum_of_surfacic_PAR > 0.:
#             additional_C_for_sink_organ = total_available_C * n.PAR_intercepted / n.surface / sum_of_surfacic_PAR
#         # TODO: Note that in either case, only the fraction of C corresponding to (light_intercepted_by_sink_organs)
#         #  / (total_light_intercepted_by all_organs) is used...
#         # We adjust the properties of the organ:
#         # n.labile_C += additional_C_for_sink_organ
#         n.C_input_rate_for_sink = additional_C_for_sink_organ / time_step_in_seconds
#         C_used_by_sink_organs += additional_C_for_sink_organ
#
#     # Eventually, we store the remaining available C:
#     shoot_system.shoot_labile_C_pool += total_available_C - C_used_by_sink_organs
#     # The remaining amount of C here may be allocated to roots!
#
#     return
# #-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
def actual_shoot_growth(shoot_MTG,
                        yield_growth = 0.8,
                        structural_mass_density = 1e5,
                        structural_C_content = 0.44,
                        structural_N_content = 0.02,
                        leaf_minimal_thickness = 1e-3,
                        leaf_shape = "diamond_shape",
                        time_step_in_seconds = 60.*60.*24.,
                        air_temperature_in_degree_Celsius = 10.):
    """
    This function computes the actual growth of the shoot system, depending on the C and N available in each organ,
    as well as the ageing of each shoot element.
    
    Parameters
    ----------
    shoot_MTG: the MTG corresponding to the shoot system
    yield_growth: the fraction of the total C used for growth that will be incorporated to the new mass (the rest is respired) (gC per gC)
    struct_mass_density: the density of the structural mass (gram of structural maqss per m3)
    structural_C_content: the C content of the structural mass (gC per gram of structural mass)
    structural_N_content: the N content of the structural mass (gN per gram of structural mass)
    leaf_minimal_thickness: the minimal thickness of a leaf (m)
    time_step_in_seconds: the time step (seconds)
    air_temperature_in_degree_Celsius: the mean temperature experienced by the shoots (degree Celsius)
    
    Returns
    -------
    total_actual_C_consumption_rate_by_shoot_growth: the rate of labile C consumption for producing the actual new
                                                     structural mass (including the respiratory cost)
    total_actual_N_consumption_rate_by_shoot_growth: the rate of labile N consumption for producing the actual new
                                                     structural mass
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

    # We initialise a C & N growth consumption rate at the scale of the shoot system (gC or gN per plant per second):
    total_actual_C_consumption_rate_by_shoot_growth = 0.
    total_actual_N_consumption_rate_by_shoot_growth = 0.

    # We cover each organ n of the shoot:
    for vid in g.vertices_iter(scale=1):

        # We define the current organ as "n":
        n = g.node(vid)

        # We reset the rates of actual growth consumption:
        n.actual_C_consumption_rate_by_shoot_growth = 0.
        n.actual_N_consumption_rate_by_shoot_growth = 0.

        # First, we increment the time since the shoot element has been formed:
        n.actual_time_since_formation +=  time_step_in_seconds
        n.thermal_time_since_formation += thermal_time_step_in_seconds

        # If the element does not have any potential growth, the we move to the next element:
        if n.potential_volume_increment <=0.:
            continue
        # If the element does not have any C or N left, it cannot grow:
        if n.labile_C <=0. or n.labile_N <=0.:
            continue

        # COMPUTING THE MAXIMAL VOLUME INCREMENT:
        #########################################

        # The maximal volume increment is defined according to the amount of labile C or N available for growth:
        max_volume_increment_by_C = n.labile_C / (structural_mass_density * structural_C_content / yield_growth)
        max_volume_increment_by_N = n.labile_N / (structural_mass_density * structural_N_content)
        # The final volume is the minimal volume between the potential volume calculated from potential growth
        # and the maximal volumes that can be created from the available C and N:
        final_volume_increment = min(n.potential_volume_increment, max_volume_increment_by_C, max_volume_increment_by_N)

        # COMPUTING THE FINAL DIMENSIONS OF THE SHOOT ELEMENT
        #####################################################

        # CASE 1: The shoot element is an internode
        #------------------------------------------
        if n.label == "Internode":
            # We compute the maximal length increment and radius increment based on the final volume increment:
            possible_length_increment = final_volume_increment / (n.radius ** 2 * pi)
            possible_radius_increment = sqrt( (n.volume + final_volume_increment) / (n.length * pi) ) - n.radius
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
                    # We give priority to elongation and elongate the internode according to the final volume increment:
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
                    new_radius = sqrt((n.volume + remaining_volume_increment) / (n.length * pi))
                    # We make sure that the new radius is not lower than the previous one:
                    if (new_radius - n.radius) < 0.:
                        # print("!!! WATCH OUT: the new computed radius of the internode", vid,
                        #       "was lower than the final radius by",
                        #       abs(new_radius - n.radius) / n.radius * 100, "%!")
                        # In such case, we don't modify the radius and we reduce the final volume increment
                        # to the volume increment by elongation only:
                        final_volume_increment = volume_increment_by_elongation
                    else:
                        n.radius = new_radius

            # Eventually, we recompute the new surface, volume and structural mass of the internode:
            n.surface = n.length * n.radius * 2 * pi
            n.volume = n.length * n.radius**2 * pi
            n.struct_mass = n.volume * structural_mass_density

        # CASE 2: The shoot element is a leaf
        #------------------------------------
        elif n.label == "Leaf":
            surface_increment = final_volume_increment / leaf_minimal_thickness
            n.surface += surface_increment
            if leaf_shape == "diamond_shape":
                # We consider the leaf as a diamond characterized by its length and its width, which is calculated as
                # the product of the length with an allometric coefficient. The surface of the diamond is then defined as:
                # surface = (length / 2.) * (param.leaf_width_to_length_ratio * length) / 2.
                # i.e. surface = length ^ 2 * param.leaf_width_to_length_ratio / 4.
                # Hence:
                length_increment = sqrt(n.surface * 4 / (param.leaf_width_to_length_ratio)) - n.length
                # The length is incremented with the corresponding result:
                n.length += length_increment
                # The radius, which is half of the width, is also calculated based on the allometric coefficient
                new_radius = n.length * param.leaf_width_to_length_ratio / 2.
            else:
                # We assume that the leaf is an ellipse characterized by its length and its width, which is calculated as
                # the product of the length with an allometric coefficient. The surface of the ellipse is then defined as:
                # surface = (length / 2.) * (param.leaf_width_to_length_ratio * length) / 2. * pi
                # i.e. surface = length ^ 2 * param.leaf_width_to_length_ratio / 4. * pi
                # Hence:
                length_increment = sqrt(n.surface * 4 / (param.leaf_width_to_length_ratio * pi)) - n.length
                # The length is incremented with the corresponding result:
                n.length += length_increment
                # The radius, which is half of the width, is also calculated based on the allometric coefficient
                new_radius = n.length * param.leaf_width_to_length_ratio / 2.
            # We make sure that the new radius is not lower than the previous one:
            if (new_radius - n.radius) / n.radius < - 1e-12:
                print("!!! WATCH OUT: the new radius of the leaf", vid,
                      "was lower than the final radius by",
                      abs(n.potential_radius - n.radius) / n.radius * 100, "%!")
                print("!!! We therefore did not modify the radius of this leaf !")
            else:
                # If the new radius is correct, we record it:
                n.radius = new_radius
            # We update the surface of the leaf based on the final length and radius:
            if leaf_shape == "diamond_shape":
                # If the leaf is considered as a diamond:
                n.surface = n.length / 2. * n.radius
            else:
                # Else, we assume that the leaf is shaped as an ellipse:
                n.surface = n.length / 2. * n.radius * pi
            # We update other characteristics based on this growth:
            n.volume += final_volume_increment
            n.struct_mass += final_volume_increment * structural_mass_density

        # In any case, the C and N used for creating the corresponding structural mass is calculated:
        n.actual_C_consumption_rate_by_shoot_growth = (final_volume_increment * structural_mass_density
                                                       * structural_C_content / yield_growth / time_step_in_seconds)
        n.actual_N_consumption_rate_by_shoot_growth = (final_volume_increment * structural_mass_density
                                                       * structural_N_content / time_step_in_seconds)

        # The C respired because of growth is also calculated:
        n.shoot_growth_respiration_rate = n.actual_C_consumption_rate_by_shoot_growth * (1 - yield_growth)
        # And we increase the total consumption rates with the contribution from this element:
        total_actual_C_consumption_rate_by_shoot_growth += n.actual_C_consumption_rate_by_shoot_growth
        total_actual_N_consumption_rate_by_shoot_growth += n.actual_N_consumption_rate_by_shoot_growth

    return total_actual_C_consumption_rate_by_shoot_growth, total_actual_N_consumption_rate_by_shoot_growth
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function for computing the ageing of each shoot organ:
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
    [Nothing is returned; the MTG has been updated with the properties 'actual_time_since_formation' (seconds)
    and 'thermal_time_since_formation' (seconds at reference temperature)]
    """

    # We consider a temperature-modifier of time, in order to calculate the equivalent of thermal time:
    temperature_adjustment = temperature_modification(
        temperature_in_Celsius=air_temperature_in_degree_Celsius,
        process_at_T_ref=1,
        T_ref=param.shoot_growth_T_ref,
        A=param.shoot_growth_A,
        B=param.shoot_growth_B,
        C=param.shoot_growth_C)

    # We cover each shoot organ:
    for vid in g.vertices_iter(scale=1):
        element = g.node(vid)
        # We increase the time since this element has been formed:
        element.actual_time_since_formation += time_step_in_seconds
        element.thermal_time_since_formation += time_step_in_seconds * temperature_adjustment

    return
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# Function for computing shoot ramification:

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
# (sachant que les entrenoeuds courts ne donnent pas de ramifications). A partir du premier entrenoeud qui ramifie,
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

def computing_ramification(shoot_MTG, initial_internode_length=0.0005, initial_leaf_length=0.001,
                           time_step_in_seconds=60.*60.):
    """

    Parameters
    ----------
    shoot_MTG
    initial_internode_length
    initial_leaf_length

    Returns
    -------

    """
    # We rename the shoot_MTG for simplification:
    g = shoot_MTG

    # We get the special property "ramification_allowed" at the base of the shoot system:
    base = g.node(0)
    ramification_allowed = base.ramification_allowed

    # To cover the MTG in a specific order, we need to define its "root":
    root_gen = g.component_roots_at_scale_iter(g.root, scale=1)
    root = next(root_gen)

    # DEFINING THE CONDITIONS OF TURNING DORMANT BUDS INTO ACTIVE BUDS:
    # We cover a first time each organ n of the shoot, from the base to the tip:
    for vid in pre_order(g, root):
        if g.node(vid).label != "Bud":
            continue
        else:
            bud = g.node(vid)
            # If no ramification has occured until now, we check if one bud somewhere in the shoot system meets the
            # conditions for triggering ramification for this bud and the bud above it:
            if not ramification_allowed:
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
                    print("> The first ramification should emerge from bud", vid)
                    # We modify the general property "ramification_allowed" for the shoot system:
                    base.ramification_allowed = True
            else:
                # We define the conditions of bud emergence according to its thermal age only:
                if bud.thermal_time_since_formation >= param.bud_minimal_dormancy_period:
                    bud.emergence_possibility = True

    # TRANSFORMING THE ACTIVE BUDS INTO A PHYTOMER AND PHYTOMER:
    # We initialize the total costs of phytomer production:
    total_C_costs_for_new_phytomers = 0.
    total_N_costs_for_new_phytomers = 0.
    # We select the vertices corresponding to the buds that are allowed to emerge:
    list_of_apices = [g.node(vid) for vid in g.vertices_iter(scale=1) if g.node(vid).label == "Bud"
                      and g.node(vid).emergence_possibility == True]
    # As long as the list of emerging buds to treat is not empty:
    while list_of_apices != []:
        # We cover each bud present in the list:
        for apex in list_of_apices:
            # print("For the bud", apex.index(),"the initial angle down is", apex.angle_down)
            # We add a new phytomer at the top of the current apex:
            [internode, leaf, bud] = ADDING_A_PHYTOMER(apex,
                                                       internode_radius=apex.radius,
                                                       internode_length=initial_internode_length,
                                                       leaf_length=initial_leaf_length,
                                                       bud_length=0.,
                                                       properties_identical_to_parent=False,
                                                       nil_properties=True
                                                       )
            print("> A new ramification has emerged from bud", apex.index())
            # We correct the angle down of the new internode, so that its direction corresponds to the one from the initial bud:
            internode.angle_down = apex.angle_down
            # print("For the new internode created from the bud", apex.index(), "the angle down is", internode.angle_down)
            # By default, we set the bud as dormant:
            bud.emergence_possibility = False
            # We remove the previous apex, as a new one will now be formed on the top of the new phytomer:
            g.remove_vertex(apex.index(), reparent_child=True)
            # And we add a new apex at the top of the new internode, of length 0:
            new_apex = ADDING_A_CHILD_ORGAN(parent_element=internode,
                                            edge_type="<",
                                            label="Apex",
                                            properties_identical_to_parent=False,
                                            nil_properties=True)
            new_apex.radius = internode.radius
            # We remove the former apex from the list, so that it is not considered next time:
            list_of_apices.remove(apex)
            # TODO: WATCH OUT BELOW, IN CASE YOU WANT TO CONTINUE ADDING PHYTOMERS
            # # Finally, we update the list of apices by adding the new apex that has been formed,
            # # in case it could itself form a new phytomer:
            # list_of_apices.append(new_apex)

            # REGISTERING THE COST OF THE FORMATION FOR THE NEW PHYTOMER:
            # We (re)initialize the structural mass produced:
            new_structural_mass = 0.
            # For each new element:
            for element in [internode, leaf, bud]:
                # We add the structural mass of each new element:
                new_structural_mass += element.struct_mass
                # We register the cost of C and N consumption within each new element:
                element.actual_C_consumption_rate_by_shoot_growth \
                    = element.struct_mass * param.structural_C_content / param.yield_growth / time_step_in_seconds
                element.actual_N_consumption_rate_by_shoot_growth \
                    = element.struct_mass * param.structural_N_content / time_step_in_seconds
            # We calculate the corresponding C (gC) and N (gN) cost for the plant:
            C_cost = new_structural_mass * param.structural_C_content / param.yield_growth
            N_cost = new_structural_mass * param.structural_N_content
            # We increase the total costs of C&N with the contribution of the current phytomer:
            total_C_costs_for_new_phytomers += C_cost
            total_N_costs_for_new_phytomers += N_cost

    return total_C_costs_for_new_phytomers, total_N_costs_for_new_phytomers
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function for computing shoot organ senescence:
def computing_senescence(shoot_MTG,
                         time_step_in_seconds=1/(60.*60.*24.),
                         air_temperature_in_degree_Celsius=10.):

    """
    This function defines whether each shoot organ is in senescence state or not.
    Parameters
    ----------
    shoot_MTG: the MTG containing each shoot organ
    time_step_in_seconds: the time step of the simulation (seconds)
    air_temperature_in_degree_Celsius: the average temperature experienced by shoot organs

    Returns
    -------
    [The MTG has been updated with the property "senescence_status"]
    """

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
            if (n.labile_C_deficit > 0.
                    and (n.thermal_time_since_formation + thermal_time_step_in_seconds) >= param.thermal_youth_duration):
                # We declare that the element should enter senescence:
                n.senescence_status = True

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function for computing leaf fall:
def computing_leaf_fall(shoot_system,
                        shoot_MTG,
                        time_step_in_seconds=1/(60.*60.*24.),
                        air_temperature_in_degree_Celsius=10.):
    """
    The function defines whether a leaf has fallen or not. The fall is activated when the leaf is senescent and when
    its labile C and N content are below a criticial threshold.

    Parameters
    ----------
    shoot_system: the object corresponding to the shoot system (not the MTG!)
    shoot_MTG: the MTG corresponding to the shoot system
    time_step_in_seconds: the time step of the simulation (seconds)
    air_temperature_in_degree_Celsius: the average temperature experienced by the shoot organs

    Returns
    -------
    [The MTG has been updated with a special label for fallen leaves; the shoot system has been updated with the rate
    of C loss caused by leaf fall and the shoot structural necromass]
    """

    # We rename the shoot MTG for simplicity:
    g=shoot_MTG

    # We get the list of all living leaves in the shoot MTG:
    list_of_leaves = [g.node(vid) for vid in g.vertices_iter(scale=1) if g.node(vid).label == "Leaf"]

    # We cover each leaf:
    for leaf in list_of_leaves:
        # If the leaf has to fall:
        if (leaf.senescence_status == True
                and (leaf.labile_C / leaf.struct_mass) <= param.labile_C_content_treshold_for_leaf_fall
                and (leaf.labile_N / leaf.struct_mass) <= param.labile_N_content_treshold_for_leaf_fall):
            # We assign a new label to the fallen leaf:
            leaf.label="Fallen_leaf"
            # And we integrate its fall in the general C and N information about the shoot system:
            shoot_system.shoot_struct_necromass += leaf.struct_mass
            shoot_system.struct_C_loss_rate_by_leaf_fall += leaf.struct_mass * param.structural_C_content
            shoot_system.labile_C_loss_rate_by_leaf_fall += (leaf.labile_C - leaf.labile_C_deficit)
            shoot_system.struct_N_loss_rate_by_leaf_fall += leaf.struct_mass * param.structural_N_content
            shoot_system.labile_N_loss_rate_by_leaf_fall += (leaf.labile_N - leaf.labile_N_deficit)

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# # Function for computing the local C & N balance in each shoot organ:
# def shoot_organ_CN_balance(shoot_MTG,
#                            time_step_in_seconds = 60.*60.*24.,
#                            resetting_rates=False,
#                            resetting_local_amounts=False):
#     """
#     The function performs the C and N balance in each shoot organ, and updates their labile C and N content from this.
#     It also computes the final amount of available C and N in the shoot system based on this balance. Finally, it can
#     reset rates and labile C and N pool in each organ, if required.
#
#     Parameters
#     ----------
#     shoot_MTG: the MTG corresponding to the shoot system
#     time_step_in_seconds: the time step of the simulation (seconds)
#     resetting_rates: if True, all rates affecting C and N in each organ will be reset to 0 at the end
#     resetting_local_amounts: if True, labile C and N in each organ will be reset to 0 at the end
#
#     Returns
#     -------
#     shoot_available_C: the resulting amount of C available in the whole shoot system when taking into account
#                        the balance in each organ (gC per plant)
#     shoot_available_N: the resulting amount of N available in the whole shoot system when taking into account
#                        the balance in each organ (gN per plant)
#     """
#
#     # We rename the MTG for simplification:
#     g = shoot_MTG
#
#     # We initialize the amount of available C & N in the shoot:
#     shoot_available_C = 0.
#     shoot_available_N = 0.
#
#     # We cover each organ of the shoot:
#     for vid in g.vertices_iter(scale=1):
#
#         # We define the current element as n:
#         n = g.node(vid)
#         # We make sure that the current element has a positive length:
#         if n.length<=0.:
#             continue
#
#         # COMPUTING THE NEW RESULTING NET INCREASE RATES IN C & N:
#         # We compute the net rate of C increase in the labile C pool associated to the organ (gC per node per second):
#         net_C_increase_rate = (n.photosynthesized_C_rate
#                                + n.shoot_C_distribution_rate
#                                + n.C_input_rate_for_sink
#                                # TODO: should C_storage_rate be included in the balance of the structural mass or here?
#                                - n.shoot_maintenance_respiration_rate
#                                - n.actual_C_consumption_rate_by_shoot_growth)
#         # We compute the net rate of N increase in the labile N pool associated to the organ:
#         net_N_increase_rate = (n.shoot_N_distribution_rate
#                                - n.actual_N_consumption_rate_by_shoot_growth
#                                - n.shoot_N_mobilization_rate)
#
#         # COMPUTING THE NEW RESULTING NET GAIN IN LABILE C:
#         n.labile_C += net_C_increase_rate * time_step_in_seconds - n.labile_C_deficit
#         if n.labile_C < 0.:
#             print("    !!! Warning: for the shoot element", n.index(), "a deficit in labile C was obtained! We therefore set this pool to 0 gC and will record its deficit for the next time!")
#             n.labile_C_deficit = - n.labile_C
#             n.labile_C = 0.
#         else:
#             n.labile_C_deficit = 0.
#
#         # COMPUTING THE NEW RESULTING NET GAIN IN LABILE N:
#         n.labile_N += net_N_increase_rate * time_step_in_seconds - n.labile_N_deficit
#         if n.labile_N < 0.:
#             print("    !!! Warning: for the shoot element", n.index(), "a deficit in labile N was obtained! We therefore set this pool to 0 gN and will record its deficit for the next time!")
#             n.labile_N_deficit = - n.labile_N
#             n.labile_N = 0.
#         else:
#             n.labile_N_deficit = 0.
#
#         # We add the local amounts (or deficits) of labile C & N to the total amount of available C & N in the shoot:
#         shoot_available_C += (n.labile_C - n.labile_C_deficit)
#         shoot_available_N += (n.labile_N - n.labile_N_deficit)
#
#         # If necessary, we reset the local amounts to 0:
#         if resetting_local_amounts:
#             n.labile_C = 0.
#             n.labile_C_deficit = 0.
#             n.labile_N = 0.
#             n.labile_N_deficit = 0.
#
#         # If necessary, we reset the local exchange rates to 0:
#         if resetting_rates:
#             n.photosynthesized_C_rate = 0.
#             n.C_input_rate= 0.
#             n.shoot_maintenance_respiration_rate = 0.
#             n.actual_C_consumption_rate_by_shoot_growth = 0.
#             n.N_input_rate = 0.
#             n.actual_N_consumption_rate_by_shoot_growth = 0.
#             n.shoot_N_mobilization_rate = 0.
#
#     return shoot_available_C, shoot_available_N
# #-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function for computing C&N balance both at the shoot organ and at the shoot system scale:
def shoot_CN_balance_at_organ_and_system_scales(shoot_system, shoot_MTG, time_step_in_seconds = 60.*60,
                                                resetting_local_amounts=False, resetting_local_rates=False):
    """
    ...

    Parameters
    ----------
    shoot_system: the shoot system object (not the MTG)
    shoot_MTG: the shoot MTG containing each shoot organ
    time_step_in_seconds: the time step of the simulation (seconds)

    Returns
    -------

    """

    # We rename the shoot system as "s" and the shoot MTG as "g":
    s = shoot_system
    g = shoot_MTG

    # INITIALIZATION OF SHOOT SYSTEM LABILE POOLS AND RATES:
    shoot_labile_C = 0.
    shoot_labile_N = 0.
    shoot_structural_mass = 0.
    total_photosynthesized_C_rate = 0.
    total_C_input_rate = 0.
    total_shoot_maintenance_respiration_rate = 0.
    total_actual_C_consumption_rate_by_shoot_growth = 0.
    total_N_input_rate = 0.
    total_actual_N_consumption_rate_by_shoot_growth = 0.
    total_shoot_N_mobilization_rate = 0.

    # C&N BALANCE IN EACH ORGAN:
    ############################

    # We cover each organ of the shoot MTG:
    for vid in g.vertices_iter(scale=1):

        # We define the current element as n:
        n = g.node(vid)
        # We make sure that the current element has a positive length:
        if n.length <= 0.:
            continue

        # COMPUTING THE NEW RESULTING NET INCREASE RATE IN C & N IN THIS ORGAN
        # (WITHOUT CONSIDERING REDISTRIBUTION PROCESSES!):
        # We compute the net rate of C increase in the labile C pool associated to the organ (gC per node per second):
        net_C_increase_rate = (
                n.photosynthesized_C_rate
                #TODO: should C_storage_rate be included in the balance of the structural mass or here?
                - n.shoot_maintenance_respiration_rate
                - n.actual_C_consumption_rate_by_shoot_growth
        )
        # We compute the net rate of N increase in the labile N pool associated to the organ:
        net_N_increase_rate = - n.actual_N_consumption_rate_by_shoot_growth

        # COMPUTING THE NEW RESULTING NET GAIN IN LABILE C:
        n.labile_C += net_C_increase_rate * time_step_in_seconds
        if n.labile_C < 0.:
            print("    !!! Warning: for the shoot element", n.index(),
                  "a deficit in labile C was obtained! We therefore set this pool to 0 gC and will record its deficit for the next time!")
            n.labile_C_deficit += - n.labile_C
            n.labile_C = 0.
        else:
            n.labile_C_deficit = 0.

        # COMPUTING THE NEW RESULTING NET GAIN IN LABILE N:
        n.labile_N += net_N_increase_rate * time_step_in_seconds
        if n.labile_N < 0.:
            print("    !!! Warning: for the shoot element", n.index(),
                  "a deficit in labile N was obtained! We therefore set this pool to 0 gN and will record its deficit for the next time!")
            n.labile_N_deficit += - n.labile_N
            n.labile_N = 0.
        else:
            n.labile_N_deficit = 0.

        # INCREMENTING THE SHOOT SYSTEM PROPERTIES:
        # We add the local amounts (or deficits) of labile C & N to the total amount of available C & N in the shoot:
        shoot_labile_C += (n.labile_C - n.labile_C_deficit)
        shoot_labile_N += (n.labile_N - n.labile_N_deficit)
        shoot_structural_mass += n.struct_mass
        # We also increment the rates at the scale of the shoot system:
        total_photosynthesized_C_rate += n.photosynthesized_C_rate
        total_shoot_maintenance_respiration_rate += n.shoot_maintenance_respiration_rate
        total_actual_C_consumption_rate_by_shoot_growth += n.actual_C_consumption_rate_by_shoot_growth
        total_actual_N_consumption_rate_by_shoot_growth += n.actual_N_consumption_rate_by_shoot_growth

        # If necessary, we reset the local amounts to 0:
        if resetting_local_amounts:
            n.labile_C = 0.
            n.labile_C_deficit = 0.
            n.labile_N = 0.
            n.labile_N_deficit = 0.

        # If necessary, we reset the local exchange rates to 0:
        if resetting_local_rates:
            n.photosynthesized_C_rate = 0.
            n.shoot_maintenance_respiration_rate = 0.
            n.actual_C_consumption_rate_by_shoot_growth = 0.
            n.actual_N_consumption_rate_by_shoot_growth = 0.

    # UPDATING LABILE POOLS AND RATES IN THE SHOOT SYSTEM:
    ######################################################

    # Considering an additional supply of C and N from the germinating mother seed:
    shoot_labile_C += s.C_supply_from_seed_rate * time_step_in_seconds
    shoot_labile_N += s.N_supply_from_seed_rate * time_step_in_seconds

    # If the net shoot labile C of the whole shoot system is positive:
    if (shoot_labile_C - s.shoot_labile_C_deficit) >= 0.:
        # Then we assign it to the shoot labile C and reset the deficit to 0:
        s.shoot_labile_C_pool = (shoot_labile_C - s.shoot_labile_C_deficit)
        s.shoot_labile_C_deficit = 0.
    else:
        # Otherwise, we record the deficit in labile C:
        s.shoot_labile_C_pool = 0.
        s.shoot_labile_C_deficit = -(shoot_labile_C - s.shoot_labile_C_deficit)

    # If the net shoot labile N of the whole shoot system is positive:
    if (shoot_labile_N - s.shoot_labile_N_deficit) >= 0.:
        # Then we assign it to the shoot labile C and reset the deficit to 0:
        s.shoot_labile_N_pool = (shoot_labile_N - s.shoot_labile_N_deficit)
        s.shoot_labile_N_deficit = 0.
    else:
        # Otherwise, we record the deficit in labile C:
        s.shoot_labile_N_pool = 0.
        s.shoot_labile_N_deficit = -(shoot_labile_N - s.shoot_labile_N_deficit)

    # We update the shoot structural mass as well:
    s.shoot_struct_mass = shoot_structural_mass

    # We also record the summed rates:
    s.photosynthesized_C_rate = total_photosynthesized_C_rate
    s.shoot_maintenance_respiration_rate = total_shoot_maintenance_respiration_rate
    # #TODO: WATCH OUT - we do not compute plants' actual growth costs since they are calculated elsewhere
    #  based on the call of each growth function:
    # s.actual_C_consumption_rate_by_shoot_growth = total_actual_C_consumption_rate_by_shoot_growth
    # s.actual_N_consumption_rate_by_shoot_growth = total_actual_N_consumption_rate_by_shoot_growth

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function for resetting rates and other temporary variables in the shoot MTG to 0:
def shoot_MTG_reinitializing(shoot_g, reinitializing_rates=True, reinitializing_amounts=False):
    """
    This function resets all temporary variables to 0 in each element of the MTG.

    Parameters
    ----------
    shoot_g: the shoot MTG to be reinitialized

    Returns
    -------
    [Nothing - the shoot MTG has been updated]
    """

    # We define an empty dictionary where 0. will be assigned to each vid of the shoot MTG:
    empty_dict = {vid: 0 for vid in shoot_g.vertices()}

    if reinitializing_rates:
        # We reset each variable that is normally modified at each time step:
        shoot_g.properties()["potential_radius"].update(empty_dict)
        shoot_g.properties()["potential_length"].update(empty_dict)
        shoot_g.properties()["photosynthesized_C_rate"].update(empty_dict)
        shoot_g.properties()["shoot_maintenance_respiration_rate"].update(empty_dict)
        shoot_g.properties()["C_input_rate_for_sink"].update(empty_dict)
        shoot_g.properties()["N_input_rate"].update(empty_dict)
        shoot_g.properties()["shoot_N_mobilization_rate"].update(empty_dict)
        # shoot_g.properties()["net_shoot_N_mobilization_rate"].update(empty_dict)
        shoot_g.properties()["actual_C_consumption_rate_by_shoot_growth"].update(empty_dict)
        shoot_g.properties()["actual_N_consumption_rate_by_shoot_growth"].update(empty_dict)
        shoot_g.properties()["shoot_C_distribution_rate"].update(empty_dict)
        shoot_g.properties()["shoot_N_distribution_rate"].update(empty_dict)
        shoot_g.properties()["PAR_intercepted"].update(empty_dict)
        shoot_g.properties()["surfacic_PAR_intercepted"].update(empty_dict)
        shoot_g.properties()["Caribu_Eabs"].update(empty_dict)
        shoot_g.properties()["Caribu_Ei"].update(empty_dict)
        shoot_g.properties()["Caribu_area"].update(empty_dict)
        shoot_g.properties()["photosynthetic_area"].update(empty_dict)
        # WATCH OUT: we don't reinitialize light intensity share as it needs to remain positive,
        # even if there is no light at the current time step!
        # shoot_g.properties()["light_intensity_share"].update(empty_dict)

    if reinitializing_amounts:
        shoot_g.properties()["labile_C"].update(empty_dict)
        shoot_g.properties()["labile_C_deficit"].update(empty_dict)
        shoot_g.properties()["labile_N"].update(empty_dict)
        shoot_g.properties()["labile_N_deficit"].update(empty_dict)

    return
#-----------------------------------------------------------------------------------------------------------------------

########################################################################################################################
########################################################################################################################

# MAIN SIMULATION:
##################

if __name__ == '__main__':

    # We create a reference MTG:
    shoot_g = create_shoot_MTG(number_of_phytomers=3,
                               initial_internode_height=0.001,
                               initial_internode_radius=0.0006,
                               initial_leaf_length=0.003)

    # PLOTTING ONE SINGLE PLANT:
    x_cam = 0.1
    y_cam = 0.1
    z_cam = 0.
    zoom_factor = 100.
    plotting_shoot_MTG(shoot_g,
                       plant_species="rapeseed",
                       # plant_species="faba_bean",
                       # leaf_shape="diamond_shape",
                       leaf_shape="rapeseed_leaf_shape",
                       # leaf_shape="faba_bean_leaf_shape",
                       petiole_showing=False,
                       width=1200, height=900,
                       upscaling_factor=zoom_factor,
                       x_center = 0. * zoom_factor,
                       y_center = 0. * zoom_factor,
                       z_center = 0. * zoom_factor,
                       x_cam=x_cam * zoom_factor, y_cam=y_cam * zoom_factor, z_cam=z_cam * zoom_factor,
                       grid=False,
                       # background_color=[0, 0, 0],
                       # single_color=[50,0,0],
                       property_name='length', cmap='jet', vmin=0, vmax=6e-3, lognorm=False,
                       # property_name='labile_C', cmap='jet', vmin=1e-5, vmax=1e1, lognorm=True,
                       # property_name='PAR_intercepted', cmap='jet', vmin=0, vmax=600, lognorm=False,
                       # property_name='Caribu_Eabs', cmap='jet', vmin=0., vmax=1., lognorm=False,
                       # property_name='z1', cmap='jet', vmin=0., vmax=0.5, lognorm=False,
                       recording_plot=True
                       # plot_name = os.path.join(plot_directory, 'plot_%.5d.png' % step)
                       )
    print("Plot done!")
