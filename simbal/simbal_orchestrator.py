#  -*- coding: utf-8 -*-

"""
    The script 'simbal_orchestrator' is used to run simulations with both shoots and roots.

    :copyright: see AUTHORS.
    :license: see LICENSE for details.
"""

import openalea.mtg
from openalea.mtg import *
from openalea.mtg import turtle as turt
from openalea.mtg.traversal import pre_order, post_order
from openalea.mtg.plantframe import color
import openalea.plantgl.all as pgl

from alinea.caribu.plantgl_adaptor import scene_to_cscene, pgl_to_triangles

import AllocF as AllocF
import RootSimple as RootSimple
import simbal_parameters as param
import openalea.rhizodep.model as rhizodep_param #TODO: Solve the problem of parameters updating for RhizoDep functions...!

from tools import (loading_parameters_from_file, buildDic, formatted_inputs, temperature_modification,
                   recording_MTG_properties, circle_coordinates, updating_coordinates, plotting_whole_plant,
                   get_shoot_visitor, get_root_visitor)
                   # create_a_scene_with_plants_popuplation)

import os
import numpy as np
import pandas as pd
from math import sqrt, pi, floor, trunc, ceil, exp, isnan, cos, sin, acos
from decimal import Decimal

########################################################################################################################
########################################################################################################################

# We define a class that will contain the properties of the interface between shoot and roots (amounts and rates):
class whole_plant(object):
    # We initiate the object:
    def __init__(self):
        self.shoot_struct_mass = 0.
        self.root_struct_mass = 0.

        self.shoot_labile_C_pool = 0. # gC per shoot system
        self.shoot_labile_C_deficit = 0. # gC per shoot system
        self.shoot_labile_N_pool = 0. # gN per shoot system
        self.shoot_labile_N_deficit = 0. # gN per shoot system
        self.root_labile_C_pool = 0.  # gC per root system
        self.root_labile_C_deficit = 0. # gC per root system
        self.root_labile_N_pool = 0.  # gN per root system
        self.root_labile_N_deficit = 0. # gN per root system

        self.net_C_transfer_from_shoot_to_root_rate = 0. # gC per second
        self.net_N_transfer_from_root_to_shoot_rate = 0. # gN per second

        self.transpiration_rate = 0. # gram of water per second
        self.total_PAR_intercepted = 0. # micromoles of photons per plant per second
        self.photosynthesized_C_rate = 0. # gC per second
        self.N_uptake_rate = 0. # gN per second
        self.N2_fixation_rate = 0.  # gN per second

        self.shoot_maintenance_respiration_rate = 0. # gC per second
        self.root_maintenance_respiration_rate = 0. # gC per second

        self.actual_C_consumption_rate_by_shoot_growth = 0. # gC per second
        self.actual_N_consumption_rate_by_shoot_growth = 0. # gN per second
        self.actual_C_consumption_rate_by_root_growth = 0. # gC per second
        self.actual_N_consumption_rate_by_root_growth = 0. # gN per second

        self.net_shoot_C_immobilization_rate = 0. # gC per second
        self.net_shoot_N_immobilization_rate = 0. # gN per second
        self.net_root_C_immobilization_rate = 0. # gC per second
        self.net_root_N_immobilization_rate = 0. # gN per second

        self.C_loss_rate_by_leaf_fall = 0. # gC per second
        self.N_loss_rate_by_leaf_fall = 0.  # gN per second
        self.C_rhizodeposition_rate = 0. # gC per second
        self.N_rhizodeposition_rate = 0. # gN per second
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# We define a function that computes C & N exchanges between shoot and root systems:
def computing_shoot_root_exchange_rates(plant, computing_C_exchange=False, computing_N_exchange=False,
                                        time_step_in_seconds=1*(60*60*24),
                                        C_option="by_diffusion", N_option="by_diffusion",
                                        water_concentration = 10.):

    """

    Parameters
    ----------
    plant: plant object containing the variables at the scale of the whole plant
    computing_C_exchange
    computing_N_exchange
    time_step_in_seconds
    C_option: option to transfer C from shoots to roots - either "by_diffusion" or "full_transfer_from_shoot_to_root"
    N_option
    water_concentration

    Returns
    -------

    """

    # COMPUTING THE EXCHANGE OF CARBON:
    ###################################

    if computing_C_exchange:

        # We first compute the net amount of labile C in shoots and roots, taking into account possible deficits:
        net_shoot_C = (plant.shoot_labile_C_pool - plant.shoot_labile_C_deficit)
        net_root_C = (plant.root_labile_C_pool - plant.root_labile_C_deficit)
        # To avoid problems, we do not allow the net amounts to be negative because of deficits:
        if net_shoot_C < 0.:
            net_shoot_C = 0.
        if net_root_C < 0.:
            net_root_C = 0.

        # If C is to be exchanged by diffusion based on the labile pools in shoots and roots:
        if C_option == "by_diffusion":
            # Then the net C exchange rate is calculated using a gradient of concentration, a diffusion coefficient,
            # and the section of exchange between shoot and root:
            plant.exchange_section = 0.01 # TODO: the exchange section between shoots and roots should be calculated somewhere else!
            if plant.shoot_struct_mass > 0. and plant.root_struct_mass > 0.:
                net_C_transfer_from_shoot_to_root_rate = ((param.C_diffusion_coefficient
                                                          * (net_shoot_C / plant.shoot_struct_mass
                                                             - net_root_C / plant.root_struct_mass))
                                                          * plant.exchange_section)
            else:
                print("!!! ERROR: the structural mass of shoot and root are", plant.shoot_struct_mass, "and", plant.root_struct_mass, "gDW, respectively!")
                print("           > We set the net flow rate of C from shoot to root to 0.")
                net_C_transfer_from_shoot_to_root_rate = 0.
        # Otherwise, if the C present in the shoot is meant to be transfer in totality to the root system:
        elif C_option == "full_transfer_from_shoot_to_root":
            net_C_transfer_from_shoot_to_root_rate = net_shoot_C / time_step_in_seconds
        else:
            print("WATCH OUT: the option for computing the transfer of C from shoot to root is not recognized! "
                  "The transfer will correspond to the full amount of labile C in the source.")
            net_C_transfer_from_shoot_to_root_rate = net_shoot_C / time_step_in_seconds

        # We now verify that the net amount transfered would not be higher than the available amount in the source:
        expected_C_exchange = net_C_transfer_from_shoot_to_root_rate * time_step_in_seconds

        # CASE 1: The expected C transfered from shoot to root is positive, but there is not enough C in the shoot:
        if expected_C_exchange > net_shoot_C:
            # Then the amount of C exchanged is limited to the amount of C available in shoots:
            net_final_C_exchange = net_shoot_C
        # CASE 2: The flow from root to shoot is positive but there is not enough C in the root to sustain it:
        elif expected_C_exchange < 0. and abs(expected_C_exchange) > net_root_C:
            # Then the amount of C exchanged is limited to the amount of C available in roots:
            net_final_C_exchange = - net_root_C
        # CASE 3: There is enough C in the source compartment to sustain the transfer, so the rate is unaltered:
        else:
            net_final_C_exchange = expected_C_exchange
        # We register the final transfer rate:
        plant.net_C_transfer_from_shoot_to_root_rate = net_final_C_exchange / time_step_in_seconds

    # COMPUTING THE EXCHANGE OF NITROGEN:
    #####################################

    if computing_N_exchange:

        # We first compute the net amount of labile N in shoots and roots, taking into account possible deficits:
        net_shoot_N = (plant.shoot_labile_N_pool - plant.shoot_labile_N_deficit)
        net_root_N = (plant.root_labile_N_pool - plant.root_labile_N_deficit)
        # To avoid problems, we do not allow the net amounts to be negative because of deficits:
        if net_shoot_N < 0.:
            net_shoot_N = 0.
        if net_root_N < 0.:
            net_root_N = 0.

        # If N is to be exchanged by diffusion based on the labile pools in shoots and roots:
        if N_option == "by_diffusion":
            # Then the net N exchange rate is calculated using a gradient of concentration, a diffusion coefficient,
            # and the section of exchange between shoot and root:
            plant.exchange_section = 0.01 # TODO: the exchange section between shoots and roots should be calculated somewhere else!
            if plant.shoot_struct_mass > 0. and plant.root_struct_mass > 0.:
                net_N_transfer_from_root_to_shoot_rate = ((param.C_diffusion_coefficient
                                                          * (net_root_N / plant.root_struct_mass
                                                             - net_shoot_N / plant.shoot_struct_mass))
                                                          * plant.exchange_section)
            else:
                print("!!! ERROR: the structural mass of shoot and root are", plant.shoot_struct_mass, "and", plant.root_struct_mass, "gDW, respectively!")
                print("           > We set the net flow rate of N from root to shoot to 0.")
                net_N_transfer_from_root_to_shoot_rate = 0.
        # Otherwise, if the nitrogen is meant to be transferred to shoots along with the transpired water,
        # a transpiration flux has to be provided, from which the N flow is calculated based on the concentration
        # of N in the root labile pool:
        elif N_option == "by_water_flux":
            if plant.transpiration_rate >= 0.:
                net_N_transfer_from_root_to_shoot_rate = (plant.transpiration_rate / water_concentration
                                                            * net_root_N / plant.root_struct_mass)
                # NOTE: N_transfer_rate in gN per plant and per second
                #       transpiration_rate in gram of water per plant and per second
                #       water_concentration in gram of water per gDW of root structural mass
                #       net_root_N in gN per plant
                #       root_struct_mass in gDW of root structural mass
            else:
                print("ERROR - the transpiration flux was negative! The flow rate of N from root to shoot was set to 0.")
                net_N_transfer_from_root_to_shoot_rate = 0.
        # Otherwise, if the N present in the root is meant to be transfered in totality to the shoot system:
        elif C_option == "full_transfer_from_root_to_shoot":
            net_N_transfer_from_root_to_shoot_rate = net_root_C / time_step_in_seconds
        else:
            print("WATCH OUT: the option for computing the transfer of N from root to shoot is not recognized! "
                  "The transfer will correspond to the full amount of labile N in the source.")
            net_N_transfer_from_root_to_shoot_rate = net_root_N / time_step_in_seconds

        # We now verify that the net amount transfered would not be higher than the available amount in the source:
        expected_N_exchange = net_N_transfer_from_root_to_shoot_rate * time_step_in_seconds

        # CASE 1: The expected N transfered from root to shoot is positive, but there is not enough N in the root:
        if expected_N_exchange > net_root_N:
            # Then the amount of C exchanged is limited to the amount of N available in roots:
            net_final_N_exchange = net_root_N
        # CASE 2: The flow from root to shoot is positive but there is not enough C in the root to sustain it:
        elif expected_N_exchange < 0. and abs(expected_N_exchange) > net_shoot_N:
            # Then the amount of N exchanged is limited to the amount of N available in shoots:
            net_final_N_exchange = - net_shoot_N
        # CASE 3: There is enough N in the source compartment to sustain the transfer, so the rate is unaltered:
        else:
            net_final_N_exchange = expected_N_exchange
        # We register the final transfer rate:
        plant.net_N_transfer_from_root_to_shoot_rate = net_final_N_exchange / time_step_in_seconds

    return

# # EXAMPLE:
#
# i = whole_plant()
# i.shoot_struct_mass = 2.
# i.root_struct_mass = 0.5
#
# i.shoot_labile_C_pool = 5.
# i.shoot_labile_N_pool = 0.01
# i.root_labile_C_pool = 0.5
# i.root_labile_N_pool = 0.02
#
# i.transpiration_rate = 0.1
#
# shoot_root_exchange(i, computing_C_exchange=True, computing_N_exchange=True,
#                     C_option="by_diffusion", N_option="by_water_flux")
# print("The new flux of C from shoot to root is", i.C_flux_to_root, "gC per second.")
# print("The new flux of N from root to shoot is", i.N_flux_to_shoot, "gN per second.")
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

def recording_main_plant_properties(plant, time_in_days):
    results_dict = {"Time (days)": time_in_days,
                    "Shoot structural mass (gDW)": plant.shoot_struct_mass,
                    "Root structural mass (gDW)": plant.root_struct_mass,
                    "Shoot labile C pool (gC per plant)": plant.shoot_labile_C_pool,
                    "Shoot labile C pool deficit (gC per plant)": plant.shoot_labile_C_deficit,
                    "Shoot labile N pool (gN per plant)": plant.shoot_labile_N_pool,
                    "Shoot labile N pool deficit (gN per plant)": plant.shoot_labile_N_deficit,
                    "Shoot reserve C pool (gC per plant)": plant.shoot_reserve_C_pool,
                    "Shoot reserve N pool (gN per plant)": plant.shoot_reserve_N_pool,
                    "Root labile C pool (gC per plant)": plant.root_labile_C_pool,
                    "Root labile C pool deficit (gC per plant)": plant.root_labile_C_deficit,
                    "Root labile N pool (gN per plant)": plant.root_labile_N_pool,
                    "Root labile N pool deficit (gN per plant)": plant.root_labile_N_deficit,
                    "Root reserve C pool (gC per plant)": plant.root_reserve_C_pool,
                    "Root reserve N pool (gN per plant)": plant.root_reserve_N_pool,
                    "Rate of net shoot-to-root C transfer (gC per plant per second)": plant.net_C_transfer_from_shoot_to_root_rate,
                    "Rate of net root-to-shoot N transfer (gN per plant per second)": plant.net_N_transfer_from_root_to_shoot_rate,
                    "Rate of water transpiration (gram per second)": plant.transpiration_rate,
                    "Rate of light interception (micromol of photons per plant per second)": plant.total_PAR_intercepted,
                    "Rate of photosynthesis (gC per plant per second)": plant.photosynthesized_C_rate,
                    "Rate of shoot maintenance respiration (gC per plant per second)": plant.shoot_maintenance_respiration_rate,
                    "Rate of C transformation into shoot structural mass (gC per plant per second)": plant.actual_C_consumption_rate_by_shoot_growth * param.yield_growth,
                    "Rate of C respiration by shoot growth (gC per plant per second)": plant.actual_C_consumption_rate_by_shoot_growth * (
                                1 - param.yield_growth),
                    "Rate of N transformation into shoot structural mass (gN per plant per second)": plant.actual_N_consumption_rate_by_shoot_growth,
                    "Rate of shoot net C immobilization (gC per plant per second)": plant.net_shoot_C_immobilization_rate,
                    "Rate of shoot net N immobilization (gN per plant per second)": plant.net_shoot_N_immobilization_rate,
                    "Rate of shoot C loss by leaf fall (gC per plant per second)": plant.C_loss_rate_by_leaf_fall,
                    "Rate of shoot N loss by leaf fall (gC per plant per second)": plant.N_loss_rate_by_leaf_fall,
                    "Rate of N uptake by roots (gN per plant per second)": plant.N_uptake_rate,
                    "Rate of N2 fixation by roots (gN per plant per second)": plant.N2_fixation_rate,
                    "Rate of root maintenance respiration (gC per plant per second)": plant.root_maintenance_respiration_rate,
                    "Rate of C rhizodeposition (gC per plant per second)": plant.C_rhizodeposition_rate,
                    "Rate of N rhizodeposition (gC per plant per second)": plant.N_rhizodeposition_rate,
                    "Rate of C transformation into root structural mass (gC per plant per second)": plant.actual_C_consumption_rate_by_root_growth * param.yield_growth,
                    "Rate of C respiration by root growth (gC per plant per second)": plant.actual_C_consumption_rate_by_root_growth * (
                                1 - param.yield_growth),
                    "Rate of N transformation into root structural mass (gN per plant per second)": plant.actual_N_consumption_rate_by_root_growth,
                    "Rate of root net C immobilization (gC per plant per second)": plant.net_root_C_immobilization_rate,
                    "Rate of root net N immobilization (gN per plant per second)": plant.net_root_N_immobilization_rate
                    }
    return results_dict


# Function for calling the main functions in a specific order and run the whole simulation:
def run_simulation(my_path=None,
                   simulation_period_in_days=20., time_step_in_days=1/24., starting_time_in_days=0.,
                   input_file=None,
                   input_file_time_step_in_days=1.,
                   forcing_constant_aboveground_environment=False,
                   constant_air_temperature_in_Celsius=10, constant_light_irradiation=400,
                   forcing_constant_belowground_environment=False,
                   constant_soil_temperature_in_Celsius=10, constant_soil_available_nitrogen=1,
                   plant_density_in_square_meter=35,
                   inter_row_distance = 0.45,
                   row_length=1.,
                   scene_width=1.,
                   C_allocation_option = "No_priority_to_leaves",
                   simulation_option="shoot_and_roots",
                   initial_shoot_labile_C_concentration = 1e-3, # in gC per gram of structural mass
                   initial_shoot_labile_N_concentration = 1e-4, # in gN per gram of structural mass
                   initial_root_labile_C_concentration = 1e-3, # in gC per gram of structural mass
                   initial_root_labile_N_concentration = 1e-2 # in gN per gram of structural mass
                   ):

    """
    This function organizes the main simulation. After plant initiation, it calls the different modules and functions
    from SIMBAL for each time step, and records the outputs.
    Parameters
    ----------
    my_path
    simulation_period_in_days
    time_step_in_days
    starting_time_in_days
    input_file
    input_file_time_step_in_days
    forcing_constant_aboveground_environment
    constant_air_temperature_in_Celsius
    constant_light_irradiation
    forcing_constant_belowground_environment
    constant_soil_temperature_in_Celsius
    constant_soil_available_nitrogen
    plant_density_in_square_meter
    inter_row_distance
    row_length
    scene_width
    C_allocation_option
    initial_shoot_labile_C_concentration
    initial_shoot_labile_N_concentration
    initial_root_labile_C_concentration
    initial_root_labile_N_concentration

    Returns
    -------

    """

    # SETTING FILE PATHS AND FOLDERS:
    #--------------------------------
    if my_path is None:
        # We define the reference directory (here, we move to the parent directory of the current directory):
        my_path = os.path.dirname(os.getcwd())

    # We consider the outputs directory:
    outputs_directory = os.path.join(my_path, "simulation", "outputs")
    # If the plot directory doesn't exist:
    if not os.path.exists(outputs_directory):
        # Then we create it:
        os.mkdir(outputs_directory)
    else:
        # Otherwise, we delete all files that are already present inside:
        for root, dirs, files in os.walk(outputs_directory):
            for file in files:
                os.remove(os.path.join(root, file))

    # We consider the plot directory:
    plot_directory = os.path.join(my_path, "simulation", "outputs","plots")
    # If the plot directory doesn't exist:
    if not os.path.exists(plot_directory):
        # Then we create it:
        os.mkdir(plot_directory)
    else:
        # Otherwise, we delete all the images that are already present inside:
        for root, dirs, files in os.walk(plot_directory):
            for file in files:
                os.remove(os.path.join(root, file))

    # We consider the directory containing the attributes of shoot properties over time:
    shoot_MTG_directory = os.path.join(my_path, "simulation", "outputs", "shoot_MTG_properties")
    # If this directory doesn't exist:
    if not os.path.exists(shoot_MTG_directory):
        # Then we create it:
        os.mkdir(shoot_MTG_directory)
    else:
        # Otherwise, we delete all the images that are already present inside:
        for root, dirs, files in os.walk(shoot_MTG_directory):
            for file in files:
                os.remove(os.path.join(root, file))

    # UPDATING PARAMETERS AND INSTRUCTIONS:
    #--------------------------------------
    # The values of the parameters in the file 'simbal_parameters.py' are updated according to the values in the file
    # 'simbal_parameters.xlsx':
    scenario = loading_parameters_from_file(parameters_filepath = os.path.join(my_path, "simulation", "inputs", "simbal_parameters.xlsx"),
                                            outputs_path = os.path.join(my_path, "simulation", "outputs"))
    # We create a dictionary containing the parameters and instructions specified in the scenario:
    scenario_parameters = buildDic(scenario)
    # We update the parameters from the original python file of SIMBAL with the parameters taken from this dictionnary:
    param.__dict__.update(scenario)

    if simulation_option != "shoot_only":
        # We do the same for the RhizoDep parameters (IS IT NEEDED? NO...):
        rhizodep_param.__dict__.update(scenario)

    # We read the other specific instructions of the scenario that do not correspond to the parameters in the list
    # (Note: The following function "get" looks if the first argument is present in the scenario parameters.
    # If not, it returns the default value indicated in the second argument, which may be the default argument
    # of the function 'run_simulation' itself):
    SIMULATION_OPTION = scenario_parameters.get('simulation_option',
                                                simulation_option)
    SIMULATION_PERIOD_IN_DAYS = scenario_parameters.get('simulation_period_in_days',
                                                        simulation_period_in_days)
    TIME_STEP_IN_DAYS = scenario_parameters.get('time_step_in_days',
                                                time_step_in_days)
    STARTING_TIME_IN_DAYS = scenario_parameters.get('starting_time_in_days',
                                                starting_time_in_days)
    INPUT_FILE = scenario_parameters.get('input_file',
                                         input_file)
    INPUT_FILE_TIME_STEP_IN_DAYS = scenario_parameters.get('input_file_time_step_in_days',
                                                           input_file_time_step_in_days)
    FORCING_CONSTANT_ABOVEGROUND_ENVIRONMENT = scenario_parameters.get('forcing_constant_aboveground_environment',
                                                                       forcing_constant_aboveground_environment)
    CONSTANT_LIGHT_IRRADIATION = scenario_parameters.get('constant_light_irradiation',
                                                         constant_light_irradiation)
    CONSTANT_AIR_TEMPERATURE = scenario_parameters.get('constant_air_temperature_in_Celsius',
                                                       constant_air_temperature_in_Celsius)
    FORCING_CONSTANT_BELOWGROUND_ENVIRONMENT = scenario_parameters.get('forcing_constant_belowground_environment',
                                                                       forcing_constant_belowground_environment)
    CONSTANT_SOIL_TEMPERATURE = scenario_parameters.get('constant_soil_temperature_in_Celsius',
                                                        constant_soil_temperature_in_Celsius)
    CONSTANT_SOIL_AVAILABLE_N = scenario_parameters.get('constant_soil_available_nitrogen',
                                                        constant_soil_available_nitrogen)

    # SETTING INITIAL CONDITIONS
    #---------------------------
    
    # INITIALIZING THE WHOLE PLANT:
    # We create a "plant", i.e. a shoot-root interface without structure that will store variables at the scale of the plant:
    plant = whole_plant()
    
    # INITIALIZING SHOOT STRUCTURE:
    # In case we only want to simulate the behavior of the roots, keeping the shoot properties constant over time:
    if SIMULATION_OPTION == "roots_only":
        # We create a MTG "shoot_g" corresponding to the shoot, without any dimensions (it will not be shown in the graph):
        shoot_g = AllocF.create_shoot_MTG(number_of_phytomers=1,
                                          initial_internode_height=0.,
                                          initial_internode_radius=0.,
                                          initial_leaf_length=0.)
        # We define an arbitrary shoot structural mass, which will remain constant over the simulation:
        plant.shoot_struct_mass = 0.1
    # Otherwise, a fully functional shoot MTG is created:
    else:
        # We create a MTG "shoot_g" corresponding to the shoot (stem and leaves):
        shoot_g = AllocF.create_shoot_MTG(number_of_phytomers=3,
                                          initial_internode_height=0.0005,
                                          initial_internode_radius=0.0008,
                                          initial_leaf_length=0.003)
        # We update the shoot system general properties based on the shoot MTG:
        for vid in shoot_g.vertices_iter(scale=1):
            n = shoot_g.node(vid)
            try:
                plant.shoot_struct_mass += n.struct_mass
            except:
                print("!!! ERROR: we could not add the properties of the shoot element", vid, "to the whole plant!")
        # We add an extra-property to the base of the shoot system, corresponding to the time since the last phytomer
        # was formed (used together with the phyllochron later on):
        shoot_g.node(0).thermal_time_since_last_phytomer_formation = 0.
        # We record the properties in the initial state of the shoot MTG in an external file:
        recording_MTG_properties(shoot_g, os.path.join(shoot_MTG_directory, 'MTG_00000.csv'))

    # INITIALIZING ROOT STRUCTURE:
    # In case we only want to simulate the behavior of the shoot, keeping the root properties constant over time:
    if SIMULATION_OPTION == "shoot_only":
        # We create another MTG corresponding to the root system, without any length (it will not appear on the graph):
        root_g = RootSimple.initiate_mtg(initial_segment_length=0.)
        # We define an arbitrary initial root mass that will remain constant over the simulation:
        plant.root_struct_mass = 0.1
    # Otherwise, a fully functional root MTG is created:
    else:
        # We create another MTG corresponding to the root system:
        root_g = RootSimple.initiate_mtg(initial_segment_length=0.003)
        # We reinitialize the variables in the root MTG:
        RootSimple.reinitializing_growth_variables(root_g)
        # We update the distance from tip for each root element in each root axis:
        RootSimple.update_distance_from_tip(root_g)
        # We update the surfaces and volumes of each root element:
        RootSimple.update_surfaces_and_volumes(root_g)
        # # And, when using RootSimple, we create the final root system "r" from it:
        # r = RootSimple.root_system(root_g)
        # We adjust the properties of the whole root system according to whole root MTG:
        root_g_properties = RootSimple.summing_and_possibly_homogenizing(root_g)
        # We record the total root structural mass within the root system:
        plant.root_struct_mass = root_g_properties["total_living_root_struct_mass"]
        # We calculate the total exchange surface of the root system:
        RootSimple.updating_exchange_surface(plant, root_g)

    # INITIALIZING PLANT'S RESOURCES:
    # We initialize the pools of labile and reserve C (gC per plant) and N (gN per plant) in the whole plant object:
    plant.shoot_labile_C_pool = initial_shoot_labile_C_concentration * plant.shoot_struct_mass
    plant.shoot_labile_N_pool = initial_shoot_labile_N_concentration * plant.shoot_struct_mass
    plant.shoot_reserve_C_pool = 0.
    plant.shoot_reserve_N_pool = 0.
    plant.root_labile_C_pool = initial_root_labile_C_concentration * plant.root_struct_mass
    plant.root_labile_N_pool = initial_root_labile_N_concentration * plant.root_struct_mass
    plant.root_reserve_C_pool = 0.
    plant.root_reserve_N_pool = 0.

    # We also initialize
    results_dict = recording_main_plant_properties(plant, time_in_days=0.)
    # plant_data_frame  = pd.DataFrame([results_dict.values()], columns=results_dict.keys())
    plant_data_frame = pd.DataFrame(results_dict, index=[0])

    # DEFINING THE PERIOD OF SIMULATION:
    # ---------------------------------

    # We define the steps of the simulation based on the desired period:
    initial_time = STARTING_TIME_IN_DAYS
    simulation_period_in_days = SIMULATION_PERIOD_IN_DAYS
    time_step_in_seconds = TIME_STEP_IN_DAYS * 60*60*24
    n_steps = trunc(SIMULATION_PERIOD_IN_DAYS / TIME_STEP_IN_DAYS)
    # We initialize a counter that will be used over successive time steps together with the phyllochron:
    elapsed_time_in_days_since_last_leaf_emergence = 0.

    # READING THE INPUT FILE:
    # -----------------------
    if INPUT_FILE and not (FORCING_CONSTANT_ABOVEGROUND_ENVIRONMENT and FORCING_CONSTANT_BELOWGROUND_ENVIRONMENT):
        # We use the function 'formatted inputs' to create a table containing the input data for each time step:
        input_frame = formatted_inputs(original_input_file=INPUT_FILE,
                                       final_input_file=os.path.join(outputs_directory, 'updated_input_file.csv'),
                                       original_time_step_in_days=INPUT_FILE_TIME_STEP_IN_DAYS,
                                       final_time_step_in_days=TIME_STEP_IN_DAYS,
                                       simulation_period_in_days=SIMULATION_PERIOD_IN_DAYS,
                                       do_not_execute_if_file_with_suitable_size_exists=False)

        # We then initialize the step and time according to the first line of the new inputs dataframe:
        initial_step_number = input_frame['step_number'].loc[0]
        initial_time_in_days = input_frame['initial_time_in_days'].loc[0]
    else:
        # We then initialize the step and time at 0:
        initial_step_number = 0
        initial_time_in_days = 0.
    # In any case, the value of 'step' is initialized:
    step = initial_step_number

    # CONSIDERING CAMERA ROTATION WHEN PLOTTING:
    #-------------------------------------------
    if param.camera_movement is True:
        # In case the camera needs to rotate around the target and/or to move away from it or to move up/down:
        if param.step_back_coefficient != 0. or param.step_up_coefficient != 0.:
            # If the required number of positions around the target is positive:
            if param.camera_rotation_n_points > 0:
                # We calculate the number of loops required around the target to go up to the end of the simulation,
                # knowing that the camera will progressively move away or down at each time step while rotating:
                number_of_loops = int(ceil(n_steps / param.camera_rotation_n_points))
            # Otherwise, the number of coordinates required will be directly the number of steps of the simulation:
            else:
                number_of_loops =  n_steps
        # Otherwise, only one loop is necessary, and the camera will reuse the same coordinates after one complete rotation:
        else:
            number_of_loops = 1
        # We calculate the x,y,z coordinates of the different positions of the camera.
        [x_coordinates, y_coordinates, z_coordinates] = circle_coordinates(x_center=param.x_center,
                                                                           y_center=param.y_center,
                                                                           z_center=param.z_camera,
                                                                           radius=sqrt((param.x_camera**2 + param.y_camera**2)),
                                                                           n_points=param.camera_rotation_n_points,
                                                                           n_loops=number_of_loops,
                                                                           step_back_coefficient=param.step_back_coefficient,
                                                                           step_up_coefficient=param.step_up_coefficient)
        # We initialize the index in the list of the coordinates:
        index_camera = 0

    ####################################################################################################################

    # SIMULATING OVER TIME:
    #######################
    #######################

    # For each time step covering the simulation period:
    for step in range(initial_step_number, initial_step_number + n_steps):

        # We calculate the initial and final time of the current time step (in days):
        initial_time_in_days = initial_time + (step) * time_step_in_days
        final_time_in_days = initial_time + (step+1) * time_step_in_days

        # DEFINING THE INPUTS AT THIS TIME STEP:
        # --------------------------------------
        # If constant aboveground conditions have been chosen:
        if FORCING_CONSTANT_ABOVEGROUND_ENVIRONMENT:
            # We set up the constant values:
            light_input = CONSTANT_LIGHT_IRRADIATION
            air_temperature = CONSTANT_AIR_TEMPERATURE
        # Otherwise, if there is an input file to be used:
        elif input_file:
            # We extract the information on incident light and air temperature from the input file:
            light_input = input_frame.loc[step - initial_step_number, 'light_irradiation_in_micromol_per_m2_per_second']
            air_temperature = input_frame.loc[step - initial_step_number, 'air_temperature_in_degree_Celsius']
        else:
            print("!!! ERROR: no input file has been provided!")

        # Similarly, if constant belowground conditions have been chosen:
        if FORCING_CONSTANT_BELOWGROUND_ENVIRONMENT:
            soil_temperature = CONSTANT_SOIL_TEMPERATURE
            soil_available_nitrogen = CONSTANT_SOIL_AVAILABLE_N
        elif input_file:
            soil_temperature = input_frame.loc[step - initial_step_number, 'soil_temperature_in_degree_Celsius']
            soil_available_nitrogen = input_frame.loc[step - initial_step_number, 'soil_available_nitrogen_in_gN_per_gDW_of_soil']

        print("")
        print("From", "{:.2f}".format(Decimal(initial_time_in_days)),"to", "{:.2f}".format(Decimal(final_time_in_days)),"days:")
        print("-------------------------")
        print("    Air temperature is", "{:.2f}".format(Decimal(air_temperature)), "degree Celsius; soil temperature is",
              "{:.2f}".format(Decimal(soil_temperature)), "degree Celsius.")
        print("    Incoming light is", "{:.0f}".format(Decimal(light_input)), "micromol of photons per square meter per second.")
        print("    Soil available N is", "{:.2E}".format(Decimal(soil_available_nitrogen)), "gN per gram of dry soil.")
        print("")

        # CALCULATING AN EQUIVALENT OF THERMAL TIME:
        # ------------------------------------------
        # We calculate a coefficient that will modify the different "ages" experienced by shoots, roots
        # (and possibly soil organisms) according to temperature:
        temperature_time_adjustment_shoot = temperature_modification(temperature_in_Celsius=air_temperature,
                                                                     process_at_T_ref=1,
                                                                     T_ref=param.shoot_growth_T_ref,
                                                                     A=param.shoot_growth_A,
                                                                     B=param.shoot_growth_B,
                                                                     C=param.shoot_growth_C)
        temperature_time_adjustment_soil = temperature_modification(temperature_in_Celsius=soil_temperature,
                                                                    process_at_T_ref=1,
                                                                    T_ref=param.root_growth_T_ref,
                                                                    A=param.root_growth_A,
                                                                    B=param.root_growth_B,
                                                                    C=param.root_growth_C)
        # We calculate the equivalent time step in thermal time:
        thermal_time_step_in_seconds_shoots = time_step_in_seconds * temperature_time_adjustment_shoot
        # TODO in the future: consider a vertical gradient of temperature, and hence a vertical modification of the adjustments

        # # TRANSFER OF C BETWEEN SHOOTS AND ROOTS:
        # #########################################
        #
        # # We first compute the rate of net C transfer from shoots to roots:
        # computing_shoot_root_exchange_rates(plant, computing_C_exchange=True, computing_N_exchange=False,
        #                                     C_option="by_diffusion")
        # # Note that only the rate has been computed here, the labile C pools have not been updated yet!

        # CALLING SUCCESSIVELY THE DIFFERENT ROOT MODULES:
        ##################################################

        if SIMULATION_OPTION != "shoot_only":

            print("    Root labile C pool is initially",
                  "{:.2E}".format(Decimal(plant.root_labile_C_pool - plant.root_labile_C_deficit)),
                  "gram of C, and labile N pool is",
                  "{:.2E}".format(Decimal(plant.root_labile_N_pool - plant.root_labile_N_deficit)), "gram of N.")

            # N INPUTS TO THE ROOT SYSTEM:
            # We compute the rate of N uptake from soil (gN per plant per second):
            RootSimple.computing_N_uptake_rate(root_system=plant, root_MTG=root_g,
                                               N_soil = soil_available_nitrogen,
                                               soil_temperature_in_degree_Celsius=soil_temperature)
            # We compute the rate of N uptake from atmosphere, i.e. the rate of N2 fixation by root nodules, if any
            # (gN per plant per second):
            RootSimple.computing_N2_fixation_rate(root_system=plant, root_MTG=root_g,
                                                  soil_temperature_in_degree_Celsius=soil_temperature)

            # NET C LOSSES FROM THE ROOT SYSTEM (WITHOUT GROWTH):
            # We compute the rate of rhizodeposition (gC per plant per second):
            RootSimple.computing_rhizodeposition_rate(root_system=plant, root_MTG=root_g,
                                                      time_step_in_seconds=time_step_in_seconds,
                                                      soil_temperature_in_degree_Celsius=soil_temperature)
            # We compute the rate of maintenance respiration (gC per plant per second):
            RootSimple.computing_maintenance_rate(root_system=plant,
                                                  soil_temperature_in_degree_Celsius=soil_temperature)
            # We compute the rate of C exchange with reserve (gC per plant per second):
            RootSimple.computing_exchange_rates_with_reserve(root_system=plant,
                                                             soil_temperature_in_degree_Celsius=soil_temperature)

            # ROOT GROWTH:
            # We reset all growth-associated variables before performing the new growth:
            RootSimple.reinitializing_growth_variables(root_g)
            # We calculate the potential growth of each root element according to ArchiSimple:
            RootSimple.potential_growth(root_g, time_step_in_seconds=time_step_in_seconds,
                                        ArchiSimple=True,
                                        soil_temperature_in_Celsius=soil_temperature)
            # We calculate the "satisfaction coefficients" for growth, depending on the amount of labile C and N left:
            satisf_for_C, satisf_for_N \
                = RootSimple.computing_growth_demand_and_satisfaction_coefficients(root_g,
                                                                                   labile_C_pool=plant.root_labile_C_pool,
                                                                                   labile_N_pool=plant.root_labile_N_pool)
            # NOTE: the labile pools of C & N have been updated with the costs of growth by this function.
            print("    The satisfaction coefficient for C is", "{:.2f}".format(Decimal(satisf_for_C)),
                  "and the one for N is", "{:.2f}".format(Decimal(satisf_for_N)))
            # We use the function ArchiSimple_growth to adapt the potential growth to the available labile C and N,
            # and we calculate the consumption rates for C & N, without updating the labile pools of C and N:
            RootSimple.ArchiSimple_growth_with_C_and_N(plant, root_g,
                                                       satisf_for_C, satisf_for_N,
                                                       time_step_in_seconds=time_step_in_seconds)
            # We proceed to the segmentation of the whole root system:
            RootSimple.segmentation_and_primordia_formation(root_g, time_step_in_seconds,
                                                            soil_temperature_in_Celsius=soil_temperature, random=True,
                                                            root_order_limitation=False)
            # TODO: Include root hairs dynamics once the cost of C & N is properly registered!
            # # We modifiy root hairs characteristics according to their specific dynamics:
            # root_hairs_dynamics(g, time_step_in_seconds=time_step_in_seconds,
            #                     soil_temperature_in_Celsius=soil_temperature)
            # We update the distance from tip for each root element in each root axis:
            RootSimple.update_distance_from_tip(root_g)
            # We update the surfaces and volumes of each root element:
            RootSimple.update_surfaces_and_volumes(root_g)
            RootSimple.updating_exchange_surface(plant, root_g)

            # FINAL ROOT C&N BALANCE:
            # We update the pools of root labile C and N at the scale of the whole root system:
            RootSimple.root_CN_balance(plant, time_step_in_seconds=time_step_in_seconds)

            print("    Root labile C pool is now",
                  "{:.2E}".format(Decimal(plant.root_labile_C_pool - plant.root_labile_C_deficit)),
                  "gC and root labile N pool is",
                  "{:.2E}".format(Decimal(plant.root_labile_N_pool - plant.root_labile_N_deficit)), "gN.")
            print("    The root system now contains", root_g.nb_vertices(), "different root elements.")

        # TRANSFER OF N FROM ROOTS TO SHOOTS:
        #####################################

        print("    > Computing root-to-shoot N transfer...")
        # We first compute the rate of net N transfer from roots to shoots:
        computing_shoot_root_exchange_rates(plant, computing_C_exchange=False, computing_N_exchange=True,
                                            N_option="by_diffusion")
        print("    > The net transfer of N is",
              "{:.2E}".format(Decimal(plant.net_N_transfer_from_root_to_shoot_rate * time_step_in_seconds)), "gN.")
        # We then limit the transfer of N from roots to shoots (or from shoots to roots) to the maximal amount of N
        # available in the source compartment.
        # We first calculate the expected amount of N that should be exchanged based on the exchange rate that has been computed:
        expected_N_exchange = plant.net_N_transfer_from_root_to_shoot_rate * time_step_in_seconds
        # CASE 1: The flow of N from root to shoot is positive, but there is not enough N in the root to sustain it.
        if expected_N_exchange > (plant.root_labile_N_pool - plant.root_labile_N_deficit):
            # Then the amount of N exchanged is limited to the amount of N available in roots:
            final_N_exchange = (plant.root_labile_N_pool - plant.root_labile_N_deficit)
            # The amount of labile N in shoots and roots is updated according to the new computed exchange:
            if SIMULATION_OPTION !="roots_only":
                plant.shoot_labile_N_pool += final_N_exchange
            if SIMULATION_OPTION != "shoot_only":
                plant.root_labile_N_pool -= final_N_exchange
        # CASE 2: The flux from shoot to root is positive but there is not enough N in the shoot to sustain it:
        elif expected_N_exchange < 0. and abs(expected_N_exchange) > (plant.shoot_labile_N_pool - plant.shoot_labile_N_deficit):
            # Then the amount of N exchanged is limited to the amount of N available in shoots:
            final_N_exchange = (plant.shoot_labile_N_pool - plant.shoot_labile_N_deficit)
            # The amount of labile N in shoots and roots is updated according to the new computed exchange:
            if SIMULATION_OPTION != "roots_only":
                plant.shoot_labile_N_pool -= final_N_exchange
            if SIMULATION_OPTION != "shoot_only":
                plant.root_labile_N_pool += final_N_exchange
        # CASE 3 (EXPECTED CASE): There is enough N in the source compartment to sustain the transfer.
        else:
            # Then the amount of labile N in shoots and roots is updated according to the computed exchange:
            if SIMULATION_OPTION != "roots_only":
                plant.shoot_labile_N_pool += expected_N_exchange
            if SIMULATION_OPTION != "shoot_only":
                plant.root_labile_N_pool -= expected_N_exchange

        # CALLING SUCCESSIVELY THE DIFFERENT SHOOT MODULES:
        ###################################################

        if SIMULATION_OPTION != "roots_only":

            # # OPTION WITH PRIORITY OF CONSUMPTION FOR SUCCESSIVE PROCESSES:
            # ###############################################################

            # print("\n    Shoot labile C pool is initially",
            #       "{:.2E}".format(Decimal(plant.shoot_labile_C_pool - plant.shoot_labile_C_deficit)),
            #       "gC and shoot labile N pool is",
            #       "{:.2E}".format(Decimal(plant.shoot_labile_N_pool - plant.shoot_labile_N_deficit)), "gN.")
            #
            # # UPTAKE OF C BY THE SHOOT ORGANS:
            # # We calculate light interception based on the shoot architecture, for each organ:
            # AllocF.computing_light_interception(shoot_g,
            #                                     incoming_light_density=light_input,
            #                                     number_of_plants_per_square_meter=plant_density_in_square_meter,
            #                                     inter_row_distance=inter_row_distance,
            #                                     scene_length=row_length,
            #                                     scene_width=scene_width,
            #                                     leaf_shape=param.leaf_shape
            #                                     )
            # # We calculate the rate of C production of each organ based on the PAR intercepted by the plant
            # # (and we also get the total PAR intercepted by the plant):
            # total_PAR_intercepted = \
            #     AllocF.computing_photosynthesis_rate(shoot_g,
            #                                          air_temperature_in_degree_Celsius=air_temperature)
            # # !!! At this stage, the pool of shoot labile C has not been updated yet!
            #
            # # MAINTENANCE IN THE SHOOT ORGANS:
            # # We calculate the C respiration cost associated to maintenance in each organ:
            # AllocF.computing_maintenance_rate(shoot_g, air_temperature_in_degree_Celsius=air_temperature)
            # # !!! At this stage, the pool of shoot labile C has not been updated yet!
            #
            # #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
            # # C&N BALANCE IN EACH ORGAN OF THE SHOOT SYSTEM:
            # # This function performs the C & N balance and updates the local labile pools of C & N IN EACH ORGAN,
            # # and resets the local exchange rates of each organ to 0, while calculating the sum of C&N available
            # # in all organs:
            # organs_available_C, organs_available_N \
            #     = AllocF.shoot_organ_CN_balance(shoot_g,
            #                                     time_step_in_seconds=time_step_in_seconds,
            #                                     resetting_rates=True)
            # # !!! The pool of shoot labile C & N has not been updated yet, but the ones in each organ have been updated!
            # #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
            #
            # # REDISTRIBUTION OF C AMONG ORGANS:
            # # We calculate the total amount of labile C available in the shoot system:
            # shoot_available_C = (plant.shoot_labile_C_pool - plant.shoot_labile_C_deficit) + organs_available_C
            # # We computes the rate of labile C redistribution in each organ according to its surfacic PAR:
            # AllocF.computing_shoot_C_redistribution(shoot_g,
            #                                         shoot_available_C = shoot_available_C,
            #                                         intercepted_light = total_PAR_intercepted,
            #                                         time_step_in_seconds=time_step_in_seconds,
            #                                         option="surfacic_PAR")
            #
            # # DISTRIBUTION OF N AMONG ORGANS:
            # # We calculate the total amount of labile N available in the shoot system:
            # shoot_available_N = (plant.shoot_labile_N_pool - plant.shoot_labile_N_deficit) + organs_available_N
            # # We computes the rate of labile N distribution in each organ according to its surfacic PAR:
            # AllocF.computing_shoot_N_distribution(shoot_g,
            #                                       shoot_available_N = shoot_available_N,
            #                                       time_step_in_seconds=time_step_in_seconds)
            #
            # #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
            # # C&N BALANCE IN EACH ORGAN OF THE SHOOT SYSTEM:
            # # This function performs the C & N balance and updates the local labile pools of C & N IN EACH ORGAN,
            # # and resets the local exchange rates of each organ to 0, while calculating the sum of C&N available
            # # in all organs:
            # AllocF.shoot_organ_CN_balance(shoot_g,
            #                               time_step_in_seconds=time_step_in_seconds,
            #                               resetting_rates=True)
            # # !!! The pool of shoot labile C & N has not been updated yet, but the ones in each organ have been updated!
            # # We reset the shoot labile C and N to 0, as it has been included in the previous calculation:
            # plant.shoot_labile_C_pool = 0.
            # plant.shoot_labile_N_pool = 0.
            # #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
            #
            # # GROWTH IN THE SHOOTS:
            # # We start considering the possibility of ramification:
            # AllocF.computing_ramification(shoot_MTG=shoot_g,
            #                               initial_internode_length=param.initial_internode_length,
            #                               initial_leaf_length=param.initial_leaf_length)
            #
            # # We calculate the possible increase in volume of each organ:
            # AllocF.potential_shoot_growth(shoot_g,
            #                               SFQN=param.SFQN, #TODO: is SFQN constant over the whole range of (0, LA_max)? Is it two distinct curves?
            #                               LA_max=param.LA_max,
            #                               leaf_minimal_thickness=param.leaf_minimal_thickness,
            #                               elong_internode=param.elong_internode,
            #                               time_step_in_seconds=time_step_in_days*60.*60.*24.,
            #                               air_temperature_in_degree_Celsius = air_temperature)
            #
            # # We calculate the fraction of C that can be redistributed by source leaves, and distribute a fraction of it
            # # in sink leaves:
            # AllocF.computing_C_redistribution_for_shoot_growth(shoot_system=plant,
            #                                                    shoot_MTG=shoot_g,
            #                                                    yield_growth=param.yield_growth,
            #                                                    structural_mass_density=param.struct_mass_density,
            #                                                    structural_C_content=param.structural_C_content,
            #                                                    time_step_in_seconds=time_step_in_days*60.*60.*24.,
            #                                                    option="surfacic_PAR")
            # # Note that the remaining labile C that has not been allocated to sinks has returned to the shoot labile C pool!
            #
            # #-----------------------------------------------------------------------------------------------------------
            # # TRANSFER OF C FROM SHOOTS TO ROOTS:
            # #####################################
            # print("    > Computing shoot-to-root C transfer...")
            # # At this stage, we compute the rate of net C transfer from shoots to roots:
            # computing_shoot_root_exchange_rates(plant, computing_C_exchange=True, computing_N_exchange=False,
            #                                     C_option="by_diffusion")
            # print("    > The net transfer of C was",
            #       "{:.2E}".format(Decimal(plant.net_C_transfer_from_shoot_to_root_rate * time_step_in_seconds)), "gC.")
            # # Note that only the rate has been computed here, the labile C pools have not been updated yet!
            # # ----------------------------------------------------------------------------------------------------------
            #
            # # We calculate the actual increase in size and the corresponding C & N costs, as well as the ageing of each
            # # organ, based on the local amounts of labile C:
            # AllocF.actual_shoot_growth(shoot_g,
            #                            yield_growth=param.yield_growth,
            #                            structural_mass_density=param.struct_mass_density,
            #                            structural_C_content=param.structural_C_content,
            #                            structural_N_content=param.structural_N_content,
            #                            leaf_minimal_thickness=param.leaf_minimal_thickness,
            #                            time_step_in_seconds=time_step_in_days*60.*60.*24.,
            #                            air_temperature_in_degree_Celsius=air_temperature
            #                            )
            # # At this stage, the shoot structure has been updated and the consumption rates have been calculated
            # # for each organ, but no actual consumption of labile C&N at the scale of organ or shoot system has occured yet.
            #
            # # We now consider the initiation of the structure of new phytomers (if any).
            # # We increment the time since the last phytomer was formed at the top of the main stem:
            # thermal_time_since_last_phytomer_formation\
            #     = shoot_g.node(0).thermal_time_since_last_phytomer_formation + thermal_time_step_in_seconds_shoots
            # # We get the value of the prescribed phyllochron (in seconds equivalent to the reference temperature):
            # phyllochron = param.phyllochron
            # # If the period since the last phytomer was formed is higher than the phyllochron:
            # if thermal_time_since_last_phytomer_formation > phyllochron:
            #     print("A new phytomer should be formed!")
            #     # Then at least one phytomer should be formed on the main stem.
            #     # We keep looping until all phytomers have been formed:
            #     while thermal_time_since_last_phytomer_formation - phyllochron > 0:
            #         # We subtract the value of the phyllochron to the time since last phytomer was formed:
            #         thermal_time_since_last_phytomer_formation -= phyllochron
            #         # We add a set of new phytomers in the shoot system:
            #         AllocF.new_phytomers_initiation(shoot_system=plant,
            #                                         shoot_MTG=shoot_g,
            #                                         new_initiation=True,
            #                                         initial_leaf_length=param.initial_leaf_length,
            #                                         initial_internode_length=param.initial_internode_length)
            #         print("> A new phytomer has been formed!")
            #         # Note: The cost of this initiation has been subtracted to the labile C and N pools in the shoot.
            # # Once all phytomers have been added, we update the time since the last one was formed,
            # # which is stored at the base of the shoot system:
            # shoot_g.node(0).thermal_time_since_last_phytomer_formation = thermal_time_since_last_phytomer_formation
            # # We finally update all dimensions of the organs:
            # updating_coordinates(shoot_g)
            #
            # # EXCHANGE WITH RESERVE:
            # # This function reallocates all the remaining C available to the different organs according to their PAR,
            # # and stores it as additional structural mass:
            # AllocF.computing_shoot_C_storage(shoot_system=plant,
            #                                  shoot_MTG=shoot_g,
            #                                  time_step_in_seconds=time_step_in_days*60.*60.*24.)
            # # This function mobilizes N from the organs and stores it at the base of the shoot system:
            # AllocF.computing_shoot_N_remobilization(shoot_g, k_N_mobilization=param.k_N_mobilization)
            # # This function immobilizes N in each shoot organ:
            # AllocF.computing_shoot_N_immobilization(shoot_g)
            #
            # # SENESCENCE AND LEAF FALL:
            # AllocF.computing_senescence(shoot_g, time_step_in_seconds=time_step_in_seconds, air_temperature_in_degree_Celsius=air_temperature)
            # AllocF.computing_leaf_fall(shoot_g, time_step_in_seconds=time_step_in_seconds, air_temperature_in_degree_Celsius=air_temperature)
            #
            # #///////////////////////////////////////////////////////////////////////////////////////////////////////////
            # # C&N BALANCE IN THE SHOOT ORGANS:
            # # This function performs the actual C & N balance, updates the labile pools of C & N in each organ and reset
            # # the exchange rates of the whole shoot system to 0:
            # organs_available_C, organs_available_N \
            #     = AllocF.shoot_organ_CN_balance(shoot_g, time_step_in_seconds=time_step_in_seconds, resetting_rates=True)
            # # //////////////////////////////////////////////////////////////////////////////////////////////////////////
            #
            # # C&N BALANCE IN THE SHOOT SYSTEM:
            # plant.shoot_labile_C_pool += organs_available_C
            # plant.shoot_labile_N_pool += organs_available_N
            #
            # print("    The shoot system now contains", shoot_g.nb_vertices(), "different shoot elements.")
            # print("    Shoot labile C pool is now",
            #       "{:.2E}".format(Decimal(plant.shoot_labile_C_pool - plant.shoot_labile_C_deficit)),
            #       "gC and labile N pool is",
            #       "{:.2E}".format(Decimal(plant.shoot_labile_N_pool - plant.shoot_labile_N_deficit)), "gN.")
            # print("")

            # # TRANSFER OF C BETWEEN SHOOTS AND ROOTS:
            # #########################################
            #
            # # print(" --- Before C transfer to root, the root labile C is", plant.root_labile_C_pool)
            # # We update the properties of the whole plant for the C transfer to shoots:
            # computing_shoot_root_exchange_rates(plant, computing_C_exchange=True, computing_N_exchange=False,
            #                                     C_option="by_diffusion")
            # # We limit the transfer of C from shoots to roots (or from roots to shoots) to the maximal amount of labile C
            # # available in the source compartment.
            # # We first calculate the expected amount of C that should be exchanged based on the flux:
            # expected_C_exchange = plant.net_C_transfer_from_shoot_to_root_rate * time_step_in_seconds
            # # CASE 1: The flux from shoot to root is positive, but there is not enough C in the shoot to sustain it:
            # if expected_C_exchange > (plant.shoot_labile_C_pool - plant.shoot_labile_C_deficit):
            #     # Then the amount of C exchanged is limited to the amount of C available in shoots:
            #     final_C_exchange = (plant.shoot_labile_C_pool - plant.shoot_labile_C_deficit)
            #     #  The amount of labile C in shoots and roots is updated according to the new computed exchange:
            #     if SIMULATION_OPTION != "roots_only":
            #         plant.shoot_labile_C_pool -= final_C_exchange
            #     if SIMULATION_OPTION != "shoot_only":
            #         plant.root_labile_C_pool += final_C_exchange
            # # CASE 2: The flux from root to shoot is positive but there is not enough C in the root to sustain it:
            # elif expected_C_exchange < 0. and abs(expected_C_exchange) > (plant.root_labile_C_pool - plant.root_labile_C_deficit):
            #     final_C_exchange = (plant.root_labile_C_pool - plant.root_labile_C_deficit)
            #     # The amount of labile C in shoots and roots is updated according to the new computed exchange:
            #     if SIMULATION_OPTION != "roots_only":
            #         plant.shoot_labile_C_pool += final_C_exchange
            #     if SIMULATION_OPTION != "shoot_only":
            #         plant.root_labile_C_pool -= final_C_exchange
            # # CASE 3 (EXPECTED CASE): There is enough C in the source compartment to sustain the transfer:
            # else:
            #     # Then the amount of labile C in shoots and roots is updated according to the initial computed exchange:
            #     if SIMULATION_OPTION != "roots_only":
            #         plant.shoot_labile_C_pool -= expected_C_exchange
            #     if SIMULATION_OPTION != "shoot_only":
            #         plant.root_labile_C_pool += expected_C_exchange

            # OPTION WITH NO PRIORITY OF CONSUMPTION:
            #########################################

            # We register the initial state of the labile C & N pools in the shoot:
            initial_shoot_labile_C = plant.shoot_labile_C_pool
            initial_shoot_labile_N = plant.shoot_labile_N_pool

            print("\n    Shoot labile C pool is initially",
                  "{:.2E}".format(Decimal(initial_shoot_labile_C)),
                  "gC and shoot labile N pool is",
                  "{:.2E}".format(Decimal(initial_shoot_labile_N)), "gN.")

            # LIGHT INTERCEPTION:
            # We calculate light interception based on the shoot architecture, for each organ:
            total_PAR_intercepted, sum_of_surfacic_PAR =\
                AllocF.computing_light_interception(shoot_g,
                                                incoming_light_density=light_input,
                                                number_of_plants_per_square_meter=plant_density_in_square_meter,
                                                inter_row_distance=inter_row_distance,
                                                scene_length=row_length,
                                                scene_width=scene_width,
                                                leaf_shape=param.leaf_shape)
            plant.total_PAR_intercepted = total_PAR_intercepted

            # DISTRIBUTION OF C AMONG ORGANS:
            # We calculate the total amount of labile C available in the shoot system:
            shoot_available_C = plant.shoot_labile_C_pool
            if shoot_available_C < 0.:
                print("WATCH OUT: the amount of shoot labile C was negative! No redistribution of C occured.")
            else:
                print("    => Redistribution of",
                      "{:.2E}".format(Decimal(shoot_available_C)),
                      "gC among the shoot organs...")
                # We computes the rate of labile C redistribution in each organ according to its surfacic PAR:
                AllocF.computing_shoot_C_redistribution(shoot_g,
                                                        shoot_available_C=shoot_available_C,
                                                        intercepted_light=None,
                                                        time_step_in_seconds=time_step_in_seconds,
                                                        option="surfacic_PAR")
                # We reset the general shoot labile C pools, as the labile C has been transfered to each organ:
                plant.shoot_labile_C_pool = 0.
                # plant.shoot_labile_C_deficit = 0.

            # DISTRIBUTION OF N AMONG ORGANS:
            # We calculate the total amount of labile N available in the shoot system:
            shoot_available_N = plant.shoot_labile_N_pool
            if shoot_available_N < 0.:
                print("WATCH OUT: the amount of shoot labile N was negative! No redistribution of N occured.")
            else:
                print("    => Redistribution of",
                      "{:.2E}".format(Decimal(shoot_available_N)),
                      "gN among the shoot organs...")
                # We computes the rate of labile N distribution in each organ according to its surfacic PAR:
                N_redistributed = AllocF.computing_shoot_N_distribution(shoot_g,
                                                                        shoot_available_N=shoot_available_N,
                                                                        time_step_in_seconds=time_step_in_seconds)
                # We reset the general shoot labile C pools, as the labile C has been transfered to each organ:
                plant.shoot_labile_N_pool = 0.
                # plant.shoot_labile_N_deficit = 0.

            #///////////////////////////////////////////////////////////////////////////////////////////////////////////
            # C&N BALANCE IN THE SHOOT ORGANS:
            # This function performs the actual C & N balance, updates the labile pools of C & N in each organ and reset
            # the exchange rates of the whole shoot system to 0:
            AllocF.shoot_organ_CN_balance(shoot_g, time_step_in_seconds=time_step_in_seconds,
                                          resetting_rates=True, resetting_local_amounts=False)
            #///////////////////////////////////////////////////////////////////////////////////////////////////////////

            # PHOTOSYNTHESIS:
            # We calculate the rate of C production of each organ based on the PAR intercepted by the plant
            # (and we also get the total PAR intercepted by the plant):
            total_PAR_intercepted, total_photosynthesis_rate = \
                AllocF.computing_photosynthesis_rate(shoot_g,
                                                     air_temperature_in_degree_Celsius=air_temperature)
            plant.photosynthesized_C_rate = total_photosynthesis_rate
            print("    Photosynthesis will generate",
                  "{:.2E}".format(Decimal(total_photosynthesis_rate*time_step_in_seconds)),
                  "gC...")
            # Note: at this stage, the labile C amounts have not been updated in the organs or at the scale of the plant!

            # MAINTENANCE IN THE SHOOT ORGANS:
            # We calculate the C respiration cost associated to maintenance in each organ:
            total_maintenance_rate = AllocF.computing_maintenance_rate(shoot_g,
                                                                       air_temperature_in_degree_Celsius=air_temperature)
            plant.shoot_maintenance_respiration_rate = total_maintenance_rate
            print("    The maintenance will consume",
                  "{:.2E}".format(Decimal(total_maintenance_rate*time_step_in_seconds)),
                  "gC over the shoot system...")

            # GROWTH IN THE SHOOT ORGANS:
            # We start considering the possibility of ramification:
            C_cost_for_new_ramification, N_cost_for_new_ramification \
                = AllocF.computing_ramification(shoot_MTG=shoot_g,
                                                initial_internode_length=param.initial_internode_length,
                                                initial_leaf_length=param.initial_leaf_length)
            # We update the pools of shoot labile C & N with the costs of ramification:
            plant.shoot_labile_C_pool -= C_cost_for_new_ramification
            plant.shoot_labile_N_pool -= N_cost_for_new_ramification
            if C_cost_for_new_ramification > 0. or N_cost_for_new_ramification > 0.:
                print("> The costs to sustain ramification were",
                      "{:.2E}".format(Decimal(C_cost_for_new_ramification)), "gC and",
                      "{:.2E}".format(Decimal(N_cost_for_new_ramification)), "gN.")

            # We calculate the possible increase in volume of each organ:
            AllocF.potential_shoot_growth(shoot_g,
                                          SFQN=param.SFQN,
                                          # TODO: is SFQN constant over the whole range of (0, LA_max)? Is it two distinct curves?
                                          LA_max=param.LA_max,
                                          leaf_minimal_thickness=param.leaf_minimal_thickness,
                                          elong_internode=param.elong_internode,
                                          time_step_in_seconds=time_step_in_days * 60. * 60. * 24.,
                                          air_temperature_in_degree_Celsius=air_temperature)
            # We calculate the actual increase in size and the corresponding C & N costs, as well as the ageing of each
            # organ, based on the local amounts of labile C:
            C_consumption_rate_for_extension, N_consumption_rate_for_extension = \
                AllocF.actual_shoot_growth(shoot_g,
                                           yield_growth=param.yield_growth,
                                           structural_mass_density=param.struct_mass_density,
                                           structural_C_content=param.structural_C_content,
                                           structural_N_content=param.structural_N_content,
                                           leaf_minimal_thickness=param.leaf_minimal_thickness,
                                           time_step_in_seconds=time_step_in_days * 60. * 60. * 24.,
                                           air_temperature_in_degree_Celsius=air_temperature)
            # At this stage, the shoot structure has been updated and the consumption rates have been calculated
            # for each organ, but no actual consumption of labile C&N at the scale of organ or shoot system has occured yet.
            print("    The extension of organs will consume",
                  "{:.2E}".format(Decimal(C_consumption_rate_for_extension * time_step_in_seconds)),
                  "gC and",
                  "{:.2E}".format(Decimal(N_consumption_rate_for_extension * time_step_in_seconds)),
                  "gN over the shoot system...")
            # We now consider the initiation of the structure of new phytomers (if any).
            # We increment the time since the last phytomer was formed at the top of the main stem:
            thermal_time_since_last_phytomer_formation \
                = shoot_g.node(0).thermal_time_since_last_phytomer_formation + thermal_time_step_in_seconds_shoots
            # We get the value of the prescribed phyllochron (in seconds equivalent to the reference temperature):
            phyllochron = param.phyllochron
            # We initialize the total costs of phytomers initiation:
            total_C_costs_for_new_phytomers = 0.
            total_N_costs_for_new_phytomers = 0.
            # If the period since the last phytomer was formed is higher than the phyllochron:
            if thermal_time_since_last_phytomer_formation > phyllochron:
                print("A new phytomer should be formed!")
                # Then at least one phytomer should be formed on the main stem.
                # We keep looping until all phytomers have been formed:
                while thermal_time_since_last_phytomer_formation - phyllochron > 0:
                    # We subtract the value of the phyllochron to the time since last phytomer was formed:
                    thermal_time_since_last_phytomer_formation -= phyllochron
                    # We add a set of new phytomers in the shoot system:
                    C_cost_for_new_phytomers, N_cost_for_new_phytomers = AllocF.new_phytomers_initiation(
                        shoot_system=plant,
                        shoot_MTG=shoot_g,
                        new_initiation=True,
                        initial_leaf_length=param.initial_leaf_length,
                        initial_internode_length=param.initial_internode_length)
                    # We include these costs in the labile C and N pools in the shoot system:
                    plant.shoot_labile_C_pool -= C_cost_for_new_phytomers
                    plant.shoot_labile_N_pool -= N_cost_for_new_phytomers
                    print("> A new phytomer has been formed! The C & N costs were",
                          "{:.2E}".format(Decimal(C_cost_for_new_phytomers)), "gC and",
                          "{:.2E}".format(Decimal(N_cost_for_new_phytomers)),"gN.")
                    # We increase the total costs with the cost of this new phytomer:
                    total_C_costs_for_new_phytomers += C_cost_for_new_phytomers
                    total_N_costs_for_new_phytomers += total_N_costs_for_new_phytomers
            # Once all phytomers have been added, we update the time since the last one was formed,
            # which is stored at the base of the shoot system:
            shoot_g.node(0).thermal_time_since_last_phytomer_formation = thermal_time_since_last_phytomer_formation
            # We finally update all dimensions of the organs:
            updating_coordinates(shoot_g)

            # And we register the total consumption rate for growth in the plant object:
            plant.actual_C_consumption_rate_by_shoot_growth\
                = (C_consumption_rate_for_extension
                   + (C_cost_for_new_ramification + total_C_costs_for_new_phytomers) / time_step_in_seconds)
            plant.actual_N_consumption_rate_by_shoot_growth \
                = (N_consumption_rate_for_extension
                   + (N_cost_for_new_ramification + total_N_costs_for_new_phytomers) / time_step_in_seconds)

            # # EXCHANGE WITH RESERVE:
            # # This function reallocates all the remaining C available to the different organs according to their PAR,
            # # and stores it as additional structural mass:
            # AllocF.computing_shoot_C_storage(shoot_system=plant,
            #                                  shoot_MTG=shoot_g,
            #                                  time_step_in_seconds=time_step_in_days * 60. * 60. * 24.)
            # # This function mobilizes N from the organs and stores it at the base of the shoot system:
            # AllocF.computing_shoot_N_remobilization(shoot_g, k_N_mobilization=param.k_N_mobilization)
            # # This function immobilizes N in each shoot organ:
            # AllocF.computing_shoot_N_immobilization(shoot_g)
            #
            # # SENESCENCE AND LEAF FALL:
            # AllocF.computing_senescence(shoot_g, time_step_in_seconds=time_step_in_seconds,
            #                             air_temperature_in_degree_Celsius=air_temperature)
            # AllocF.computing_leaf_fall(shoot_g, time_step_in_seconds=time_step_in_seconds,
            #                            air_temperature_in_degree_Celsius=air_temperature)

            #///////////////////////////////////////////////////////////////////////////////////////////////////////////
            # C&N BALANCE IN THE SHOOT ORGANS:
            # This function performs the actual C & N balance, updates the labile pools of C & N in each organ and reset
            # the labile C and N exchange rates locally to 0:
            organs_available_C, organs_available_N \
                = AllocF.shoot_organ_CN_balance(shoot_g, time_step_in_seconds=time_step_in_seconds,
                                                resetting_rates=True, resetting_local_amounts=True)
            #///////////////////////////////////////////////////////////////////////////////////////////////////////////
            # We move the available labile C & N from the organ scale to shoot system scale:
            plant.shoot_labile_C_pool += organs_available_C
            plant.shoot_labile_N_pool += organs_available_N
            # At this stage, all the labile C and N has been transfered at the scale of the shoot system, organs are "empty"!

            print("    The shoot system now contains", shoot_g.nb_vertices(), "different shoot elements.")
            print("    Shoot labile C pool is now",
                  "{:.2E}".format(Decimal(plant.shoot_labile_C_pool)),
                  "gC and labile N pool is",
                  "{:.2E}".format(Decimal(plant.shoot_labile_N_pool)), "gN.")

            # TRANSFER OF C FROM SHOOTS TO ROOTS:
            #####################################
            print("    > Computing shoot-to-root C transfer...")
            # At this stage, we compute the rate of net C transfer from shoots to roots:
            computing_shoot_root_exchange_rates(plant, computing_C_exchange=True, computing_N_exchange=False,
                                                C_option="by_diffusion")
            print("    > The net transfer of C was",
                  "{:.2E}".format(Decimal(plant.net_C_transfer_from_shoot_to_root_rate * time_step_in_seconds)), "gC.")
            # Note that only the rate has been computed here, the labile C pools have not been updated yet!
            # We limit the transfer of C from shoots to roots (or from roots to shoots) to the maximal amount of labile C
            # available in the source compartment.
            # We first calculate the expected amount of C that should be exchanged based on the flux:
            expected_C_exchange = plant.net_C_transfer_from_shoot_to_root_rate * time_step_in_seconds
            # CASE 1: The flux from shoot to root is positive, but there is not enough C in the shoot to sustain it:
            if expected_C_exchange > (plant.shoot_labile_C_pool - plant.shoot_labile_C_deficit):
                # Then the amount of C exchanged is limited to the amount of C available in shoots:
                final_C_exchange = (plant.shoot_labile_C_pool - plant.shoot_labile_C_deficit)
                #  The amount of labile C in shoots and roots is updated according to the new computed exchange:
                if SIMULATION_OPTION != "roots_only":
                    plant.shoot_labile_C_pool -= final_C_exchange
                if SIMULATION_OPTION != "shoot_only":
                    plant.root_labile_C_pool += final_C_exchange
            # CASE 2: The flux from root to shoot is positive but there is not enough C in the root to sustain it:
            elif expected_C_exchange < 0. and abs(expected_C_exchange) > (plant.root_labile_C_pool - plant.root_labile_C_deficit):
                final_C_exchange = (plant.root_labile_C_pool - plant.root_labile_C_deficit)
                # The amount of labile C in shoots and roots is updated according to the new computed exchange:
                if SIMULATION_OPTION != "roots_only":
                    plant.shoot_labile_C_pool += final_C_exchange
                if SIMULATION_OPTION != "shoot_only":
                    plant.root_labile_C_pool -= final_C_exchange
            # CASE 3 (EXPECTED CASE): There is enough C in the source compartment to sustain the transfer:
            else:
                # Then the amount of labile C in shoots and roots is updated according to the initial computed exchange:
                if SIMULATION_OPTION != "roots_only":
                    plant.shoot_labile_C_pool -= expected_C_exchange
                if SIMULATION_OPTION != "shoot_only":
                    plant.root_labile_C_pool += expected_C_exchange

        # CHECKING C & N BALANCE AND RECORDING:
        #######################################

        net_shoot_C_variation_rate = (total_photosynthesis_rate
                                      - total_maintenance_rate
                                      - C_consumption_rate_for_extension
                                      - total_C_costs_for_new_phytomers / time_step_in_seconds
                                      - C_cost_for_new_ramification / time_step_in_seconds)
        expected_shoot_labile_C_pool = initial_shoot_labile_C + net_shoot_C_variation_rate * time_step_in_seconds
        print("    The expected shoot labile C pool is",
              "{:.2E}".format(Decimal(expected_shoot_labile_C_pool)), "gC per plant.")
        print("")

        # We record the main properties of the plant and add it to a general table:
        # if results_table == []:
        #     # We create a list containing all properties of the plant, which will be used to record an output table:
        #     list_of_properties = [property for property in plant.__dict__.keys()]
        #     column_names = ["Time (years)"] + list_of_properties
        #     # We initialize the table with the headers:
        #     results_table = [column_names]
        # plant_properties = [time_step_in_days * 365]
        # for property in plant.__dict__.keys():
        #     # We add the value of this property to the list:
        #     plant_properties.append(getattr(plant, property, "NA"))
        # results_table.append(plant_properties)

        # PLOTTING ONE SINGLE PLANT OVER TIME:
        ######################################

        # If the camera is to be moved around the scene, we modify the position of the camera:
        if param.camera_movement is True:
            x_cam = x_coordinates[index_camera]
            y_cam = y_coordinates[index_camera]
            z_cam = z_coordinates[index_camera]
            # We define the index of the coordinates to read at the next step:
            index_camera = index_camera + 1
            # If this index is higher than the number of coordinates in each vector and if no other movement
            # of the camera (away or down) has been set over time:
            if (index_camera >= param.camera_rotation_n_points
                    and param.step_back_coefficient == 0. and param.step_up_coefficient == 0.):
                # Then we reset the index to 0, so that the camera goes back to the same locations as previously stored:
                index_camera = 0
            # If "distancing" is required, the camera will continue to read the further coordinates already calculated.
        # Otherwise, we simply use the parameters specified for the camera:
        else:
            x_cam = param.x_camera
            y_cam = param.y_camera
            z_cam = param.z_camera

        # We display the MTG in its new state:
        plotting_whole_plant(shoot_MTG=shoot_g, root_MTG=root_g,
                             width=1200, height=1200,
                             x_center=param.x_center, y_center=param.y_center, z_center=param.z_center,
                             x_cam=x_cam, y_cam=y_cam, z_cam=z_cam,
                             grid=False,
                             leaf_shape=param.leaf_shape,
                             background_color=[94, 76, 64],
                             # single_color_shoot=[100,200,100],
                             # single_color_root=[200,200,100],
                             single_color_shoot=None,
                             # single_color_root=None,
                             property_name_shoot='PAR_intercepted', cmap_shoot='inferno', vmin_shoot=0, vmax_shoot=400, lognorm_shoot=False,
                             # property_name_root='radius', cmap_root='jet', vmin_root=1e-4, vmax_root=1e-2, lognorm_root=True,
                             recording_plot=True,
                             plot_name=os.path.join(plot_directory, 'plot_%.5d.png' % step))

        # RECORDING THE PROPERTIES of THE SHOOT MTG FOR THIS TIME STEP:
        # We can record the main properties of the plant:
        dataframe = recording_MTG_properties(shoot_g, file_name = os.path.join(shoot_MTG_directory, 'MTG_%.5d.csv' % step))

        # RECORDING THE FINAL RESULTS AT PLANT SCALE:
        results_dict = recording_main_plant_properties(plant, final_time_in_days)
        new_df = pd.DataFrame([results_dict.values()], columns=results_dict.keys())
        plant_data_frame = pd.concat([plant_data_frame, new_df])

        # new_row = pd.Series(results_dict)
        # plant_data_frame = pd.concat([plant_data_frame, new_row.to_frame().T], ignore_index=False)


        # We record the dataframe as a csv file:
        try:
            plant_data_frame.to_csv(os.path.join(outputs_directory,'plant_results.csv'), na_rep='NA', index=False, header=True)
        except:
            plant_data_frame.to_csv(os.path.join(outputs_directory,'plant_results_BACKUP.csv'), na_rep='NA', index=False, header=True)

    print("")
    print("Simulation is done!")

    return
#-----------------------------------------------------------------------------------------------------------------------
########################################################################################################################
########################################################################################################################

# MAIN SIMULATION:
##################

if __name__ == '__main__':

    # We run the whole simulation according to the instructions and parameters defined in the 'simbal_parameters.xlsx'
    # file (otherwise, default parameters from the 'simbal_parameters.py' file or from the function will be used):
    run_simulation()