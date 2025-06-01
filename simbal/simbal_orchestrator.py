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

import AllocF as AllocF
import RootSimple as RootSimple
import simbal_parameters as param
from tools import (loading_parameters_from_file, formatted_inputs, temperature_modification, recording_MTG_properties,
                   circle_coordinates, plotting_whole_plant)

import os
import numpy as np
import pandas as pd
from math import sqrt, pi, floor, trunc, exp, isnan, cos, sin, acos
from decimal import Decimal

########################################################################################################################
########################################################################################################################

# We define a class that will be used to contain the properties of the interface between shoot and roots:
class whole_plant(object):
    # We initiate the object:
    def __init__(self):
        self.shoot_struct_mass = 0.
        self.root_struct_mass = 0.

        self.shoot_labile_C_pool = 0. # mol of C per shoot system
        self.shoot_labile_C_deficit = 0. # mol of C per shoot system
        self.shoot_labile_N_pool = 0. # mol of N per shoot system
        self.shoot_labile_N_deficit = 0. # mol of N per shoot system
        self.root_labile_C_pool = 0.  # mol of C per root system
        self.root_labile_C_deficit = 0. # mol of C per root system
        self.root_labile_N_pool = 0.  # mol of N per root system
        self.root_labile_N_deficit = 0. # mol of N per root system

        self.C_flux_to_shoot = 0.  # mol of C per second
        self.N_flux_to_shoot = 0.  # mol of N per second
        self.C_flux_to_root = 0.  # mol of C per second
        self.N_flux_to_root = 0.  # mol of N per second

        self.net_C_transfer_from_shoot_to_root_rate = 0. # mol of C per second
        self.net_N_transfer_from_root_to_shoot_rate = 0. # mol of N per second
        self.transpiration_flux = 0. # mol of water per second

        self.photosynthesized_C_rate = 0. # mol of C per second
        self.shoot_maintenance_respiration_rate = 0. # mol of C per second
        self.actual_C_consumption_rate_by_shoot_growth = 0. # mol of C per second
        self.actual_N_consumption_rate_by_shoot_growth = 0. # mol of N per second
        self.net_shoot_N_immobilization_rate = 0. # mol of N per second
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# We define a function that computes C & N exchanges between shoot and root systems:
def shoot_root_exchange(i, computing_C_exchange=False, computing_N_exchange=False, time_step_in_seconds=1*(60*60*24),
                        C_option="by_diffusion", N_option="by_diffusion",
                        water_concentration = 10.):

    """

    Parameters
    ----------
    i: interface object between shoot and root, containing the variables at the scale of the whole plant
    computing_C_exchange
    computing_N_exchange
    time_step_in_seconds
    C_option: option to transfer C from shoots to roots - either "by_diffusion" or "full_transfer_from_shoot_to_root"
    N_option
    resistance_to_C
    resistance_to_N
    water_concentration

    Returns
    -------

    """

    # COMPUTING THE EXCHANGE OF CARBON:
    if computing_C_exchange:
        # If C is to be exchanged by diffusion based on the labile pools in shoots and roots:
        if C_option == "by_diffusion":
            # Then the net C flux is calculated using a gradient of concentration, a diffusion coefficient,
            # and the section of exchange:
            i.exchange_section = 0.01
            if i.shoot_struct_mass > 0. and i.root_struct_mass > 0.:
                shoot_to_root_C_flux = (param.C_diffusion_coefficient
                                        * (i.shoot_labile_C_pool / i.shoot_struct_mass
                                           - i.root_labile_C_pool / i.root_struct_mass)) * i.exchange_section
            else:
                print("!!! ERROR: the structural mass of shoot and root are", i.shoot_struct_mass, "and", i.root_struct_mass, "gDW, respectively!")
                print("           > We set the flux of C from shoot to root to 0.")
                shoot_to_root_C_flux = 0.
            # We register the transfer rate:
            i.net_C_transfer_from_shoot_to_root_rate = shoot_to_root_C_flux

        # Otherwise, if the C present in the shoot is meant to be transfer in totality to the root system:
        elif C_option == "full_transfer_from_shoot_to_root":
            i.net_C_transfer_from_shoot_to_root_rate = i.shoot_labile_C_pool / time_step_in_seconds
        else:
            print("WATCH OUT: the option for computing the transfer of C from shoot to root is not recognized! The flux has been set to zero.")
            i.net_C_transfer_from_shoot_to_root_rate = 0.

    # COMPUTING THE EXCHANGE OF NITROGEN:
    if computing_N_exchange:
        # If N is to be exchanged by diffusion based on the labile pools in shoots and roots:
        if N_option == "by_diffusion":
            # Then the net N flux is calculated using a gradient of concentration, a diffusion coefficient,
            # and the section of exchange:
            i.exchange_section = 0.01
            if i.shoot_struct_mass > 0. and i.root_struct_mass > 0.:
                root_to_shoot_N_flux = (param.N_diffusion_coefficient
                                        * (i.root_labile_N_pool / i.root_struct_mass
                                           - i.shoot_labile_N_pool / i.shoot_struct_mass)) * i.exchange_section
            else:
                print("!!! ERROR: the structural mass of shoot and root are", i.shoot_struct_mass, "and",
                      i.root_struct_mass, "gDW, respectively!")
                print("           > We set the flux of N from root to shoot to 0.")
                root_to_shoot_N_flux = 0.
            # We register the transfer rate:
            i.net_N_transfer_from_root_to_shoot_rate = root_to_shoot_N_flux

        # Otherwise, if the nitrogen is meant to be transferred to shoots along with the transpired water,
        # a transpiration flux has to be provided, from which the N flux is calculated based on the concentration
        # of N in the root labile pool:
        elif N_option == "by_water_flux":
            if i.transpiration_flux >= 0.:
                i.net_N_transfer_from_root_to_shoot_rate = (i.transpiration_flux / water_concentration
                                                            * i.root_labile_N_pool / i.root_struct_mass)
                # NOTE: N_flux in mol of N per plant and per second
                #       transpiration_flux in mol of water per plant and per second
                #       water_concentration in mol of water per gDW of root structural mass
                #       root_labile_N_pool in mol of N per plant
                #       root_struct_mass in gDW of root structural mass
            else:
                print("ERROR - the transpiration flux was negative! The flux of N from root to shoot was set to 0.")
                i.net_N_transfer_from_root_to_shoot_rate = 0.
        # Otherwise, if the N present in the root is meant to be transfer in totality to the shoot system:
        elif N_option == "full_transfer_from_root_to_shoot":
            i.net_N_transfer_from_root_to_shoot_rate = i.root_labile_N_pool / time_step_in_seconds
        else:
            print("WATCH OUT: the option for computing the transfer of N from root to shoot is not recognized! The flux has been set to zero.")
            i.net_N_transfer_from_root_to_shoot_rate = 0.

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
# i.transpiration_flux = 0.1
#
# shoot_root_exchange(i, computing_C_exchange=True, computing_N_exchange=True,
#                     C_option="by_diffusion", N_option="by_water_flux")
# print("The new flux of C from shoot to root is", i.C_flux_to_root, "mol of C per second.")
# print("The new flux of N from root to shoot is", i.N_flux_to_shoot, "mol of N per second.")
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
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
                   initial_shoot_labile_C_concentration = 1e-3, # in gC per gram of structural mass
                   initial_shoot_labile_N_concentration = 1e-4, # in gN per gram of structural mass
                   initial_root_labile_C_concentration = 1e-3, # in gC per gram of structural mass
                   initial_root_labile_N_concentration = 1e-4 # in gN per gram of structural mass
                   ):


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

    # UPDATING PARAMETERS:
    #---------------------
    # The values of the parameters in the file 'simbal_parameters.py' are updated according to the values in the file
    # 'simbal_parameters.xlsx':
    loading_parameters_from_file(parameters_filepath = os.path.join(my_path, "simulation", "inputs", "simbal_parameters.xlsx"),
                                 outputs_path = os.path.join(my_path, "simulation", "outputs"))

    # SETTING INITIAL CONDITIONS
    #---------------------------
    
    # INITIALIZING WHOLE PLANT:
    # We create a "plant", i.e. a shoot-root interface that will store variables at the scale of the plant:
    plant = whole_plant()
    
    # INITIALIZING SHOOT STRUCTURE:
    # We create a MTG "shoot_g" corresponding to the shoot (stem and leaves):
    shoot_g = AllocF.create_shoot_MTG(number_of_phytomers=5)
    # We update the shoot system general properties based on the shoot MTG:
    for vid in shoot_g.vertices_iter(scale=1):
        n = shoot_g.node(vid)
        try:
            plant.shoot_struct_mass += n.struct_mass
        except:
            print("!!! ERROR: we could not add the properties of the shoot element", vid, "to the whole plant!")
    # We add an extra-property to the base of the shoot system, corresponding to the time since the last phytomer
    # was formed (used together with the phyllochrone later on):
    shoot_g.node(0).thermal_time_since_last_phytomer_formation = 0.
    # We record the properties in the initial state of the shoot MTG in an external file:
    recording_MTG_properties(shoot_g, os.path.join(shoot_MTG_directory, 'MTG_00000.csv'))

    # INITIALIZING ROOT STRUCTURE:
    # We create another MTG corresponding to the root system:
    root_g = RootSimple.initiate_mtg(initial_segment_length=5e-3)
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
    # We initialize the pools of labile and reserve C (gC per plant) and N (gN per plant):
    plant.root_labile_C_pool = initial_root_labile_C_concentration * plant.root_struct_mass
    plant.root_labile_N_pool = initial_root_labile_N_concentration * plant.root_struct_mass
    plant.root_reserve_C_pool = 0.
    plant.root_reserve_N_pool = 0.
    plant.shoot_labile_C_pool = initial_shoot_labile_C_concentration * plant.shoot_struct_mass
    plant.shoot_labile_N_pool = initial_shoot_labile_N_concentration * plant.shoot_struct_mass

    # DEFINING THE PERIOD OF SIMULATION:
    # ---------------------------------

    # We define the steps of the simulation based on the desired period:
    initial_time = starting_time_in_days
    simulation_period_in_days = simulation_period_in_days
    time_step_in_seconds = time_step_in_days * 60*60*24
    n_steps = trunc(simulation_period_in_days / time_step_in_days)
    # We initialize a counter that will be used over successive time steps together with the phyllochrone:
    elapsed_time_in_days_since_last_leaf_emergence = 0.

    # READING THE INPUT FILE:
    # -----------------------
    if input_file and not (forcing_constant_aboveground_environment and forcing_constant_belowground_environment):
        # We use the function 'formatted inputs' to create a table containing the input data for each required step,
        # depending on the chosen time step:
        input_frame = formatted_inputs(original_input_file=input_file,
                                       final_input_file=os.path.join(outputs_directory, 'updated_input_file.csv'),
                                       original_time_step_in_days=input_file_time_step_in_days,
                                       final_time_step_in_days=time_step_in_days,
                                       simulation_period_in_days=simulation_period_in_days,
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
    if param.camera_rotation is True:
        # We calculate the x,y,z coordinates of the different positions of the camera:
        [x_coordinates, y_coordinates, z_coordinates] = circle_coordinates(x_center=param.x_center,
                                                                           y_center=param.y_center,
                                                                           z_center=param.z_camera,
                                                                           radius=sqrt((param.x_camera**2 + param.y_camera**2)),
                                                                           n_points=param.camera_rotation_n_points)
        # We initialize the index in the list of the coordinates:
        index_camera = 0

    # LOOPING OVER TIME:
    ####################

    # For each time step covering the period:
    for step in range(initial_step_number, initial_step_number + n_steps):

        # We calculate the current time in hours:
        initial_time_in_days = initial_time + (step) * time_step_in_days
        final_time_in_days = initial_time + (step+1) * time_step_in_days

        # DEFINING THE INPUTS AT THIS TIME STEP:
        # --------------------------------------
        # If constant aboveground conditions have been chosen:
        if forcing_constant_aboveground_environment:
            # We set up the constant values:
            light_input = constant_light_irradiation
            air_temperature = constant_air_temperature_in_Celsius
        # Otherwise, if there is an input file to be used:
        elif input_file:
            light_input = input_frame.loc[step - initial_step_number, 'light_irradiation_in_micromol_per_m2_per_second']
            air_temperature = input_frame.loc[step - initial_step_number, 'air_temperature_in_degree_Celsius']
        else:
            print("!!! ERROR: no input file has been provided!")
        # Similarly, if constant belowground conditions have been chosen:
        if forcing_constant_belowground_environment:
            soil_temperature = constant_soil_temperature_in_Celsius
            soil_available_nitrogen = constant_soil_available_nitrogen
        elif input_file:
            soil_temperature = input_frame.loc[step - initial_step_number, 'soil_temperature_in_degree_Celsius']
            soil_available_nitrogen = input_frame.loc[step - initial_step_number, 'soil_available_nitrogen_in_gN_per_gDW_of_soil']

        print("")
        print("From", "{:.2f}".format(Decimal(initial_time_in_days)),"to", "{:.2f}".format(Decimal(final_time_in_days)),"days:")
        print("-------------------------")
        print("    Air temperature is", "{:.2f}".format(Decimal(air_temperature)), "degree Celsius; soil temperature is",
              "{:.2f}".format(Decimal(soil_temperature)), "degree Celsius.")
        print("    Incoming light is", "{:.2f}".format(Decimal(light_input)), "micromol of photons per square meter per second.")
        print("    Soil available N is", "{:.2f}".format(Decimal(soil_available_nitrogen)), "gN per gram of dry soil.")
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

        # CALLING SUCCESSIVELY THE DIFFERENT ROOT MODULES:
        ##################################################

        # We update the amount of labile C in roots after considering the C input to the roots previously stored:
        RootSimple.updating_root_labile_C_with_input(root_system = plant, time_step_in_seconds=time_step_in_seconds)
        # We reset the input of C to 0:
        plant.net_C_transfer_from_shoot_to_root_rate = 0.
        print("    Labile C pool in the root system after C input is initially", "{:.2E}".format(Decimal(plant.root_labile_C_pool)),
              "gram of C and labile N pool is", "{:.2E}".format(Decimal(plant.root_labile_N_pool)), "gram of N.")

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
        # We reset all growth-associated variabes before performing the new growth:
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
        # and we update the labile pools of C and N:
        RootSimple.ArchiSimple_growth_with_C_and_N(plant, root_g,
                                                   satisf_for_C, satisf_for_N,
                                                   time_step_in_seconds=time_step_in_seconds)
        # We proceed to the segmentation of the whole root system
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
        # We update the pools of root labile C and N:
        RootSimple.root_CN_balance(plant, time_step_in_seconds=time_step_in_seconds)

        print("    Labile C pool in the root system is now", "{:.2E}".format(Decimal(plant.root_labile_C_pool)),
              "gC and labile N pool is", "{:.2E}".format(Decimal(plant.root_labile_N_pool)), "gN.")
        print("    The root system contains", root_g.nb_vertices(), "different root elements.")

        # TRANSFER OF N BETWEEN ROOTS AND SHOOTS:
        #########################################

        # print(" --- Before N transfer to shoot, the shoot labile N is", plant.shoot_labile_N_pool)
        # We first compute the possible rate of net N transfer from roots to shoots:
        shoot_root_exchange(plant, computing_C_exchange=False, computing_N_exchange=True,
                            N_option="by_diffusion")
        # We then limit the transfer of N from roots to shoots (or from shoots to roots) to the maximal amount of N
        # available in the source compartment.
        # We first calculate the expected amount of N that should be exchanged based on the flux:
        expected_N_exchange = plant.net_N_transfer_from_root_to_shoot_rate * time_step_in_seconds
        # CASE 1: The flow of N from root to shoot is positive, but there is not enough N in the root to sustain it.
        if expected_N_exchange > (plant.root_labile_N_pool - plant.root_labile_N_deficit):
            # Then the amount of N exchanged is limited to the amount of N available in roots:
            final_N_exchange = (plant.root_labile_N_pool - plant.root_labile_N_deficit)
            # The amount of labile N in shoots and roots is updated according to the new computed exchange:
            plant.shoot_labile_N_pool += final_N_exchange
            plant.root_labile_N_pool -= final_N_exchange
        # CASE 2: The flux from shoot to root is positive but there is not enough N in the shoot to sustain it:
        elif expected_N_exchange < 0. and abs(expected_N_exchange) > (plant.shoot_labile_N_pool - plant.shoot_labile_N_deficit):
            # Then the amount of N exchanged is limited to the amount of N available in shoots:
            final_N_exchange = (plant.shoot_labile_N_pool - plant.shoot_labile_N_deficit)
            # The amount of labile N in shoots and roots is updated according to the new computed exchange:
            plant.shoot_labile_N_pool -= final_N_exchange
            plant.root_labile_N_pool += final_N_exchange
        # CASE 3 (EXPECTED CASE): There is enough N in the source compartment to sustain the transfer.
        else:
            # Then the amount of labile N in shoots and roots is updated according to the computed exchange:
            plant.shoot_labile_N_pool += expected_N_exchange
            plant.root_labile_N_pool -= expected_N_exchange
        # print(" --- After N transfer to shoot, the shoot labile N is", plant.shoot_labile_N_pool)

        # CALLING SUCCESSIVELY THE DIFFERENT SHOOT MODULES:
        ###################################################

        print("    Labile C pool in the shoot system is initially", "{:.2E}".format(Decimal(plant.shoot_labile_C_pool)),
              "gC and labile N pool is", "{:.2E}".format(Decimal(plant.shoot_labile_N_pool)), "gN.")

        # # We first examine whether new phytomers will be produced at this time step.
        # # We reset the authorization to add a new phytomer:
        # phytomer_formation = False
        # # We increment the time elapsed since the last phytomer formation with the current time step:
        # elapsed_time_in_days_since_last_leaf_emergence += time_step_in_days
        # # If the elapsed time since last phytomer formation is higher than the prescribed phyllochrone:
        # if elapsed_time_in_days_since_last_leaf_emergence >= param.phyllochrone_in_days:
        #     # We authorize the formation of a new phytomer:
        #     phytomer_formation = True
        #     # And we reset the elapsed time:
        #     elapsed_time_in_days_since_last_leaf_emergence = elapsed_time_in_days_since_last_leaf_emergence - param.phyllochrone_in_days

        # UPTAKE OF C BY THE SHOOTS:
        # We calculate light interception based on the architecture:
        AllocF.computing_light_interception(shoot_g,
                                            number_of_plants_per_square_meter=plant_density_in_square_meter,
                                            inter_row_distance=inter_row_distance,
                                            scene_length=row_length,
                                            scene_width=scene_width,
                                            leaf_shape=param.leaf_shape
                                            )
        # We calculate the rate of C production based on the PAR intercepted by the plant:
        AllocF.computing_photosynthesis_rate(shoot_g, air_temperature_in_degree_Celsius=air_temperature)

        # MAINTENANCE IN THE SHOOTS:
        # We calculate the C respiration cost associated to maintenance:
        AllocF.computing_maintenance_rate(shoot_g, air_temperature_in_degree_Celsius=air_temperature)

        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if C_allocation_option=="No_priority_to_leaves":
            x=1.
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        # At this stage, we distribute the residual C in each organ according to PAR:
        AllocF.computing_shoot_C_redistribution(shoot_g,
                                                shoot_available_C = plant.shoot_labile_C_pool,
                                                intercepted_light = None, time_step_in_seconds=time_step_in_seconds)

        # We allocate the N available from the base of the shoot up to each leaf, based on PAR:
        AllocF.computing_shoot_N_distribution(shoot_g,
                                              shoot_available_N =  plant.shoot_labile_N_pool,
                                              # shoot_available_N=None,
                                              time_step_in_seconds=time_step_in_seconds)

        # AllocF.shoot_organ_CN_balance(shoot_g, time_step_in_seconds=time_step_in_seconds, resetting_rates=True)

        # GROWTH IN THE SHOOTS:
        # We calculate the possible increase in volume of each organ:
        AllocF.potential_shoot_growth(shoot_g,
                                      SFQN=param.SFQN, #TODO: is SFQN constant over the whole range of (0, LA_max)? Is it two distinct curves?
                                      LA_max=param.LA_max,
                                      leaf_minimal_thickness=param.leaf_minimal_thickness,
                                      elong_internode=param.elong_internode,
                                      time_step_in_seconds=time_step_in_days*60.*60.*24.,
                                      air_temperature_in_degree_Celsius = air_temperature
                                      )
        # We calculate the fraction of C that can be redistributed by source leaves, and distribute a fraction of it
        # in sink leaves:
        AllocF.computing_C_redistribution_for_shoot_growth(shoot_g,
                                                           yield_growth=param.yield_growth,
                                                           structural_mass_density=param.struct_mass_density,
                                                           structural_C_content=param.structural_C_content,
                                                           time_step_in_seconds=time_step_in_days*60.*60.*24.
                                                           )
        # We calculate the actual increase in size and the corresponding C & N costs, as well as the ageing of each organ:
        AllocF.actual_shoot_growth(shoot_g,
                                   yield_growth=param.yield_growth,
                                   structural_mass_density=param.struct_mass_density,
                                   structural_C_content=param.structural_C_content,
                                   structural_N_content=param.structural_N_content,
                                   leaf_minimal_thickness=param.leaf_minimal_thickness,
                                   time_step_in_seconds=time_step_in_days*60.*60.*24.,
                                   air_temperature_in_degree_Celsius=air_temperature
                                   )
        # We now consider the initiation of the structure of new phytomers (if any).
        # We increment the time since the last phytomer was formed at the top of the main stem:
        thermal_time_since_last_phytomer_formation\
            = shoot_g.node(0).thermal_time_since_last_phytomer_formation + thermal_time_step_in_seconds_shoots
        # We get the value of the prescribed phyllochrone (in seconds equivalent to the reference temperature):
        phyllochron = param.phyllochron_in_days * (60 * 60 * 24)
        # print("Time since the last phytomer was formed is now", thermal_time_since_last_phytomer_formation/(60.*60.*24.),
        #       "days, and the phyllochrone is", param.phyllochron_in_days," days.")
        # If the period since the last phytomer was formed is higher than the phyllochrone:
        if thermal_time_since_last_phytomer_formation > phyllochron:
            print("A new phytomer should be formed!")
            # Then at least one phytomer should be formed on the main stem.
            # We keep looping until all phytomers have been formed:
            while thermal_time_since_last_phytomer_formation - phyllochron > 0:
                # We subtract the value of the phyllochrone to the time since last phytomer was formed:
                thermal_time_since_last_phytomer_formation -= phyllochron
                # We add a set of new phytomers in the shoot system:
                AllocF.new_phytomers_initiation(shoot_g, new_initiation=True)
                print("> A new phytomer has been formed!")
                # TODO: Check how resources should be allocated to the initiation of new phytomers, if the dimensions are positive
        # Once all phytomers have been added, we update the time since the last one was formed, which is stored at the base of the shoot system:
        shoot_g.node(0).thermal_time_since_last_phytomer_formation = thermal_time_since_last_phytomer_formation
        # We finally update all dimensions of the organs:
        AllocF.updating_coordinates(shoot_g)

        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if C_allocation_option=="C_allocation_to_roots_before_C_storage":
            # We process the C allocation to roots after photosynthesis and final growth, but before C redistribution:
            x = 1
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        # EXCHANGE WITH RESERVE:
        # This function reallocates all the remaining C available to the different organs according to their PAR,
        # and stores it as additional structural mass:
        AllocF.computing_shoot_C_storage(shoot_g,time_step_in_seconds=time_step_in_days*60.*60.*24.)

        # # This function mobilizes N from the organs and stores it at the base of the shoot system:
        # AllocF.computing_shoot_N_remobilization(shoot_g, k_N_mobilization=param.k_N_mobilization)
        #
        # # This function immobilizes N in each shoot organ:
        # AllocF.computing_shoot_N_immobilization(shoot_g)

        # C&N BALANCE IN THE SHOOTS:
        # This function performs the actual C & N balance and updates the labile pools of C & N in each organ:
        AllocF.shoot_organ_CN_balance(shoot_g, time_step_in_seconds=time_step_in_seconds)

        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if C_allocation_option=="Full_priority_to_leaves":
            # We process the C allocation to roots only at the end, after all shoots have photosynthesized, grown and exchanged C:
            x = 1
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        # SENESCENCE AND LEAF FALL:
        AllocF.computing_senescence(shoot_g)

        # # Updating the shoot system general properties:
        # for vid in shoot_g.vertices_iter(scale=1):
        #     n = shoot_g.node(vid)
        #     try:
        #         plant.shoot_struct_mass += n.struct_mass
        #         plant.shoot_labile_C_pool += n.C_labile
        #         plant.shoot_labile_N_pool += n.N_labile
        #     except:
        #         print("!!! ERROR: we could not add the properties of the shoot element", vid)

        print("    The shoot system now contains", shoot_g.nb_vertices(), "different shoot elements.")
        print("    Labile C pool in the shoot system is now", "{:.2E}".format(Decimal(plant.shoot_labile_C_pool)),
              "gC and labile N pool is", "{:.2E}".format(Decimal(plant.shoot_labile_N_pool)), "gN.")
        print("")

        # TRANSFER OF C BETWEEN SHOOTS AND ROOTS:
        #########################################

        # print(" --- Before C transfer to root, the root labile C is", plant.root_labile_C_pool)
        # We update the properties of the whole plant for the C transfer to shoots:
        shoot_root_exchange(plant, computing_C_exchange=True, computing_N_exchange=False,
                            C_option="by_diffusion")
        # We limit the transfer of C from shoots to roots (or from roots to shoots) to the maximal amount of labile C
        # available in the source compartment.
        # We first calculate the expected amount of C that should be exchanged based on the flux:
        expected_C_exchange = plant.net_C_transfer_from_shoot_to_root_rate * time_step_in_seconds
        # CASE 1: The flux from shoot to root is positive, but there is not enough C in the shoot to sustain it:
        if expected_C_exchange > (plant.shoot_labile_C_pool - plant.shoot_labile_C_deficit):
            # Then the amount of C exchanged is limited to the amount of C available in shoots:
            final_C_exchange = (plant.shoot_labile_C_pool - plant.shoot_labile_C_deficit)
            #  The amount of labile C in shoots and roots is updated according to the new computed exchange:
            plant.shoot_labile_C_pool -= final_C_exchange
            plant.root_labile_C_pool += final_C_exchange
        # CASE 2: The flux from root to shoot is positive but there is not enough C in the root to sustain it:
        elif expected_C_exchange < 0. and abs(expected_C_exchange) > (plant.root_labile_C_pool - plant.root_labile_C_deficit):
            final_C_exchange = (plant.root_labile_C_pool - plant.root_labile_C_deficit)
            # The amount of labile C in shoots and roots is updated according to the new computed exchange:
            plant.shoot_labile_C_pool += final_C_exchange
            plant.root_labile_C_pool -= final_C_exchange
        # CASE 3 (EXPECTED CASE): There is enough C in the source compartment to sustain the transfer:
        else:
            # Then the amount of labile C in shoots and roots is updated according to the initial computed exchange:
            plant.shoot_labile_C_pool -= expected_C_exchange
            plant.root_labile_C_pool += expected_N_exchange

        # PLOTTING ONE SINGLE PLANT OVER TIME:
        ######################################

        # If the camera is to be moved around the scene, we modify the position of the camera:
        if param.camera_rotation is True:
            x_cam = x_coordinates[index_camera]
            y_cam = y_coordinates[index_camera]
            z_cam = z_coordinates[index_camera]
            # We define the index of the coordinates to read at the next step:
            index_camera = index_camera + 1
            # If this index is higher than the number of coordinates in each vector:
            if index_camera >= param.camera_rotation_n_points:
                # Then we reset the index to 0:
                index_camera = 0
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
                             single_color_shoot=[100,200,100],
                             single_color_root=[200,200,100],
                             # single_color_shoot=None,
                             # single_color_root=None,
                             # property_name_shoot='PAR_intercepted', cmap_shoot='inferno', vmin_shoot=0, vmax_shoot=400, lognorm_shoot=False,
                             # property_name_root='radius', cmap_root='jet', vmin_root=1e-4, vmax_root=1e-2, lognorm_root=True,
                             recording_plot=True,
                             plot_name=os.path.join(plot_directory, 'plot_%.5d.png' % step))

        # RECORDING THE PROPERTIES:
        # We can record the main properties of the plant:
        dataframe = recording_MTG_properties(shoot_g, file_name = os.path.join(shoot_MTG_directory, 'MTG_%.5d.csv' % step))

    print("")
    print("Simulation is done!")

    return
#-----------------------------------------------------------------------------------------------------------------------
########################################################################################################################
########################################################################################################################

# MAIN SIMULATION:
##################

if __name__ == '__main__':
    run_simulation(my_path="C:/Users/frees/SIMBAL",
                   simulation_period_in_days=20., time_step_in_days=1/24., starting_time_in_days=0.,
                   input_file="C:/Users/frees/SIMBAL/simulation/inputs/inputs_over_time.csv",
                   input_file_time_step_in_days=1.,
                   forcing_constant_aboveground_environment=False)
