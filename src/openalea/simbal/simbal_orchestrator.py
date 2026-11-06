#  -*- coding: utf-8 -*-

"""
    The script 'simbal_orchestrator' is used to run simulations with both shoots and roots.

    :copyright: see AUTHORS.
    :license: see LICENSE for details.
"""

import openalea.simbal.AllocF as AllocF
import openalea.simbal.RootSimple as RootSimple
import openalea.simbal.simbal_parameters as param
import openalea.rhizodep.parameters as rhizodep_param
import openalea.plantgl.all as pgl

from openalea.simbal.tools import (loading_parameters_from_file, buildDic, formatted_inputs, temperature_modification,
                                   recording_MTG_properties, circle_coordinates, updating_coordinates,
                                   prepare_plants_population, prepare_heterogeneous_plants_population,
                                   create_plants_population_scene, plotting_whole_plant, plotting_whole_scene)

import os
import shutil
import pandas as pd
import numpy as np
from math import sqrt, trunc, ceil
from decimal import Decimal
from sympy import isolate
import ast
import pickle

from src.openalea.simbal.simbal_parameters import plastochron


########################################################################################################################
########################################################################################################################

# We define a class that will contain the properties of the interface between shoot and roots (amounts and rates):
class whole_plant(object):
    """
    This class defines a whole plant, defined with general properties for shoots and roots.
    """
    # We initiate the object:
    def __init__(self):
        self.shoot_struct_mass = 0. # gram of structural mass
        self.shoot_struct_necromass = 0. # gram of structural mass
        self.root_struct_mass = 0. # gram of structural mass
        self.root_struct_necromass = 0. # gram of structural mass

        self.shoot_leaf_surface = 0. # m2

        self.shoot_labile_C_pool = 0. # gC per shoot system
        self.shoot_labile_C_deficit = 0. # gC per shoot system
        self.shoot_reserve_C_pool = 0. # gC per shoot system
        self.shoot_labile_N_pool = 0. # gN per shoot system
        self.shoot_labile_N_deficit = 0. # gN per shoot system
        self.shoot_reserve_N_pool = 0.  # gN per shoot system
        self.root_labile_C_pool = 0.  # gC per root system
        self.root_labile_C_deficit = 0. # gC per root system
        self.root_reserve_C_pool = 0.  # gC per root system
        self.root_labile_N_pool = 0.  # gN per root system
        self.root_labile_N_deficit = 0. # gN per root system
        self.root_reserve_N_pool = 0.  # gN per root system

        self.net_C_transfer_from_shoot_to_root_rate = 0. # gC per plant per second
        self.net_N_transfer_from_root_to_shoot_rate = 0. # gN per plant per second
        self.exchange_section = 0. # m2 of interface between shoots and roots

        self.transpiration_rate = 0. # gram of water per plant per second
        self.total_PAR_intercepted = 0. # micromoles of photons per plant per second
        self.photosynthesized_C_rate = 0. # gC per plant per second
        self.C_supply_from_seed_rate = 0. # gC per plant per second
        self.N_supply_from_seed_rate = 0. # gN per plant per second
        self.N_uptake_rate = 0. # gN per plant per second
        self.N2_fixation_rate = 0.  # gN per plant per second
        self.N2_fixation_rate_C_cost = 0.  # gC per plant per second

        self.shoot_maintenance_respiration_rate = 0. # gC per plant per second
        self.root_maintenance_respiration_rate = 0. # gC per plant per second

        self.actual_C_consumption_rate_by_shoot_growth = 0. # gC per plant per second
        self.actual_N_consumption_rate_by_shoot_growth = 0. # gN per plant per second
        self.actual_C_consumption_rate_by_root_growth = 0. # gC per plant per second
        self.actual_N_consumption_rate_by_root_growth = 0. # gN per plant per second

        self.net_shoot_C_immobilization_rate = 0. # gC per plant per second
        self.net_shoot_N_immobilization_rate = 0. # gN per plant per second
        self.net_root_C_immobilization_rate = 0. # gC per plant per second
        self.net_root_N_immobilization_rate = 0. # gN per plant per second

        self.struct_C_loss_rate_by_leaf_fall = 0. # gC per plant per second
        self.labile_C_loss_rate_by_leaf_fall = 0. # gC per plant per second
        self.struct_N_loss_rate_by_leaf_fall = 0. # gN per plant per second
        self.labile_N_loss_rate_by_leaf_fall = 0. # per plant gN per second
        self.C_rhizodeposition_rate = 0. # gC per plant per second
        self.N_rhizodeposition_rate = 0. # gN per plant per second

        self.root_mortality_rate = 0. # gram of structural mass per second
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# We define a function that computes C & N exchanges between shoot and root systems:
def computing_shoot_root_exchange_rates(plant,
                                        computing_C_exchange=False,
                                        computing_N_exchange=False,
                                        time_step_in_seconds=1*(60*60*24),
                                        C_option="by_diffusion",
                                        N_option="by_diffusion"):

    """
    This function computes the rates of C and/or N exchange between shoots and roots according to different possible
    rules. It updates the rates of exchange within the plant object, but does not alter the pools of labile C and N.
    The rate of C exchange is either controlled by a diffusion process based on the gradient of concentration of labile C
    (gC per gram of dry structural mass) between shoots and roots and a surface of exchange, or by the maximal achievable
    transport (i.e. all remaining labile C in the shoots will be transfered to the roots). The rate of N exchange is
    controlled by similar options, and can also be controlled by the transpiration rate from root xylem to shoot xylem.

    Parameters
    ----------
    plant: plant object containing the variables at the scale of the whole plant
    computing_C_exchange: if True, the function will compute the net C exchange from shoots to roots
    computing_N_exchange: if True, the function will compute the net N exchange from roots to shoots
    time_step_in_seconds: the time step of the simulation (seconds)
    C_option: option to transfer C from shoots to roots - either "by_diffusion" or "full_transfer_from_shoot_to_root"
    N_option: option to transfer N from roots to shoots - either "by_diffusion", "full_transfer_from_shoot_to_root",
                     or "by_water_flux" (in this latter case, the transpiration is driving the transfer of N)

    Returns
    -------
    [The plant object has been updated with the properties "net_C_transfer_from_shoot_to_root_rate"
    and/or "net_N_transfer_from_root_to_shoot_rate" - but the labile C and N pools or their deficits have not been updated]
    """

    # COMPUTING THE EXCHANGE OF CARBON:
    ###################################

    if computing_C_exchange:

        # We first compute the net amount of labile C in shoots and roots, taking into account possible deficits:
        net_shoot_C = (plant.shoot_labile_C_pool - plant.shoot_labile_C_deficit)
        net_root_C = (plant.root_labile_C_pool - plant.root_labile_C_deficit)
        # TODO: WATCH OUT BELOW:
        # To avoid problems linked to deficits, we do not allow the net amounts to be negative because of deficits:
        if net_shoot_C < 0.:
            net_shoot_C = 0.
        if net_root_C < 0.:
            net_root_C = 0.

        # If C is to be exchanged by diffusion based on the labile pools in shoots and roots:
        if C_option == "by_diffusion":
            # Then the net C exchange rate is calculated using a gradient of concentration, a diffusion coefficient,
            # and the section of exchange between shoot and root. We first make sure that there is an exchange surface:
            if plant.exchange_section is None or plant.exchange_section == 0.:
                print("!!! WARNING: The exchange surface between shoots and roots was not correct! It has been set to 1 cm2.")
                plant.exchange_section = 1E-4
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

        # We now verify that the net amount transferred would not be higher than the available amount in the source:
        expected_C_exchange = net_C_transfer_from_shoot_to_root_rate * time_step_in_seconds

        # CASE 1: The expected C transfered from shoot to root is positive, but there is not enough C in the shoot:
        if expected_C_exchange > net_shoot_C:
            # Then the amount of C exchanged is limited to a fraction of the amount of C available in shoots:
            fraction = param.fraction_of_shoot_C_leftover_remaining_in_shoots
            net_final_C_exchange = net_shoot_C * fraction
            print("WATCH OUT: all the C from shoots should have been allocated to roots, i.e.", "{:.2E}".format(Decimal(net_shoot_C)), "gC!")
            print("We therefore allocated only", fraction*100,"% of it to the roots." )
        # CASE 2: The flow from root to shoot is positive but there is not enough C in the root to sustain it:
        elif expected_C_exchange < 0. and abs(expected_C_exchange) > net_root_C:
            # Then the amount of C exchanged is limited to the amount of C available in roots:
            net_final_C_exchange = - net_root_C
            print("WATCH OUT: a net transfer of C from roots to shoots of", "{:.2E}".format(Decimal(-net_root_C)), "gC has occured!")
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
            if plant.exchange_section is None or plant.exchange_section == 0.:
                print("!!! WARNING: The exchange surface between shoots and roots could not be found! It has been set to 1 cm2.")
                plant.exchange_section = 1E-4
            if plant.shoot_struct_mass > 0. and plant.root_struct_mass > 0.:
                net_N_transfer_from_root_to_shoot_rate = ((param.N_diffusion_coefficient
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
            if plant.transpiration_rate > 0.:
                net_N_transfer_from_root_to_shoot_rate = (plant.transpiration_rate / plant.root_water_concentration
                                                            * net_root_N / plant.root_struct_mass)
                # NOTE: N_transfer_rate in gN per plant and per second
                #       transpiration_rate in gram of water per plant and per second
                #       root_water_concentration in gram of water per gDW of root structural mass
                #       net_root_N in gN per plant
                #       root_struct_mass in gDW of root structural mass
            else:
                print("ERROR - the transpiration flux was negative or nil! The flow rate of N from root to shoot was set to 0.")
                net_N_transfer_from_root_to_shoot_rate = 0.
        # Otherwise, if the N present in the root is meant to be transfered in totality to the shoot system:
        elif N_option == "full_transfer_from_root_to_shoot":
            net_N_transfer_from_root_to_shoot_rate = net_root_N / time_step_in_seconds
        else:
            print("WATCH OUT: the option for computing the transfer of N from root to shoot is not recognized! "
                  "The transfer will correspond to the full amount of labile N in the source.")
            net_N_transfer_from_root_to_shoot_rate = net_root_N / time_step_in_seconds

        # We now verify that the net amount transferred would not be higher than the available amount in the source:
        expected_N_exchange = net_N_transfer_from_root_to_shoot_rate * time_step_in_seconds

        # CASE 1: The expected N transferred from root to shoot is positive, but there is not enough N in the root:
        if expected_N_exchange > net_root_N:
            # Then the amount of N exchanged is limited to a fraction of the amount of N available in roots:
            fraction = param.fraction_of_root_N_leftover_remaining_in_roots
            net_final_N_exchange = net_root_N * fraction
            print("WATCH OUT: all the N from roots should have been allocated to shoots!")
            print("We therefore allocated only", fraction*100,"% of it to the roots." )
        # CASE 2: The flow from shoot to root is positive but there is not enough N in the shoot to sustain it:
        elif expected_N_exchange < 0. and -expected_N_exchange > -net_shoot_N:
            # Then the amount of N exchanged is limited to the amount of N available in shoots:
            net_final_N_exchange = - net_shoot_N
            print("WATCH OUT: a net transfer of N from shoots to roots of", "{:.2E}".format(Decimal(-net_shoot_N)), "gN has occured!")
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
# Function for recording the properties of the plant:
def recording_main_plant_properties(plant, shoot_MTG, root_MTG, final_time_in_days, thermal_time_in_degree_days):
    """
    This function records the current properties of the plant object and returns them in a dictionary.

    Parameters
    ----------
    plant: the object corresponding to the whole plant, containing shoot and root-related variables
    final_time_in_days: the time step of the simulation (seconds)

    Returns
    -------
    A dictionary containing the main properties of the plant object (e.g. labile C and N pools, structural mass
    and rates in shoots and roots
    """
    results_dict = {"Final time (days)": final_time_in_days,
                    "Final thermal time (degree-days)": thermal_time_in_degree_days,
                    "Shoot structural mass (gDW)":
                        plant.shoot_struct_mass,
                    "Root structural mass (gDW)":
                        plant.root_struct_mass,
                    "Shoot labile C pool (gC per plant)":
                        plant.shoot_labile_C_pool,
                    "Shoot labile C pool deficit (gC per plant)":
                        plant.shoot_labile_C_deficit,
                    "Shoot labile N pool (gN per plant)":
                        plant.shoot_labile_N_pool,
                    "Shoot labile N pool deficit (gN per plant)":
                        plant.shoot_labile_N_deficit,
                    "Shoot reserve C pool (gC per plant)":
                        plant.shoot_reserve_C_pool,
                    "Shoot reserve N pool (gN per plant)":
                        plant.shoot_reserve_N_pool,
                    "Root labile C pool (gC per plant)":
                        plant.root_labile_C_pool,
                    "Root labile C pool deficit (gC per plant)":
                        plant.root_labile_C_deficit,
                    "Root labile N pool (gN per plant)":
                        plant.root_labile_N_pool,
                    "Root labile N pool deficit (gN per plant)":
                        plant.root_labile_N_deficit,
                    "Root reserve C pool (gC per plant)":
                        plant.root_reserve_C_pool,
                    "Root reserve N pool (gN per plant)":
                        plant.root_reserve_N_pool,
                    "Rate of net shoot-to-root C transfer (gC per plant per second)":
                        plant.net_C_transfer_from_shoot_to_root_rate,
                    "Rate of net root-to-shoot N transfer (gN per plant per second)":
                        plant.net_N_transfer_from_root_to_shoot_rate,
                    "Rate of water transpiration (gram per second)":
                        plant.transpiration_rate,
                    "Rate of light interception (micromol of photons per plant per second)":
                        plant.total_PAR_intercepted,
                    "Rate of photosynthesis (gC per plant per second)":
                        plant.photosynthesized_C_rate,
                    "Rate of C supply from the mother seed (gC per plant per second)":
                        plant.C_supply_from_seed_rate,
                    "Rate of N supply from the mother seed (gN per plant per second)":
                        plant.N_supply_from_seed_rate,
                    "Rate of shoot maintenance respiration (gC per plant per second)":
                        plant.shoot_maintenance_respiration_rate,
                    "Rate of C transformation into shoot structural mass (gC per plant per second)":
                        plant.actual_C_consumption_rate_by_shoot_growth * param.yield_growth,
                    "Rate of C respiration by shoot growth (gC per plant per second)":
                        plant.actual_C_consumption_rate_by_shoot_growth * (1 - param.yield_growth),
                    "Rate of N transformation into shoot structural mass (gN per plant per second)":
                        plant.actual_N_consumption_rate_by_shoot_growth,
                    "Rate of shoot net C immobilization (gC per plant per second)":
                        plant.net_shoot_C_immobilization_rate,
                    "Rate of shoot net N immobilization (gN per plant per second)":
                        plant.net_shoot_N_immobilization_rate,
                    "Rate of shoot structural C loss by leaf fall (gC per plant per second)":
                        plant.struct_C_loss_rate_by_leaf_fall,
                    "Rate of shoot labile C loss by leaf fall (gC per plant per second)":
                        plant.labile_C_loss_rate_by_leaf_fall,
                    "Rate of shoot structural N loss by leaf fall (gC per plant per second)":
                        plant.struct_N_loss_rate_by_leaf_fall,
                    "Rate of shoot labile N loss by leaf fall (gC per plant per second)":
                        plant.labile_N_loss_rate_by_leaf_fall,
                    "Rate of N uptake by roots (gN per plant per second)":
                        plant.N_uptake_rate,
                    "Rate of N2 fixation by roots (gN per plant per second)":
                        plant.N2_fixation_rate,
                    "Rate of C consumption by N2 fixation by roots (gC per plant per second)":
                        plant.N2_fixation_rate_C_cost,
                    "Rate of root maintenance respiration (gC per plant per second)":
                        plant.root_maintenance_respiration_rate,
                    "Rate of C rhizodeposition (gC per plant per second)":
                        plant.C_rhizodeposition_rate,
                    "Rate of N rhizodeposition (gC per plant per second)":
                        plant.N_rhizodeposition_rate,
                    "Rate of C transformation into root structural mass (gC per plant per second)":
                        plant.actual_C_consumption_rate_by_root_growth * param.yield_growth,
                    "Rate of C respiration by root growth (gC per plant per second)":
                        plant.actual_C_consumption_rate_by_root_growth * (1 - param.yield_growth),
                    "Rate of N transformation into root structural mass (gN per plant per second)":
                        plant.actual_N_consumption_rate_by_root_growth,
                    "Rate of root net C immobilization (gC per plant per second)":
                        plant.net_root_C_immobilization_rate,
                    "Rate of root net N immobilization (gN per plant per second)":
                        plant.net_root_N_immobilization_rate,
                    "Number of leaves":
                        len([vid for vid in shoot_MTG.vertices() if shoot_MTG.node(vid).label == 'Leaf']),
                    "Number of internodes":
                        len([vid for vid in shoot_MTG.vertices() if shoot_MTG.node(vid).label == 'Internode']),
                    "Number of root apices":
                        len([vid for vid in root_MTG.vertices() if root_MTG.node(vid).label == 'Apex']),
                    "Number of root segments":
                        len([vid for vid in root_MTG.vertices() if root_MTG.node(vid).label == 'Segment'])
                    }
    return results_dict
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
def balance_check(plant, time_step_in_seconds = 60.*60., genotype_name="",
                  previous_shoot_C=0., previous_root_C=0., previous_shoot_N=0., previous_root_N=0.):
    """
    This function controls whether the information stored in the plant object are coherent in terms of C and N balance,
    by comparing the current amounts of C or N in the shoots and the roots with the amounts expected given the previous
    recorded state and the new inputs and outputs for each compartment.

    Parameters
    ----------
    plant: the plant object where the information is stored.
    time_step_in_seconds: the time step of the simulation (seconds)
    previous_shoot_C: the total amount of C in the shoots determined at the previous time step (gC per plant)
    previous_root_C: the total amount of C in the roots determined at the previous time step (gC per plant)
    previous_shoot_N: the total amount of N in the shoots determined at the previous time step (gN per plant)
    previous_root_N: the total amount of N in the roots determined at the previous time step (gN per plant)

    Returns
    -------
    current_C_and_N_levels: a list containing the current amounts of shoot C, root C, shoot N and root N

    """

    # CHECKING C BALANCE:
    #--------------------
    # We calculate the current amount of C present in the shoots and in the roots:
    shoot_total_C = (plant.shoot_struct_mass * param.structural_C_content
                    + (plant.shoot_labile_C_pool - plant.shoot_labile_C_deficit)
                    + plant.shoot_reserve_C_pool)
    root_total_C = (plant.root_struct_mass * param.structural_C_content
                    + (plant.root_labile_C_pool - plant.root_labile_C_deficit)
                    + plant.root_reserve_C_pool)
    # We calculate all inputs and outputs of C within the shoot system and the root system:
    shoot_C_inputs = (plant.photosynthesized_C_rate
                      + plant.C_supply_from_seed_rate) * time_step_in_seconds
    shoot_C_outputs = (plant.shoot_maintenance_respiration_rate
                       + plant.actual_C_consumption_rate_by_shoot_growth * (1 - param.yield_growth) # i.e. growth respiration rate
                       + plant.struct_C_loss_rate_by_leaf_fall
                       + plant.labile_C_loss_rate_by_leaf_fall
                       + plant.net_C_transfer_from_shoot_to_root_rate) * time_step_in_seconds
    root_C_inputs = (plant.net_C_transfer_from_shoot_to_root_rate) * time_step_in_seconds
    root_C_outputs = (plant.N2_fixation_rate_C_cost
                      + plant.root_maintenance_respiration_rate
                      + plant.actual_C_consumption_rate_by_root_growth * (1 - param.yield_growth)
                      + plant.C_rhizodeposition_rate) * time_step_in_seconds

    # We finally compare the current amount of C with the previous state and the new inputs and outputs:
    difference = ((shoot_total_C + root_total_C)
                  - ((previous_shoot_C + previous_root_C)
                     + (shoot_C_inputs + root_C_inputs)
                     - (shoot_C_outputs + root_C_outputs)))
    if abs(difference) > 1e-16:
        if genotype_name == "":
            print("!!! ERROR: the C balance in the whole plant is wrong!",
                  "The net difference between the current state and the expected state is",
                  "{:.2E}".format(Decimal(difference)), "gC!")
        else:
            print("!!! ERROR: the C balance in the whole plant of", genotype_name,
                  "is wrong! The net difference between the current state and the expected state is",
                  "{:.2E}".format(Decimal(difference)), "gC!")
    difference = shoot_total_C - (previous_shoot_C + shoot_C_inputs - shoot_C_outputs)
    if abs(difference) > 1e-16:
        print("    > The C balance in the shoots is wrong! The net difference between the current state and the expected state is", "{:.2E}".format(Decimal(difference)), "gC!")
    if abs(difference) > 1e-8:
        print("PROBLEM!")
    difference = root_total_C - (previous_root_C + root_C_inputs - root_C_outputs)
    if abs(difference) > 1e-16:
        print("    > The C balance in the roots is wrong! The net difference between the current state and the expected state is", "{:.2E}".format(Decimal(difference)), "gC!")

    # CHECKING N BALANCE:
    #--------------------
    # We calculate the current amount of N present in the shoots and in the roots:
    shoot_total_N = (plant.shoot_struct_mass * param.structural_N_content
                    + (plant.shoot_labile_N_pool - plant.shoot_labile_N_deficit)
                    + plant.shoot_reserve_N_pool)
    root_total_N = (plant.root_struct_mass * param.structural_N_content
                    + (plant.root_labile_N_pool - plant.root_labile_N_deficit)
                    + plant.root_reserve_N_pool)
    # We calculate all inputs and outputs of N within the shoot system and the root system:
    shoot_N_inputs = (plant.N_supply_from_seed_rate
                      + plant.net_N_transfer_from_root_to_shoot_rate) * time_step_in_seconds
    shoot_N_outputs = (plant.struct_N_loss_rate_by_leaf_fall
                       + plant.labile_N_loss_rate_by_leaf_fall) * time_step_in_seconds
    root_N_inputs = (plant.N_uptake_rate + plant.N2_fixation_rate) * time_step_in_seconds
    root_N_outputs = (plant.N_rhizodeposition_rate
                      + plant.net_N_transfer_from_root_to_shoot_rate) * time_step_in_seconds

    # We finally compare the current amount of N with the previous state and the new inputs and outputs:
    difference = ((shoot_total_N + root_total_N)
                  - ((previous_shoot_N + previous_root_N)
                     + (shoot_N_inputs + root_N_inputs)
                     - (shoot_N_outputs + root_N_outputs)))
    if abs(difference) > 1e-15:
        if genotype_name == "":
            print("!!! ERROR: the N balance in the whole plant is wrong!",
                  "The net difference between the current state and the expected state is",
                  "{:.2E}".format(Decimal(difference)), "gN!")
        else:
            print("!!! ERROR: the N balance in the whole plant of", genotype_name,
                  "is wrong! The net difference between the current state and the expected state is",
                  "{:.2E}".format(Decimal(difference)), "gN!")
    difference = shoot_total_N - (previous_shoot_N + shoot_N_inputs - shoot_N_outputs)
    if abs(difference) > 1e-15:
        print("    > The N balance in the shoots is wrong! The net difference between the current state and the expected state is", "{:.2E}".format(Decimal(difference)), "gN!")
    difference = root_total_N - (previous_root_N + root_N_inputs - root_N_outputs)
    if abs(difference) > 1e-15:
        print("    > The N balance in the roots is wrong! The net difference between the current state and the expected state is", "{:.2E}".format(Decimal(difference)), "gN!")

    # We finally return the current amounts of C and N in the shoots and the roots:
    return shoot_total_C, root_total_C, shoot_total_N, root_total_N
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# # Function for calling the main functions in a specific order and run the whole simulation:
# def run_simulation_for_one_plant(my_path=None,
#                                  parameters_file="simbal_parameters.xlsx",
#                                  simulation_period_in_days=20., time_step_in_days=1/24.,
#                                  starting_time_in_days=0., starting_thermal_time_in_degree_days=0.,
#                                  input_file=None,
#                                  input_file_time_step_in_days=1.,
#                                  forcing_constant_aboveground_environment=False,
#                                  constant_air_temperature_in_Celsius=10, constant_light_irradiation=400,
#                                  forcing_constant_belowground_environment=False,
#                                  constant_soil_temperature_in_Celsius=10, constant_soil_available_nitrogen=1,
#                                  plant_density_in_square_meter=35,
#                                  inter_row_distance = 0.45,
#                                  row_length=1.,
#                                  scene_width=1.,
#                                  zoom_factor=100.,
#                                  simulation_option="shoot_and_roots",
#                                  initial_shoot_labile_C_concentration = 1e-2, # in gC per gram of structural mass
#                                  initial_shoot_labile_N_concentration = 1e-3, # in gN per gram of structural mass
#                                  ):
#
#     """
#     This function organizes the main simulation. After plant initiation, it calls the different modules and functions
#     from SIMBAL for each time step, and records the outputs.
#     Parameters
#     ----------
#     my_path
#     simulation_period_in_days
#     time_step_in_days
#     starting_time_in_days
#     input_file
#     input_file_time_step_in_days
#     forcing_constant_aboveground_environment
#     constant_air_temperature_in_Celsius
#     constant_light_irradiation
#     forcing_constant_belowground_environment
#     constant_soil_temperature_in_Celsius
#     constant_soil_available_nitrogen
#     plant_density_in_square_meter
#     inter_row_distance
#     row_length
#     scene_width
#     C_allocation_option
#     initial_shoot_labile_C_concentration
#     initial_shoot_labile_N_concentration
#     initial_root_labile_C_concentration
#     initial_root_labile_N_concentration
#
#     Returns
#     -------
#
#     """
#
#     # SETTING FILE PATHS AND FOLDERS:
#     #--------------------------------
#     if my_path is None:
#         # We define the reference directory (here, we move to the parent directory of the current directory):
#         # my_path = os.path.dirname(os.getcwd())
#         # We first define the path and the file to read as a .csv:
#         my_path = os.path.join('../../..', os.getcwd())
#
#     # We consider the outputs directory:
#     outputs_directory = os.path.join(my_path, "simulation", "outputs")
#     # If the plot directory doesn't exist:
#     if not os.path.exists(outputs_directory):
#         # Then we create it:
#         os.mkdir(outputs_directory)
#     else:
#         # Otherwise, we delete all files that are already present inside:
#         for root, dirs, files in os.walk(outputs_directory):
#             for file in files:
#                 os.remove(os.path.join(root, file))
#
#     # We consider the plot directory:
#     plot_directory = os.path.join(my_path, "simulation", "outputs","plots")
#     # If the plot directory doesn't exist:
#     if not os.path.exists(plot_directory):
#         # Then we create it:
#         os.mkdir(plot_directory)
#     else:
#         # Otherwise, we delete all the images that are already present inside:
#         for root, dirs, files in os.walk(plot_directory):
#             for file in files:
#                 os.remove(os.path.join(root, file))
#
#     # We consider the directory containing the attributes of shoot properties over time:
#     shoot_MTG_properties_directory = os.path.join(my_path, "simulation", "outputs", "shoot_MTG_properties")
#     # If this directory doesn't exist:
#     if not os.path.exists(shoot_MTG_properties_directory):
#         # Then we create it:
#         os.mkdir(shoot_MTG_properties_directory)
#     else:
#         # Otherwise, we delete all the images that are already present inside:
#         for root, dirs, files in os.walk(shoot_MTG_properties_directory):
#             for file in files:
#                 os.remove(os.path.join(root, file))
#
#     # We consider the directory containing the attributes of root properties over time:
#     root_MTG_properties_directory = os.path.join(my_path, "simulation", "outputs", "root_MTG_properties")
#     # If this directory doesn't exist:
#     if not os.path.exists(root_MTG_properties_directory):
#         # Then we create it:
#         os.mkdir(root_MTG_properties_directory)
#     else:
#         # Otherwise, we delete all the images that are already present inside:
#         for root, dirs, files in os.walk(root_MTG_properties_directory):
#             for file in files:
#                 os.remove(os.path.join(root, file))
#
#     # UPDATING PARAMETERS AND INSTRUCTIONS:
#     #--------------------------------------
#     # The values of the parameters in the file 'simbal_parameters.py' are updated according to the values in the file
#     # 'simbal_parameters.xlsx':
#     scenario = loading_parameters_from_file(parameters_filepath = os.path.join(my_path, "simulation", "inputs", parameters_file),
#                                             outputs_path = os.path.join(my_path, "simulation", "outputs"))
#     # We create a dictionary containing the parameters and instructions specified in the scenario:
#     scenario_parameters = buildDic(scenario)
#     # We update the parameters from the original python file of SIMBAL with the parameters taken from this dictionnary:
#     param.__dict__.update(scenario)
#
#     if simulation_option != "shoot_only":
#         # We update the parameters from the original python file of RhizoDep with the parameters taken from this dictionnary:
#         rhizodep_param.__dict__.update(scenario)
#
#     # We read the other specific instructions of the scenario that do not correspond to the parameters in the list
#     # (Note: The following function "get" looks if the first argument is present in the scenario parameters.
#     # If not, it returns the default value indicated in the second argument, which may be the default argument
#     # of the function 'run_simulation' itself):
#     SIMULATION_OPTION = scenario_parameters.get('simulation_option',
#                                                 simulation_option)
#     SIMULATION_PERIOD_IN_DAYS = scenario_parameters.get('simulation_period_in_days',
#                                                         simulation_period_in_days)
#     TIME_STEP_IN_DAYS = scenario_parameters.get('time_step_in_days',
#                                                 time_step_in_days)
#     STARTING_TIME_IN_DAYS = scenario_parameters.get('starting_time_in_days',
#                                                 starting_time_in_days)
#     STARTING_THERMAL_TIME_IN_DAYS = scenario_parameters.get('starting_thermal_time_in_degree_days',
#                                                 starting_thermal_time_in_degree_days)
#     if input_file is None:
#         INPUT_FILE = scenario_parameters.get('input_file',
#                                              input_file)
#     else:
#         INPUT_FILE = input_file
#     INPUT_FILE_TIME_STEP_IN_DAYS = scenario_parameters.get('input_file_time_step_in_days',
#                                                            input_file_time_step_in_days)
#     FORCING_CONSTANT_ABOVEGROUND_ENVIRONMENT = scenario_parameters.get('forcing_constant_aboveground_environment',
#                                                                        forcing_constant_aboveground_environment)
#     CONSTANT_LIGHT_IRRADIATION = scenario_parameters.get('constant_light_irradiation',
#                                                          constant_light_irradiation)
#     CONSTANT_AIR_TEMPERATURE = scenario_parameters.get('constant_air_temperature_in_Celsius',
#                                                        constant_air_temperature_in_Celsius)
#     FORCING_CONSTANT_BELOWGROUND_ENVIRONMENT = scenario_parameters.get('forcing_constant_belowground_environment',
#                                                                        forcing_constant_belowground_environment)
#     CONSTANT_SOIL_TEMPERATURE = scenario_parameters.get('constant_soil_temperature_in_Celsius',
#                                                         constant_soil_temperature_in_Celsius)
#     CONSTANT_SOIL_AVAILABLE_N = scenario_parameters.get('constant_soil_available_nitrogen',
#                                                         constant_soil_available_nitrogen)
#
#     # SETTING INITIAL CONDITIONS
#     #---------------------------
#
#     # INITIALIZING THE WHOLE PLANT:
#     # We create a "plant", i.e. a shoot-root interface without structure that will store variables at the scale of the plant:
#     plant = whole_plant()
#
#     # INITIALIZING SHOOT STRUCTURE:
#     # In case we only want to simulate the behavior of the roots, keeping the shoot properties constant over time:
#     if SIMULATION_OPTION == "roots_only":
#         # We create a MTG "shoot_g" corresponding to the shoot, without any dimensions (it will not be shown in the graph):
#         shoot_g = AllocF.create_shoot_MTG(number_of_phytomers=1,
#                                           initial_internode_height=0.,
#                                           initial_internode_radius=0.,
#                                           initial_leaf_length=0.,
#                                           leaf_surface_max=0.,
#                                           initial_labile_C_concentration=0,
#                                           initial_labile_N_concentration=0,
#                                           sowing_depth=0.
#                                           )
#         # We define an arbitrary shoot structural mass (gram), which will remain constant over the simulation:
#         plant.shoot_struct_mass = 0.1
#         plant.shoot_labile_C_pool = 0.
#         plant.shoot_labile_N_pool = 0.
#     # Otherwise, a fully functional shoot MTG is created:
#     else:
#         # We create a MTG "shoot_g" corresponding to the shoot (stem and leaves):
#         shoot_g = AllocF.create_shoot_MTG(number_of_cotyledons=2, # We define two cotyledons
#                                           number_of_phytomers=0,
#                                           initial_internode_radius=rhizodep_param.D_ini / 2.,
#                                           initial_hypocotyl_length=param.initial_hypocotyl_length,
#                                           initial_internode_height=param.initial_internode_length,
#                                           initial_cotyledon_length=param.initial_cotyledon_length,
#                                           initial_leaf_length=param.initial_leaf_length,
#                                           initial_labile_C_concentration=initial_shoot_labile_C_concentration,
#                                           initial_labile_N_concentration=initial_shoot_labile_N_concentration,
#                                           sowing_depth=param.sowing_depth,
#                                           )
#
#         # We update the shoot system general properties based on the shoot MTG:
#         for vid in shoot_g.vertices_iter(scale=1):
#             n = shoot_g.node(vid)
#             try:
#                 plant.shoot_struct_mass += n.struct_mass
#                 plant.shoot_labile_C_pool += n.labile_C
#                 plant.shoot_labile_C_deficit += n.labile_C_deficit
#                 plant.shoot_labile_N_pool += n.labile_N
#                 plant.shoot_labile_N_deficit += n.labile_N_deficit
#             except:
#                 print("!!! ERROR: we could not compute the initial properties of the shoot element", vid)
#         # We add an extra-property to the base of the shoot system, corresponding to the time since the last phytomer
#         # was formed (used together with the phyllochron later on):
#         shoot_g.node(0).thermal_time_since_last_phytomer_formation = 0.
#         # We record the properties in the initial state of the shoot MTG in an external file:
#         recording_MTG_properties(shoot_g, os.path.join(shoot_MTG_properties_directory, 'MTG_00000.csv'))
#
#     # INITIALIZING ROOT STRUCTURE:
#     # In case we only want to simulate the behavior of the shoot, keeping the root properties constant over time:
#     if SIMULATION_OPTION == "shoot_only":
#         # We create another MTG corresponding to the root system, without any length (it will not appear on the graph):
#         root_g = RootSimple.initiate_mtg(initial_segment_length=0.)
#         # We define an arbitrary initial root mass (gram) that will remain constant over the simulation:
#         plant.root_struct_mass = 0.01
#     # Otherwise, a fully functional root MTG is created:
#     else:
#         # We create another MTG corresponding to the root system:
#         root_g = RootSimple.initiate_mtg(initial_segment_length=0.,
#                                          initial_apex_length=0.001)
#         # We reinitialize the variables in the root MTG:
#         RootSimple.reinitializing_growth_variables(root_g)
#         # We update the distance from tip for each root element in each root axis:
#         RootSimple.update_distance_from_tip(root_g)
#         # We update the surfaces and volumes of each root element:
#         RootSimple.update_surfaces_and_volumes(root_g)
#         # We adjust the properties of the whole root system according to whole root MTG:
#         root_g_properties = RootSimple.summing_and_possibly_homogenizing(root_g)
#         # We record the total root structural mass within the root system:
#         plant.root_struct_mass = root_g_properties["total_living_root_struct_mass"]
#         # We calculate the total exchange surface of the root system:
#         RootSimple.updating_exchange_surface(plant, root_g)
#         # We record the properties in the initial state of the root MTG in an external file:
#         recording_MTG_properties(root_g, os.path.join(root_MTG_properties_directory, 'MTG_00000.csv'))
#
#     # We also initialize the dictionary containing the general properties of the plant:
#     results_dict = recording_main_plant_properties(plant, shoot_g, root_g, final_time_in_days=0., thermal_time_in_degree_days=0.)
#     # We add a few more entries:
#     results_dict["Air temperature (degree_Celsius)"] = "NA"
#     results_dict["Soil temperature (degree_Celsius)"] = "NA"
#     results_dict["Light irradiation (micromol of photons per m2 per second)"] = "NA"
#     results_dict["Soil available N (gN per gDW of soil)"] = "NA"
#     results_dict["Step"] = 0
#     # We initialize a dataframe containing these properties for each time step:
#     plant_data_frame = pd.DataFrame(results_dict, index=[0])
#
#     # And we initialize the C and N balance by computing the initial amounts of C and N in the shoots and roots:
#     shoot_total_C = (plant.shoot_struct_mass * param.structural_C_content
#                      + (plant.shoot_labile_C_pool - plant.shoot_labile_C_deficit)
#                      + plant.shoot_reserve_C_pool)
#     root_total_C = (plant.root_struct_mass * param.structural_C_content
#                     + (plant.root_labile_C_pool - plant.root_labile_C_deficit)
#                     + plant.root_reserve_C_pool)
#     shoot_total_N = (plant.shoot_struct_mass * param.structural_N_content
#                      + (plant.shoot_labile_N_pool - plant.shoot_labile_N_deficit)
#                      + plant.shoot_reserve_N_pool)
#     root_total_N = (plant.root_struct_mass * param.structural_N_content
#                     + (plant.root_labile_N_pool - plant.root_labile_N_deficit)
#                     + plant.root_reserve_N_pool)
#
#     # DEFINING THE PERIOD OF SIMULATION:
#     # ----------------------------------
#
#     # We define the steps of the simulation based on the desired period:
#     initial_time = STARTING_TIME_IN_DAYS
#     simulation_period_in_days = SIMULATION_PERIOD_IN_DAYS
#     time_step_in_seconds = TIME_STEP_IN_DAYS * 60*60*24
#     n_steps = trunc(SIMULATION_PERIOD_IN_DAYS / TIME_STEP_IN_DAYS)
#     # We initialize a counter that will be used over successive time steps together with the phyllochron:
#     elapsed_time_in_days_since_last_leaf_emergence = 0.
#
#     # We initiate actual and thermal time at the end of each step:
#     actual_time_in_days = STARTING_TIME_IN_DAYS
#     thermal_time_in_degree_days = STARTING_THERMAL_TIME_IN_DAYS
#
#     # READING THE INPUT FILE:
#     # -----------------------
#     if INPUT_FILE and not (FORCING_CONSTANT_ABOVEGROUND_ENVIRONMENT and FORCING_CONSTANT_BELOWGROUND_ENVIRONMENT):
#         # We use the function 'formatted inputs' to create a table containing the input data for each time step:
#         input_frame = formatted_inputs(original_input_file=INPUT_FILE,
#                                        final_input_file=os.path.join(outputs_directory, 'updated_input_file.csv'),
#                                        original_time_step_in_days=INPUT_FILE_TIME_STEP_IN_DAYS,
#                                        final_time_step_in_days=TIME_STEP_IN_DAYS,
#                                        simulation_period_in_days=SIMULATION_PERIOD_IN_DAYS,
#                                        do_not_execute_if_file_with_suitable_size_exists=False)
#
#         # We then initialize the step and time according to the first line of the new inputs dataframe:
#         initial_step_number = input_frame['step_number'].loc[0]
#         initial_time_in_days = input_frame['initial_time_in_days'].loc[0]
#     else:
#         # We then initialize the step and time at 1:
#         initial_step_number = 1
#         initial_time_in_days = 0.
#     # In any case, the value of 'step' is initialized:
#     step = initial_step_number
#
#     # CONSIDERING CAMERA ROTATION WHEN PLOTTING:
#     #-------------------------------------------
#     if param.camera_movement is True:
#         # In case the camera needs to rotate around the target and/or to move away from it or to move up/down:
#         if param.step_back_coefficient != 0. or param.step_up_coefficient != 0.:
#             # If the required number of positions around the target is positive:
#             if param.camera_rotation_n_points > 0:
#                 # We calculate the number of loops required around the target to go up to the end of the simulation,
#                 # knowing that the camera will progressively move away or down at each time step while rotating:
#                 number_of_loops = int(ceil(n_steps / param.camera_rotation_n_points))
#             # Otherwise, the number of coordinates required will be directly the number of steps of the simulation:
#             else:
#                 number_of_loops =  n_steps
#         # Otherwise, only one loop is necessary, and the camera will reuse the same coordinates after one complete rotation:
#         else:
#             number_of_loops = 1
#         # We calculate the x,y,z coordinates of the different positions of the camera.
#         [x_coordinates, y_coordinates, z_coordinates] = circle_coordinates(x_center=param.x_center,
#                                                                            y_center=param.y_center,
#                                                                            z_center=param.z_camera,
#                                                                            radius=sqrt((param.x_camera-param.x_center)**2 + (param.y_camera-param.y_center)**2),
#                                                                            n_points=param.camera_rotation_n_points,
#                                                                            n_loops=number_of_loops,
#                                                                            step_back_coefficient=param.step_back_coefficient,
#                                                                            step_up_coefficient=param.step_up_coefficient)
#         # We initialize the index in the list of the coordinates:
#         index_camera = 0
#
#     ####################################################################################################################
#
#     # SIMULATING OVER TIME:
#     #######################
#     #######################
#
#     # For each time step covering the simulation period:
#     for step in range(initial_step_number, initial_step_number + n_steps):
#
#         # We calculate the initial and final time of the current time step (in days):
#         initial_time_in_days = initial_time + (step - 1) * TIME_STEP_IN_DAYS
#         final_time_in_days = initial_time + (step) * TIME_STEP_IN_DAYS
#
#         # DEFINING THE INPUTS AT THIS TIME STEP:
#         # --------------------------------------
#         # If constant aboveground conditions have been chosen:
#         if FORCING_CONSTANT_ABOVEGROUND_ENVIRONMENT:
#             # We set up the constant values:
#             light_input = CONSTANT_LIGHT_IRRADIATION
#             air_temperature = CONSTANT_AIR_TEMPERATURE
#         # Otherwise, if there is an input file to be used:
#         elif input_file:
#             # We extract the information on incident light and air temperature from the input file:
#             light_input = input_frame.loc[step - initial_step_number, 'light_irradiation_in_micromol_per_m2_per_second']
#             air_temperature = input_frame.loc[step - initial_step_number, 'air_temperature_in_degree_Celsius']
#         else:
#             print("!!! ERROR: no input file has been provided!")
#
#         # Similarly, if constant belowground conditions have been chosen:
#         if FORCING_CONSTANT_BELOWGROUND_ENVIRONMENT:
#             soil_temperature = CONSTANT_SOIL_TEMPERATURE
#             soil_available_nitrogen = CONSTANT_SOIL_AVAILABLE_N
#         elif input_file:
#             soil_temperature = input_frame.loc[step - initial_step_number, 'soil_temperature_in_degree_Celsius']
#             soil_available_nitrogen = input_frame.loc[step - initial_step_number, 'soil_available_nitrogen_in_gN_per_gDW_of_soil']
#
#         print("")
#         print("From", "{:.2f}".format(Decimal(starting_time_in_days)),"to", "{:.2f}".format(Decimal(final_time_in_days)),"days:")
#         print("-------------------------")
#         print("    Air temperature is", "{:.2f}".format(Decimal(air_temperature)), "degree Celsius; soil temperature is",
#               "{:.2f}".format(Decimal(soil_temperature)), "degree Celsius.")
#         print("    Incoming light is", "{:.0f}".format(Decimal(light_input)), "micromol of photons per square meter per second.")
#         print("    Soil available N is", "{:.2E}".format(Decimal(soil_available_nitrogen)), "gN per gram of dry soil.")
#         print("")
#
#         # CALCULATING AN EQUIVALENT OF THERMAL TIME:
#         # ------------------------------------------
#         # We calculate a coefficient that will modify the different "ages" experienced by shoots, roots
#         # (and possibly soil organisms) according to temperature:
#         temperature_time_adjustment_shoot = temperature_modification(temperature_in_Celsius=air_temperature,
#                                                                      process_at_T_ref=1,
#                                                                      T_ref=param.shoot_growth_T_ref,
#                                                                      A=param.shoot_growth_A,
#                                                                      B=param.shoot_growth_B,
#                                                                      C=param.shoot_growth_C)
#         temperature_time_adjustment_soil = temperature_modification(temperature_in_Celsius=soil_temperature,
#                                                                     process_at_T_ref=1,
#                                                                     T_ref=param.root_growth_T_ref,
#                                                                     A=param.root_growth_A,
#                                                                     B=param.root_growth_B,
#                                                                     C=param.root_growth_C)
#         # We calculate the equivalent time step in thermal time:
#         thermal_time_step_in_seconds_shoots = time_step_in_seconds * temperature_time_adjustment_shoot
#         # TODO in the future: consider a vertical gradient of temperature, and hence a vertical modification of the adjustments
#
#         # We eventually compute the final and actual thermal time for this time step:
#         actual_time_in_days = final_time_in_days
#         thermal_time_in_degree_days += (thermal_time_step_in_seconds_shoots)/(24.*60.*60.) * param.T_ref_growth
#
#         # RESETTING RATES AT PLANT SCALE:
#         #################################
#
#         # We reset the rates of exchange within the plant to 0:
#         plant.net_C_transfer_from_shoot_to_root_rate = 0.  # gC per second
#         plant.net_N_transfer_from_root_to_shoot_rate = 0.  # gN per second
#         plant.transpiration_rate = 0.  # gram of water per second
#         plant.total_PAR_intercepted = 0.  # micromoles of photons per plant per second
#         plant.photosynthesized_C_rate = 0.  # gC per second
#         plant.C_supply_from_seed_rate = 0. # gC per second
#         plant.N_supply_from_seed_rate = 0. # gN per second
#         plant.N_uptake_rate = 0.  # gN per second
#         plant.N2_fixation_rate = 0.  # gN per second
#         plant.shoot_maintenance_respiration_rate = 0.  # gC per second
#         plant.root_maintenance_respiration_rate = 0.  # gC per second
#         plant.actual_C_consumption_rate_by_shoot_growth = 0.  # gC per second
#         plant.actual_N_consumption_rate_by_shoot_growth = 0.  # gN per second
#         plant.actual_C_consumption_rate_by_root_growth = 0.  # gC per second
#         plant.actual_N_consumption_rate_by_root_growth = 0.  # gN per second
#
#         # TRANSFER OF C FROM SHOOTS TO ROOTS:
#         #####################################
#
#         print("    > Computing shoot-to-root C transfer...")
#         # We compute the rate of net C transfer from shoots to roots:
#         computing_shoot_root_exchange_rates(plant, computing_C_exchange=True, computing_N_exchange=False,
#                                             time_step_in_seconds=time_step_in_seconds,
#                                             C_option="by_diffusion")
#         print("    > The net transfer of C should be", "{:.2E}".format(Decimal(plant.net_C_transfer_from_shoot_to_root_rate * time_step_in_seconds)), "gC.")
#         # Note that only the rate has been computed here, the labile C pools have not been updated yet!
#         # We limit the transfer of C from shoots to roots (or from roots to shoots) to the maximal amount of labile C
#         # available in the source compartment.
#         # We first calculate the expected amount of C that should be exchanged based on the flux:
#         expected_C_exchange = plant.net_C_transfer_from_shoot_to_root_rate * time_step_in_seconds
#
#         # CASE 1: The flux from shoot to root is positive, but there is not enough C in the shoot to sustain it:
#         if expected_C_exchange > 0. and expected_C_exchange > (plant.shoot_labile_C_pool - plant.shoot_labile_C_deficit):
#             # Then the amount of C exchanged is limited to the amount of C available in shoots:
#             final_C_exchange = plant.shoot_labile_C_pool
#             #  The amount of labile C in shoots and roots is updated according to the new computed exchange:
#             plant.shoot_labile_C_pool -= final_C_exchange
#             plant.root_labile_C_pool += final_C_exchange
#         # CASE 2: The flux from root to shoot is positive but there is not enough C in the root to sustain it:
#         elif expected_C_exchange < 0. and abs(expected_C_exchange) > (plant.root_labile_C_pool - plant.root_labile_C_deficit):
#             final_C_exchange = plant.root_labile_C_pool
#             # The amount of labile C in shoots and roots is updated according to the new computed exchange:
#             plant.shoot_labile_C_pool += final_C_exchange
#             plant.root_labile_C_pool -= final_C_exchange
#         # CASE 3 (EXPECTED CASE): There is enough C in the source compartment to sustain the transfer:
#         else:
#             final_C_exchange = expected_C_exchange
#             # Then the amount of labile C in shoots and roots is updated according to the initial computed exchange:
#             plant.shoot_labile_C_pool -= expected_C_exchange
#             plant.root_labile_C_pool += expected_C_exchange
#         print("    ... The net transfer of C was eventually", "{:.2E}".format(Decimal(final_C_exchange)), "gC.")
#
#         # # CASE 1: The flux from shoot to root is positive, but there is not enough C in the shoot to sustain it:
#         # if expected_C_exchange > 0. and expected_C_exchange > (plant.shoot_labile_C_pool - plant.shoot_labile_C_deficit):
#         #     # Then the amount of C exchanged is limited to the amount of C available in shoots:
#         #     final_C_exchange = plant.shoot_labile_C_pool
#         #     #  The amount of labile C in shoots and roots is updated according to the new computed exchange:
#         #     if SIMULATION_OPTION != "roots_only":
#         #         plant.shoot_labile_C_pool -= final_C_exchange
#         #     if SIMULATION_OPTION != "shoot_only":
#         #         plant.root_labile_C_pool += final_C_exchange
#         # # CASE 2: The flux from root to shoot is positive but there is not enough C in the root to sustain it:
#         # elif expected_C_exchange < 0. and abs(expected_C_exchange) > (plant.root_labile_C_pool - plant.root_labile_C_deficit):
#         #     final_C_exchange = plant.root_labile_C_pool
#         #     # The amount of labile C in shoots and roots is updated according to the new computed exchange:
#         #     if SIMULATION_OPTION != "roots_only":
#         #         plant.shoot_labile_C_pool += final_C_exchange
#         #     if SIMULATION_OPTION != "shoot_only":
#         #         plant.root_labile_C_pool -= final_C_exchange
#         # # CASE 3 (EXPECTED CASE): There is enough C in the source compartment to sustain the transfer:
#         # else:
#         #     final_C_exchange = expected_C_exchange
#         #     # Then the amount of labile C in shoots and roots is updated according to the initial computed exchange:
#         #     if SIMULATION_OPTION != "roots_only":
#         #         plant.shoot_labile_C_pool -= expected_C_exchange
#         #     if SIMULATION_OPTION != "shoot_only":
#         #         plant.root_labile_C_pool += expected_C_exchange
#         # print("    ... The net transfer of C was eventually", "{:.2E}".format(Decimal(final_C_exchange)), "gC.")
#
#         # TRANSFER OF N FROM ROOTS TO SHOOTS:
#         #####################################
#
#         print("    > Computing root-to-shoot N transfer...")
#         # We compute the rate of net N transfer from roots to shoots:
#         computing_shoot_root_exchange_rates(plant, computing_C_exchange=False, computing_N_exchange=True,
#                                             time_step_in_seconds=time_step_in_seconds,
#                                             N_option="by_diffusion")
#         print("    > The net transfer of N should be",
#               "{:.2E}".format(Decimal(plant.net_N_transfer_from_root_to_shoot_rate * time_step_in_seconds)),
#               "gN.")
#         # We then limit the transfer of N from roots to shoots (or from shoots to roots) to the maximal amount of N
#         # available in the source compartment.
#         # We first calculate the expected amount of N that should be exchanged based on the exchange rate that has been computed:
#         expected_N_exchange = plant.net_N_transfer_from_root_to_shoot_rate * time_step_in_seconds
#
#         # CASE 1: The flow of N from root to shoot is positive, but there is not enough N in the root to sustain it.
#         if expected_N_exchange > 0. and expected_N_exchange > (plant.root_labile_N_pool - plant.root_labile_N_deficit):
#             # Then the amount of N exchanged is limited to the amount of N available in roots:
#             final_N_exchange = plant.root_labile_N_pool
#             # The amount of labile N in shoots and roots is updated according to the new computed exchange:
#             plant.shoot_labile_N_pool += final_N_exchange
#             plant.root_labile_N_pool -= final_N_exchange
#         # CASE 2: The flux from shoot to root is positive but there is not enough N in the shoot to sustain it:
#         elif expected_N_exchange < 0. and abs(expected_N_exchange) > (plant.shoot_labile_N_pool - plant.shoot_labile_N_deficit):
#             # Then the amount of N exchanged is limited to the amount of N available in shoots:
#             final_N_exchange = plant.shoot_labile_N_pool
#             # The amount of labile N in shoots and roots is updated according to the new computed exchange:
#             plant.shoot_labile_N_pool -= final_N_exchange
#             plant.root_labile_N_pool += final_N_exchange
#         # CASE 3 (EXPECTED CASE): There is enough N in the source compartment to sustain the transfer.
#         else:
#             # Then the amount of labile N in shoots and roots is updated according to the computed exchange:
#             final_N_exchange = expected_N_exchange
#             plant.shoot_labile_N_pool += expected_N_exchange
#             plant.root_labile_N_pool -= expected_N_exchange
#         print("    ... The net transfer of N was eventually", "{:.2E}".format(Decimal(final_N_exchange)), "gC.")
#
#         # # CASE 1: The flow of N from root to shoot is positive, but there is not enough N in the root to sustain it.
#         # if expected_N_exchange > 0. and expected_N_exchange > (plant.root_labile_N_pool - plant.root_labile_N_deficit):
#         #     # Then the amount of N exchanged is limited to the amount of N available in roots:
#         #     final_N_exchange = plant.root_labile_N_pool
#         #     # The amount of labile N in shoots and roots is updated according to the new computed exchange:
#         #     if SIMULATION_OPTION != "roots_only":
#         #         plant.shoot_labile_N_pool += final_N_exchange
#         #     if SIMULATION_OPTION != "shoot_only":
#         #         plant.root_labile_N_pool -= final_N_exchange
#         # # CASE 2: The flux from shoot to root is positive but there is not enough N in the shoot to sustain it:
#         # elif expected_N_exchange < 0. and abs(expected_N_exchange) > (
#         #         plant.shoot_labile_N_pool - plant.shoot_labile_N_deficit):
#         #     # Then the amount of N exchanged is limited to the amount of N available in shoots:
#         #     final_N_exchange = plant.shoot_labile_N_pool
#         #     # The amount of labile N in shoots and roots is updated according to the new computed exchange:
#         #     if SIMULATION_OPTION != "roots_only":
#         #         plant.shoot_labile_N_pool -= final_N_exchange
#         #     if SIMULATION_OPTION != "shoot_only":
#         #         plant.root_labile_N_pool += final_N_exchange
#         # # CASE 3 (EXPECTED CASE): There is enough N in the source compartment to sustain the transfer.
#         # else:
#         #     # Then the amount of labile N in shoots and roots is updated according to the computed exchange:
#         #     if SIMULATION_OPTION != "roots_only":
#         #         final_N_exchange = expected_N_exchange
#         #         plant.shoot_labile_N_pool += expected_N_exchange
#         #     if SIMULATION_OPTION != "shoot_only":
#         #         plant.root_labile_N_pool -= expected_N_exchange
#         # print("    ... The net transfer of N was eventually", "{:.2E}".format(Decimal(final_N_exchange)), "gC.")
#
#         # CALLING SUCCESSIVELY THE DIFFERENT ROOT MODULES:
#         ##################################################
#
#         if SIMULATION_OPTION != "shoot_only":
#
#             print("    Root labile C pool is initially",
#                   "{:.2E}".format(Decimal(plant.root_labile_C_pool - plant.root_labile_C_deficit)),
#                   "gram of C, and labile N pool is",
#                   "{:.2E}".format(Decimal(plant.root_labile_N_pool - plant.root_labile_N_deficit)), "gram of N.")
#
#             # N INPUTS TO THE ROOT SYSTEM:
#             # We compute the rate of N uptake from soil (gN per plant per second):
#             RootSimple.computing_N_uptake_rate(root_system=plant, root_MTG=root_g,
#                                                N_soil = soil_available_nitrogen,
#                                                soil_temperature_in_degree_Celsius=soil_temperature)
#             # We compute the rate of N uptake from atmosphere, i.e. the rate of N2 fixation by root nodules, if any
#             # (gN per plant per second):
#             RootSimple.computing_N2_fixation_rate(root_system=plant, root_MTG=root_g,
#                                                   soil_temperature_in_degree_Celsius=soil_temperature)
#
#             # NET C & N LOSSES FROM THE ROOT SYSTEM (WITHOUT GROWTH):
#             # We compute the rate of rhizodeposition (gC or gN per plant per second):
#             RootSimple.computing_rhizodeposition_rate(root_system=plant, root_MTG=root_g,
#                                                       option="loss_as_input_rate_fraction",
#                                                       C_rhizodeposition_fraction=0.30,
#                                                       N_rhizodeposition_fraction=0.10,
#                                                       time_step_in_seconds=time_step_in_seconds,
#                                                       soil_temperature_in_degree_Celsius=soil_temperature)
#             # We compute the rate of maintenance respiration (gC per plant per second):
#             RootSimple.computing_maintenance_rate(root_system=plant,
#                                                   soil_temperature_in_degree_Celsius=soil_temperature)
#             # We compute the rate of C exchange with reserve (gC per plant per second):
#             RootSimple.computing_exchange_rates_with_reserve(root_system=plant,
#                                                              soil_temperature_in_degree_Celsius=soil_temperature)
#
#             # ROOT GROWTH:
#             # We reset all growth-associated variables before performing the new growth:
#             RootSimple.reinitializing_growth_variables(root_g)
#             # We calculate the potential growth of each root element according to ArchiSimple:
#             RootSimple.potential_growth(root_g, time_step_in_seconds=time_step_in_seconds,
#                                         ArchiSimple=True,
#                                         soil_temperature_in_Celsius=soil_temperature)
#             # We calculate the "satisfaction coefficients" for growth, depending on the amount of labile C and N left:
#             satisf_for_C, satisf_for_N \
#                 = RootSimple.computing_growth_demand_and_satisfaction_coefficients(root_g,
#                                                                                    labile_C_pool=plant.root_labile_C_pool - plant.root_labile_C_deficit,
#                                                                                    labile_N_pool=plant.root_labile_N_pool - plant.root_labile_N_deficit)
#             print("    The satisfaction coefficient for C is", "{:.2f}".format(Decimal(satisf_for_C)),
#                   "and the one for N is", "{:.2f}".format(Decimal(satisf_for_N)))
#             # We use the function ArchiSimple_growth to adapt the potential growth to the available labile C and N,
#             # and we calculate the consumption rates for C & N, without updating the labile pools of C and N:
#             RootSimple.ArchiSimple_growth_with_C_and_N(plant, root_g,
#                                                        satisf_for_C, satisf_for_N,
#                                                        time_step_in_seconds=time_step_in_seconds)
#             # We proceed to the segmentation of the whole root system:
#             RootSimple.segmentation_and_primordia_formation(root_g, time_step_in_seconds,
#                                                             soil_temperature_in_Celsius=soil_temperature, random=True,
#                                                             root_order_limitation=False)
#             # TODO: Include root hairs dynamics once the cost of C & N is properly registered!
#             # # We modifiy root hairs characteristics according to their specific dynamics:
#             # root_hairs_dynamics(g, time_step_in_seconds=time_step_in_seconds,
#             #                     soil_temperature_in_Celsius=soil_temperature)
#             # We update the distance from tip for each root element in each root axis:
#             RootSimple.update_distance_from_tip(root_g)
#             # We update the surfaces and volumes of each root element:
#             RootSimple.update_surfaces_and_volumes(root_g)
#             RootSimple.updating_exchange_surface(plant, root_g)
#
#             print("    Before C&N root balance, root labile C pool is",
#                   "{:.2E}".format(Decimal(plant.root_labile_C_pool - plant.root_labile_C_deficit)),
#                   "gC and root labile N pool is",
#                   "{:.2E}".format(Decimal(plant.root_labile_N_pool - plant.root_labile_N_deficit)), "gN.")
#
#             # FINAL ROOT C&N BALANCE:
#             # We update the pools of root labile C and N at the scale of the whole root system:
#             RootSimple.root_CN_balance(plant, time_step_in_seconds=time_step_in_seconds)
#
#             print("    Root labile C pool is now",
#                   "{:.2E}".format(Decimal(plant.root_labile_C_pool - plant.root_labile_C_deficit)),
#                   "gC and root labile N pool is",
#                   "{:.2E}".format(Decimal(plant.root_labile_N_pool - plant.root_labile_N_deficit)), "gN.")
#             print("    The root system now contains", root_g.nb_vertices()-1, "different root elements.")
#
#         # CALLING SUCCESSIVELY THE DIFFERENT SHOOT MODULES:
#         ###################################################
#
#         if SIMULATION_OPTION != "roots_only":
#
#             # We register the initial state of the labile C & N pools in the shoot:
#             initial_shoot_labile_C = plant.shoot_labile_C_pool - plant.shoot_labile_C_deficit
#             initial_shoot_labile_N = plant.shoot_labile_N_pool - plant.shoot_labile_N_deficit
#
#             print("\n    The shoot system of initially contains", shoot_g.nb_vertices()-1, "elements.")
#             print("    Shoot labile C pool is initially",
#                   "{:.2E}".format(Decimal(initial_shoot_labile_C)),
#                   "gC and shoot labile N pool is",
#                   "{:.2E}".format(Decimal(initial_shoot_labile_N)), "gN.")
#
#             # LIGHT INTERCEPTION:
#             #####################
#
#             # We create an updated list containing one shoot MTG of each type copied and distributed over the scene:
#             scene_information, list_of_scene_shoot_MTG, list_of_scene_root_MTG \
#                 = prepare_heterogeneous_plants_population(list_of_distinct_plant_names=[param.plant_genotype],
#                                                           list_of_distinct_shoot_MTG=[shoot_g],
#                                                           list_of_distinct_root_MTG=[],
#                                                           list_of_numbers_of_plants_per_square_meter=[plant_density_in_square_meter],
#                                                           list_of_inter_row_distances=[inter_row_distance],
#                                                           list_of_number_of_rows_to_show=[5],
#                                                           list_of_scene_lengths=[row_length],
#                                                           list_of_scene_widths=[scene_width],
#                                                           upscaling_factor=zoom_factor)
#
#             # We create an updated scene, with reindexed shapes compared to the original indices in the list of MTG:
#             scene, scene_information = create_plants_population_scene(scene_information=scene_information,
#                                                                       list_of_shoot_MTG=list_of_scene_shoot_MTG,
#                                                                       list_of_root_MTG=[],
#                                                                       leaf_shape="diamond_shape",
#                                                                       simplified_leaf_shape_for_surface_calculation="diamond_shape",
#                                                                       upscaling_factor=zoom_factor)
#
#             # We compute the light interception on the whole scene, and collect the intercepted PAR for each distinct plant type:
#             list_of_total_PAR_intercepted, list_of_sum_of_surfacic_PAR = \
#                 AllocF.computing_light_interception_on_plant_population(scene=scene,
#                                                                         scene_information=scene_information,
#                                                                         list_of_distinct_shoot_MTG=[shoot_g],
#                                                                         list_of_scene_MTG=list_of_scene_shoot_MTG,
#                                                                         incoming_light_density=light_input,
#                                                                         zenithal_light=False,
#                                                                         infinite=False,
#                                                                         averaging_light_intensity_share=True,
#                                                                         period_of_light_intensity_averaging_in_days=0.5,
#                                                                         time_step_in_days=time_step_in_days
#                                                                         )
#
#             # # We calculate light interception based on the shoot architecture, for each organ:
#             # total_PAR_intercepted, sum_of_surfacic_PAR = \
#             #     AllocF.computing_light_interception(shoot_g,
#             #                                         incoming_light_density=light_input,
#             #                                         number_of_plants_per_square_meter=plant_density_in_square_meter,
#             #                                         inter_row_distance=inter_row_distance,
#             #                                         scene_length=row_length,
#             #                                         scene_width=scene_width,
#             #                                         # leaf_shape=param.plotting_leaf_shape,
#             #                                         #  WATCH OUT: here we simplify the shape to reduce computation time!
#             #                                         leaf_shape=param.light_interception_leaf_shape,
#             #                                         averaging_light_intensity_share=True,
#             #                                         period_of_light_intensity_averaging_in_days=0.5,
#             #                                         time_step_in_days=time_step_in_days
#             #                                         )
#             # We record the PAR intercepted in the properties of the plant:
#             plant.total_PAR_intercepted = list_of_total_PAR_intercepted[0]
#
#             # DISTRIBUTING C & N IN EACH SHOOT ORGAN:
#             #########################################
#
#             # Reallocating labile N within each shoot organ:
#             #-----------------------------------------------
#             # The total amount of labile N to redistribute into the shoot organs is defined as:
#             shoot_available_N = plant.shoot_labile_N_pool - plant.shoot_labile_N_deficit
#             if shoot_available_N <= 0.:
#                 print("WATCH OUT: there is no labile N left in the shoot to be redistributed!")
#                 # We cover each organ of the shoot:
#                 for vid in shoot_g.vertices_iter(scale=1):
#                     # We access the current organ "n":
#                     n = shoot_g.node(vid)
#                     if n.surface <= 0.:
#                         continue
#                     labile_N_target = 0.
#                     # We calculate the rate of net N input to this organ corresponding from this redistribution:
#                     n.shoot_N_distribution_rate = (labile_N_target - (
#                             n.labile_N - n.labile_N_deficit)) / time_step_in_seconds
#                     # We now redefine the local pool of labile N according to this target level:
#                     n.labile_N = 0.
#                     n.labile_N_deficit = -labile_N_target
#             else:
#                 # We cover each organ of the shoot:
#                 for vid in shoot_g.vertices_iter(scale=1):
#                     # We access the current organ "n":
#                     n = shoot_g.node(vid)
#                     if n.surface <= 0.:
#                         continue
#                     # We define the fraction of shoot available N allocated to this organ based on the previously
#                     # computed light intensity share:
#                     fraction = n.light_intensity_share
#                     labile_N_target = shoot_available_N * fraction
#                     # We calculate the rate of net N input to this organ corresponding from this redistribution:
#                     n.shoot_N_distribution_rate = (labile_N_target - (
#                             n.labile_N - n.labile_N_deficit)) / time_step_in_seconds
#                     # We now redefine the local pool of labile N according to this target level:
#                     n.labile_N = labile_N_target
#                     n.labile_N_deficit = 0.
#
#             # Redistributing labile C into each shoot organ:
#             #-----------------------------------------------
#             # The total amount of labile C to redistribute into the shoot organs is defined as:
#             shoot_available_C = plant.shoot_labile_C_pool - plant.shoot_labile_C_deficit
#             if shoot_available_C <= 0.:
#                 # We cover each organ n of the shoot:
#                 for vid in shoot_g.vertices_iter(scale=1):
#                     # We access the current organ "n":
#                     n = shoot_g.node(vid)
#                     if n.surface <= 0.:
#                         continue
#                     labile_C_target = 0.
#                     # We calculate the rate of net N input to this organ corresponding from this redistribution:
#                     n.shoot_C_distribution_rate = (labile_C_target - (
#                             n.labile_C - n.labile_C_deficit)) / time_step_in_seconds
#                     # We now redefine the local pool of labile C according to this target level:
#                     n.labile_C_deficit = -labile_C_target
#                     n.labile_C = 0.
#             else:
#                 # We cover each organ n of the shoot:
#                 for vid in shoot_g.vertices_iter(scale=1):
#                     # We access the current organ "n":
#                     n = shoot_g.node(vid)
#                     if n.surface <= 0.:
#                         continue
#                     # We define the fraction of shoot available N allocated to this organ based on the previously
#                     # computed light intensity share:
#                     fraction = n.light_intensity_share
#                     labile_C_target = shoot_available_C * fraction
#                     # We calculate the rate of net N input to this organ corresponding from this redistribution:
#                     n.shoot_C_distribution_rate = (labile_C_target - (
#                             n.labile_C - n.labile_C_deficit)) / time_step_in_seconds
#                     # We now redefine the local pool of labile C according to this target level:
#                     n.labile_C = labile_C_target
#                     n.labile_C_deficit = 0.
#
#             # PHOTOSYNTHESIS:
#             # We calculate the rate of C production of each organ based on the PAR intercepted by the plant:
#             total_PAR_intercepted, total_photosynthesis_rate = \
#                 AllocF.computing_photosynthesis_rate(plant, shoot_g,
#                                                      air_temperature_in_degree_Celsius=air_temperature)
#             # We record this total photosynthesis rate into the plant object:
#             plant.photosynthesized_C_rate = total_photosynthesis_rate
#             print("    Photosynthesis will generate",
#                   "{:.2E}".format(Decimal(total_photosynthesis_rate*time_step_in_seconds)),
#                   "gC over this time step...")
#             # Note: at this stage, the labile C amounts have not been updated in the organs or at the scale of the plant!
#
#             # POSSIBILE SUPPLY OF RESOURCES FROM THE GERMINATING SEED:
#             # If time is lower than the prescribed duration of supply from the seed:
#             if thermal_time_in_degree_days * (60.*60.*24.) <= param.seed_supply_duration * param.T_ref_growth:
#                 # Then we allow to compute the supply of C and N resources from the germinating seed to the new seedling:
#                 AllocF.computing_supply_from_seed_rate(plant, air_temperature_in_degree_Celsius = air_temperature)
#
#             # MAINTENANCE IN THE SHOOT ORGANS:
#             # We calculate the C respiration cost associated to maintenance in each organ:
#             total_maintenance_rate = AllocF.computing_maintenance_rate(shoot_g,
#                                                                        air_temperature_in_degree_Celsius=air_temperature)
#             plant.shoot_maintenance_respiration_rate = total_maintenance_rate
#             # print("    The maintenance will consume",
#             #       "{:.2E}".format(Decimal(total_maintenance_rate*time_step_in_seconds)),
#             #       "gC over the shoot system...")
#
#             # GROWTH IN THE SHOOT ORGANS:
#             #############################
#
#             # ORGAN EXTENSION; We calculate the possible increase in volume of each organ:
#             AllocF.potential_shoot_growth(shoot_g,
#                                           SFQN=param.SFQN, # TODO: is SFQN constant over the whole range of (0, LA_max)? Is it two distinct curves?
#                                           leaf_minimal_thickness=param.leaf_minimal_thickness,
#                                           leaf_shape=param.simplified_leaf_shape_for_surface_calculation,
#                                           elong_hypocotyl=param.elong_internode * 10.,
#                                           elong_internode_slow=param.elong_internode / 10.,
#                                           elong_internode_fast=param.elong_internode,
#                                           time_step_in_seconds=TIME_STEP_IN_DAYS * 60. * 60. * 24.,
#                                           air_temperature_in_degree_Celsius=air_temperature)
#             # We calculate the actual increase in size and the corresponding C & N costs, as well as the ageing of each
#             # organ, based on the local amounts of labile C:
#             C_consumption_rate_for_extension, N_consumption_rate_for_extension = \
#                 AllocF.actual_shoot_growth(shoot_g,
#                                            yield_growth=param.yield_growth,
#                                            structural_mass_density=param.struct_mass_density,
#                                            structural_C_content=param.structural_C_content,
#                                            structural_N_content=param.structural_N_content,
#                                            leaf_minimal_thickness=param.leaf_minimal_thickness,
#                                            leaf_shape=param.simplified_leaf_shape_for_surface_calculation,
#                                            time_step_in_seconds=TIME_STEP_IN_DAYS * 60. * 60. * 24.,
#                                            air_temperature_in_degree_Celsius=air_temperature)
#             # At this stage, the shoot structure has been updated and the consumption rates have been calculated
#             # for each organ, but no actual consumption of labile C&N at the scale of organ or shoot system has occured yet.
#             print("    The extension of organs will consume",
#                   "{:.2E}".format(Decimal(C_consumption_rate_for_extension * time_step_in_seconds)),
#                   "gC and",
#                   "{:.2E}".format(Decimal(N_consumption_rate_for_extension * time_step_in_seconds)),
#                   "gN over the shoot system...")
#
#             # NEW PHYTOMERS: We now consider the initiation of the structure of new phytomers, based on the phyllochrone.
#             # We increment the time since the last phytomer was formed at the top of the main stem:
#             thermal_time_since_last_phytomer_formation \
#                 = shoot_g.node(0).thermal_time_since_last_phytomer_formation + thermal_time_step_in_seconds_shoots
#             # We get the value of the prescribed plastochron and phyllochron (in seconds equivalent to the reference temperature):
#             plastochron = param.plastochron
#             phyllochron = param.phyllochron
#             # We initialize the total costs of phytomers initiation:
#             total_C_costs_for_new_phytomers = 0.
#             total_N_costs_for_new_phytomers = 0.
#             # If the period since the last phytomer was formed is higher than the phyllochron:
#             if thermal_time_since_last_phytomer_formation > phyllochron:
#                 print("A new phytomer should be formed!")
#                 # Then at least one phytomer should be formed on the main stem.
#                 # We keep looping until all phytomers have been formed:
#                 while thermal_time_since_last_phytomer_formation - phyllochron > 0:
#                     # We calculate the exact thermal time elapsed since the phytomer has actually emerged:
#                     thermal_time_since_emergence\
#                         = (thermal_time_since_last_phytomer_formation
#                            - np.floor(thermal_time_since_last_phytomer_formation / phyllochron) * phyllochron)
#                     # We compute then the exact thermal time elapsed since the phytomer was initiated in the apex,
#                     # which is linearily correlated to the time since emergence by the following relationship:
#                     # t_initiation = plastochron + t_emergence * plastochron  / phyllochron
#                     thermal_time_since_initiation = plastochron + thermal_time_since_emergence * plastochron / phyllochron
#                     # The actual times are then calculated:
#                     actual_time_since_emergence = thermal_time_since_emergence / temperature_time_adjustment_shoot
#                     #TODO: WATCH OUT BELOW - may not be correct!
#                     actual_time_since_initiation = thermal_time_since_initiation / temperature_time_adjustment_shoot
#                     # We add a set of new phytomers in the shoot system:
#                     C_cost_for_new_phytomers, N_cost_for_new_phytomers = AllocF.new_phytomers_emergence(
#                         shoot_system=plant,
#                         shoot_MTG=shoot_g,
#                         thermal_time_since_initiation=thermal_time_since_initiation,
#                         actual_time_since_initiation=actual_time_since_initiation,
#                         thermal_time_since_emergence=thermal_time_since_emergence,
#                         actual_time_since_emergence=actual_time_since_emergence,
#                         initial_leaf_length=param.initial_leaf_length,
#                         initial_internode_length=param.initial_internode_length,
#                         time_step_in_seconds=time_step_in_seconds)
#                     print("> A new phytomer has been formed! The C & N costs were",
#                           "{:.2E}".format(Decimal(C_cost_for_new_phytomers)), "gC and",
#                           "{:.2E}".format(Decimal(N_cost_for_new_phytomers)),"gN.")
#                     # We increase the total costs with the cost of this new phytomer:
#                     total_C_costs_for_new_phytomers += C_cost_for_new_phytomers
#                     total_N_costs_for_new_phytomers += total_N_costs_for_new_phytomers
#                     # We now subtract the value of the phyllochron to the time since last phytomer was formed, and loop:
#                     thermal_time_since_last_phytomer_formation -= phyllochron
#             # Once all phytomers have been added, we update the time since the last one was formed,
#             # which is stored at the base of the shoot system:
#             shoot_g.node(0).thermal_time_since_last_phytomer_formation = thermal_time_since_last_phytomer_formation
#
#             # RAMIFICATION - We consider the possibility of ramification:
#             C_cost_for_new_ramification, N_cost_for_new_ramification \
#                 = AllocF.computing_ramification(shoot_MTG=shoot_g,
#                                                 initial_internode_length=param.initial_internode_length,
#                                                 initial_leaf_length=param.initial_leaf_length,
#                                                 time_step_in_seconds=time_step_in_seconds)
#             if C_cost_for_new_ramification > 0. or N_cost_for_new_ramification > 0.:
#                 print("> The costs to sustain ramification with the creation of new phytomers were",
#                       "{:.2E}".format(Decimal(C_cost_for_new_ramification)), "gC and",
#                       "{:.2E}".format(Decimal(N_cost_for_new_ramification)), "gN.")
#
#             # We finally update the coordinates of the new organs:
#             updating_coordinates(shoot_g)
#             # And we register the total consumption rate for growth in the plant object:
#             plant.actual_C_consumption_rate_by_shoot_growth \
#                 = (C_consumption_rate_for_extension
#                    + (C_cost_for_new_ramification + total_C_costs_for_new_phytomers) / time_step_in_seconds)
#             plant.actual_N_consumption_rate_by_shoot_growth \
#                 = (N_consumption_rate_for_extension
#                    + (N_cost_for_new_ramification + total_N_costs_for_new_phytomers) / time_step_in_seconds)
#
#             if plant.actual_C_consumption_rate_by_shoot_growth > 0.:
#                 print("...")
#
#
#             # We perform the C&N balance both within each organ and at the shoot system scale:
#             AllocF.shoot_CN_balance_at_organ_and_system_scales(shoot_system=plant,
#                                                                shoot_MTG=shoot_g,
#                                                                time_step_in_seconds=time_step_in_seconds)
#
#             print("    The shoot system now contains", shoot_g.nb_vertices()-1, "different shoot elements.")
#             print("    Shoot labile C pool is now",
#                   "{:.2E}".format(Decimal(plant.shoot_labile_C_pool - plant.shoot_labile_C_deficit)),
#                   "gC and labile N pool is",
#                   "{:.2E}".format(Decimal(plant.shoot_labile_N_pool - plant.shoot_labile_N_deficit)), "gN.")
#
#         # AGEING AND CONSEQUENCES:
#         ##########################
#
#         AllocF.geometry_evolution(shoot_MTG = shoot_g,
#                                   plant_thermal_age_in_degree_days = thermal_time_in_degree_days)
#
#         # CHECKING C & N BALANCE:
#         #########################
#
#         # We compute the current amounts of total C and N in the shoots and roots, and compare it with the previous ones
#         # when altered by the new inputs and outputs computed at this time step:
#         shoot_total_C, root_total_C, shoot_total_N, root_total_N \
#             = balance_check(plant,
#                             time_step_in_seconds = time_step_in_seconds,
#                             previous_shoot_C = shoot_total_C,
#                             previous_root_C = root_total_C,
#                             previous_shoot_N = shoot_total_N,
#                             previous_root_N = root_total_N)
#         print("")
#
#         # PLOTTING ONE SINGLE PLANT OVER TIME:
#         ######################################
#
#         # If the camera is to be moved around the scene, we modify the position of the camera:
#         if param.camera_movement is True:
#             x_cam = x_coordinates[index_camera]
#             y_cam = y_coordinates[index_camera]
#             z_cam = z_coordinates[index_camera]
#             # We define the index of the coordinates to read at the next step:
#             index_camera = index_camera + 1
#             # If this index is higher than the number of coordinates in each vector and if no other movement
#             # of the camera (away or down) has been set over time:
#             if (index_camera >= param.camera_rotation_n_points
#                     and param.step_back_coefficient == 0. and param.step_up_coefficient == 0.):
#                 # Then we reset the index to 0, so that the camera goes back to the same locations as previously stored:
#                 index_camera = 0
#             # If "distancing" is required, the camera will continue to read the further coordinates already calculated.
#         # Otherwise, we simply use the parameters specified for the camera:
#         else:
#             x_cam = param.x_camera
#             y_cam = param.y_camera
#             z_cam = param.z_camera
#
#         # We display the MTG in its new state:
#         plotting_whole_plant(shoot_MTG=shoot_g, root_MTG=root_g,
#                              width=1200, height=1200,
#                              upscaling_factor=zoom_factor,
#                              origin_coordinates=[0.,0.,-param.sowing_depth],
#                              x_center=param.x_center*zoom_factor,
#                              y_center=param.y_center*zoom_factor,
#                              z_center=param.z_center*zoom_factor,
#                              x_cam=x_cam*zoom_factor, y_cam=y_cam*zoom_factor, z_cam=z_cam*zoom_factor,
#                              grid=False,
#                              plant_species=param.plant_genotype,
#                              leaf_shape=param.plotting_leaf_shape,
#                              simplified_leaf_shape_for_surface_calculation=param.simplified_leaf_shape_for_surface_calculation,
#                              # background_color=[94, 76, 64],
#                              background_color=[255,255,255],
#                              single_color_shoot=[60,100,40],
#                              single_color_root=[250,150,100],
#                              # single_color_shoot=None,
#                              # single_color_root=None,
#                              # property_name_shoot='surfacic_PAR_intercepted', cmap_shoot='cividis', vmin_shoot=0., vmax_shoot=1000, lognorm_shoot=False,
#                              # property_name_shoot='PAR_intercepted', cmap_shoot='inferno', vmin_shoot=0, vmax_shoot=400, lognorm_shoot=False,
#                              # property_name_root='radius', cmap_root='jet', vmin_root=1e-4, vmax_root=1e-2, lognorm_root=True,
#                              recording_plot=True,
#                              plot_name=os.path.join(plot_directory, 'plot_%.5d.png' % step),
#                              shoots_overlapping_roots=True)
#
#         # RECORDING THE PROPERTIES OF THE SHOOT MTG FOR THIS TIME STEP:
#         # We can record the main properties of the plant:
#         recording_MTG_properties(shoot_g, file_name = os.path.join(shoot_MTG_properties_directory, 'MTG_%.5d.csv' % step))
#
#         # RECORDING THE PROPERTIES OF THE ROOT MTG FOR THIS TIME STEP:
#         # We can record the main properties of the plant:
#         recording_MTG_properties(root_g, file_name = os.path.join(root_MTG_properties_directory, 'MTG_%.5d.csv' % step))
#
#         # RECORDING THE FINAL RESULTS AT PLANT SCALE:
#         results_dict = recording_main_plant_properties(plant, shoot_g, root_g,
#                                                        final_time_in_days=actual_time_in_days,
#                                                        thermal_time_in_degree_days=thermal_time_in_degree_days)
#         results_dict["Air temperature (degree_Celsius)"] = air_temperature
#         results_dict["Soil temperature (degree_Celsius)"] = soil_temperature
#         results_dict["Light irradiation (micromol of photons per m2 per second)"] = light_input
#         results_dict["Soil available N (gN per gDW of soil)"] = soil_available_nitrogen
#         results_dict["Step"] = step
#         new_df = pd.DataFrame([results_dict.values()], columns=results_dict.keys())
#         plant_data_frame = pd.concat([plant_data_frame, new_df])
#         # We record the dataframe as a csv file:
#         try:
#             plant_data_frame.to_csv(os.path.join(outputs_directory,'plant_results.csv'), na_rep='NA', index=False, header=True)
#         except:
#             plant_data_frame.to_csv(os.path.join(outputs_directory,'plant_results_BACKUP.csv'), na_rep='NA', index=False, header=True)
#
#     print("")
#     print("Simulation is done!")
#
#     return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function for calling the main functions in a specific order and run the whole simulation:
def run_simulation_for_a_plant_cover(my_path=None):

    """
    Returns
    -------
    """

    # We make sure that there is a path to the main folder:
    if my_path is None:
        # We define the reference directory (here, we move to the parent directory of the current directory):
        # my_path = os.path.dirname(os.getcwd())
        # We first define the path and the file to read as a .csv:
        my_path = os.path.join('../../..', os.getcwd())

    # We initialize an upscaling factor that allows to have better resolution on PlantGL graphs:
    zoom_factor = 100

    # We initiate dictionaries that will contain items for each of the distinct genotype to be treated:
    parameters_dict = {}
    shoot_MTG_dict = {}
    root_MTG_dict = {}
    growth_stages_dict = {}
    outputs_directories = {}
    plot_directories = {}
    shoot_MTG_properties_directories = {}
    root_MTG_properties_directories = {}
    shoot_MTG_directories = {}
    root_MTG_directories = {}
    plant_objects = {}
    plant_data_frame_dict = {}
    shoot_total_C_dict = {}
    root_total_C_dict = {}
    shoot_total_N_dict = {}
    root_total_N_dict = {}
    cam_x_dict = {}
    cam_y_dict = {}
    cam_z_dict = {}
    cam_index_dict = {}
    recording_steps_list_dict = {}

    # We initially consider the simulation as possible, and will change this if a condition is not met:
    simulation_allowed = True
    # We initially consider no light computation, and will change this if the instructions are not "roots_only"
    # for both the general population and for each genotype:
    computing_light_interception = False

    # LOADING THE MAIN PARAMETERS AT THE SCALE OF THE POPULATION:
    #------------------------------------------------------------

    # We get the general parameters related to the population of plants:
    parameters = loading_parameters_from_file(
        parameters_filepath=os.path.join(my_path, "simulation", "inputs", "simbal_parameters_population.xlsx"),
        outputs_path=os.path.join(my_path, "simulation", "outputs"))
    # We store this set of parameters in the parameters dictionary:
    parameters_dict["population"] = parameters
    # We update the general parameters with this set of parameters:
    param.__dict__.update(parameters)
    # If the shoots are to be computed:
    if parameters["simulation_option"] != "roots_only":
        # Then we will need to compute light interception:
        computing_light_interception = True

    # We get scene-related specific parameters:
    try:
        # We get the list of plant genotype names:
        genotype_names_list = ast.literal_eval(parameters["plant_genotype"])
        # (NOTE: ast.literal_eval allows to convert string into a Python list or dictionary)
        # We get information about the way plants are to be positionned in the scene:
        numbers_of_plants_per_square_meter_dict=ast.literal_eval(parameters["plant_density"])
        sowing_depths_dict = ast.literal_eval(parameters["sowing_depth"])
        inter_row_distances_dict=ast.literal_eval(parameters["inter_row_distance"])
        number_of_rows_to_show_dict=ast.literal_eval(parameters["number_of_rows_to_show"])
        scene_lengths_dict=ast.literal_eval(parameters["scene_length"])
        scene_widths_dict=ast.literal_eval(parameters["scene_width"])
    except:
        print("ERROR: failure to import lists and dictionaries from the population parameter file!")
        simulation_allowed = False

    # READING INSTRUCTIONS FROM THE MAIN PARAMETER FILE:
    simulation_period_in_days = param.simulation_period_in_days
    time_step_in_days = param.time_step_in_days
    starting_time_in_days = param.starting_time_in_days
    starting_thermal_time_in_degree_days = param.starting_thermal_time_in_degree_days
    input_file = param.input_file
    input_file_time_step_in_days = param.input_file_time_step_in_days
    forcing_constant_aboveground_environment = param.forcing_constant_aboveground_environment
    constant_air_temperature_in_Celsius = param.constant_air_temperature_in_Celsius
    constant_light_irradiation = param.constant_light_irradiation
    forcing_constant_belowground_environment = param.forcing_constant_belowground_environment
    constant_soil_temperature_in_Celsius = param.constant_soil_temperature_in_Celsius
    constant_soil_available_nitrogen = param.constant_soil_available_nitrogen
    plotting_population = param.plotting_population
    plotting_single_plant = param.plotting_single_plant

    # We create a list of names containing first the population, then the names of each genotype:
    population_and_genotypes_names_list = ["population"]
    population_and_genotypes_names_list.extend(genotype_names_list)

    # SETTING FILE PATHS AND FOLDERS:
    #--------------------------------

    # We consider a general output directory:
    outputs_directory = os.path.join(my_path, "simulation", "outputs")
    # We add this output directory to the dictionary of output directories:
    outputs_directories["population"] = outputs_directory
    # If the directory doesn't exist:
    if not os.path.exists(outputs_directory):
        # Then we create it:
        os.mkdir(outputs_directory)
    else:
        # We remove the folder and all its components:
        shutil.rmtree(outputs_directory)
        # Then we recreate it:
        os.mkdir(outputs_directory)

    # We consider a plot directory for showing the whole population:
    plot_directory = os.path.join(my_path, "simulation", "outputs", "population_plots")
    # We add this plot directory to the dictionary of plot directories:
    plot_directories["population"] = plot_directory
    # If the directory doesn't exist:
    if not os.path.exists(plot_directory):
        # Then we create it:
        os.mkdir(plot_directory)
    else:
        # Otherwise, we delete all files that are already present inside:
        for root, dirs, files in os.walk(plot_directory):
            for file in files:
                os.remove(os.path.join(root, file))

    #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # FOR EACH DISTINCT PLANT GENOTYPE TO BE SIMULATED:
    for plant_genotype_name in genotype_names_list:

        # We consider the outputs directory for the current MTG type:
        sub_outputs_directory = os.path.join(outputs_directory, "outputs_" + plant_genotype_name)
        outputs_directories[plant_genotype_name]=sub_outputs_directory
        # If the directory doesn't exist:
        if not os.path.exists(sub_outputs_directory):
            # Then we create it:
            os.mkdir(sub_outputs_directory)
        else:
            # Otherwise, we delete all files that are already present inside:
            for root, dirs, files in os.walk(sub_outputs_directory):
                for file in files:
                    os.remove(os.path.join(root, file))

        # We consider the plot subdirectory:
        plot_directory = os.path.join(sub_outputs_directory, "plots")
        plot_directories[plant_genotype_name]=plot_directory
        # If the directory doesn't exist:
        if not os.path.exists(plot_directory):
            # Then we create it:
            os.mkdir(plot_directory)
        else:
            # Otherwise, we delete all the images that are already present inside:
            for root, dirs, files in os.walk(plot_directory):
                for file in files:
                    os.remove(os.path.join(root, file))

        # We consider the directory containing the attributes of shoot properties over time:
        shoot_MTG_properties_directory = os.path.join(sub_outputs_directory, "shoot_MTG_properties")
        shoot_MTG_properties_directories[plant_genotype_name] = shoot_MTG_properties_directory
        # If this directory doesn't exist:
        if not os.path.exists(shoot_MTG_properties_directory):
            # Then we create it:
            os.mkdir(shoot_MTG_properties_directory)
        else:
            # Otherwise, we delete all the images that are already present inside:
            for root, dirs, files in os.walk(shoot_MTG_properties_directory):
                for file in files:
                    os.remove(os.path.join(root, file))

        # We consider the directory containing the attributes of root properties over time:
        root_MTG_properties_directory = os.path.join(sub_outputs_directory, "root_MTG_properties")
        root_MTG_properties_directories[plant_genotype_name] = root_MTG_properties_directory
        # If this directory doesn't exist:
        if not os.path.exists(root_MTG_properties_directory):
            # Then we create it:
            os.mkdir(root_MTG_properties_directory)
        else:
            # Otherwise, we delete all the images that are already present inside:
            for root, dirs, files in os.walk(root_MTG_properties_directory):
                for file in files:
                    os.remove(os.path.join(root, file))

        # We consider the directory containing the shoot MTG files recorded over time:
        shoot_MTG_files_saving_directory = os.path.join(sub_outputs_directory, "shoot_MTG_files")
        shoot_MTG_directories[plant_genotype_name] = shoot_MTG_files_saving_directory
        # If this directory doesn't exist:
        if not os.path.exists(shoot_MTG_files_saving_directory):
            # Then we create it:
            os.mkdir(shoot_MTG_files_saving_directory)
        else:
            # Otherwise, we delete all the images that are already present inside:
            for root, dirs, files in os.walk(shoot_MTG_files_saving_directory):
                for file in files:
                    os.remove(os.path.join(root, file))

        # We consider the directory containing the root MTG files recorded over time:
        root_MTG_files_saving_directory = os.path.join(sub_outputs_directory, "root_MTG_files")
        root_MTG_directories[plant_genotype_name] = root_MTG_files_saving_directory
        # If this directory doesn't exist:
        if not os.path.exists(root_MTG_files_saving_directory):
            # Then we create it:
            os.mkdir(root_MTG_files_saving_directory)
        else:
            # Otherwise, we delete all the images that are already present inside:
            for root, dirs, files in os.walk(root_MTG_files_saving_directory):
                for file in files:
                    os.remove(os.path.join(root, file))

        # UPDATING PARAMETER FOR THE CURRENT GENOTYPE:
        #---------------------------------------------

        # The values of the parameters in the simulation for the given MTG are updated:
        parameters = loading_parameters_from_file(
            parameters_filepath=os.path.join(main_simbal_folder, "simulation", "inputs",
                                             "simbal_parameters_" + plant_genotype_name + ".xlsx"),
            outputs_path=os.path.join(main_simbal_folder, "simulation", "outputs",
                                      "outputs_" + plant_genotype_name))
        param.__dict__.update(parameters)
        rhizodep_param.__dict__.update(parameters)
        # We store this set of parameters in the parameters dictionary under the corresponding genotype name:
        parameters_dict[plant_genotype_name] = parameters
        # If the shoots are to be computed:
        if parameters["simulation_option"] != "roots_only":
            # Then we will need to compute light interception:
            computing_light_interception = True

        # CREATING OR LOADING MTG FILES:
        #-------------------------------

        # We check whether existing shoot or root MTG files should be loaded or not:
        start_from_a_known_shoot_MTG = param.start_from_a_known_shoot_MTG
        start_from_a_known_root_MTG = param.start_from_a_known_root_MTG

        # SHOOT MTG: -If we decide to start from a given shoot MTG to be loaded:
        if start_from_a_known_shoot_MTG and param.simulation_option != "roots_only":
            filename = os.path.join(my_path, "simulation", "inputs", 'initial_shoot_MTG_' + plant_genotype_name + '.pckl')
            if os.path.exists(filename):
                # We load the MTG file and name it "g":
                f = open(filename, 'rb')
                shoot_MTG_dict[plant_genotype_name] = pickle.load(f)
                f.close()
                print("The MTG", filename, "has been loaded!")
            else:
                print("!!! ERROR: the file", filename, "could not be found!")
                simulation_allowed=False
        # Otherwise we initiate the properties of the MTG "shoot_g":
        else:
            # We create the shoot MTG and store it in a dictionary under the corresponding genotype name:
            shoot_g = AllocF.create_shoot_MTG(
                number_of_cotyledons=2,  # We start with two cotyledones by default
                number_of_phytomers=0,
                initial_internode_radius=param.initial_internode_radius,
                initial_hypocotyl_length=0.001,
                initial_internode_height=0.,
                initial_cotyledon_length=param.initial_cotyledon_length,
                initial_leaf_length=param.initial_leaf_length,
                initial_labile_C_concentration=0.,
                initial_labile_N_concentration=0.,
                sowing_depth=sowing_depths_dict[plant_genotype_name]
            )
            shoot_MTG_dict[plant_genotype_name] = shoot_g

        # ROOT MTG - If we decide to start from a given root MTG to be loaded:
        if start_from_a_known_root_MTG and param.simulation_option != "shoots_only":
            filename = os.path.join(my_path, "simulation", "inputs", 'initial_root_MTG_' + plant_genotype_name + '.pckl')
            if os.path.exists(filename):
                # We load the MTG file and name it "g":
                f_root = open(filename, 'rb')
                root_MTG_dict[plant_genotype_name] = pickle.load(f_root)
                f_root.close()
                print("The MTG", filename, "has been loaded!")
            else:
                print("!!! ERROR: the file", filename, "could not be found!")
                simulation_allowed = False
        # Otherwise we initiate the properties of the MTG "root_g":
        else:
            # We create the root MTG and store it in a dictionary under the corresponding genotype name:
            root_g = RootSimple.initiate_mtg(initial_segment_length=0., initial_apex_length=0.001)
            root_MTG_dict[plant_genotype_name] = root_g

        # SETTING INITIAL CONDITIONS
        #---------------------------

        # INITIALIZING THE WHOLE PLANT:
        # We create a "plant", i.e. a shoot-root interface without structure that will store variables at the scale of the plant:
        plant = whole_plant()
        plant_objects[plant_genotype_name] = plant

        # INITIALIZING WHOLE SHOOT PROPERTIES:
        # In case we only want to simulate the behavior of the roots, keeping the shoot properties constant over time:
        if param.simulation_option == "roots_only":
            # We create a corresponding simple shoot MTG, without any dimensions (it will not be shown in the graph):
            shoot_g = AllocF.create_shoot_MTG(number_of_cotyledons=2,
                                                                    number_of_phytomers=0,
                                                                    initial_internode_height=0.,
                                                                    initial_internode_radius=0.,
                                                                    initial_leaf_length=0.,
                                                                    initial_labile_C_concentration=0.,
                                                                    initial_labile_N_concentration=0.,
                                                                    sowing_depth=sowing_depths_dict[plant_genotype_name]
                                                                    )
            shoot_MTG_dict[plant_genotype_name] = shoot_g
            # We define an arbitrary shoot structural mass and labile C and N pools,
            # which will remain constant over the simulation:
            plant.shoot_struct_mass = 0.01
            plant.shoot_labile_C_pool = 0.1
            plant.shoot_labile_N_pool = 0.01
        else:
            # Otherwise, the shoot MTG has already been defined/initiated, so we just have to update the properties
            # of the corresponding whole plant object:
            shoot_MTG_properties_directory = shoot_MTG_properties_directories[plant_genotype_name]
            shoot_MTG_files_directory = shoot_MTG_directories[plant_genotype_name]
            # We update the shoot system general properties based on the shoot MTG:
            for vid in shoot_g.vertices_iter(scale=1):
                n = shoot_g.node(vid)
                try:
                    plant.shoot_struct_mass += n.struct_mass
                    plant.shoot_labile_C_pool += n.labile_C
                    plant.shoot_labile_C_deficit += n.labile_C_deficit
                    plant.shoot_labile_N_pool += n.labile_N
                    plant.shoot_labile_N_deficit += n.labile_N_deficit
                except:
                    print("!!! ERROR: we could not compute the initial properties of the shoot element", vid, "in the shoot MTG", i)
            # We add an extra-property to the apex of the main stem, corresponding to the time since the last phytomer
            # was formed (used together with the phyllochron later on):
            base = shoot_g.node(0)
            base.thermal_time_since_last_phytomer_formation = 0.
            # We record the properties in the initial state of the shoot MTG in an external file:
            recording_MTG_properties(shoot_g, os.path.join(shoot_MTG_properties_directory, 'MTG_00000.csv'))
            # We record the initial shoot MTG file:
            g_file_name = os.path.join(shoot_MTG_files_directory, 'shoot%.5d.pckl')
            with open(g_file_name % 0, 'wb') as output:
                pickle.dump(shoot_g, output, protocol=2)

        # INITIALIZING WHOLE ROOT PROPERTIES:
        # In case we only want to simulate the behavior of the shoot, keeping the root properties constant over time:
        if param.simulation_option == "shoot_only":
            # We create a corresponding simple root MTG, without any dimensions (it will not be shown in the graph):
            root_g = root_MTG_dict[plant_genotype_name]
            # We define an arbitrary initial root mass and root labile C&N that will remain constant over the simulation:
            plant.root_struct_mass = 0.01
            plant.root_labile_C_pool = 0.01
            plant.root_labile_N_pool = 0.1
        else:
            # Otherwise, the root MTG has already been defined/initiated, so we just have to update the properties
            # of the corresponding whole plant object:
            root_MTG_properties_directory = root_MTG_properties_directories[plant_genotype_name]
            root_MTG_files_directory = root_MTG_directories[plant_genotype_name]
            # We update the distance from tip for each root element in each root axis:
            RootSimple.update_distance_from_tip(root_g)
            # We update the surfaces and volumes of each root element:
            RootSimple.update_surfaces_and_volumes(root_g)
            # We adjust the properties of the whole root system according to whole root MTG:
            root_g_properties = RootSimple.summing_and_possibly_homogenizing(root_g)
            # We record the total root structural mass within the root system:
            plant.root_struct_mass = root_g_properties["total_living_root_struct_mass"]
            # We calculate the total exchange surface of the root system:
            RootSimple.updating_exchange_surface(plant, root_g)
            # We record the properties in the initial state of the root MTG in an external file:
            recording_MTG_properties(root_g, os.path.join(root_MTG_properties_directory, 'MTG_00000.csv'))
            # We record the initial root MTG file:
            g_file_name = os.path.join(root_MTG_files_directory, 'root%.5d.pckl')
            with open(g_file_name % 0, 'wb') as output:
                pickle.dump(root_g, output, protocol=2)

        # RECORDING THE INITIAL PROPERTIES OF THE WHOLE PLANT:
        # We initialize the dictionary containing the general properties of the plant:
        results_dict = recording_main_plant_properties(plant, shoot_g, root_g,
                                                       final_time_in_days=starting_time_in_days,
                                                       thermal_time_in_degree_days=starting_thermal_time_in_degree_days)
        # We add a few more entries:
        results_dict["Air temperature (degree_Celsius)"] = "NA"
        results_dict["Soil temperature (degree_Celsius)"] = "NA"
        results_dict["Light irradiation (micromol of photons per m2 per second)"] = "NA"
        results_dict["Soil available N (gN per gDW of soil)"] = "NA"
        results_dict["Step"] = 0

        # We initialize a dataframe containing these properties for each time step:
        plant_data_frame_dict[plant_genotype_name] = pd.DataFrame(results_dict, index=[0])

        # And we initialize the C and N balance by computing the initial amounts of C and N in the shoots and roots:
        shoot_total_C_dict[plant_genotype_name] = (plant.shoot_struct_mass * param.structural_C_content
                                             + (plant.shoot_labile_C_pool - plant.shoot_labile_C_deficit)
                                             + plant.shoot_reserve_C_pool)
        root_total_C_dict[plant_genotype_name] = (plant.root_struct_mass * param.structural_C_content
                                            + (plant.root_labile_C_pool - plant.root_labile_C_deficit)
                                            + plant.root_reserve_C_pool)
        shoot_total_N_dict[plant_genotype_name] = (plant.shoot_struct_mass * param.structural_N_content
                                             + (plant.shoot_labile_N_pool - plant.shoot_labile_N_deficit)
                                             + plant.shoot_reserve_N_pool)
        root_total_N_dict[plant_genotype_name] = (plant.root_struct_mass * param.structural_N_content
                                            + (plant.root_labile_N_pool - plant.root_labile_N_deficit)
                                            + plant.root_reserve_N_pool)

    #///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    # DEFINING THE PERIOD OF SIMULATION:
    # ----------------------------------

    # We define the steps of the simulation based on the desired period:
    time_step_in_seconds = time_step_in_days * 60*60*24
    n_steps = trunc(simulation_period_in_days / time_step_in_days)
    # We initialize a counter that will be used over successive time steps together with the phyllochron:
    elapsed_time_in_days_since_last_leaf_emergence = 0.

    # We initiate actual and thermal time at the end of each step:
    actual_time_in_days = starting_time_in_days + time_step_in_days
    thermal_time_in_degree_days = starting_thermal_time_in_degree_days
    # We initialize the initial step to a default value:
    initial_step_number = 1

    # READING THE INPUT FILE:
    # -----------------------
    # If there is an input file to be read and at least one input is not forced by a constant:
    if input_file and not (forcing_constant_aboveground_environment and forcing_constant_belowground_environment):
        # We use the function 'formatted inputs' to create a table containing the input data for each time step:
        input_frame = formatted_inputs(original_input_file=os.path.join(my_path, "simulation", "inputs", input_file),
                                       final_input_file=os.path.join(outputs_directory, 'updated_input_file.csv'),
                                       original_time_step_in_days=input_file_time_step_in_days,
                                       final_time_step_in_days=time_step_in_days,
                                       simulation_period_in_days=simulation_period_in_days,
                                       do_not_execute_if_file_with_suitable_size_exists=False)
        # We then modify the initial step according to the match between the "step_number" column and the "time_in_days"
        # column, based on the starting_time_in_days prescribed in the parameters:
        initial_step_number \
            = input_frame[input_frame.initial_time_in_days >= starting_time_in_days]['step_number'].iloc[0]

    #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # FOR EACH DISTINCT PLANT GENOTYPE TO BE SIMULATED:
    for plant_genotype_name in genotype_names_list:
        # We update the parameters specific to this genotype:
        param.__dict__.update(parameters_dict[plant_genotype_name])
        # For saving all results at specified time intervals, we create a list with corresponding time steps:
        recording_steps_list_dict[plant_genotype_name] = list(
            range(initial_step_number - 1,
                  initial_step_number + n_steps,
                  round(param.recording_interval_in_days / time_step_in_days)))
    #///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    # FOR THE POPULATION AND EACH DISTINCT PLANT GENOTYPE TO BE SIMULATED:
    for name in population_and_genotypes_names_list:
        # We update the parameters from the original python file with the parameters read in the Excel file:
        param.__dict__.update(parameters_dict[name])

        # CONSIDERING CAMERA ROTATION WHEN PLOTTING:
        # -------------------------------------------

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
                    number_of_loops = n_steps
            # Otherwise, only one loop is necessary, and the camera will reuse the same coordinates after one complete rotation:
            else:
                number_of_loops = 1
            # We calculate the x,y,z coordinates of the different positions of the camera:
            [x,y,z] = circle_coordinates(x_center=param.x_center,
                                         y_center=param.y_center,
                                         z_center=param.z_camera,
                                         radius=sqrt((param.x_camera-param.x_center)**2 + (param.y_camera-param.y_center)**2),
                                         n_points=param.camera_rotation_n_points,
                                         n_loops=number_of_loops,
                                         step_back_coefficient=param.step_back_coefficient,
                                         step_up_coefficient=param.step_up_coefficient)
            cam_x_dict[name] = x
            cam_y_dict[name] = y
            cam_z_dict[name] = z
            # We initialize the index of the camera, so that it will get the first coordinate of x, y or z vectors:
            cam_index_dict[name] = 0
    #///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    # SIMULATING OVER TIME:
    #######################
    #######################

    #ttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt
    # FOR EACH TIME STEP COVERING THE SIMULATION PERIOD:
    for step in range(initial_step_number, initial_step_number + n_steps):

        # We calculate the initial and final time of the current time step (in days):
        initial_time_in_days = starting_time_in_days + (step - 1) * time_step_in_days
        final_time_in_days = starting_time_in_days + (step) * time_step_in_days

        # DEFINING THE INPUTS AT THIS TIME STEP:
        # --------------------------------------
        # If constant aboveground conditions have been chosen:
        if forcing_constant_aboveground_environment:
            # We set up the constant values:
            light_input = constant_light_irradiation
            air_temperature = constant_air_temperature_in_Celsius
            # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            # FOR EACH DISTINCT PLANT GENOTYPE TO BE SIMULATED:
            for plant_genotype_name in genotype_names_list:
                # We assume that germination of each genotype has already started:
                growth_stages_dict[plant_genotype_name] = 1
            # ///////////////////////////////////////////////////////////////////////////////////////////////////////////
        # Otherwise, if there is an input file to be used:
        elif input_file:
            # We extract the information on incident light and air temperature from the input file:
            light_input = input_frame.loc[step - initial_step_number, 'light_irradiation_in_micromol_per_m2_per_second']
            air_temperature = input_frame.loc[step - initial_step_number, 'air_temperature_in_degree_Celsius']

            #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            # FOR EACH DISTINCT PLANT GENOTYPE TO BE SIMULATED:
            # We get its "growth stage" from the input file, e.g. 0 = absent:
            for plant_genotype_name in genotype_names_list:
                # We reconstruct the name of the variable in the input file:
                variable_name = plant_genotype_name + "_stage"
                # We get the corresponding value and store it in a dictionary:
                growth_stages_dict[plant_genotype_name] = input_frame.loc[step - initial_step_number, variable_name]
            #///////////////////////////////////////////////////////////////////////////////////////////////////////////
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
        print("    Incoming light is", "{:.0f}".format(Decimal(light_input)), "micromol of photons per square meter per second.")
        print("    Soil available N is", "{:.2E}".format(Decimal(soil_available_nitrogen)), "gN per gram of dry soil.")

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
        thermal_time_step_in_seconds_roots = time_step_in_seconds * temperature_time_adjustment_soil
        # TODO in the future: consider a vertical gradient of temperature, and hence a vertical modification of the adjustments

        # We eventually compute the final and actual thermal time for this time step:
        actual_time_in_days = final_time_in_days
        thermal_time_in_degree_days += (thermal_time_step_in_seconds_shoots)/(24.*60.*60.) * param.T_ref_growth

        #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        # FOR EACH DISTINCT PLANT GENOTYPE TO BE SIMULATED:
        for plant_genotype_name in genotype_names_list:

            # VERIFYING THE PRESENCE OF THE PLANT:
            if growth_stages_dict[plant_genotype_name] == 0:
                # If the plant is absent for this genotype, we move to the next genotype:
                continue

            # ACCESSING THE CURRENT PLANT:
            plant = plant_objects[plant_genotype_name]

            # UPDATING PARAMETERS FOR THE CURRENT PLANT
            scenario = parameters_dict[plant_genotype_name]
            # We update the parameters from the original python file of SIMBAL with the parameters taken from this dictionnary:
            param.__dict__.update(scenario)
            # We update the parameters from the original python file of RhizoDep with the parameters taken from this dictionnary:
            rhizodep_param.__dict__.update(scenario)

            # RESETTING RATES AT PLANT SCALE:
            #################################

            # We reset the rates of exchange within the plant to 0:
            plant.net_C_transfer_from_shoot_to_root_rate = 0.  # gC per second
            plant.net_N_transfer_from_root_to_shoot_rate = 0.  # gN per second
            plant.transpiration_rate = 0.  # gram of water per second
            plant.total_PAR_intercepted = 0.  # micromoles of photons per plant per second
            plant.photosynthesized_C_rate = 0.  # gC per second
            plant.C_supply_from_seed_rate = 0. # gC per second
            plant.N_supply_from_seed_rate = 0. # gN per second
            plant.N_uptake_rate = 0.  # gN per second
            plant.N2_fixation_rate = 0.  # gN per second
            plant.shoot_maintenance_respiration_rate = 0.  # gC per second
            plant.root_maintenance_respiration_rate = 0.  # gC per second
            plant.actual_C_consumption_rate_by_shoot_growth = 0.  # gC per second
            plant.actual_N_consumption_rate_by_shoot_growth = 0.  # gN per second
            plant.actual_C_consumption_rate_by_root_growth = 0.  # gC per second
            plant.actual_N_consumption_rate_by_root_growth = 0.  # gN per second

            # TRANSFER OF C FROM SHOOTS TO ROOTS:
            #####################################

            print("\n    > Computing shoot-to-root C transfer for " + plant_genotype_name + "...")
            # We compute the rate of net C transfer from shoots to roots:
            computing_shoot_root_exchange_rates(plant, computing_C_exchange=True, computing_N_exchange=False,
                                                time_step_in_seconds=time_step_in_seconds,
                                                C_option="by_diffusion")
            print("    > The net transfer of C should be",
                  "{:.2E}".format(Decimal(plant.net_C_transfer_from_shoot_to_root_rate * time_step_in_seconds)), "gC.")
            # Note that only the rate has been computed here, the labile C pools have not been updated yet!
            # We limit the transfer of C from shoots to roots (or from roots to shoots) to the maximal amount of labile C
            # available in the source compartment.
            # We first calculate the expected amount of C that should be exchanged based on the flux:
            expected_C_exchange = plant.net_C_transfer_from_shoot_to_root_rate * time_step_in_seconds

            # CASE 1: The flux from shoot to root is positive, but there is not enough C in the shoot to sustain it:
            if expected_C_exchange > 0. and expected_C_exchange > (
                    plant.shoot_labile_C_pool - plant.shoot_labile_C_deficit):
                # Then the amount of C exchanged is limited to the amount of C available in shoots:
                final_C_exchange = plant.shoot_labile_C_pool
                #  The amount of labile C in shoots and roots is updated according to the new computed exchange:
                plant.shoot_labile_C_pool -= final_C_exchange
                plant.root_labile_C_pool += final_C_exchange
            # CASE 2: The flux from root to shoot is positive but there is not enough C in the root to sustain it:
            elif expected_C_exchange < 0. and abs(expected_C_exchange) > (
                    plant.root_labile_C_pool - plant.root_labile_C_deficit):
                final_C_exchange = plant.root_labile_C_pool
                # The amount of labile C in shoots and roots is updated according to the new computed exchange:
                plant.shoot_labile_C_pool += final_C_exchange
                plant.root_labile_C_pool -= final_C_exchange
            # CASE 3 (EXPECTED CASE): There is enough C in the source compartment to sustain the transfer:
            else:
                final_C_exchange = expected_C_exchange
                # Then the amount of labile C in shoots and roots is updated according to the initial computed exchange:
                plant.shoot_labile_C_pool -= expected_C_exchange
                plant.root_labile_C_pool += expected_C_exchange
            print("    ... The net transfer of C was eventually", "{:.2E}".format(Decimal(final_C_exchange)), "gC.")

            # TRANSFER OF N FROM ROOTS TO SHOOTS:
            #####################################

            print("    > Computing root-to-shoot N transfer...")
            # We compute the rate of net N transfer from roots to shoots:
            computing_shoot_root_exchange_rates(plant, computing_C_exchange=False, computing_N_exchange=True,
                                                time_step_in_seconds=time_step_in_seconds,
                                                N_option="by_diffusion")
            print("    > The net transfer of N should be",
                  "{:.2E}".format(Decimal(plant.net_N_transfer_from_root_to_shoot_rate * time_step_in_seconds)),
                  "gN.")
            # We then limit the transfer of N from roots to shoots (or from shoots to roots) to the maximal amount of N
            # available in the source compartment.
            # We first calculate the expected amount of N that should be exchanged based on the exchange rate that has been computed:
            expected_N_exchange = plant.net_N_transfer_from_root_to_shoot_rate * time_step_in_seconds

            # CASE 1: The flow of N from root to shoot is positive, but there is not enough N in the root to sustain it.
            if expected_N_exchange > 0. and expected_N_exchange > (
                    plant.root_labile_N_pool - plant.root_labile_N_deficit):
                # Then the amount of N exchanged is limited to the amount of N available in roots:
                final_N_exchange = plant.root_labile_N_pool
                # The amount of labile N in shoots and roots is updated according to the new computed exchange:
                plant.shoot_labile_N_pool += final_N_exchange
                plant.root_labile_N_pool -= final_N_exchange
            # CASE 2: The flux from shoot to root is positive but there is not enough N in the shoot to sustain it:
            elif expected_N_exchange < 0. and abs(expected_N_exchange) > (
                    plant.shoot_labile_N_pool - plant.shoot_labile_N_deficit):
                # Then the amount of N exchanged is limited to the amount of N available in shoots:
                final_N_exchange = plant.shoot_labile_N_pool
                # The amount of labile N in shoots and roots is updated according to the new computed exchange:
                plant.shoot_labile_N_pool -= final_N_exchange
                plant.root_labile_N_pool += final_N_exchange
            # CASE 3 (EXPECTED CASE): There is enough N in the source compartment to sustain the transfer.
            else:
                # Then the amount of labile N in shoots and roots is updated according to the computed exchange:
                final_N_exchange = expected_N_exchange
                plant.shoot_labile_N_pool += expected_N_exchange
                plant.root_labile_N_pool -= expected_N_exchange
            print("    ... The net transfer of N was eventually", "{:.2E}".format(Decimal(final_N_exchange)), "gC.")

            # CALLING SUCCESSIVELY THE DIFFERENT ROOT MODULES:
            ##################################################

            if param.simulation_option == "shoot_only":

                # We set a constant rate of N uptake by the root system:
                plant.N_uptake_rate = 1e-5
                plant.root_labile_N_pool += plant.N_uptake_rate * time_step_in_seconds

            else:

                # We access the root MTG:
                root_g = root_MTG_dict[plant_genotype_name]

                print("\n    For", plant_genotype_name, "root labile C pool is initially",
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

                # NET C & N LOSSES FROM THE ROOT SYSTEM (WITHOUT GROWTH):
                # We compute the rate of rhizodeposition (gC or gN per plant per second):
                RootSimple.computing_rhizodeposition_rate(root_system=plant, root_MTG=root_g,
                                                          option="loss_as_input_rate_fraction",
                                                          C_rhizodeposition_fraction=0.30,
                                                          N_rhizodeposition_fraction=0.10,
                                                          time_step_in_seconds=time_step_in_seconds,
                                                          soil_temperature_in_degree_Celsius=soil_temperature)
                # We compute the rate of maintenance respiration (gC per plant per second):
                RootSimple.computing_maintenance_rate(root_system=plant,
                                                      soil_temperature_in_degree_Celsius=soil_temperature)
                #TODO: WATCH OUT BELOW
                # # We compute the rate of C exchange with reserve (gC per plant per second):
                # RootSimple.computing_exchange_rates_with_reserve(root_system=plant,
                #                                                  soil_temperature_in_degree_Celsius=soil_temperature)

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
                                                                                       labile_C_pool=plant.root_labile_C_pool - plant.root_labile_C_deficit,
                                                                                       labile_N_pool=plant.root_labile_N_pool - plant.root_labile_N_deficit)
                print("    The satisfaction coefficient for C is", "{:.2f}".format(Decimal(satisf_for_C)),
                      "and the one for N is", "{:.2f}".format(Decimal(satisf_for_N)))
                # We use the function ArchiSimple_growth to adapt the potential growth to the available labile C and N,
                # update the new root structure, and calculate the consumption rates for C & N,
                # without updating the labile pools of C and N:
                RootSimple.ArchiSimple_growth_with_C_and_N(plant, root_g,
                                                           satisf_for_C, satisf_for_N,
                                                           time_step_in_seconds=time_step_in_seconds)
                # We proceed to the segmentation of the whole root system:
                RootSimple.segmentation_and_primordia_formation(root_g, time_step_in_seconds,
                                                                soil_temperature_in_Celsius=soil_temperature, random=True,
                                                                root_order_limitation=param.root_order_limitation,
                                                                root_order_treshold=param.root_order_treshold)
                # TODO: Include root hairs dynamics once the cost of C & N is properly registered!
                # # We modifiy root hairs characteristics according to their specific dynamics:
                # root_hairs_dynamics(g, time_step_in_seconds=time_step_in_seconds,
                #                     soil_temperature_in_Celsius=soil_temperature)
                # We update the distance from tip for each root element in each root axis:
                RootSimple.update_distance_from_tip(root_g)
                # We update the surfaces and volumes of each root element:
                RootSimple.update_surfaces_and_volumes(root_g)
                RootSimple.updating_exchange_surface(plant, root_g)

                print("    Before C&N root balance, root labile C pool is",
                      "{:.2E}".format(Decimal(plant.root_labile_C_pool - plant.root_labile_C_deficit)),
                      "gC and root labile N pool is",
                      "{:.2E}".format(Decimal(plant.root_labile_N_pool - plant.root_labile_N_deficit)), "gN.")

                # FINAL ROOT C&N BALANCE:
                # We update the pools of root labile C and N at the scale of the whole root system:
                RootSimple.root_CN_balance(plant, time_step_in_seconds=time_step_in_seconds)

                print("    Root labile C pool is now",
                      "{:.2E}".format(Decimal(plant.root_labile_C_pool - plant.root_labile_C_deficit)),
                      "gC and root labile N pool is",
                      "{:.2E}".format(Decimal(plant.root_labile_N_pool - plant.root_labile_N_deficit)), "gN.")
                print("    The root system now contains", root_g.nb_vertices()-1, "different root elements.")

            # CALLING SUCCESSIVELY THE DIFFERENT SHOOT MODULES:
            ###################################################

            # REINITIALIZING SHOOT VARIABLES:
            # If no shoot system is simulated we simply set constant labile C and photosynthesis rate:
            if param.simulation_option == "roots_only":
                # We set a constant rate of C uptake by the shoot system:
                plant.photosynthesized_C_rate = 1e-4
                plant.shoot_labile_C_pool += plant.photosynthesized_C_rate * time_step_in_seconds
            else:
                # We access the shoot MTG:
                shoot_g = shoot_MTG_dict[plant_genotype_name]
                # We reset all local properties to 0:
                AllocF.shoot_MTG_reinitializing(shoot_g, reinitializing_rates=True, reinitializing_amounts=True)
        #///////////////////////////////////////////////////////////////////////////////////////////////////////////////

        # LIGHT INTERCEPTION:
        #********************

        # If light interception is necessary for the population and/or at least one genotype, we compute it:
        if computing_light_interception:
            print("\n    Preparing the scene for light interception...")
            # We create a dictionary of scene details and an updated list containing one shoot MTG of each genotype
            # copied and distributed over the scene:
            scene_information, list_of_all_shoot_MTG, list_of_all_root_MTG = prepare_heterogeneous_plants_population(
                genotype_names_list=genotype_names_list,
                parameters_dict=parameters_dict,
                shoot_MTG_dict=shoot_MTG_dict,
                root_MTG_dict=root_MTG_dict,
                option="shoots_only",
                numbers_of_plants_per_square_meter_dict=numbers_of_plants_per_square_meter_dict,
                sowing_depths_dict=sowing_depths_dict,
                inter_row_distances_dict=inter_row_distances_dict,
                number_of_rows_to_show_dict=number_of_rows_to_show_dict,
                scene_lengths_dict=scene_lengths_dict,
                scene_widths_dict=scene_widths_dict,
                upscaling_factor=zoom_factor)
            # We create an updated scene, with reindexed shapes compared to the original indices in the list of MTG:
            scene, scene_information = create_plants_population_scene(scene_information=scene_information,
                                                                      genotype_names_list=genotype_names_list,
                                                                      parameters_dict=parameters_dict,
                                                                      list_of_shoot_MTG=list_of_all_shoot_MTG,
                                                                      list_of_root_MTG=[],
                                                                      option="shoots_only",
                                                                      leaf_shape=param.light_interception_leaf_shape,
                                                                      simplified_leaf_shape_for_surface_calculation=param.simplified_leaf_shape_for_surface_calculation,
                                                                      upscaling_factor=zoom_factor)
            print("    Intercepting light...")
            # We compute the light interception on the whole scene, and collect the intercepted PAR for each distinct plant type:
            total_PAR_intercepted_dict, sum_of_surfacic_PAR_dict = \
                AllocF.computing_light_interception_on_plant_population(scene=scene,
                                                                        scene_information=scene_information,
                                                                        genotype_names_list=genotype_names_list,
                                                                        shoot_MTG_dict=shoot_MTG_dict,
                                                                        list_of_scene_MTG=list_of_all_shoot_MTG,
                                                                        incoming_light_density=light_input,
                                                                        zenithal_light=False,
                                                                        infinite=False,
                                                                        averaging_light_intensity_share=True,
                                                                        period_of_light_intensity_averaging_in_days=1.,
                                                                        time_step_in_days=time_step_in_days
                                                                        )

        # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        # FOR EACH DISTINCT PLANT GENOTYPE TO BE SIMULATED:
        for plant_genotype_name in genotype_names_list:

            # VERIFYING THE PRESENCE OF THE PLANT:
            if growth_stages_dict[plant_genotype_name] == 0:
                # If the plant is absent for this genotype, we move to the next genotype:
                continue

            # ACCESSING THE CURRENT PLANT:
            plant = plant_objects[plant_genotype_name]

            # UPDATING PARAMETERS FOR THE CURRENT PLANT:
            # We get the set of parameters for this genotype:
            scenario = parameters_dict[plant_genotype_name]
            # We update the parameters from the original python file of SIMBAL with this set of parameters:
            param.__dict__.update(scenario)

            # If the shoot system is to be computed:
            if param.simulation_option != "roots_only":

                # We access the shoot MTG:
                shoot_g = shoot_MTG_dict[plant_genotype_name]

                # We register the initial state of the labile C & N pools in the shoot:
                initial_shoot_labile_C = plant.shoot_labile_C_pool - plant.shoot_labile_C_deficit
                initial_shoot_labile_N = plant.shoot_labile_N_pool - plant.shoot_labile_N_deficit

                print("\n    The shoot system of", plant_genotype_name, "initially contains", shoot_g.nb_vertices() - 1,
                      "elements.")
                print("    Its shoot labile C pool is initially",
                      "{:.2E}".format(Decimal(initial_shoot_labile_C)),
                      "gC and its shoot labile N pool is",
                      "{:.2E}".format(Decimal(initial_shoot_labile_N)), "gN.")

                total_PAR_intercepted = total_PAR_intercepted_dict[plant_genotype_name]
                sum_of_surfacic_PAR = sum_of_surfacic_PAR_dict[plant_genotype_name]
                # We record the PAR intercepted as a property of the whole plant:
                plant.total_PAR_intercepted = total_PAR_intercepted

                # DISTRIBUTING C & N IN EACH SHOOT ORGAN:
                #****************************************

                # Redistributing labile C into each shoot organ:
                #-----------------------------------------------
                # The total amount of labile C to redistribute into the shoot organs is defined as:
                shoot_available_C = plant.shoot_labile_C_pool - plant.shoot_labile_C_deficit
                # SPECIAL CASE: if the tip of hypocothyl is belowground, it can grow and receives all available labile C
                hypocotyl = shoot_g.node(1)
                if hypocotyl.z2 < 0.:
                    # We calculate the rate of net N input to this organ corresponding to N distribution:
                    hypocotyl.shoot_C_distribution_rate = (shoot_available_C - (
                                hypocotyl.labile_C - hypocotyl.labile_C_deficit)) / time_step_in_seconds
                    hypocotyl.labile_C = shoot_available_C
                    hypocotyl.labile_C_deficit = 0.
                else:
                    # We cover each organ n of the shoot:
                    for vid in shoot_g.vertices_iter(scale=1):
                        # We access the current organ "n":
                        n = shoot_g.node(vid)
                        if n.surface <= 0. or n.length <= 0.:
                            continue
                        # If there is no labile C to be distributed, the labile C in the shoot element will remain 0:
                        if shoot_available_C <= 0.:
                            labile_C_target = 0.
                            # We calculate the rate of net C input to this organ corresponding to this redistribution:
                            n.shoot_C_distribution_rate = (labile_C_target - (
                                    n.labile_C - n.labile_C_deficit)) / time_step_in_seconds
                            # We now redefine the local pool of labile C according to this target level:
                            n.labile_C_deficit = -labile_C_target
                            n.labile_C = 0.
                        # Otherwise, if there is labile C in the shoot to be redistributed:
                        else:
                            #TODO: WATCH OUT BELOW
                            # We define the fraction of shoot available C allocated to this organ based
                            # on the previously computed light intensity share.
                            # If no light interception was computed so far:
                            if len(n.previous_light_intensity_shares) == 0:
                                # Then the fraction of labile C allocated to this organ will be based on its structural mass:
                                fraction = n.struct_mass / plant.shoot_struct_mass
                            else:
                                fraction = n.light_intensity_share
                            # The eventual shoot labile C to be targeted after redistribution is:
                            labile_C_target = shoot_available_C * fraction
                            # We calculate the rate of net N input to this organ corresponding from this redistribution:
                            n.shoot_C_distribution_rate = (labile_C_target - (n.labile_C - n.labile_C_deficit)) / time_step_in_seconds
                            # We now redefine the local pool of labile C according to this target level:
                            n.labile_C = labile_C_target
                            n.labile_C_deficit = 0.

                # Redistributing labile N within each shoot organ:
                #-------------------------------------------------
                # The total amount of labile N to redistribute into the shoot organs is defined as:
                shoot_available_N = plant.shoot_labile_N_pool - plant.shoot_labile_N_deficit
                # SPECIAL CASE: if the hypocothyl can grow, then it receives all available labile N
                hypocotyl = shoot_g.node(1)
                if hypocotyl.z2 < 0.:
                    # We calculate the rate of net N input to this organ corresponding to N distribution:
                    hypocotyl.shoot_N_distribution_rate = (shoot_available_N - (
                                hypocotyl.labile_N - hypocotyl.labile_N_deficit)) / time_step_in_seconds
                    hypocotyl.labile_N = shoot_available_N
                    hypocotyl.labile_N_deficit = 0.
                else:
                    # We cover each organ of the shoot:
                    for vid in shoot_g.vertices_iter(scale=1):
                        # We access the current organ "n":
                        n = shoot_g.node(vid)
                        if n.surface <= 0.:
                            continue
                        # If there is no labile N to be distributed, the shoot labile N will remain 0:
                        if shoot_available_N <= 0.:
                            labile_N_target = 0.
                            # We calculate the rate of net N input to this organ corresponding from this redistribution:
                            n.shoot_N_distribution_rate = (labile_N_target - (
                                    n.labile_N - n.labile_N_deficit)) / time_step_in_seconds
                            # We now redefine the local pool of labile N according to this target level:
                            n.labile_N = 0.
                            n.labile_N_deficit = -labile_N_target
                        # Otherwise, if there is labile N in the shoot:
                        else:
                            #TODO: WATCH OUT BELOW
                            # We define the fraction of shoot available N allocated to this organ based on the previously
                            # computed light intensity share.
                            # If no light interception was computed so far:
                            if len(n.previous_light_intensity_shares) == 0:
                                # Then the fraction of labile N allocated to this organ will be based on its structural mass:
                                fraction =  n.struct_mass / plant.shoot_struct_mass
                            else:
                                fraction = n.light_intensity_share
                            # The eventual shoot labile N to be targeted after redistribution is:
                            labile_N_target = shoot_available_N * fraction
                            # We calculate the rate of net N input to this organ corresponding from this redistribution:
                            n.shoot_N_distribution_rate = (labile_N_target - (n.labile_N - n.labile_N_deficit)) / time_step_in_seconds
                            # We now redefine the local pool of labile N according to this target level:
                            # n.labile_N = 0.
                            # n.labile_N_deficit = -labile_N_target
                            n.labile_N = labile_N_target
                            n.labile_N_deficit = 0.

                # POSSIBLE SUPPLY OF RESOURCES FROM THE GERMINATING SEED:
                # If the elapsed thermal time is lower than the prescribed duration of supply from the seed:
                if thermal_time_in_degree_days * (60.*60.*24.) <= param.seed_supply_duration * param.T_ref_growth:
                    # Then we allow to compute the supply of C and N resources from the germinating seed to the new seedling:
                    AllocF.computing_supply_from_seed_rate(plant, air_temperature_in_degree_Celsius = air_temperature)
                    print("    The seed will generate",
                          "{:.2E}".format(Decimal(plant.C_supply_from_seed_rate*time_step_in_seconds)), "gC and",
                          "{:.2E}".format(Decimal(plant.N_supply_from_seed_rate*time_step_in_seconds)),"gN over this time step...")
                    # Note: at this stage, the labile C and N amounts have not been updated in the organs or at the scale of the plant!

                # PHOTOSYNTHESIS:
                # We calculate the rate of C production of each organ based on the PAR intercepted by the plant:
                total_PAR_intercepted, total_photosynthesis_rate = \
                    AllocF.computing_photosynthesis_rate(plant, shoot_g,
                                                         air_temperature_in_degree_Celsius=air_temperature)
                # We record this total photosynthesis rate into the plant object:
                plant.photosynthesized_C_rate = total_photosynthesis_rate
                print("    Photosynthesis will generate",
                      "{:.2E}".format(Decimal(total_photosynthesis_rate*time_step_in_seconds)),
                      "gC over this time step...")
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
                # EXTENSION: We calculate the possible increase in volume for each organ:
                AllocF.potential_shoot_growth(shoot_g,
                                              SFQN=param.SFQN, # TODO: is SFQN constant over the whole range of (0, LA_max)? Is it two distinct curves?
                                              leaf_minimal_thickness=param.leaf_minimal_thickness,
                                              leaf_shape=param.simplified_leaf_shape_for_surface_calculation,
                                              elong_hypocotyl=param.elong_internode,
                                              elong_internode_slow=param.elong_internode/10.,
                                              elong_internode_fast=param.elong_internode,
                                              time_step_in_seconds=time_step_in_seconds,
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
                                               leaf_shape=param.simplified_leaf_shape_for_surface_calculation,
                                               time_step_in_seconds=time_step_in_seconds,
                                               air_temperature_in_degree_Celsius=air_temperature)
                # At this stage, the shoot structure has been updated and the consumption rates have been calculated
                # for each organ, but no actual consumption of labile C&N at the scale of organ or shoot system has occured yet.
                print("    The extension of organs will consume",
                      "{:.2E}".format(Decimal(C_consumption_rate_for_extension * time_step_in_seconds)),
                      "gC and",
                      "{:.2E}".format(Decimal(N_consumption_rate_for_extension * time_step_in_seconds)),
                      "gN over the shoot system...")

                # NEW PHYTOMERS: We now consider the initiation of the structure of new phytomers, based on the phyllochron.
                # We increment the time since the last phytomer was formed at the top of the main stem:
                thermal_time_since_last_phytomer_formation \
                    = shoot_g.node(0).thermal_time_since_last_phytomer_formation + thermal_time_step_in_seconds_shoots
                # We get the value of the prescribed plastochron and phyllochron (in seconds equivalent to the reference temperature):
                plastochron = param.plastochron
                phyllochron = param.phyllochron
                # We initialize the total costs of phytomers initiation:
                total_C_costs_for_new_phytomers = 0.
                total_N_costs_for_new_phytomers = 0.
                # If the period since the last phytomer was formed is higher than the phyllochron:
                if (thermal_time_since_last_phytomer_formation > phyllochron
                        and shoot_g.node(1).thermal_time_since_aboveground_emergence > 0.):
                    print("A new phytomer should be formed!")
                    # Then at least one phytomer should be formed on the main stem, and additional one on ramifications,
                    # if any. We keep looping until all phytomers have been formed:
                    while thermal_time_since_last_phytomer_formation - phyllochron > 0:
                        # We calculate the exact thermal time elapsed since the phytomer has actually emerged:
                        thermal_time_since_emergence \
                            = (thermal_time_since_last_phytomer_formation
                               - np.floor(thermal_time_since_last_phytomer_formation / phyllochron) * phyllochron)
                        # We compute then the exact thermal time elapsed since the phytomer was initiated in the apex,
                        # which is linearily correlated to the time since emergence by the following relationship:
                        # t_initiation = plastochron + t_emergence * plastochron  / phyllochron
                        thermal_time_since_initiation = plastochron + thermal_time_since_emergence * plastochron / phyllochron
                        # The actual times are then calculated:
                        actual_time_since_emergence = thermal_time_since_emergence / temperature_time_adjustment_shoot
                        # TODO: WATCH OUT BELOW - may not be correct!
                        actual_time_since_initiation = thermal_time_since_initiation / temperature_time_adjustment_shoot
                        # We add a set of new phytomers in the shoot system:
                        C_cost_for_new_phytomers, N_cost_for_new_phytomers = AllocF.new_phytomers_emergence(
                            shoot_system=plant,
                            shoot_MTG=shoot_g,
                            thermal_time_since_initiation=thermal_time_since_initiation,
                            actual_time_since_initiation=actual_time_since_initiation,
                            thermal_time_since_emergence=thermal_time_since_emergence,
                            actual_time_since_emergence=actual_time_since_emergence,
                            initial_leaf_length=param.initial_leaf_length,
                            initial_internode_length=param.initial_internode_length,
                            time_step_in_seconds=time_step_in_seconds)
                        print("> New phytomers have been formed! The C & N costs were",
                              "{:.2E}".format(Decimal(C_cost_for_new_phytomers)), "gC and",
                              "{:.2E}".format(Decimal(N_cost_for_new_phytomers)), "gN.")
                        # We increase the total costs with the cost of this new phytomer:
                        total_C_costs_for_new_phytomers += C_cost_for_new_phytomers
                        total_N_costs_for_new_phytomers += total_N_costs_for_new_phytomers
                        # We now subtract the value of the phyllochron to the time since last phytomer was formed, and loop:
                        thermal_time_since_last_phytomer_formation -= phyllochron
                # Once all phytomers have been added, we update the time since the last one was formed,
                # which is stored at the base of the shoot system:
                shoot_g.node(0).thermal_time_since_last_phytomer_formation = thermal_time_since_last_phytomer_formation

                # RAMIFICATION - We consider the possibility of ramification:
                C_cost_for_new_ramification, N_cost_for_new_ramification \
                    = AllocF.computing_ramification(shoot_MTG=shoot_g,
                                                    initial_internode_length=param.initial_internode_length,
                                                    initial_leaf_length=param.initial_leaf_length,
                                                    time_step_in_seconds=time_step_in_seconds)
                if C_cost_for_new_ramification > 0. or N_cost_for_new_ramification > 0.:
                    print("> The costs to sustain ramification with the creation of new phytomers were",
                          "{:.2E}".format(Decimal(C_cost_for_new_ramification)), "gC and",
                          "{:.2E}".format(Decimal(N_cost_for_new_ramification)), "gN.")

                # We finally update the coordinates of the new organs:
                updating_coordinates(shoot_g)
                # And we register the total consumption rate for growth in the plant object:
                plant.actual_C_consumption_rate_by_shoot_growth \
                    = (C_consumption_rate_for_extension
                       + (C_cost_for_new_ramification + total_C_costs_for_new_phytomers) / time_step_in_seconds)
                plant.actual_N_consumption_rate_by_shoot_growth \
                    = (N_consumption_rate_for_extension
                       + (N_cost_for_new_ramification + total_N_costs_for_new_phytomers) / time_step_in_seconds)


                # C&N BALANCE both within each organ and at the shoot system scale:
                AllocF.shoot_CN_balance_at_organ_and_system_scales(shoot_system=plant,
                                                                   shoot_MTG=shoot_g,
                                                                   time_step_in_seconds=time_step_in_seconds)
                print("    The shoot system now contains", shoot_g.nb_vertices()-1, "different shoot elements.")
                print("    Shoot labile C pool is now",
                      "{:.2E}".format(Decimal(plant.shoot_labile_C_pool - plant.shoot_labile_C_deficit)),
                      "gC and labile N pool is",
                      "{:.2E}".format(Decimal(plant.shoot_labile_N_pool - plant.shoot_labile_N_deficit)), "gN.")

                # AGEING AND CONSEQUENCES:
                # We increase the age of each organ:
                AllocF.ageing(shoot_MTG=shoot_g,
                              time_step_in_seconds=time_step_in_seconds,
                              air_temperature_in_degree_Celsius=air_temperature,
                              zoom_factor=zoom_factor)
                # We update the geometry with plant / organ age, e.g. insertion angles:
                AllocF.geometry_evolution(shoot_MTG = shoot_g,
                                          plant_genotype_name = plant_genotype_name,
                                          plant_thermal_age_in_degree_days = thermal_time_in_degree_days)

            # CHECKING WHOLE PLANT C & N BALANCE:
            #####################################

            print("    C and N balance check on", plant_genotype_name + "...")

            # We compute the current amounts of total C and N in the shoots and roots, and compare it with the previous ones
            # when altered by the new inputs and outputs computed at this time step:
            (shoot_total_C_dict[plant_genotype_name], root_total_C_dict[plant_genotype_name],
             shoot_total_N_dict[plant_genotype_name], root_total_N_dict[plant_genotype_name]) \
                = balance_check(plant,
                                genotype_name=plant_genotype_name,
                                time_step_in_seconds = time_step_in_seconds,
                                previous_shoot_C = shoot_total_C_dict[plant_genotype_name],
                                previous_root_C = root_total_C_dict[plant_genotype_name],
                                previous_shoot_N = shoot_total_N_dict[plant_genotype_name],
                                previous_root_N = root_total_N_dict[plant_genotype_name])

            # RECORDING THE CURRENT STATE:
            ##############################

            # If the current iteration correspond to the time where one full time interval for recording has been reached:
            if step in recording_steps_list_dict[plant_genotype_name]:

                print("    Recording the results for", plant_genotype_name + "...")

                # ACCESSING THE CURRENT PLANT AND DIRECTORIES:
                shoot_g = shoot_MTG_dict[plant_genotype_name]
                root_g = root_MTG_dict[plant_genotype_name]
                sub_outputs_directory = outputs_directories[plant_genotype_name]
                shoot_MTG_properties_directory = shoot_MTG_properties_directories[plant_genotype_name]
                root_MTG_properties_directory = root_MTG_properties_directories[plant_genotype_name]
                shoot_MTG_files_directory = shoot_MTG_directories[plant_genotype_name]
                root_MTG_files_directory = root_MTG_directories[plant_genotype_name]

                # RECORDING THE PROPERTIES OF THE SHOOT MTG FOR THIS TIME STEP:
                if param.simulation_option != "roots_only":
                    # We can record the main properties of the plant:
                    recording_MTG_properties(shoot_g, file_name=os.path.join(shoot_MTG_properties_directory, 'MTG_%.5d.csv' % step))

                # RECORDING THE PROPERTIES OF THE ROOT MTG FOR THIS TIME STEP:
                if param.simulation_option != "shoot_only":
                    # We can record the main properties of the plant:
                    recording_MTG_properties(root_g, file_name=os.path.join(root_MTG_properties_directory, 'MTG_%.5d.csv' % step))

                # RECORDING THE SHOOT MTG FILE FOR THIS TIME STEP:
                if param.simulation_option != "roots_only":
                    g_file_name = os.path.join(shoot_MTG_files_directory, 'shoot%.5d.pckl')
                    with open(g_file_name % step, 'wb') as output:
                        pickle.dump(shoot_g, output, protocol=2)

                # RECORDING THE ROOT MTG FILE FOR THIS TIME STEP:
                if param.simulation_option != "shoots_only":
                    g_file_name = os.path.join(root_MTG_files_directory, 'root%.5d.pckl')
                    with open(g_file_name % step, 'wb') as output:
                        pickle.dump(root_g, output, protocol=2)

            # RECORDING THE FINAL RESULTS AT PLANT SCALE:
            results_dict = recording_main_plant_properties(plant, shoot_g, root_g,
                                                           final_time_in_days=actual_time_in_days,
                                                           thermal_time_in_degree_days=thermal_time_in_degree_days)
            results_dict["Air temperature (degree_Celsius)"] = air_temperature
            results_dict["Soil temperature (degree_Celsius)"] = soil_temperature
            results_dict["Light irradiation (micromol of photons per m2 per second)"] = light_input
            results_dict["Soil available N (gN per gDW of soil)"] = soil_available_nitrogen
            results_dict["Step"] = step
            new_df = pd.DataFrame([results_dict.values()], columns=results_dict.keys())

            # We get the dataframe where previous results have been recorded:
            plant_data_frame = plant_data_frame_dict[plant_genotype_name]
            # We add the new results:
            plant_data_frame = pd.concat([plant_data_frame, new_df])
            plant_data_frame_dict[plant_genotype_name] = plant_data_frame
            # We record the dataframe as a csv file:
            try:
                plant_data_frame.to_csv(os.path.join(sub_outputs_directory, 'plant_results.csv'), na_rep='NA',
                                        index=False, header=True)
            except:
                plant_data_frame.to_csv(os.path.join(sub_outputs_directory, 'plant_results_BACKUP.csv'),
                                        na_rep='NA', index=False, header=True)

        #///////////////////////////////////////////////////////////////////////////////////////////////////////////////

        # PLOTTING:
        ###########

        print("\n    Plotting...")

        if plotting_single_plant or plotting_population:
            print("")

            # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
            # FOR THE POPULATION AND FOR EACH DISTINCT PLANT GENOTYPE TO BE SIMULATED:
            for name in population_and_genotypes_names_list:

                # ACCESSING THE PLOT DIRECTORY AND UPDATING PARAMETERS:
                plot_directory = plot_directories[name]
                scenario = parameters_dict[name]
                # We update the parameters from the original python file of SIMBAL with the parameters taken from this dictionnary:
                param.__dict__.update(scenario)

                # If the camera is to be moved around the scene, we modify the position of the camera:
                if param.camera_movement is True:
                    x_center = param.x_center
                    y_center = param.y_center
                    z_center = param.z_center
                    index_camera = cam_index_dict[name]
                    x_cam = cam_x_dict[name][index_camera]
                    y_cam = cam_y_dict[name][index_camera]
                    z_cam = cam_z_dict[name][index_camera]
                    # We define the index of the coordinates to read at the next step:
                    cam_index_dict[name] = index_camera + 1
                    # If this index is higher than the number of coordinates in each vector and if no other movement
                    # of the camera (away or down) has been set over time:
                    if (cam_index_dict[name] >= param.camera_rotation_n_points
                            and param.step_back_coefficient == 0. and param.step_up_coefficient == 0.):
                        # Then we reset the index to 0, so that the camera goes back to the same locations as previously stored:
                        cam_index_dict[name] = 0
                    # If "distancing" is required, the camera will continue to read the further coordinates already calculated.
                # Otherwise, we simply use the parameters specified for the camera:
                else:
                    x_center = param.x_center
                    y_center = param.y_center
                    z_center = param.z_center
                    x_cam = param.x_camera
                    y_cam = param.y_camera
                    z_cam = param.z_camera

                # We consider the possibility of single colors to show for background, shoots and roots:
                if param.background_color is None or pd.isna(param.background_color):
                    plot_background_color = None
                else:
                    plot_background_color = ast.literal_eval(param.background_color)
                if param.single_color_shoot is None or pd.isna(param.single_color_shoot):
                    plot_single_color_shoot = None
                else:
                    plot_single_color_shoot = ast.literal_eval(param.single_color_shoot)
                if param.single_color_root is None or pd.isna(param.single_color_root):
                    plot_single_color_root = None
                else:
                    plot_single_color_root = ast.literal_eval(param.single_color_root)

                # PLOTTING THE POPULATION:
                if name == "population" and plotting_population:
                    plotting_whole_scene(genotype_names_list=genotype_names_list,
                                         parameters_dict=parameters_dict,
                                         shoot_MTG_dict=shoot_MTG_dict,
                                         root_MTG_dict=root_MTG_dict,
                                         numbers_of_plants_per_square_meter_dict=numbers_of_plants_per_square_meter_dict,
                                         inter_row_distances_dict=inter_row_distances_dict,
                                         number_of_rows_to_show_dict=number_of_rows_to_show_dict,
                                         scene_lengths_dict=scene_lengths_dict,
                                         scene_widths_dict=scene_widths_dict,
                                         sowing_depths_dict=sowing_depths_dict,
                                         width=1200, height=1200,
                                         upscaling_factor=zoom_factor,
                                         x_center=x_center * zoom_factor,
                                         y_center=y_center * zoom_factor,
                                         z_center=z_center * zoom_factor,
                                         x_cam=x_cam * zoom_factor,
                                         y_cam=y_cam * zoom_factor,
                                         z_cam=z_cam * zoom_factor,
                                         grid=False,
                                         background_color=plot_background_color,
                                         single_color_shoot=plot_single_color_shoot,
                                         single_color_root=plot_single_color_root,
                                         property_name_shoot=param.displayed_property_shoot,
                                         cmap_shoot=param.color_map_shoot,
                                         vmin_shoot=param.displayed_min_value_shoot,
                                         vmax_shoot=param.displayed_max_value_shoot,
                                         lognorm_shoot=param.log_scale_shoot,
                                         property_name_root=param.displayed_property_root,
                                         cmap_root=param.color_map_root,
                                         vmin_root=param.displayed_min_value_root,
                                         vmax_root=param.displayed_max_value_root,
                                         lognorm_root=param.log_scale_root,
                                         plotting_leaf_shape=param.plotting_leaf_shape,
                                         simplified_leaf_shape_for_surface_calculation=param.simplified_leaf_shape_for_surface_calculation,
                                         displaying_soil_surface=True,
                                         recording_plot=param.recording_images,
                                         plot_name=os.path.join(plot_directory, 'plot_%.5d.png' % step)
                                         )

                # PLOTTING ONE SINGLE PLANT:
                if name != "population" and param.plotting_single_plant and growth_stages_dict[name] != 0:
                    shoot_g = shoot_MTG_dict[name]
                    root_g = root_MTG_dict[name]

                    # We display the MTG in its new state:
                    plotting_whole_plant(shoot_MTG=shoot_g, root_MTG=root_g,
                                         width=1200, height=1200,
                                         upscaling_factor=zoom_factor,
                                         origin_coordinates=[0.,0.,-sowing_depths_dict[name]*zoom_factor],
                                         x_center=x_center*zoom_factor, y_center=y_center*zoom_factor, z_center=z_center*zoom_factor,
                                         x_cam=x_cam*zoom_factor, y_cam=y_cam*zoom_factor, z_cam=z_cam*zoom_factor,
                                         grid=False,
                                         plant_species=param.plant_genotype,
                                         leaf_shape=param.plotting_leaf_shape,
                                         simplified_leaf_shape_for_surface_calculation=param.simplified_leaf_shape_for_surface_calculation,
                                         background_color=plot_background_color,
                                         single_color_shoot=plot_single_color_shoot,
                                         single_color_root=plot_single_color_root,
                                         property_name_shoot=param.displayed_property_shoot,
                                         cmap_shoot=param.color_map_shoot,
                                         vmin_shoot=param.displayed_min_value_shoot,
                                         vmax_shoot=param.displayed_max_value_shoot,
                                         lognorm_shoot=param.log_scale_shoot,
                                         property_name_root=param.displayed_property_root,
                                         cmap_root=param.color_map_root,
                                         vmin_root=param.displayed_min_value_root,
                                         vmax_root=param.displayed_max_value_root,
                                         lognorm_root=param.log_scale_root,
                                         recording_plot=param.recording_images,
                                         plot_name=os.path.join(plot_directory, 'plot_%.5d.png' % step),
                                         shoots_overlapping_roots=True)
            #///////////////////////////////////////////////////////////////////////////////////////////////////////////

    #ttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt

    print("")
    print("Simulation is done!")

    return
#-----------------------------------------------------------------------------------------------------------------------
########################################################################################################################
########################################################################################################################

# MAIN SIMULATION:
##################

if __name__ == '__main__':

    current_wd = os.getcwd()
    main_simbal_folder = os.path.join(os.getcwd(), os.pardir, os.pardir, os.pardir)
    # main_simbal_folder = 'C:\\Users\\frees\\fspm\\simbal'

    # # RUNNING SIMULATION ON ONE MONOCROP:
    # #####################################
    # # We run the whole simulation according to the instructions and parameters defined in the 'simbal_parameters.xlsx'
    # # file (otherwise, default parameters from the 'simbal_parameters.py' file or from the function will be used):
    # run_simulation_for_one_plant(my_path=main_simbal_folder,
    #                              parameters_file="simbal_parameters.xlsx",
    #                              input_file=os.path.join(main_simbal_folder, "simulation", "inputs", "inputs_over_time.csv"),
    #                              input_file_time_step_in_days=1/24.)

    # RUNNING SIMULATION ON AN HETEROGENEOUS PLANT COVER:
    #####################################################

    # We run the simulation for the heterogeneous plant cover:
    run_simulation_for_a_plant_cover(my_path=main_simbal_folder)
    print("Simulation done!")

#TODO: finish ramification