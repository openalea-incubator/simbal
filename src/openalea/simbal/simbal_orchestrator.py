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

from openalea.simbal.tools import (loading_parameters_from_file, buildDic, formatted_inputs, temperature_modification,
                   recording_MTG_properties, circle_coordinates, updating_coordinates, plotting_whole_plant)
                   # create_a_scene_with_plants_population)

import os
import pandas as pd
from math import sqrt, trunc, ceil
from decimal import Decimal

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
            print("WATCH OUT: all the C from shoots should have been allocated to roots, i.e.", "{:.2E}".format(Decimal(net_final_C_exchange)), "gC!")
            print("We therefore allocated only", fraction*100,"% of it to the roots." )
        # CASE 2: The flow from root to shoot is positive but there is not enough C in the root to sustain it:
        elif expected_C_exchange < 0. and abs(expected_C_exchange) > net_root_C:
            # Then the amount of C exchanged is limited to the amount of C available in roots:
            net_final_C_exchange = - net_root_C
            print("WATCH OUT: a net transfer of C from roots to shoots of", -net_root_C, "gC has occured!")
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
        elif expected_N_exchange < 0. and -expected_N_exchange > net_shoot_N:
            # Then the amount of N exchanged is limited to the amount of N available in shoots:
            net_final_N_exchange = - net_shoot_N
            print("WATCH OUT: a net transfer of N from shoots to roots of", -net_shoot_N, "gN has occured!")
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
def recording_main_plant_properties(plant, shoot_MTG, root_MTG, time_in_days):
    """
    This function records the current properties of the plant object and returns them in a dictionnary.

    Parameters
    ----------
    plant: the object corresponding to the whole plant, containing shoot and root-related variables
    time_in_days: the time step of the simulation (seconds)

    Returns
    -------
    A dictionnary containing the main properties of the plant object (e.g. labile C and N pools, structural mass
    and rates in shoots and roots
    """
    results_dict = {"Final time (days)": time_in_days,
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
def balance_check(plant, time_step_in_seconds = 60.*60.,
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
                       + plant.actual_C_consumption_rate_by_shoot_growth * (1 - param.yield_growth)
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
        print("!!! ERROR: the C balance in the whole plant is wrong! The net difference between the current state and the expected state is", "{:.2E}".format(Decimal(difference)), "gC!")
    difference = shoot_total_C - (previous_shoot_C + shoot_C_inputs - shoot_C_outputs)
    if abs(difference) > 1e-16:
        print("!!! ERROR: the C balance in the shoots is wrong! The net difference between the current state and the expected state is", "{:.2E}".format(Decimal(difference)), "gC!")
    difference = root_total_C - (previous_root_C + root_C_inputs - root_C_outputs)
    if abs(difference) > 1e-16:
        print("!!! ERROR: the C balance in the roots is wrong! The net difference between the current state and the expected state is", "{:.2E}".format(Decimal(difference)), "gC!")

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
    if abs(difference) > 1e-16:
        print("!!! ERROR: the N balance in the whole plant is wrong! The net difference between the current state and the expected state is", "{:.2E}".format(Decimal(difference)), "gN!")
    difference = shoot_total_N - (previous_shoot_N + shoot_N_inputs - shoot_N_outputs)
    if abs(difference) > 1e-16:
        print("!!! ERROR: the N balance in the shoots is wrong! The net difference between the current state and the expected state is", "{:.2E}".format(Decimal(difference)), "gN!")
    difference = root_total_N - (previous_root_N + root_N_inputs - root_N_outputs)
    if abs(difference) > 1e-16:
        print("!!! ERROR: the N balance in the roots is wrong! The net difference between the current state and the expected state is", "{:.2E}".format(Decimal(difference)), "gN!")

    # We finally return the current amounts of C and N in the shoots and the roots:
    return shoot_total_C, root_total_C, shoot_total_N, root_total_N
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
                   simulation_option="shoot_and_roots",
                   initial_shoot_labile_C_concentration = 1e-2, # in gC per gram of structural mass
                   initial_shoot_labile_N_concentration = 1e-3, # in gN per gram of structural mass
                   initial_root_labile_C_concentration = 1e-2, # in gC per gram of structural mass
                   initial_root_labile_N_concentration = 1e-3 # in gN per gram of structural mass
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
        # my_path = os.path.dirname(os.getcwd())
        # We first define the path and the file to read as a .csv:
        my_path = os.path.join('../../..', os.getcwd())

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

    # We consider the directory containing the attributes of root properties over time:
    root_MTG_directory = os.path.join(my_path, "simulation", "outputs", "root_MTG_properties")
    # If this directory doesn't exist:
    if not os.path.exists(root_MTG_directory):
        # Then we create it:
        os.mkdir(root_MTG_directory)
    else:
        # Otherwise, we delete all the images that are already present inside:
        for root, dirs, files in os.walk(root_MTG_directory):
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
        # We update the parameters from the original python file of RhizoDep with the parameters taken from this dictionnary:
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
    if input_file is None:
        INPUT_FILE = scenario_parameters.get('input_file',
                                             input_file)
    else:
        INPUT_FILE = input_file
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
                                          initial_leaf_length=0.,
                                          leaf_surface_max=0.,
                                          initial_labile_C_content=0,
                                          initial_labile_N_content=0
                                          )
        # We define an arbitrary shoot structural mass (gram), which will remain constant over the simulation:
        plant.shoot_struct_mass = 0.1
        plant.shoot_labile_C_pool = 0.
        plant.shoot_labile_N_pool = 0.
    # Otherwise, a fully functional shoot MTG is created:
    else:
        # We create a MTG "shoot_g" corresponding to the shoot (stem and leaves):
        shoot_g = AllocF.create_shoot_MTG(number_of_cotyledons=2, # We define two cotyledons
                                          number_of_phytomers=0,
                                          initial_internode_radius=rhizodep_param.D_ini / 2.,
                                          initial_hypocotyle_length=0.002,
                                          initial_internode_height=0.,
                                          initial_leaf_length=0.0005,
                                          initial_labile_C_content=initial_shoot_labile_C_concentration,
                                          initial_labile_N_content=initial_shoot_labile_N_concentration
                                          )

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
                print("!!! ERROR: we could not compute the initial properties of the shoot element", vid)
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
        # We define an arbitrary initial root mass (gram) that will remain constant over the simulation:
        plant.root_struct_mass = 0.01
    # Otherwise, a fully functional root MTG is created:
    else:
        # We create another MTG corresponding to the root system:
        root_g = RootSimple.initiate_mtg(initial_segment_length=0.,
                                         initial_apex_length=0.001)
        # We reinitialize the variables in the root MTG:
        RootSimple.reinitializing_growth_variables(root_g)
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
        recording_MTG_properties(root_g, os.path.join(root_MTG_directory, 'MTG_00000.csv'))

    # We also initialize the dictionary containing the general properties of the plant:
    results_dict = recording_main_plant_properties(plant, shoot_g, root_g, time_in_days=0.)
    # We add a few more entries:
    results_dict["Air temperature (degree_Celsius)"] = "NA"
    results_dict["Soil temperature (degree_Celsius)"] = "NA"
    results_dict["Light irradiation (micromol of photons per m2 per second)"] = "NA"
    results_dict["Soil available N (gN per gDW of soil)"] = "NA"
    results_dict["Step"] = 0
    # We initialize a dataframe containing these properties for each time step:
    plant_data_frame = pd.DataFrame(results_dict, index=[0])

    # And we initialize the C and N balance by computing the initial amounts of C and N in the shoots and roots:
    shoot_total_C = (plant.shoot_struct_mass * param.structural_C_content
                     + (plant.shoot_labile_C_pool - plant.shoot_labile_C_deficit)
                     + plant.shoot_reserve_C_pool)
    root_total_C = (plant.root_struct_mass * param.structural_C_content
                    + (plant.root_labile_C_pool - plant.root_labile_C_deficit)
                    + plant.root_reserve_C_pool)
    shoot_total_N = (plant.shoot_struct_mass * param.structural_N_content
                     + (plant.shoot_labile_N_pool - plant.shoot_labile_N_deficit)
                     + plant.shoot_reserve_N_pool)
    root_total_N = (plant.root_struct_mass * param.structural_N_content
                    + (plant.root_labile_N_pool - plant.root_labile_N_deficit)
                    + plant.root_reserve_N_pool)

    # DEFINING THE PERIOD OF SIMULATION:
    # ----------------------------------

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
        # We then initialize the step and time at 1:
        initial_step_number = 1
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
        initial_time_in_days = initial_time + (step - 1) * TIME_STEP_IN_DAYS
        final_time_in_days = initial_time + (step) * TIME_STEP_IN_DAYS

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

        print("    > Computing shoot-to-root C transfer...")
        # We compute the rate of net C transfer from shoots to roots:
        computing_shoot_root_exchange_rates(plant, computing_C_exchange=True, computing_N_exchange=False,
                                            C_option="by_diffusion")
        print("    > The net transfer of C is",
              "{:.2E}".format(Decimal(plant.net_C_transfer_from_shoot_to_root_rate * time_step_in_seconds)), "gC.")
        # Note that only the rate has been computed here, the labile C pools have not been updated yet!
        # We limit the transfer of C from shoots to roots (or from roots to shoots) to the maximal amount of labile C
        # available in the source compartment.
        # We first calculate the expected amount of C that should be exchanged based on the flux:
        expected_C_exchange = plant.net_C_transfer_from_shoot_to_root_rate * time_step_in_seconds

        # CASE 1: The flux from shoot to root is positive, but there is not enough C in the shoot to sustain it:
        if expected_C_exchange > 0. and expected_C_exchange > (plant.shoot_labile_C_pool - plant.shoot_labile_C_deficit):
            # Then the amount of C exchanged is limited to the amount of C available in shoots:
            final_C_exchange = plant.shoot_labile_C_pool
            #  The amount of labile C in shoots and roots is updated according to the new computed exchange:
            if SIMULATION_OPTION != "roots_only":
                plant.shoot_labile_C_pool -= final_C_exchange
            if SIMULATION_OPTION != "shoot_only":
                plant.root_labile_C_pool += final_C_exchange
        # CASE 2: The flux from root to shoot is positive but there is not enough C in the root to sustain it:
        elif expected_C_exchange < 0. and abs(expected_C_exchange) > (plant.root_labile_C_pool - plant.root_labile_C_deficit):
            final_C_exchange = plant.root_labile_C_pool
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

        # TRANSFER OF N FROM ROOTS TO SHOOTS:
        #####################################

        print("    > Computing root-to-shoot N transfer...")
        # We compute the rate of net N transfer from roots to shoots:
        computing_shoot_root_exchange_rates(plant, computing_C_exchange=False, computing_N_exchange=True,
                                            N_option="by_diffusion")
        print("    > The net transfer of N is",
              "{:.2E}".format(Decimal(plant.net_N_transfer_from_root_to_shoot_rate * time_step_in_seconds)),
              "gN.")
        # We then limit the transfer of N from roots to shoots (or from shoots to roots) to the maximal amount of N
        # available in the source compartment.
        # We first calculate the expected amount of N that should be exchanged based on the exchange rate that has been computed:
        expected_N_exchange = plant.net_N_transfer_from_root_to_shoot_rate * time_step_in_seconds

        # CASE 1: The flow of N from root to shoot is positive, but there is not enough N in the root to sustain it.
        if expected_N_exchange > 0. and expected_N_exchange > (plant.root_labile_N_pool - plant.root_labile_N_deficit):
            # Then the amount of N exchanged is limited to the amount of N available in roots:
            final_N_exchange = plant.root_labile_N_pool
            # The amount of labile N in shoots and roots is updated according to the new computed exchange:
            if SIMULATION_OPTION != "roots_only":
                plant.shoot_labile_N_pool += final_N_exchange
            if SIMULATION_OPTION != "shoot_only":
                plant.root_labile_N_pool -= final_N_exchange
        # CASE 2: The flux from shoot to root is positive but there is not enough N in the shoot to sustain it:
        elif expected_N_exchange < 0. and abs(expected_N_exchange) > (
                plant.shoot_labile_N_pool - plant.shoot_labile_N_deficit):
            # Then the amount of N exchanged is limited to the amount of N available in shoots:
            final_N_exchange = plant.shoot_labile_N_pool
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

        if SIMULATION_OPTION != "roots_only":

            # We register the initial state of the labile C & N pools in the shoot:
            initial_shoot_labile_C = plant.shoot_labile_C_pool - plant.shoot_labile_C_deficit
            initial_shoot_labile_N = plant.shoot_labile_N_pool - plant.shoot_labile_N_deficit

            print("    The shoot system initially contains", shoot_g.nb_vertices()-1, "elements.")
            print("\n    Shoot labile C pool is initially",
                  "{:.2E}".format(Decimal(initial_shoot_labile_C)),
                  "gC and shoot labile N pool is",
                  "{:.2E}".format(Decimal(initial_shoot_labile_N)), "gN.")

            # LIGHT INTERCEPTION:
            #####################
            # We calculate light interception based on the shoot architecture, for each organ:
            total_PAR_intercepted, sum_of_surfacic_PAR = \
                AllocF.computing_light_interception(shoot_g,
                                                    incoming_light_density=light_input,
                                                    number_of_plants_per_square_meter=plant_density_in_square_meter,
                                                    inter_row_distance=inter_row_distance,
                                                    scene_length=row_length,
                                                    scene_width=scene_width,
                                                    leaf_shape=param.leaf_shape,
                                                    averaging_light_intensity_share=True,
                                                    period_of_light_intensity_averaging_in_days=0.5,
                                                    time_step_in_days=time_step_in_days
                                                    )
            # We record the PAR intercepted in the properties of the plant:
            plant.total_PAR_intercepted = total_PAR_intercepted

            # DISTRIBUTING C & N IN EACH SHOOT ORGAN:
            #########################################

            # Reallocating labile N within each shoot organ:
            #-----------------------------------------------
            # The total amount of labile N to redistribute into the shoot organs is defined as:
            shoot_available_N = plant.shoot_labile_N_pool - plant.shoot_labile_N_deficit
            if shoot_available_N <= 0.:
                print("WATCH OUT: there is no labile N left in the shoot to be redistributed!")
                # We cover each organ of the shoot:
                for vid in shoot_g.vertices_iter(scale=1):
                    # We access the current organ "n":
                    n = shoot_g.node(vid)
                    if n.surface <= 0.:
                        continue
                    labile_N_target = 0.
                    # We calculate the rate of net N input to this organ corresponding from this redistribution:
                    n.shoot_N_distribution_rate = (labile_N_target - (
                            n.labile_N - n.labile_N_deficit)) / time_step_in_seconds
                    # We now redefine the local pool of labile N according to this target level:
                    n.labile_N = 0.
                    n.labile_N_deficit = -labile_N_target
            else:
                # We cover each organ of the shoot:
                for vid in shoot_g.vertices_iter(scale=1):
                    # We access the current organ "n":
                    n = shoot_g.node(vid)
                    if n.surface <= 0.:
                        continue
                    # We define the fraction of shoot available N allocated to this organ based on the previously
                    # computed light intensity share:
                    fraction = n.light_intensity_share
                    labile_N_target = shoot_available_N * fraction
                    # We calculate the rate of net N input to this organ corresponding from this redistribution:
                    n.shoot_N_distribution_rate = (labile_N_target - (
                            n.labile_N - n.labile_N_deficit)) / time_step_in_seconds
                    # We now redefine the local pool of labile N according to this target level:
                    n.labile_N = labile_N_target
                    n.labile_N_deficit = 0.

            # Redistributing labile C into each shoot organ:
            #-----------------------------------------------
            # The total amount of labile C to redistribute into the shoot organs is defined as:
            shoot_available_C = plant.shoot_labile_C_pool - plant.shoot_labile_C_deficit
            if shoot_available_C <= 0.:
                print("WATCH OUT: there is no labile C left in the shoot to be redistributed!")
                # We cover each organ n of the shoot:
                for vid in shoot_g.vertices_iter(scale=1):
                    # We access the current organ "n":
                    n = shoot_g.node(vid)
                    if n.surface <= 0.:
                        continue
                    labile_C_target = 0.
                    # We calculate the rate of net N input to this organ corresponding from this redistribution:
                    n.shoot_C_distribution_rate = (labile_C_target - (
                            n.labile_C - n.labile_C_deficit)) / time_step_in_seconds
                    # We now redefine the local pool of labile C according to this target level:
                    n.labile_C_deficit = -labile_C_target
                    n.labile_C = 0.
            else:
                # We cover each organ n of the shoot:
                for vid in shoot_g.vertices_iter(scale=1):
                    # We access the current organ "n":
                    n = shoot_g.node(vid)
                    if n.surface <= 0.:
                        continue
                    # We define the fraction of shoot available N allocated to this organ based on the previously
                    # computed light intensity share:
                    fraction = n.light_intensity_share
                    labile_C_target = shoot_available_C * fraction
                    # We calculate the rate of net N input to this organ corresponding from this redistribution:
                    n.shoot_C_distribution_rate = (labile_C_target - (
                            n.labile_C - n.labile_C_deficit)) / time_step_in_seconds
                    # We now redefine the local pool of labile C according to this target level:
                    n.labile_C = labile_C_target
                    n.labile_C_deficit = 0.

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

            # POSSIBILE SUPPLY OF RESOURCES FROM THE GERMINATING SEED:
            # If time is lower than the prescribed duration of supply from the seed:
            if final_time_in_days * (60.*60.*24.) <= param.seed_supply_duration * temperature_time_adjustment_shoot:
                # Then we allow to compute the supply of C and N resources from the germinating seed to the new seedling:
                AllocF.computing_supply_from_seed_rate(plant, air_temperature_in_degree_Celsius = air_temperature)

            # MAINTENANCE IN THE SHOOT ORGANS:
            # We calculate the C respiration cost associated to maintenance in each organ:
            total_maintenance_rate = AllocF.computing_maintenance_rate(shoot_g,
                                                                       air_temperature_in_degree_Celsius=air_temperature)
            plant.shoot_maintenance_respiration_rate = total_maintenance_rate
            # print("    The maintenance will consume",
            #       "{:.2E}".format(Decimal(total_maintenance_rate*time_step_in_seconds)),
            #       "gC over the shoot system...")

            # GROWTH IN THE SHOOT ORGANS:
            #############################

            # ORGAN EXTENSION; We calculate the possible increase in volume of each organ:
            AllocF.potential_shoot_growth(shoot_g,
                                          SFQN=param.SFQN, # TODO: is SFQN constant over the whole range of (0, LA_max)? Is it two distinct curves?
                                          leaf_minimal_thickness=param.leaf_minimal_thickness,
                                          elong_internode=param.elong_internode,
                                          time_step_in_seconds=TIME_STEP_IN_DAYS * 60. * 60. * 24.,
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
                                           time_step_in_seconds=TIME_STEP_IN_DAYS * 60. * 60. * 24.,
                                           air_temperature_in_degree_Celsius=air_temperature)
            # At this stage, the shoot structure has been updated and the consumption rates have been calculated
            # for each organ, but no actual consumption of labile C&N at the scale of organ or shoot system has occured yet.
            print("    The extension of organs will consume",
                  "{:.2E}".format(Decimal(C_consumption_rate_for_extension * time_step_in_seconds)),
                  "gC and",
                  "{:.2E}".format(Decimal(N_consumption_rate_for_extension * time_step_in_seconds)),
                  "gN over the shoot system...")

            # NEW PHYTOMERS: We now consider the initiation of the structure of new phytomers, based on the phyllochrone.
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
                        initial_internode_length=param.initial_internode_length,
                        time_step_in_seconds=time_step_in_seconds)
                    print("> A new phytomer has been formed! The C & N costs were",
                          "{:.2E}".format(Decimal(C_cost_for_new_phytomers)), "gC and",
                          "{:.2E}".format(Decimal(N_cost_for_new_phytomers)),"gN.")
                    # We increase the total costs with the cost of this new phytomer:
                    total_C_costs_for_new_phytomers += C_cost_for_new_phytomers
                    total_N_costs_for_new_phytomers += total_N_costs_for_new_phytomers
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

            # We perform the C&N balance both within each organ and at the shoot system scale:
            AllocF.shoot_CN_balance_at_organ_and_system_scales(shoot_system=plant,
                                                               shoot_MTG=shoot_g,
                                                               time_step_in_seconds=time_step_in_seconds)

            print("    The shoot system now contains", shoot_g.nb_vertices()-1, "different shoot elements.")
            print("    Shoot labile C pool is now",
                  "{:.2E}".format(Decimal(plant.shoot_labile_C_pool - plant.shoot_labile_C_deficit)),
                  "gC and labile N pool is",
                  "{:.2E}".format(Decimal(plant.shoot_labile_N_pool - plant.shoot_labile_N_deficit)), "gN.")

        # CHECKING C & N BALANCE:
        #########################

        # We compute the current amounts of total C and N in the shoots and roots, and compare it with the previous ones
        # when altered by the new inputs and outputs computed at this time step:
        shoot_total_C, root_total_C, shoot_total_N, root_total_N \
            = balance_check(plant,
                            time_step_in_seconds = time_step_in_seconds,
                            previous_shoot_C = shoot_total_C,
                            previous_root_C = root_total_C,
                            previous_shoot_N = shoot_total_N,
                            previous_root_N = root_total_N)
        print("")

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

        zoom_factor = 100.
        # We display the MTG in its new state:
        plotting_whole_plant(shoot_MTG=shoot_g, root_MTG=root_g,
                             width=1200, height=1200,
                             upscaling_factor=zoom_factor,
                             x_center=param.x_center*zoom_factor,
                             y_center=param.y_center*zoom_factor,
                             z_center=param.z_center*zoom_factor,
                             x_cam=x_cam*zoom_factor, y_cam=y_cam*zoom_factor, z_cam=z_cam*zoom_factor,
                             grid=False,
                             leaf_shape=param.leaf_shape,
                             background_color=[94, 76, 64],
                             single_color_shoot=[30,120,60],
                             single_color_root=[200,200,100],
                             # single_color_shoot=None,
                             # single_color_root=None,
                             # property_name_shoot='PAR_intercepted', cmap_shoot='inferno', vmin_shoot=0, vmax_shoot=400, lognorm_shoot=False,
                             # property_name_root='radius', cmap_root='jet', vmin_root=1e-4, vmax_root=1e-2, lognorm_root=True,
                             recording_plot=True,
                             plot_name=os.path.join(plot_directory, 'plot_%.5d.png' % step),
                             shoots_overlapping_roots=True)

        # RECORDING THE PROPERTIES OF THE SHOOT MTG FOR THIS TIME STEP:
        # We can record the main properties of the plant:
        recording_MTG_properties(shoot_g, file_name = os.path.join(shoot_MTG_directory, 'MTG_%.5d.csv' % step))

        # RECORDING THE PROPERTIES OF THE ROOT MTG FOR THIS TIME STEP:
        # We can record the main properties of the plant:
        recording_MTG_properties(root_g, file_name = os.path.join(root_MTG_directory, 'MTG_%.5d.csv' % step))

        # RECORDING THE FINAL RESULTS AT PLANT SCALE:
        results_dict = recording_main_plant_properties(plant, shoot_g, root_g, final_time_in_days)
        results_dict["Air temperature (degree_Celsius)"] = air_temperature
        results_dict["Soil temperature (degree_Celsius)"] = soil_temperature
        results_dict["Light irradiation (micromol of photons per m2 per second)"] = light_input
        results_dict["Soil available N (gN per gDW of soil)"] = soil_available_nitrogen
        results_dict["Step"] = step
        new_df = pd.DataFrame([results_dict.values()], columns=results_dict.keys())
        plant_data_frame = pd.concat([plant_data_frame, new_df])
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

    current_wd = os.getcwd()
    main_simbal_folder = os.path.join(os.getcwd(), os.pardir, os.pardir, os.pardir)
    # main_simbal_folder = 'C:\\Users\\frees\\fspm\\simbal'
    # We run the whole simulation according to the instructions and parameters defined in the 'simbal_parameters.xlsx'
    # file (otherwise, default parameters from the 'simbal_parameters.py' file or from the function will be used):
    run_simulation(my_path=main_simbal_folder,
                   input_file=os.path.join(main_simbal_folder, "simulation", "inputs", "inputs_over_time.csv"),
                   input_file_time_step_in_days=1/24.)

#TODO: vérifier l'apparition des phytomères au cours du temps
