#  -*- coding: utf-8 -*-

"""
    The script 'RootSimple' contains all functions for running the development anf functionning of the roots.

    :copyright: see AUTHORS.
    :license: see LICENSE for details.
"""

from openalea.mtg import *
from openalea.mtg.traversal import pre_order, post_order
from openalea.rhizodep.model import (volume_and_external_surface_from_radius_and_length,
                                     summing_and_possibly_homogenizing, reinitializing_growth_variables,
                                     update_distance_from_tip, update_surfaces_and_volumes,
                                     potential_growth, segmentation_and_primordia_formation, root_hairs_dynamics,
                                     initiate_mtg)

import simbal_parameters as param
import openalea.rhizodep.model as rhizodep_param

from tools import temperature_modification

from math import sqrt, pi, floor, exp, isnan
from decimal import Decimal

########################################################################################################################
########################################################################################################################

# We define a class "root_system" that will be used to contain the properties of the root system:
class root_system(object):
    # We initiate the object "root_system":
    def __init__(self, g):
        # The object 'root_system' will in particular contain an attribute corresponding to the full root MTG,
        # as well as general variables for the whole root system:
        self.root_MTG = g
        self.root_labile_C_pool = 0. # gC per root system
        self.root_labile_C_deficit = 0. # gC per root system
        self.root_reserve_C_pool = 0. # gC per root system
        self.root_labile_N_pool = 0. # gN per plant
        self.root_labile_N_deficit = 0.  # gN per plant
        self.root_reserve_N_pool = 0. # gN per plant
        self.root_struct_mass = 0. # gDW
        self.root_external_surface = 0. # m2
        self.root_exchange_surface = 0. # m2
        self.root_labile_C_input = 0. # gC per root system
        self.root_N_uptake = 0. # gN per root system

#-----------------------------------------------------------------------------------------------------------------------

def updating_root_labile_C_with_input(root_system, time_step_in_seconds = 60*60):
    """
    This simple function updates the pool of labile C in the root system (gC) with the rate of net C exchange
    from shoots to roots.
    Parameters
    ----------
    root_system: the object containing the general properties of the root system
    time_step_in_seconds: the time step to be used, in seconds

    Returns
    -------
    [Nothing is returned, the root system has been updated with the new labile C pool]

    """

    # We update the labile pool:
    root_system.root_labile_C_pool += root_system.net_C_transfer_from_shoot_to_root_rate * time_step_in_seconds # gC per root system
    if root_system.root_labile_C_pool < 0.:
        print("    WATCH OUT: after exchange with shoots, the labile C pool is negative - we reset it to 0!")
        root_system.root_labile_C_deficit = - root_system.root_labile_C_pool
        root_system.root_labile_C_pool = 0.
    return

#-----------------------------------------------------------------------------------------------------------------------

def updating_exchange_surface(root_system, root_MTG):
    """
    This function covers each element of the root MTG and computes a total surface of exchange for the whole root system.
    Parameters
    ----------
    root_system: the object containing the general properties of the root system
    root_MTG: the root MTG detailing each root element

    Returns
    -------
    [Nothing is returned, the root system has been updated with the new 'root_exchange_surface']
    """

    g = root_MTG

    # We initialize the total exchange surface of the root system with the soil:
    total_exchange_surface = 0.

    # We cover all the vertices in the MTG:
    for vid in g.vertices_iter(scale=1):
        # n represents the vertex:
        n = g.node(vid)
        # If the root element can take up N:
        if n.length > 0.:
            # We calculate the surface of exchange between soil solution and parenchyma:
            S_epid = n.epidermis_surface_without_hairs
            S_hairs = n.living_root_hairs_external_surface
            S_cortex = n.cortical_parenchyma_surface
            S_stele = n.stelar_parenchyma_surface
            cond_walls = n.relative_conductance_walls
            cond_exo = n.relative_conductance_exodermis
            cond_endo = n.relative_conductance_endodermis
            exchange_surface = (S_epid + S_hairs) + cond_walls * (cond_exo * S_cortex + cond_endo * S_stele)
            n.exchange_surface = exchange_surface
            total_exchange_surface += exchange_surface

    root_system.root_exchange_surface = total_exchange_surface

    return
#-----------------------------------------------------------------------------------------------------------------------

def computing_N_uptake_rate(root_system, root_MTG, N_soil = 1e-6, soil_temperature_in_degree_Celsius=10.):

    """
    This function computes the rate of N uptake by roots (gN per root system per second), according to the concentration
    of soil available N and to soil temperature. The pool of labile N in roots is not updated by the function.

    Parameters
    ----------
    root_system: the object containing the general properties of the root system
    root_MTG: the root MTG detailing each root element
    N_soil: the concentration of available N in the soil (gN per gDW of soil)
    soil_temperature_in_degree_Celsius: the mean temperature perceived by the root system (degree Celsius)

    Returns
    -------
    [Nothing is returned, the root system has been updated with the new 'N_uptake_rate']
    """

    # We initialize the total N uptake rate of the whole root system (gN per plant per second):
    N_uptake_rate = 0.

    # We correct the maximal rate of N uptake according to soil temperature:
    corrected_Vmax = param.Vmax_N_uptake * temperature_modification(
        temperature_in_Celsius=soil_temperature_in_degree_Celsius,
        process_at_T_ref=1,
        T_ref=param.N_uptake_T_ref,
        A=param.N_uptake_A,
        B=param.N_uptake_B,
        C=param.N_uptake_C)

    # We cover all the vertices in the MTG:
    for vid in root_MTG.vertices_iter(scale=1):
        # n represents the vertex:
        n = root_MTG.node(vid)
        # If the root element can take up N:
        if n.length > 0. and N_soil > 0.:
            # Then we increase the total rate of N uptake of the root system with the N uptake by the current element:
            N_uptake_rate += n.exchange_surface * corrected_Vmax * N_soil / (param.Km_N_uptake + N_soil)

    # We assign the final result to the whole root system:
    if N_uptake_rate < 0.:
        print("ERROR !!! The computed N uptake rate is negative!")
    else:
        root_system.N_uptake_rate = N_uptake_rate
        # print("N uptake rate is", N_uptake, "gN per plant per second.")
    return

#-----------------------------------------------------------------------------------------------------------------------
def computing_N2_fixation_rate(root_system, root_MTG, soil_temperature_in_degree_Celsius=10.):
    """
    This function computes the rate of nitrogen fixation by nodules that have developed on the root system. The nodules
     are not explicitly represented. Instead, a rate for the whole root system is calculated, taking into account:
      1) the mass and the age of each root segment (only mature root element will develop their full potential),
      2) the concentration of labile N in the root system (fixation rate decreases with increasing labile root N),
      3) the concentration of labile C in the root system (fixation rate only starts above a certain treshold).
    The pool of labile N in roots is not updated by the function.

    Parameters
    ----------
    root_system: the root system containing general variables
    root_MTG: the MTG corresponding to the root system
    soil_temperature_in_degree_Celsius: the mean temperature perceived by the root system (degree Celsius)

    Returns
    -------
    [Nothing is returned, the root system has been updated with the new 'N2_fixation_rate']
    """

    # We initialize the total N2 fixation rate of the whole root system:
    N2_fixation_rate = 0.

    # We correct the maximal rate of N uptake according to soil temperature:
    temperature_corrected_rate = param.N2_fixation_rate_max * temperature_modification(
        temperature_in_Celsius=soil_temperature_in_degree_Celsius,
        process_at_T_ref=1,
        T_ref=param.N2_fixation_T_ref,
        A=param.N2_fixation_A,
        B=param.N2_fixation_B,
        C=param.N2_fixation_C)

    # We make sure that there is enough C to sustain the fixation, and not too much labile N already in the root:
    if (root_system.root_labile_C_pool > param.minimal_labile_C_concentration_for_N2_fixation
            and root_system.root_labile_N_pool <= param.maximal_labile_N_concentration_for_N2_fixation):

        # We calculate a factor of modification based on the level of root labile N
        # (i.e., the more the labile N, the less fixation there will be):
        N_modulator = 1 - root_system.root_labile_N_pool / param.maximal_labile_N_concentration_for_N2_fixation

        # We cover all the vertices in the MTG:
        for vid in root_MTG.vertices_iter(scale=1):
            # n represents the vertex:
            n = root_MTG.node(vid)
            # If the root element has emerged and is old enough to sustain N2 fixation:
            if n.thermal_time_since_emergence > param.minimal_age_for_nodulation:
                # Then we increase the total N2 fixation of the root system with the contribution of the current element,
                # knowing that a young root element will contribute less than a more mature one:
                N2_fixation_rate += (n.struct_mass * temperature_corrected_rate
                                     * max(1, (n.thermal_time_since_emergence/param.maturation_period_for_nodulation))
                                     * N_modulator)

    # We finally assign the computed rate as the corresponding property of the whole root system:
    root_system.N2_fixation_rate = N2_fixation_rate
    # print("N2 fixation rate is", N2_fixation_rate, "gN per plant per second.")

    return

#-----------------------------------------------------------------------------------------------------------------------
def updating_root_labile_N_with_input(root_system, time_step_in_seconds = 1 * (60*60*24)):
    """
    This function simply updates the pool of labile N in roots based on the rates of N uptake and N2 fixation.
    Parameters
    ----------
    root_system: the root system containing general variables
    time_step_in_seconds: the actual time step (seconds)

    Returns
    -------
    [Nothing is returned, the root system has been updated with the new 'root_labile_N_pool']
    """

    # We update the labile pool:
    root_system.root_labile_N_pool += (root_system.N_uptake_rate + root_system.N2_fixation_rate) * time_step_in_seconds # gN per root system
    if root_system.root_labile_N_pool < 0.:
        print("WATCH OUT: after N uptake from soil, the labile N pool is negative - we reset it to 0!")
        root_system.root_labile_N_pool = 0.
    # print("After N uptake, the root system has now", root_system.root_labile_N_pool, "gN available.")
    return

#-----------------------------------------------------------------------------------------------------------------------

def computing_rhizodeposition_rate(root_system, root_MTG,
                                   option="loss_as_input_rate_fraction",
                                   C_rhizodeposition_fraction=0.20,
                                   N_rhizodeposition_fraction=0.10,
                                   time_step_in_seconds=60*60,
                                   soil_temperature_in_degree_Celsius=10.):
    """
    This function computes a rate of rhizodeposition at the scale of the whole root system, based either on a fraction
    of the total input rate of C (or N) that is lost over time, or on a fraction of the total C (or N) labile pools.
    The rate of rhizodeposition can be adjusted by the temperature experienced by the root system.
    The pools of labile C and N are not updated by this function.
    Parameters
    ----------
    root_system: the root system containing general variables
    root_MTG: the MTG corresponding to the root system
    option: if "loss_as_input_rate_fraction", a fraction of C (or N) input rate is lost, else a fraction of root labile C or N is lost
    C_rhizodeposition_fraction: the fraction of C to be lost (gC per gC)
    N_rhizodeposition_fraction: the fraction of N to be lost (gN per gN)
    time_step_in_seconds: the actual time step (seconds)
    soil_temperature_in_degree_Celsius: the mean temperature perceived by the root system (degree Celsius)

    Returns
    -------
    [Nothing is returned, the root system and the root MTG have been updated with the new 'C_rhizodeposition_rate' and 'N_rhizodeposition_rate']
    """

    # We rename the root system:
    r = root_system

    # We compute a modifying factor based on temperature:
    temperature_modifier = temperature_modification(temperature_in_Celsius=soil_temperature_in_degree_Celsius,
                                                    process_at_T_ref=1,
                                                    T_ref=param.permeability_coeff_T_ref,
                                                    A=param.permeability_coeff_A,
                                                    B=param.permeability_coeff_B,
                                                    C=param.permeability_coeff_C)

    # 1) We compute the rhizodeposition at the scale of the whole root system:
    if option == "loss_as_input_rate_fraction":
        # Option 1: rhizodeposition as a fraction of the input.
        r.C_rhizodeposition_rate = (C_rhizodeposition_fraction * temperature_modifier
                                    * r.net_C_transfer_from_shoot_to_root_rate) # gC per root system per second
        r.N_rhizodeposition_rate = (N_rhizodeposition_fraction * temperature_modifier
                                    * r.N_uptake_rate) # gN per root system per second
    else:
        # Option 2: rhizodeposition as a fraction of the labile pool.
        r.C_rhizodeposition_rate = (C_rhizodeposition_fraction * temperature_modifier
                                    * r.root_labile_C_pool / time_step_in_seconds) # gC per root system per second
        r.N_rhizodeposition_rate = (N_rhizodeposition_fraction * temperature_modifier
                                    * r.root_labile_N_pool / time_step_in_seconds) # gN per root system per second

    # 2) We now distribute the rhizodeposition along the roots.
    # We cover all the vertices in the MTG:
    for vid in root_MTG.vertices_iter(scale=1):
        # n represents the vertex:
        n = root_MTG.node(vid)
        # If the root element has a positive length:
        if n.length > 0.:
            # We allocate a part of the rate of rhizodeposition to this element, based on its contribution
            # to total exchange surface:
            relative_fraction = n.exchange_surface / r.root_exchange_surface
            n.C_rhizodeposition_rate = relative_fraction * r.C_rhizodeposition_rate
            n.N_rhizodeposition_rate = relative_fraction * r.N_rhizodeposition_rate
    
    return

#-----------------------------------------------------------------------------------------------------------------------

def computing_maintenance_rate(root_system, soil_temperature_in_degree_Celsius=10.):

    """
    This function computes the rate of maintenance respiration, based on the amount of labile C in the root system
    and on soil temperature.

    Parameters
    ----------
    root_system: the root system containing general variables
    soil_temperature_in_degree_Celsius: the mean temperature perceived by the root system (degree Celsius)

    Returns
    -------
    [Nothing is returned, the root system has been updated with the new 'C_maintenance_rate']

    """

    # We compute a modifying factor based on temperature:
    temperature_modifier = temperature_modification(temperature_in_Celsius=soil_temperature_in_degree_Celsius,
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

    # We compute the rate of maintenance respiration at the scale of the whole root system:
    root_system.root_maintenance_respiration_rate = (corrected_Vmax * root_system.root_labile_C_pool
                                                     / (Km_maintenance + root_system.root_labile_C_pool))

    return

#-----------------------------------------------------------------------------------------------------------------------

def computing_exchange_rates_with_reserve(root_system, soil_temperature_in_degree_Celsius=10.):
    """
    This function computes the rate of net C immobilization into the reserve of the root system (negative if the balance
     results in a net mobilization). It also computes net N immobilisation. The labile or reserve pools are not updated.

    Parameters
    ----------
    root_system: the root system containing general variables
    soil_temperature_in_degree_Celsius: the mean temperature perceived by the root system (degree Celsius)

    Returns
    -------
    [Nothing is returned, the root system has been updated with the new 'net_C_immobilization_rate' and 'net_N_immobilization_rate']
    """

    # We rename the root system:
    r = root_system

    # We initialize temporary rates:
    C_immobilization_as_reserve_rate=0.
    C_mobilization_from_reserve_rate=0.
    N_immobilization_as_reserve_rate=0.
    N_mobilization_from_reserve_rate=0.

    # NET C IMMOBILIZATION:
    #######################

    # We adjust the units of the parameters inherited from RhizoDep (originally expressed in moles of equivalent-hexose):
    max_mobilization_rate = param.max_mobilization_rate * 6 * 12.01 # now in gC per gram of structural mass per second
    Km_mobilization = param.Km_mobilization * 6 * 12.01 # now in gC per gram of structural mass
    max_immobilization_rate = param.max_immobilization_rate * 6 * 12.01 # now in gC per gram of structural mass per second
    Km_immobilization = param.Km_immobilization * 6 * 12.01 # now in gC per gram of structural mass
    root_reserve_C_min = param.C_hexose_reserve_min * 6 * 12.01 # now in gC per gram of structural mass
    root_reserve_C_max = param.C_hexose_reserve_max * 6 * 12.01 # now in gC per gram of structural mass
    root_labile_C_min_for_reserve = param.C_hexose_root_min_for_reserve * 6 * 12.01 # now in gC per gram of structural mass

    # We correct the maximum rates according to soil temperature:
    corrected_max_mobilization_rate = max_mobilization_rate * temperature_modification(
        temperature_in_Celsius=soil_temperature_in_degree_Celsius,
        process_at_T_ref=1,
        T_ref=param.max_mobilization_rate_T_ref,
        A=param.max_mobilization_rate_A,
        B=param.max_mobilization_rate_B,
        C=param.max_mobilization_rate_C)
    corrected_max_immobilization_rate = max_immobilization_rate * temperature_modification(
        temperature_in_Celsius=soil_temperature_in_degree_Celsius,
        process_at_T_ref=1,
        T_ref=param.max_immobilization_rate_T_ref,
        A=param.max_immobilization_rate_A,
        B=param.max_immobilization_rate_B,
        C=param.max_immobilization_rate_C)

    # CALCULATIONS OF THEORETICAL MOBILIZATION RATE:
    # -----------------------------------------------
    if r.root_reserve_C_pool <= root_reserve_C_min:
        C_mobilization_from_reserve_rate = 0.
    else:
        # We calculate the potential mobilization of hexose from reserve (in gC) according to the Michaelis-Menten function:
        C_mobilization_from_reserve_rate = (corrected_max_mobilization_rate * r.root_struct_mass
                                            * r.root_reserve_C_pool / (Km_mobilization + r.root_reserve_C_pool))

    # CALCULATIONS OF THEORETICAL IMMOBILIZATION RATE:
    # -----------------------------------------------
    if r.root_labile_C_pool <= root_labile_C_min_for_reserve or r.root_reserve_C_pool >= root_reserve_C_max:
        # Then there is no possible immobilization in the reserve:
        C_immobilization_as_reserve_rate = 0.
    else:
        # We calculate the potential immobilization:
        C_immobilization_as_reserve_rate = (corrected_max_immobilization_rate * r.root_struct_mass
                                            * r.root_labile_C_pool / (Km_immobilization + r.root_labile_C_pool))

    # NET N IMMOBILIZATION:
    #######################
    #TODO: To be completed with N mobilization/immobilization?

    # RECORDING THE RESULTS:
    # -----------------------
    root_system.net_C_immobilization_rate = C_immobilization_as_reserve_rate - C_mobilization_from_reserve_rate
    root_system.net_N_immobilization_rate = N_immobilization_as_reserve_rate - N_mobilization_from_reserve_rate

    return

#-----------------------------------------------------------------------------------------------------------------------

# Function calculating a satisfaction coefficient for the growth of the whole root system:
# -----------------------------------------------------------------------------------------
def computing_growth_demand_and_satisfaction_coefficients(root_MTG, labile_C_pool, labile_N_pool,
                                                          soil_temperature_in_degree_Celsius=10.):

    # We rename the root MTG:
    g = root_MTG

    # We initialize the sum of individual demands for struct_mass:
    sum_struct_mass_demand = 0.

    # We have to cover each vertex from the apices up to the base one time:
    root_gen = g.component_roots_at_scale_iter(g.root, scale=1)
    root = next(root_gen)

    # We cover all the vertices in the MTG:
    for vid in post_order(g, root):
        # n represents the current root element:
        n = g.node(vid)

        # We calculate the initial volume of the element:
        initial_volume = volume_and_external_surface_from_radius_and_length(g, n, n.initial_radius, n.initial_length)["volume"]
        # We calculate the potential volume of the element based on the potential radius and potential length:
        potential_volume = volume_and_external_surface_from_radius_and_length(g, n, n.potential_radius, n.potential_length)["volume"]

        # The growth demand of the element in struct_mass is calculated (no cost for respiration are included here):
        n.growth_demand_in_struct_mass = (potential_volume - initial_volume) * param.root_tissue_density
        sum_struct_mass_demand += n.growth_demand_in_struct_mass

    # We finally calculate the satisfaction coefficients.
    # If there is no demand for C and N:
    if sum_struct_mass_demand <= 0:
        # Then the satisfaction coefficients are 1 by definition:
        satisf_for_C = 1.
        satisf_for_N = 1.
        return satisf_for_C, satisf_for_N

    # If there is no labile C to produce growth:
    if labile_C_pool <= 0.:
        satisf_for_C = 0.
    else:
        # Otherwise, the satisfaction coefficient correspond to the ratio between the maximal mass obtained from the
        # labile C and the C demand, INCLUDING THE CARBON COST FOR GROWTH RESPIRATION:
        C_demand = sum_struct_mass_demand * param.struct_mass_C_content / param.yield_growth
        satisf_for_C = labile_C_pool / C_demand
        # print("    > The C cost for potential growth should be", C_demand)
        # print("    > The final C cost will be", min(labile_C_pool, C_demand))
    if satisf_for_C > 1.:
        satisf_for_C = 1.

    # If there is no labile N to produce growth:
    if labile_N_pool <= 0.:
        satisf_for_N = 0.
    else:
        # Otherwise, the satisfaction coefficient correspond to the ratio between the maximal mass obtained from the
        # labile N and the N demand:
        satisf_for_N = labile_N_pool / (sum_struct_mass_demand * param.struct_mass_N_content)
    if satisf_for_N > 1.:
        satisf_for_N =1.

    return satisf_for_C, satisf_for_N

#-----------------------------------------------------------------------------------------------------------------------
# Function performing the growth of each element based on potential growth and limitation by labile C & N:
# --------------------------------------------------------------------------------------------------------
def ArchiSimple_growth_with_C_and_N(root_system, root_MTG,
                                    satisf_for_C, satisf_for_N,
                                    time_step_in_seconds= 1 * (60*60*24),
                                    soil_temperature_in_degrees_Celsius=10.):
    """
    This function computes the growth of the root system according to ArchiSimple's rules, but limits the growth
    by the amount of available labile C and N.
    """

    # We initialize the root biomass that will be produced from the labile pool of C & N:
    total_mass_produced = 0.
    # We select the MTG included within the root system r:
    r = root_system
    g = root_MTG
    # We have to cover each vertex from the apices up to the base one time:
    root_gen = g.component_roots_at_scale_iter(g.root, scale=1)
    root = next(root_gen)

    # We consider temperature adjustment:
    temperature_time_adjustment = temperature_modification(temperature_in_Celsius=soil_temperature_in_degrees_Celsius,
                                                           process_at_T_ref=param.root_growth_T_ref,
                                                           A=param.root_growth_A,
                                                           B=param.root_growth_B,
                                                           C=param.root_growth_C)

    # PERFORMING ARCHISIMPLE GROWTH:
    # -------------------------------

    # We cover all the vertices in the MTG:
    for vid in post_order(g, root):

        # n represents the current root element:
        n = g.node(vid)

        # We make sure that the element is not dead and has not already been stopped at the previous time step:
        if n.type == "Dead" or n.type == "Just_dead" or n.type == "Stopped":
            # Then we pass to the next element in the iteration:
            continue
        # We make sure that the root elements at the basis that support adventitious root are not considered:
        if n.type == "Support_for_seminal_root" or n.type == "Support_for_adventitious_root":
            # Then we pass to the next element in the iteration:
            continue

        # We initialize the initial volume of the root element:
        initial_volume = n.volume
        # We calculate the potential volume of the element based on the potential radius and potential length:
        potential_volume = volume_and_external_surface_from_radius_and_length(g, n, n.potential_radius, n.potential_length)["volume"]
        # The final volume is calculated from the potential one using the lower coefficient of satisfaction between C and N:
        limiting_coefficient = min(satisf_for_C, satisf_for_N)
        # NOTE: the satisfaction coefficient for C already includes growth respiration cost.
        # If C is the limiting factor:
        if limiting_coefficient == satisf_for_C:
            # We calculate the volume by further excluding the cost for growth respiration from the amount of available C:
            final_volume = initial_volume + (potential_volume - initial_volume) * satisf_for_C
            # NOTE: this is truly the potential volume, knowing that an additional fraction of C will be used for growth respiration.
        else:
            final_volume = initial_volume + (potential_volume - initial_volume) * satisf_for_N

        # We now calculate the maximal length of the element based on the final volume, while keeping the radius unchanged:
        length_max = final_volume / (pi * n.radius ** 2)

        # ELONGATION IS CONSIDERED FIRST:
        # If the element can elongate:
        if n.potential_length > n.initial_length:
            # If elongation is possible but is limited by the availability of C or N:
            if n.potential_length >= length_max:
                # Then elongation is done up to the limited final volume:
                n.length = length_max
            # Otherwise, elongation can be done up to the full potential:
            else:
                # Elongation is done up to the full potential:
                n.length = n.potential_length
            # The corresponding new volume is calculated:
            volume_after_elongation = (pi * n.radius ** 2) * n.length

        # ACTUAL RADIAL GROWTH IS THEN CONSIDERED:
        # If the radius of the element can increase:
        if n.potential_radius > n.initial_radius and final_volume > volume_after_elongation:

            # We then calculate the corresponding new possible radius corresponding to this maximum volume:
            if n.type == "Root_nodule":
                # If the element corresponds to a nodule, then it we calculate the radius of a theoretical sphere:
                possible_radius = ((final_volume - volume_after_elongation) * 3. / (4. * pi)) ** (1. / 3.)
            else:
                # Otherwise, we calculate the radius of a cylinder. The final volume is calculated as:
                # final volume = final_radius^2 * pi * new_length, i.e.:
                possible_radius = sqrt( final_volume / (n.length * pi) )
            if possible_radius < 0.9999 * n.initial_radius:  # We authorize a difference of 0.01% due to calculation
                # errors!
                print("!!! ERROR: the calculated new radius of element", n.index(), "is lower than the initial one!")
                print("The possible radius was", possible_radius, "and the initial radius was", n.initial_radius)

            # If the maximal radius that can be obtained is lower than the potential radius suggested by the
            # potential growth module:
            if possible_radius <= n.potential_radius:
                # Then radial growth is limited:
                n.radius = possible_radius
            else:
                # Otherwise, radial growth is done up to the full potential and the remaining hexose is calculated:
                n.radius = n.potential_radius

        # RECORDING THE FINAL VOLUME AND MASS:
        n.actual_elongation = n.length - n.initial_length
        n.volume = final_volume
        n.struct_mass = n.volume * param.root_tissue_density
        n.struct_mass_produced = (final_volume - initial_volume) * param.root_tissue_density
        # In case where the root element corresponds to an apex, the distance to the last ramification is increased:
        if n.label == "Apex":
            n.dist_to_ramif += n.actual_elongation
        if n.length < 0 or n.struct_mass < 0:
            print("!!! ERROR: the element", n.index(), "of class", n.label, "has a length of", n.length,
                  "and a mass of", n.struct_mass)

        # INCREMENTING THE TOTAL ROOT MASS THAT HAS BEEN PRODUCED
        total_mass_produced += n.struct_mass_produced

        # MODIFYING SPECIFIC PROPERTIES AFTER ELONGATION:
        # If there has been an actual elongation:
        if n.length > n.initial_length:
            # If the elongated apex corresponded to any primordium that has been allowed to emerge:
            if n.type == "Seminal_root_before_emergence" \
                    or n.type == "Adventitious_root_before_emergence" \
                    or n.type == "Normal_root_before_emergence":
                # We now consider the apex to have emerged:
                n.type = "Normal_root_after_emergence"
                # The exact time since emergence is recorded:
                n.thermal_time_since_emergence = n.thermal_potential_time_since_emergence
                n.actual_time_since_emergence = n.thermal_time_since_emergence / temperature_time_adjustment
                # The actual elongation rate is calculated:
                n.actual_elongation = n.length - n.initial_length
                n.actual_elongation_rate = n.actual_elongation / n.actual_time_since_emergence
                # Note: at this stage, no sugar has been allocated to the emerging primordium itself!
                # if n.type == "Adventitious_root_before_emergence":
                #     print("> A new adventitious root has emerged, starting from element", n.index(), "!")
            elif n.type == "Normal_root_after_emergence":
                # The actual elongation rate is calculated:
                n.actual_elongation = n.length - n.initial_length
                n.actual_elongation_rate = n.actual_elongation / time_step_in_seconds

    # EVENTUALLY, WE UPDATE THE TOTAL MASS OF THE ROOT SYSTEM:
    r.root_struct_mass += total_mass_produced
    # And we compute the consumption rate of C & N to sustain this growth:
    r.C_consumption_by_growth_rate = total_mass_produced * param.struct_mass_C_content / param.yield_growth / time_step_in_seconds
    r.N_consumption_by_growth_rate = total_mass_produced * param.struct_mass_N_content / time_step_in_seconds

    # # AND THE LABILE POOLS OF C AND N:
    # r.root_labile_C_pool -= total_mass_produced * param.struct_mass_C_content / param.yield_growth
    # r.root_labile_N_pool -= total_mass_produced * param.struct_mass_N_content

    return

def root_CN_balance(root_system, time_step_in_seconds = 60.*60):

    # We rename the root system as "r":
    r = root_system

    # We update the pool of labile C, taking into account the possible deficit in C previously recorded,
    # and all processes affecting the C pool:
    r.root_labile_C_pool += (-r.root_labile_C_deficit
                             + (r.net_C_transfer_from_shoot_to_root_rate
                                - r.C_rhizodeposition_rate
                                - r.root_maintenance_respiration_rate
                                - r.net_C_immobilization_rate
                                - r.C_consumption_by_growth_rate) * time_step_in_seconds)
    # We reset the root labile C deficit to 0:
    r.root_labile_C_deficit = 0.

    # We update the pool of labile N, taking into account the possible deficit in N previously recorded,
    # and all processes affecting the N pool:
    # print("Before N balance, the labile N in roots is", r.root_labile_N_pool)
    # print("The deficit in N was", r.root_labile_N_deficit, "and here are the fluxes involved:", r.N_uptake_rate, r.N_rhizodeposition_rate, r.net_N_immobilization_rate, r.net_N_transfer_from_root_to_shoot_rate, r.N_consumption_by_growth_rate)
    net_variation = (- r.root_labile_N_deficit +
                             (r.N_uptake_rate + r.N2_fixation_rate
                              - r.N_rhizodeposition_rate
                              - r.net_N_immobilization_rate
                              - r.net_N_transfer_from_root_to_shoot_rate
                              - r.N_consumption_by_growth_rate) * time_step_in_seconds)
    # print("The net variation of N was", net_variation, "gN")
    r.root_labile_N_pool += net_variation
    # print("After N balance, the labile N in roots is", r.root_labile_N_pool)
    # We reset the root labile N deficit to 0:
    r.root_labile_N_deficit = 0.

    # If the pools of C & N are negative, we reset them to 0 and adjust the deficits that will be taken into account
    # at the next step:
    if r.root_labile_C_pool < 0.:
        r.root_labile_C_deficit = - r.root_labile_C_pool
        r.root_labile_C_pool = 0.
    if r.root_labile_N_pool < 0.:
        r.root_labile_N_deficit = - r.root_labile_N_pool
        r.root_labile_N_pool = 0.

    return

########################################################################################################################
########################################################################################################################

# MAIN SIMULATION:
##################

if __name__ == '__main__':

    print("Performing a simulation by RootSimple alone...")

    # # We define the MTG g:
    # g = initiate_mtg()
    # reinitializing_growth_variables(g)
    # # We update the distance from tip for each root element in each root axis:
    # update_distance_from_tip(g)
    # # We update the surfaces and volumes of each root element:
    # update_surfaces_and_volumes(g)
    #
    # # And we create a root system r from it:
    # r = root_system(g)
    # # We adjust the properties of the whole root system according to g:
    # g_properties = summing_and_possibly_homogenizing(g)
    # r.root_struct_mass = g_properties["total_living_root_struct_mass"]
    # updating_exchange_surface(r, g)
    #
    # r.root_labile_C_pool = 1e-3 * r.root_struct_mass
    # r.root_labile_N_pool = 1e-4 * r.root_struct_mass
    #
    # time_step_in_seconds = 3600.
    # n_steps = 24*10
    # soil_temperature=20
    #
    # # We perform the simulation:
    # for time_step in range(1, n_steps+1):
    #
    #     # TODO: Change this artifical exchange with the shoot!
    #     r.net_C_transfer_from_shoot_to_root_rate = 1e-9
    #     r.net_N_transfer_from_root_to_shoot_rate = 1e-10
    #
    #     print("For step", time_step,":")
    #
    #     # Setting N uptake by the root system:
    #     computing_N_uptake_rate(r, g,
    #                             N_soil = 1e-3,
    #                             soil_temperature_in_degree_Celsius=soil_temperature)
    #     computing_N2_fixation_rate(r, g,
    #                                soil_temperature_in_degree_Celsius=soil_temperature)
    #
    #     # Rhizodeposition and update of labile C:
    #     computing_rhizodeposition_rate(r, g,
    #                                    time_step_in_seconds=time_step_in_seconds,
    #                                    soil_temperature_in_degree_Celsius=soil_temperature)
    #
    #     # Maintenance respiration and update of labile C:
    #     computing_maintenance_rate(r,
    #                                soil_temperature_in_degree_Celsius=soil_temperature)
    #
    #     # Exchange with reserve and update of labile C:
    #     computing_exchange_rates_with_reserve(r,
    #                                           soil_temperature_in_degree_Celsius=soil_temperature)
    #
    #     # Growth:
    #     #--------
    #     # We reset all growth-associated variabes before performing the new growth:
    #     reinitializing_growth_variables(g)
    #     # We calculate the potential growth of each root element according to ArchiSimple:
    #     potential_growth(g, time_step_in_seconds=time_step_in_seconds,
    #                      ArchiSimple=True,
    #                      soil_temperature_in_Celsius=soil_temperature)
    #     # We calculate the "satisfaction coefficients" for growth, depending on the amount of labile C and N left:
    #     satisf_for_C, satisf_for_N = computing_growth_demand_and_satisfaction_coefficients(g,
    #                                                                                        labile_C_pool = r.root_labile_C_pool,
    #                                                                                        labile_N_pool = r.root_labile_N_pool)
    #     print("    The satisfaction coefficient for C is", "{:.2f}".format(Decimal(satisf_for_C)),
    #           "and the one for N is", "{:.2f}".format(Decimal(satisf_for_N)))
    #     # We use the function ArchiSimple_growth to adapt the potential growth to the available labile C and N,
    #     # and we update the labile pools of C and N:
    #     ArchiSimple_growth_with_C_and_N(r, g, satisf_for_C, satisf_for_N, time_step_in_seconds=time_step_in_seconds)
    #     # We proceed to the segmentation of the whole root system
    #     segmentation_and_primordia_formation(g, time_step_in_seconds,
    #                                          soil_temperature_in_Celsius=soil_temperature, random=True,
    #                                          root_order_limitation=False)
    #
    #     #TODO: Include root hairs dynamics once the cost of C & N is properly registered!
    #
    #     # # We modifiy root hairs characteristics according to their specific dynamics:
    #     # root_hairs_dynamics(g, time_step_in_seconds=time_step_in_seconds,
    #     #                     soil_temperature_in_Celsius=soil_temperature)
    #
    #     # We update the distance from tip for each root element in each root axis:
    #     update_distance_from_tip(g)
    #     # We update the surfaces and volumes of each root element:
    #     update_surfaces_and_volumes(g)
    #     updating_exchange_surface(r, g)
    #
    #     # We update the pools of labile C and N:
    #     root_CN_balance(r, time_step_in_seconds=time_step_in_seconds)
    #
    #     print("    Labile C pool in the root system is now", "{:.2E}".format(Decimal(r.root_labile_C_pool)),
    #           "gram of C and labile N pool is", "{:.2E}".format(Decimal(r.root_labile_N_pool)), "gram of N.")
    #     print("    The root system contains", g.nb_vertices(), "different root elements.")
    #     print("")
