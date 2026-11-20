#  -*- coding: utf-8 -*-

"""
    The script 'RootSimple' contains all functions for running the development and functionning of the roots.

    :copyright: see AUTHORS.
    :license: see LICENSE for details.
"""
# General imports
from math import sqrt, pi, floor, exp, isnan
from decimal import Decimal

from openalea.mtg import *
from openalea.mtg.traversal import pre_order, post_order

from openalea.rhizodep.model import volume_and_external_surface_from_radius_and_length

# Loacal or global imports?
from openalea.simbal import simbal_parameters as param
from openalea.simbal.tools import temperature_modification


########################################################################################################################
########################################################################################################################

# We define a class "root_system" that will be used to contain the properties of the root system:
class root_system:
    """
    This class defines a whole root system, defined with general root properties. It does not correspond to
    a MTG, and can be used to store generic root properties at the scale of the whole root system.
    """
    # We initiate the object "root_system":
    def __init__(self, g):
        # The object 'root_system' will in particular contain an attribute corresponding to the full root MTG,
        # as well as general variables for the whole root system:
        self.root_labile_C_pool = 0. # gC per root system
        self.root_labile_C_deficit = 0. # gC per root system
        self.root_reserve_C_pool = 0. # gC per root system
        self.root_labile_N_pool = 0. # gN per plant
        self.root_labile_N_deficit = 0.  # gN per plant
        self.root_reserve_N_pool = 0. # gN per plant
        self.root_struct_mass = 0. # gDW
        self.root_external_surface = 0. # m2
        self.root_exchange_surface = 0. # m2
        
        self.net_C_transfer_from_shoot_to_root_rate = 0. # gC per plant per second
        self.net_N_transfer_from_root_to_shoot_rate = 0. # gN per plant per second
        self.N_uptake_rate = 0. # gN per plant per second
        self.N2_fixation_rate = 0. # gN per plant per second
        self.N2_fixation_rate_C_cost = 0. # gC per plant per second
        self.C_rhizodeposition_rate = 0. # gC per plant per second
        self.N_rhizodeposition_rate = 0. # gN per plant per second
        self.root_maintenance_respiration_rate = 0. # gC per plant per second
        self.net_root_C_immobilization_rate = 0. # gC per plant per second
        self.net_root_N_immobilization_rate = 0. # gN per plant per second
        self.actual_C_consumption_rate_by_root_growth = 0. # gC per plant per second
        self.actual_N_consumption_rate_by_root_growth = 0. # gN per plant per second
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function for calculating root exchange surface in each root segment and at the scale of the whole root system:
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
# Function for calculating N uptake by each root segment and the resulting uptake at the scale of the whole root system:
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

    # We also correct the maximal rate of N uptake according to root labile N (possible inhibition):
    root_N_conc = (root_system.root_labile_N_pool - root_system.root_labile_N_deficit) / root_system.root_struct_mass
    inhibiting_labile_N_concentration = 0.1
    corrected_Vmax = corrected_Vmax * (1 - root_N_conc / inhibiting_labile_N_concentration)
    if corrected_Vmax < 0.:
        corrected_Vmax = 0.

    # We cover all the vertices in the MTG:
    for vid in root_MTG.vertices_iter(scale=1):
        # n represents the vertex:
        n = root_MTG.node(vid)
        # We reset the local rate of N uptake:
        n.N_uptake_rate = 0.
        # If the root element can take up N:
        if n.length > 0. and N_soil > 0.:
            # Then we increase the total rate of N uptake of the root system with the N uptake by the current element:
            N_uptake_rate += n.exchange_surface * corrected_Vmax * N_soil / (param.Km_N_uptake + N_soil)
            # And we record locally the N uptake by this element:
            n.N_uptake_rate = N_uptake_rate

    # We assign the final result to the whole root system:
    if N_uptake_rate < 0.:
        print("ERROR !!! The computed N uptake rate is negative!")
    else:
        root_system.N_uptake_rate = N_uptake_rate
        # print("N uptake rate is", N_uptake, "gN per plant per second.")
    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function for calculating N2 fixation at the scale of the whole root system:
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

    # We initialize the total N2 fixation rate of the whole root system (gN per plant per second) and the corresponding
    # C cost (gC per plant per second):
    total_N2_fixation_rate = 0.
    total_N2_fixation_rate_C_cost = 0.

    # We make sure that there is enough C to sustain the fixation, and not too much labile N already in the root:
    if (root_system.root_labile_C_pool > param.minimal_labile_C_concentration_for_N2_fixation
            and root_system.root_labile_N_pool <= param.maximal_labile_N_concentration_for_N2_fixation):

        # We correct the maximal rate of N uptake according to soil temperature:
        temperature_corrected_rate = param.N2_fixation_rate_max * temperature_modification(
            temperature_in_Celsius=soil_temperature_in_degree_Celsius,
            process_at_T_ref=1,
            T_ref=param.N2_fixation_T_ref,
            A=param.N2_fixation_A,
            B=param.N2_fixation_B,
            C=param.N2_fixation_C)

        # We calculate a factor of modification based on the level of root labile N
        # (i.e., the more the labile N, the less fixation there will be):
        N_modulator = 1 - root_system.root_labile_N_pool / param.maximal_labile_N_concentration_for_N2_fixation

        # We cover all the vertices in the MTG:
        for vid in root_MTG.vertices_iter(scale=1):
            # n represents the vertex:
            n = root_MTG.node(vid)
            # We reset the local rate of N2 fixation and the corresponding cost:
            n.N2_fixation_rate = 0.
            n.N2_fixation_rate_C_cost = 0.
            # If the root element has emerged and is old enough to sustain N2 fixation:
            if n.thermal_time_since_emergence > param.minimal_age_for_nodulation:
                # Then we calculate the contribution of the current element, knowing that a young root element will
                # contribute less than a more mature one:
                N2_fixation_rate = (n.struct_mass * temperature_corrected_rate
                                     * max(1, (n.thermal_time_since_emergence/param.maturation_period_for_nodulation))
                                     * N_modulator)
                # We calculate the corresponding cost in C (gC per plant per second):
                N2_fixation_rate_C_cost = N2_fixation_rate * param.N2_fixation_cost
                # We increment the total fixation rate and cost:
                total_N2_fixation_rate += N2_fixation_rate
                total_N2_fixation_rate_C_cost += N2_fixation_rate_C_cost
                # We record locally the rate of N2 fixation and the corresponding cost:
                n.N2_fixation_rate = N2_fixation_rate
                n.N2_fixation_rate_C_cost = N2_fixation_rate_C_cost

    # We finally assign the computed rate as the corresponding property of the whole root system:
    root_system.N2_fixation_rate = total_N2_fixation_rate
    root_system.N2_fixation_rate_C_cost = total_N2_fixation_rate_C_cost
    # print("N2 fixation rate is", N2_fixation_rate, "gN per plant per second.")

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function for calculating rhizodeposition at the scale of the whole root system:
def computing_rhizodeposition_rate(root_system, root_MTG,
                                   option="loss_as_input_rate_fraction",
                                   C_rhizodeposition_fraction=0.30,
                                   N_rhizodeposition_fraction=0.10,
                                   time_step_in_seconds=60*60,
                                   soil_temperature_in_degree_Celsius=10.,
                                   rhizodeposits_distribution=True):
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
    rhizodeposits_distribution: if True, the function will distribute the rate of rhizodeposition among the different
    root elements of the root MTG, based on their exchange surface

    Returns
    -------
    [Nothing is returned, the root system and the root MTG have been updated with the new 'C_rhizodeposition_rate' and 'N_rhizodeposition_rate']
    """

    # We rename the root system and the root MTG:
    r = root_system
    g = root_MTG

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

    if rhizodeposits_distribution:
        # 2) We now distribute the rhizodeposition along the roots, for allowing a spatialization of rhizodeposit
        # distribution in the soil.
        # We cover all the vertices in the MTG:
        for vid in g.vertices_iter(scale=1):
            # n represents the vertex:
            n = g.node(vid)
            # If the root element has a positive length:
            if n.length > 0.:
                # We allocate a part of the rate of rhizodeposition to this element, based on its contribution
                # to total exchange surface:
                relative_fraction = n.exchange_surface / r.root_exchange_surface
                n.C_rhizodeposition_rate = relative_fraction * r.C_rhizodeposition_rate
                n.N_rhizodeposition_rate = relative_fraction * r.N_rhizodeposition_rate
    
    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function for calculating maintenance at the scale of the whole root system:
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
#-----------------------------------------------------------------------------------------------------------------------
# Function for computing the rate of exchange between labile C & N pools and reserve C & N pools:
def computing_exchange_rates_with_reserve(root_system, soil_temperature_in_degree_Celsius=10.):
    """
    This function computes the rate of net C immobilization into the reserve of the root system (negative if the balance
     results in a net mobilization). It also computes net N immobilization. The labile or reserve pools are not updated.

    Parameters
    ----------
    root_system: the root system object containing general variables
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

    # We adjust the units of the parameters inherited from RhizoDep (originally expressed in moles of equivalent-hexose):
    max_mobilization_rate = param.max_mobilization_rate * 6 * 12.01 * param.reserve_CN_ratio # now in gN per gram of structural mass per second
    Km_mobilization = param.Km_mobilization * 6 * 12.01 * param.reserve_CN_ratio # now in gN per gram of structural mass
    max_immobilization_rate = param.max_immobilization_rate * 6 * 12.01 * param.reserve_CN_ratio # now in gN per gram of structural mass per second
    Km_immobilization = param.Km_immobilization * 6 * 12.01 * param.reserve_CN_ratio # now in gN per gram of structural mass
    root_reserve_N_min = param.C_hexose_reserve_min * 6 * 12.01 * param.reserve_CN_ratio # now in gN per gram of structural mass
    root_reserve_N_max = param.C_hexose_reserve_max * 6 * 12.01 * param.reserve_CN_ratio # now in gN per gram of structural mass
    root_labile_N_min_for_reserve = param.C_hexose_root_min_for_reserve * 6 * 12.01 * param.reserve_CN_ratio # now in gN per gram of structural mass

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
    if r.root_reserve_N_pool <= root_reserve_N_min:
        N_mobilization_from_reserve_rate = 0.
    else:
        # We calculate the potential mobilization of N from reserve according to the Michaelis-Menten function:
        N_mobilization_from_reserve_rate = (corrected_max_mobilization_rate * r.root_struct_mass
                                            * r.root_reserve_N_pool / (Km_mobilization + r.root_reserve_N_pool))

    # CALCULATIONS OF THEORETICAL IMMOBILIZATION RATE:
    # -----------------------------------------------
    if r.root_labile_N_pool <= root_labile_N_min_for_reserve or r.root_reserve_N_pool >= root_reserve_N_max:
        # Then there is no possible immobilization in the reserve:
        N_immobilization_as_reserve_rate = 0.
    else:
        # We calculate the potential immobilization:
        N_immobilization_as_reserve_rate = (corrected_max_immobilization_rate * r.root_struct_mass
                                            * r.root_labile_N_pool / (Km_immobilization + r.root_labile_N_pool))

    # RECORDING THE RESULTS:
    # -----------------------
    root_system.net_root_C_immobilization_rate = C_immobilization_as_reserve_rate - C_mobilization_from_reserve_rate
    root_system.net_root_N_immobilization_rate = N_immobilization_as_reserve_rate - N_mobilization_from_reserve_rate

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function calculating a satisfaction coefficient for the growth of the whole root system:
def computing_growth_demand_and_satisfaction_coefficients(root_MTG, labile_C_pool, labile_N_pool):
    """
    This function computes the root growth demand and the satisfaction coefficients for labile C and N.

    Parameters
    ----------
    root_MTG: the root MTG containing the root segments
    labile_C_pool: the available labile C that can be used for growth (gC per root system)
    labile_N_pool: the available labile N that can be used for growth (gN per root system)

    Returns
    -------
    satisf_for_C: the coefficient of growth satisfaction for labile C (gC used per gC required)
    satisf_for_N: the coefficient of growth satisfaction for labile N (gN used per gN required)
    """

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
        C_demand = sum_struct_mass_demand * param.structural_C_content / param.yield_growth
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
        satisf_for_N = labile_N_pool / (sum_struct_mass_demand * param.structural_N_content)
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

    Parameters
    ----------
    root_system: the root system object containing its general properties (not the root MTG!)
    root_MTG: the root MTG containing the root segments
    satisf_for_C: the coefficient of growth satisfaction for labile C (gC used per gC required)
    satisf_for_N: the coefficient of growth satisfaction for labile N (gN used per gN required)
    time_step_in_seconds: the actual time step of the simulation (seconds)
    soil_temperature_in_degrees_Celsius: the temperature in degrees Celsius of the soil

    Returns
    -------
    [Nothing; the properties of the root MTG and the root system have been updated]
    """

    # We make sure that both C & N resources are available for growth:
    if satisf_for_C <= 0. or satisf_for_N <= 0.:
        return

    # We initialize the root biomass that will be produced from the labile pool of C & N:
    total_mass_produced = 0.
    # We do the same for the necromass produced at this time step by root mortality:
    total_necromass_produced = 0.
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
        if n.type == "Dead" or n.type == "Stopped": #TODO: What about the cost for growth of a root element that has stopped elongating?
            # Then we pass to the next element in the iteration:
            continue
        # If the root element has just died, we record its structural mass in the total necromass produced:
        if n.type == "Just_dead":
            total_necromass_produced += n.struct_mass
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
            # We calculate the volume according to the amount of available C:
            final_volume = initial_volume + (potential_volume - initial_volume) * satisf_for_C
            # NOTE: this is truly the potential volume, knowing that an additional fraction of C will be used for growth respiration.
        else:
            # Otherwise, the final volume is computed according to the availability of N:
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
    r.root_struct_necromass += total_necromass_produced
    r.root_mortality_rate = total_necromass_produced / time_step_in_seconds
    # And we compute the consumption rate of C & N to sustain this growth:
    r.actual_C_consumption_rate_by_root_growth = total_mass_produced * param.structural_C_content / param.yield_growth / time_step_in_seconds
    r.actual_N_consumption_rate_by_root_growth = total_mass_produced * param.structural_N_content / time_step_in_seconds

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function for reinitializing all growth-related variables at the beginning or end of a time step:
def reinitializing_growth_variables(g):
    """
    This function re-initializes different growth-related variables (e.g. potential growth variables) of the root MTG "g".
    """

    # We cover all the vertices in the MTG:
    for vid in g.vertices_iter(scale=1):
        # n represents the vertex:
        n = g.node(vid)

        # We set to 0 the growth-related variables:
        n.resp_growth = 0.
        n.struct_mass_produced = 0.
        n.root_hairs_struct_mass_produced = 0.
        n.actual_elongation = 0.
        n.actual_elongation_rate = 0.
        n.growth_demand_in_struct_mass = 0.

        # We make sure that the initial values of length, radius and struct_mass are correctly initialized:
        n.initial_length = n.length
        n.potential_length = n.length
        n.initial_radius = n.radius
        n.potential_radius = n.radius
        n.theoretical_radius = n.radius
        n.initial_struct_mass = n.struct_mass
        n.initial_living_root_hairs_struct_mass = n.living_root_hairs_struct_mass
        n.initial_external_surface = n.external_surface
        n.initial_living_root_hairs_external_surface = n.living_root_hairs_external_surface

    return
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
# Function for computing new amounts of labile C & N in the whole root system:
def root_CN_balance(root_system, time_step_in_seconds = 60.*60.):
    """
    This function computes the new amounts of labile C and labile N in the root system object, based on the initial
    amounts and the rates of C&N exchanges already stored in the root system object.

    Parameters
    ----------
    root_system: the object corresponding to the whole root system
    time_step_in_seconds: the time step (seconds)

    Returns
    -------
    [Nothing is returned]
    """

    # We rename the root system as "r":
    r = root_system

    # We update the pool of labile C, taking into account the possible deficit in C previously recorded,
    # and all processes affecting the C pool EXCEPT THE SHOOT-ROOT EXCHANGE (#TODO: WATCH OUT!):
    r.root_labile_C_pool += (-r.root_labile_C_deficit + (
            # r.net_C_transfer_from_shoot_to_root_rate
            - r.C_rhizodeposition_rate
            - r.root_maintenance_respiration_rate
            - r.N2_fixation_rate_C_cost
            - r.net_root_C_immobilization_rate
            - r.actual_C_consumption_rate_by_root_growth) * time_step_in_seconds)
    # We reset the root labile C deficit to 0:
    r.root_labile_C_deficit = 0.

    # We update the pool of labile N, taking into account the possible deficit in N previously recorded,
    # and all processes affecting the N pool EXCEPT THE SHOOT-ROOT EXCHANGE (#TODO: WATCH OUT!):
    r.root_labile_N_pool += (- r.root_labile_N_deficit + (
        # - r.net_N_transfer_from_root_to_shoot_rate
        + r.N_uptake_rate + r.N2_fixation_rate
        - r.N_rhizodeposition_rate
        - r.net_root_N_immobilization_rate
        - r.actual_N_consumption_rate_by_root_growth) * time_step_in_seconds)
    # We reset the root labile N deficit to 0:
    r.root_labile_N_deficit = 0.

    # We update the pool of reserve C:
    r.root_labile_C_pool += r.net_root_C_immobilization_rate * time_step_in_seconds

    # We update the pool of reserve N:
    r.root_labile_N_pool += r.net_root_N_immobilization_rate * time_step_in_seconds

    # If the pools of C & N are negative, we reset them to 0 and adjust the deficits that will be taken into account
    # at the next step:
    if r.root_labile_C_pool < 0.:
        r.root_labile_C_deficit = - r.root_labile_C_pool
        r.root_labile_C_pool = 0.
    if r.root_labile_N_pool < 0.:
        r.root_labile_N_deficit = - r.root_labile_N_pool
        r.root_labile_N_pool = 0.

    return
#-----------------------------------------------------------------------------------------------------------------------