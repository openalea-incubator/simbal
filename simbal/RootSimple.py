from math import sqrt, pi, floor, exp, isnan
from decimal import Decimal

from openalea.mtg import *
from openalea.mtg.traversal import pre_order, post_order

from rhizodep.model import initiate_mtg, summing_and_possibly_homogenizing, reinitializing_growth_variables, \
    update_distance_from_tip, update_surfaces_and_volumes, potential_growth, \
    volume_and_external_surface_from_radius_and_length, segmentation_and_primordia_formation, root_hairs_dynamics
import rhizodep.parameters as param

#-----------------------------------------------------------------------------------------------------------------------

# We define a class "root_system" that will be used to contain the properties of the root system:
class root_system(object):
    # We initiate the object "root_system":
    def __init__(self, g):
        # The object 'root_system' will in particular contain an attribute corresponding to a full MTG:
        self.MTG = g
        self.labile_C_pool = 0. # mol of C per root system
        self.labile_N_pool = 0. # mol of N per plant
        self.struct_mass = 0. # gDW
        self.surface = 0. # m2
        self.labile_C_input = 0. # mol of C per root system
        self.N_uptake = 0. # mol of N per root system

#-----------------------------------------------------------------------------------------------------------------------

def update_root_labile_C_with_input(r):

    # We update the labile pool:
    r.labile_C_pool += r.labile_C_input # mol of C per root system
    if r.labile_C_pool < 0.:
        print("WATCH OUT: after exchange with shoots, the labile C pool is negative - we reset it to 0!")
        r.labile_C_pool = 0.
    return

#-----------------------------------------------------------------------------------------------------------------------

def calculating_exchange_surface(r):

    g = r.MTG

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

    r.total_exchange_surface = total_exchange_surface

    return
#-----------------------------------------------------------------------------------------------------------------------

def process_N_uptake(r):

    g = r.MTG
    # We initialize the total N uptake of the whole root system:
    N_uptake = 0.

    # TODO: Change the following initialization with proper external parameters and element-specific soil N concentration
    Vmax_N_uptake = 2e-4 # mol of N per m2 of exchange surface
    Km_N_uptake = 1e-6 # mol of N per gDW of root structural mass
    N_soil = 1e-6 # mol of N per gDW of root structural mass

    # We cover all the vertices in the MTG:
    for vid in g.vertices_iter(scale=1):
        # n represents the vertex:
        n = g.node(vid)
        # If the root element can take up N:
        if n.length > 0. and N_soil > 0.:
            # Then we increase the total N uptake of the root system with the N uptake by the current element:
            N_uptake += n.exchange_surface * Vmax_N_uptake * N_soil / (Km_N_uptake + N_soil)

    if N_uptake < 0.:
        print("ERROR !!! The computed N uptake is negative!")
    else:
        r.N_uptake = N_uptake
    return

#-----------------------------------------------------------------------------------------------------------------------

def update_root_labile_N_with_input(r):

    # We update the labile pool:
    r.labile_N_pool += r.N_uptake # mol of N per root system
    if r.labile_N_pool < 0.:
        print("WATCH OUT: after N uptake from soil, the labile N pool is negative - we reset it to 0!")
        r.labile_N_pool = 0.
    return

#-----------------------------------------------------------------------------------------------------------------------

def process_rhizodeposition(r,
                            input_fraction=True,
                            C_rhizodeposition_fraction=0.20,
                            N_rhizodeposition_fraction=0.10):

    if input_fraction:
        # Option 1: rhizodeposition as a fraction of the input.
        r.rhizodeposited_C = C_rhizodeposition_fraction * r.labile_C_input # mol of C per root system
        r.rhizodeposited_N = N_rhizodeposition_fraction * r.N_uptake # mol of N per root system
    else:
        # Option 2: rhizodeposition as a fraction of the labile pool.
        r.rhizodeposited_C = C_rhizodeposition_fraction * r.labile_C_pool # mol of C per root system
        r.rhizodeposited_N = N_rhizodeposition_fraction * r.labile_N_pool # mol of N per root system

    # We update the labile pools:
    r.labile_C_pool -= r.rhizodeposited_C
    if r.labile_C_pool < 0.:
        print("WATCH OUT: rhizodeposition has taken too much C, the labile C pool has been set to 0!")
        r.labile_C_pool = 0.

    r.labile_N_pool -= r.rhizodeposited_N
    if r.labile_N_pool < 0.:
        print("WATCH OUT: rhizodeposition has taken too much N, the labile N pool has been set to 0!")
        r.labile_N_pool = 0.

    # We now distribute the rhizodeposition along the roots.
    # We cover all the vertices in the MTG:
    for vid in g.vertices_iter(scale=1):
        # n represents the vertex:
        n = g.node(vid)
        # If the root element can take up N:
        if n.length > 0.:
            relative_fraction = n.exchange_surface / r.total_exchange_surface
            n.rhizodeposited_C = relative_fraction * r.rhizodeposited_C
            n.rhizodeposited_N = relative_fraction * r.rhizodeposited_N
    
    return

#-----------------------------------------------------------------------------------------------------------------------

def process_maintenance(r, Vmax, Km):
    r.maintenance_respiration = Vmax * r.labile_C_pool / (Km + r.labile_C_pool) # mol of C per root system
    # We update the labile pool:
    r.labile_C_pool -= r.maintenance_respiration # mol of C per root system
    if r.labile_C_pool < 0.:
        print("WATCH OUT: maintenance has taken too much C, the labile C pool has been set to 0!")
        r.labile_C_pool = 0.
    return

#-----------------------------------------------------------------------------------------------------------------------

def process_exchange_with_reserve(r):

    r.C_input_to_reserve = 1 # mol of C per root system
    r.N_input_to_reserve = 1 # mol of N per root system
    r.C_output_from_reserve = 1 # mol of C per root system
    r.N_output_from_reserve = 1 # mol of N per root system

    # We update the labile pool:
    r.labile_C_pool += r.C_output_from_reserve - r.C_input_to_reserve
    if r.labile_C_pool < 0.:
        print("WATCH OUT: reserve has taken too much C, the labile C pool has been set to 0!")
        r.labile_C_pool = 0.

    r.labile_N_pool += r.N_output_from_reserve - r.N_input_to_reserve
    if r.labile_N_pool < 0.:
        print("WATCH OUT: reserve has taken too much N, the labile N pool has been set to 0!")
        r.labile_N_pool = 0.

    return

#-----------------------------------------------------------------------------------------------------------------------

# Function calculating a satisfaction coefficient for the growth of the whole root system:
# -----------------------------------------------------------------------------------------
def satisfaction_coefficients(g,labile_C_pool, labile_N_pool):

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
        satisf_for_C =1.

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
def ArchiSimple_growth_with_C_and_N(r,
                                    satisf_for_C, satisf_for_N):
    """
    This function computes the growth of the root system according to ArchiSimple's rules, but limits the growth
    by the amount of available labile C and N.
    """

    # We initialize the root biomass that will be produced from the labile pool of C & N:
    total_mass_produced = 0.
    # We select the MTG included within the root system r:
    g = r.MTG
    # We have to cover each vertex from the apices up to the base one time:
    root_gen = g.component_roots_at_scale_iter(g.root, scale=1)
    root = next(root_gen)

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
        if n.potential_radius > n.initial_radius:

            # We then calculate the corresponding new possible radius corresponding to this maximum volume:
            if n.type == "Root_nodule":
                # If the element corresponds to a nodule, then it we calculate the radius of a theoretical sphere:
                possible_radius = ((final_volume - volume_after_elongation) * 3. / (4. * pi)) ** (1. / 3.)
            else:
                # Otherwise, we calculate the radius of a cylinder:
                possible_radius = sqrt((final_volume - volume_after_elongation) / (n.length * pi))
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

    # EVENTUALLY, WE UPDATE THE TOTAL MASS OF THE ROOT SYSTEM AND THE LABILE C AND N POOLS:
    r.struct_mass += total_mass_produced
    r.labile_C_pool -= total_mass_produced * param.struct_mass_C_content / param.yield_growth
    r.labile_N_pool -= total_mass_produced * param.struct_mass_N_content

    # if r.labile_C_pool < 1e-16 or r.labile_N_pool < 1e-16:
    #     print("ERROR !!! After growth, labile C pool is ", r.labile_C_pool, "and labile N pool is", r.labile_N_pool)
    return

########################################################################################################################
########################################################################################################################

# MAIN SIMULATION:
##################

# We define the MTG g:
g = initiate_mtg()
reinitializing_growth_variables(g)
# We update the distance from tip for each root element in each root axis:
update_distance_from_tip(g)
# We update the surfaces and volumes of each root element:
update_surfaces_and_volumes(g)

# And we create a root system r from it:
r = root_system(g)
# We adjust the properties of the whole root system according to g:
g_properties = summing_and_possibly_homogenizing(g)
r.struct_mass = g_properties["total_living_root_struct_mass"]
calculating_exchange_surface(r)

time_step_in_seconds = 3600.
n_steps = 24*10
soil_temperature=20

# We perform the simulation:
for time_step in range(1, n_steps+1):

    print("For step", time_step,":")

    # Setting the new input of C to the root system:
    r.labile_C_input = 0.002
    # Increment of labile C due to C input from shoots:
    update_root_labile_C_with_input(r)

    # Setting N uptake by the root system:
    process_N_uptake(r)
    # Increment of labile N due to N inputs:
    update_root_labile_N_with_input(r)

    # Rhizodeposition and update of labile C:
    process_rhizodeposition(r)

    # Maintenance respiration and update of labile C:
    process_maintenance(r, Vmax=1, Km=1)

    # Exchange with reserve and update of labile C:
    process_exchange_with_reserve(r)

    # Growth:
    #--------
    # We reset all growth-associated variabes before performing the new growth:
    reinitializing_growth_variables(g)
    # We calculate the potential growth of each root element according to ArchiSimple:
    potential_growth(g, time_step_in_seconds=time_step_in_seconds,
                     ArchiSimple=True,
                     soil_temperature_in_Celsius=soil_temperature)
    # We calculate the "satisfaction coefficients" for growth, depending on the amount of labile C and N left:
    satisf_for_C, satisf_for_N = satisfaction_coefficients(g,
                                                           labile_C_pool = r.labile_C_pool,
                                                           labile_N_pool = r.labile_N_pool)
    print("    The satisfaction coefficient for C is", "{:.2f}".format(Decimal(satisf_for_C)),
          "and the one for N is", "{:.2f}".format(Decimal(satisf_for_N)))
    # We use the function ArchiSimple_growth to adapt the potential growth to the available labile C and N,
    # and we update the labile pools of C and N:
    ArchiSimple_growth_with_C_and_N(r, satisf_for_C, satisf_for_N)
    # We proceed to the segmentation of the whole root system
    segmentation_and_primordia_formation(g, time_step_in_seconds,
                                         soil_temperature_in_Celsius=soil_temperature, random=True,
                                         ArchiSimple=True,
                                         root_order_limitation=False)

    #TODO: Include root hairs dynamics once the cost of C & N is properly registered!

    # # We modifiy root hairs characteristics according to their specific dynamics:
    # root_hairs_dynamics(g, time_step_in_seconds=time_step_in_seconds,
    #                     soil_temperature_in_Celsius=soil_temperature)

    # We update the distance from tip for each root element in each root axis:
    update_distance_from_tip(g)
    # We update the surfaces and volumes of each root element:
    update_surfaces_and_volumes(g)
    calculating_exchange_surface(r)

    print("    Labile C pool in the root system is now", "{:.2E}".format(Decimal(r.labile_C_pool)),
          "mol of C and labile N pool is ", "{:.2E}".format(Decimal(r.labile_N_pool)), "mol of N.")
    print("    The root system contains", g.nb_vertices(), "different root elements.")
    print("")
