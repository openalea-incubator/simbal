#  -*- coding: utf-8 -*-

"""
    simbal.parameters
    ~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`simbal.parameters` defines the constant parameters of the model.

"""

# Setting the seed for generating the same random effects every time:
random_choice = 2

# Specific Carbon Assimilation (gC per micromol of photon):
SCA = 0.000001 # Arbitrary...

# Maximal rate of maintenance respiration (gC per gram of dry structural mass per second):
Vmax_maintenance = 5.55092978073827E-10 # Taken from RhizoDep 2024 (4e-08 mol of hexose per gDW per second)

# Michaelis_Menten constant for maintenance, i.e. the concentration of labile C for which the respiration rate is half
# of the maximal rate (gC per gram of dry structural mass):
Km_maintenance = 3.8578961976131E-06 # Taken from RhizoDep 2024 and Barillot et al. (2016) (2.78e-04 mol of hexose per gDW)

# Time constant corresponding to a first-order decay for N mobilization (second-1):
k_N_mobilization = 1/(24.*60.*60.) # We assume that the half life of the mobile N in the leaf is 1 day.

# Slope of the relationship of Leaf Area = f(N_labile) (gN per m2 of leaf):
SFQN = 10 # We assume that 1 gram of N leads to 10 m2 of leaves

# Phyllotactic angle, i.e. the angle formed by the directions of two consecutive leaves on the same stem (degrees):
phyllotactic_angle	= 137.5

# Angle formed by the directions of the petiole of a leaf and the direction of the supporting stem (degrees):
leaf_insertion_angle = 75

# Angle formed by the direction of a bud (possibly a future branch) and the direction of its supporting stem (degrees):
bud_insertion_angle	= 60

# Maximal area of a leaf (square meter):
LA_max = 0.1

# Minimal thickness of a leaf (meter):
leaf_minimal_thickness = 0.001

# Allometric coefficient defining the width of a leaf from its length (meter of width per meter of length):
leaf_width_to_length_ratio = 0.5

# Potential rate of internode elongation (meter per second):
elong_internode = 2.89351851851852E-07

# Growth yield as defined by Thornley and Cannell (2000), which corresponds to the fraction of C consumed by growth
# that is actually used for biomass production (gC per gC):
yield_growth = 0.8

# Density of the structural mass of shoot tissue (gram of root dry structural mass per cubic meter)
structural_mass_density = 0.1 * 1e6

# Carbon content of the dry structural mass (gC per gram of root dry structural mass):
structural_C_content = 0.44

# Nitrogen content of the dry structural mass (gN per gram of root dry structural mass):
structural_N_content = 0.01

# We define the phyllochrone, i.e. the period separating two leaf emergences (days):
phyllochrone_in_days = 5.

