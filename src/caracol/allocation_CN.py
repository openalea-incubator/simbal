""" License Cecill-C

This module compute the allocation of C and N on a shoot and root model.
allocation_CN.py is part of caracol model

:TODO: Modify the documentation
A module that 
    (1) Check available N: store, acquired, remobilized
    (2) Distribute N among organs according to rules and hypothesis
        (H1) lux --> N --> dS
        (H2) dvt --> dS --> N
    (3) Check available C : store, acquired, remobilized
    (4) Distribute C among organs according to rules

#TODO : define the roadmap
* Separate parameters from functions (C/N) from script.
"""



###This script allows to calculate the allocation of C and N in a plant constituted of 3 class: 
# 1 leaf, 1 stem and 1 root.
# Each function is defined for a given class.

#import numpy #for loading the entire package
from math import (pi, exp) #for loading only pi function from package numpy

### Parameters:
k = 1 # speed of N uptake by roots in gN.cm-2.s-1
N_remob_rate = 1 # N remobilization rate in gN.gN-1
coef_N_leaf_area = 0.046 # leaf area = f(leaf N quantity) coefficient in cm2.g-1
coef_N_stem_area = 0.005 # stem area = f(stem N quantity) coefficient in cm2.g-1
alpha = 0.002 # increase of leaf area per unit of time in cm2.day-1
e_amax = 0.95 # maximum of radiation interception
k = -0.88 # extinction coefficient in Beer Lambert law
RUE = 1 # Radiation use efficiency in gC.MJ-1.cm-2
C_remob_rate = 1 # C remobilization rate in gC.gC-1

### Definition of architecture

class root(object):
    # we define initial properties of standard root as a cylinder:
    def __init__(self):
        self.radius = 10 # in cm
        self.length = 1 #in cm


class leaf(object):
    # leaf surface is provided by the architectural model:
    def __init__(self):
        self.area = 1  # in cm2


class stem(object):
    # we define initial properties of standard stem as a cylinder:
    def __init__(self):
        self.radius = 12  # in cm
        self.length = 7  # in cm

### Functions:
# Calculation of potential C and N demand for the root


def f_N_potential_root_demand(delta_volume=1, Ndensity=1):
    """
    Calculation of N potential root demand (in gN) in function of the root volume and root N density
    :param delta_volume: root volume in cm3
    :param Ndensity: N root density in gN.cm-3
    :return: N root potential demand in gN
    N root potential demand is printed
    """
    _N_potential_root_demand = delta_volume * Ndensity
    print("N potential root demand is", _N_potential_root_demand, "gN")
    return _N_potential_root_demand


def f_C_potential_root_demand(delta_volume=1, Cdensity=1):
    """
    Calculation of C potential root demand (in gC) in function of the root segment volume and root C density
    :param delta_volume: root volume in cm3
    :param Cdensity: C root density in gC.cm-3
    :return: C root potential demand in gC
    C root potential demand is printed
    """
    C_potential_root_demand = delta_volume * Cdensity
    print("C potential root demand is", C_potential_root_demand, "gC")
    return C_potential_root_demand

# Calculation of N offer for the plant


def f_N_available(element, Nsoil, Nstock, time_step, N_remob_rate=N_remob_rate):
    """
    Calculation of N available to the plant (in gN), in function of N uptake by roots (in g) and N remobilized
    from stock (in g)
    :param element: root object
    :param Nsoil: N quantity in soil in gN
    :param Nstock: N stock in the plant in gN
    :param time_step: time step of the simulation in day
    :return: N available to the plant in gN
    Root surface is printed
    Root N uptake is printed
    N remobilized is printed
    N available to the plant is printed
    """
     # calculation of root surface on which depends N uptake by roots:
    element.surface = 2 * pi * element.radius * element.length
    print("Root surface is", element.surface, "cm²")
    # calculation of root N uptake in function of soil N, root surface and speed of N uptake per unit of root surface
    N_root_uptake = k*Nsoil*element.surface*time_step
    print("N root uptake is", N_root_uptake, "gN")
    N_plant_remobilized = N_remob_rate * Nstock
    print("N remobilized is", N_plant_remobilized, "gN")
    N_available = N_root_uptake + N_plant_remobilized
    print("N available to the plant is", N_available, "gN")
    return N_available

# Calculation of potential N demand for shoot, leaf and stem

def f_N_shoot_potential_demand(N_available, N_potential_root_demand):
    """
    Calculation of the potential shoot demand (gN) in function of N available to the plant (gN)
    and potential root N demand (gN)
    :param N_available: N available to the plant in gN
    :param N_potential_root_demand: potential root N demand in gN
    :return: N shoot potential demand in gN
    N shoot potential demand is printed
    """
    N_shoot_potential_demand = N_available - N_potential_root_demand
    print("N potential shoot demand is", N_shoot_potential_demand, "gN")
    return N_shoot_potential_demand

def leaf_stem_surface_ratio(element1,element2):
    """
    Calculation of the surface ratio between leaf and stem (in cm².cm-²)
    :param element1: leaf object
    :param element2: stem object
    :return: surface ratio between leaf and stem in cm².cm-²
    Leaf and stem surface are printed
    """
    # calculation of leaf surface:
    print("Leaf surface is", element1.area, "cm²")
    # calculation of stem surface:
    element2.area = 2 * pi * element2.radius * element2.length
    print("Stem surface is", element2.area, "cm²")
    tau = element1.area / (element1.area + element2.area)
    return tau

def f_N_leaf_potential_demand(tau, N_shoot_potential_demand):
    """
    Calculation of leaf N potential demand (in gN) in function of the ratio of leaf/stem surface (in cm².cm-²)
    and potential shoot N demand (in gN)
    :param tau: leaf stem surface ratio in cm².cm-²
    :param N_shoot_potential_demand: potential shoot N demand in gN
    :return: leaf potential N demand in gN
    Leaf surface is printed
    Potential leaf N demand is printed
    """
    N_leaf_potential_demand = tau * N_shoot_potential_demand
    print("N potential leaf demand is", N_leaf_potential_demand, "gN")
    return N_leaf_potential_demand

def f_N_stem_potential_demand(tau, N_shoot_potential_demand):
    """
    Calculation of stem N potential demand (in gN) in function of the ratio of leaf/stem surface (in cm².cm-²)
    and potential shoot N demand (in gN)
    :param tau: leaf stem surface ratio in cm².cm-²
    :param N_shoot_potential_demand: potential shoot N demand in gN
    :return: stem potential N demand in gN
    Stem surface is printed
    Potential stem N demand is printed
    """
    N_stem_potential_demand=(1-tau)*N_shoot_potential_demand
    print("N potential stem demand is", N_stem_potential_demand, "gN")
    return N_stem_potential_demand

# Potential C demand of leaf, stem and shoot

def f_potential_leaf_area(N_leaf_potential_demand, coef_N_leaf_area):
    """
    Calculation of potential leaf area (cm2) in function of potential N leaf demand (gN)
    and the coefficient of the relationship of leaf area - leaf N quantity (cm2.g-1)
    :param N_leaf_potential_demand: leaf potential N demand in gN
    :param coef_N_leaf_area: coefficient between leaf area and leaf N quantity in cm2.g-1
    :return: potential leaf area in cm2
    Potential leaf area is printed
    """
    potential_leaf_area = coef_N_leaf_area * N_leaf_potential_demand
    print("potential leaf area is", potential_leaf_area, "cm2")
    return potential_leaf_area

def f_potential_stem_area(N_stem_potential_demand, coef_N_stem_area):
    """
    Calculation of potential stem area (cm2) in function of potential N stem demand (gN)
    and the coefficient of the relationship of stem area - stem N quantity (cm2.g-1)
    :param N_stem_potential_demand: stem potential N demand in gN
    :param coef_N_stem_area: coefficient between stem area and leaf N quantity in cm2.g-1
    :return: potential stem area in cm2
    Potential stem area is printed
    """
    potential_stem_area = coef_N_stem_area * N_stem_potential_demand
    print("potential stem area is", potential_stem_area, "cm2")
    return potential_stem_area

def f_C_potential_shoot_demand(potential_photosynthetic_area, Cdensity=1):
    """
    Calculation of C potential shoot demand (in gC) in function of the potential photosynthetic area and shoot N density
    :param potential_photosynthetic_area: potential photosynthetic area in cm2
    :param Cdensity: C shoot density in gC.cm-2
    :return: C root potential demand in gC
    C shoot potential demand is printed
    """
    C_potential_shoot_demand = potential_photosynthetic_area * Cdensity
    print("C potential shoot demand is", C_potential_shoot_demand, "gC")
    return C_potential_shoot_demand

# Calculation of C offer for the plant

def f_intercepted_PAR(initial_leaf_area, k, e_amax):
    """
    Calculation of PAR intercepted by the plant with Beer Lambert law
    :param initial_leaf_area: leaf area the previous day in cm2
    :param k: extinction coefficient
    :param e_amax: maximum intercepted radiation
    :return: intercepted PAR in MJ.cm-2
    Intercepted PAR is printed
    """
    PAR_intercepted = e_amax * exp(1 - (k * initial_leaf_area)) #exp is defined in math
    print("Intercepted PAR is", PAR_intercepted, "MJ.cm-2")
    return PAR_intercepted

def f_C_available(PAR_intercepted, RUE, Cstock, time_step):
    """
    Calculation of C available to the plant (in gC), in function of C uptake (in g) and C remobilized
    from stock (in g)
    :param PAR_intercepted: PAR intercepted by the plant
    :param RUE: Radiation use efficiency
    :param Cstock: C stock in the plant in gC
    :param time_step: time step of the simulation in day
    :return: C available to the plant in gC
    C uptake is printed
    C remobilized is printed
    C available to the plant is printed
    """
    # calculation of C uptake by the plant
    C_uptake = RUE * PAR_intercepted * time_step
    print("C uptake is", C_uptake, "gC")
    C_plant_remobilized = C_remob_rate * Cstock
    print("C remobilized is", C_plant_remobilized, "gC")
    C_available = C_uptake + C_plant_remobilized
    print("C available to the plant is", C_available, "gC")
    return C_available


############################################################################
#MAIN CODE:

# initializing

root_element=root()
leaf_element=leaf()
stem_element=stem()
#root_element1.volume=pi*(root_element1.radius)**2*root_element1.length
root_element.name="Jean_Paul"
#print(root_element1.radius)
#print(root_element1.volume)

#root_element2=root_segment()
#root_element2.volume=pi*(root_element2.radius)**2*root_element2.length
#root_element2.name="Horace"

# simulation
#Calculate the C and N potential demand of root (gC and gN)
#C_pot_demand, N_pot_demand=CN_potential_root_demand(delta_volume=2)
#print("The potential demand in C is", C_pot_demand, "gC.")

#list_root=[]
#list_root.append(root_element1)
#list_root_segments.append(root_element2)
#C_potential_root_total_demand(list_elements=list_root_segments)

starting_time=0
end_time=10
time_step=2

for t in range(starting_time,end_time,time_step):
    print("Time is", t)
    # Elements are growing:
    root_element.length = root_element.length+time_step/10
    print("Root length is", root_element.length,"cm")
    leaf_area0 = leaf_element.area
    print("Previous leaf area is", leaf_area0, "cm2")
    leaf_element.area = leaf_element.area + alpha * time_step
    print("Leaf area is", leaf_element.area, "cm2")
    stem_element.length = stem_element.length+time_step/10
    print("Stem length is", stem_element.length, "cm")
    root_element.volume = pi * (root_element.radius)**2 * root_element.length
    Root_Volume=root_element.volume
    #element.volume=v
    # Calculate potential C demand for total root
    N_potential_root_demand = f_N_potential_root_demand(delta_volume=Root_Volume, Ndensity=1)
    f_C_potential_root_demand(delta_volume=Root_Volume, Cdensity=1)
    # Calculate N available to shoot and root
    N_available=f_N_available(element=root_element, Nsoil=20, Nstock=5, time_step=time_step)
    # Calculate N potential demand of shoot, leaf and stem
    N_shoot_potential_demand = f_N_shoot_potential_demand(N_available=N_available, N_potential_root_demand=N_potential_root_demand)
    tau = leaf_stem_surface_ratio(element1=leaf_element, element2=stem_element)
    N_leaf_potential_demand = f_N_leaf_potential_demand(tau=tau, N_shoot_potential_demand=N_shoot_potential_demand)
    N_stem_potential_demand = f_N_stem_potential_demand(tau=tau, N_shoot_potential_demand=N_shoot_potential_demand)
    # Calculate potiential leaf and stem area
    potential_leaf_area = f_potential_leaf_area(N_leaf_potential_demand=N_leaf_potential_demand, coef_N_leaf_area=coef_N_leaf_area)
    potential_stem_area = f_potential_stem_area(N_stem_potential_demand=N_stem_potential_demand, coef_N_stem_area=coef_N_stem_area)
    # Calculate potential total photosynthetic area
    potential_photosynthetic_area = potential_leaf_area + potential_stem_area
    # Calculate C potential shoot demand
    C_potential_shoot_demand = f_C_potential_shoot_demand(potential_photosynthetic_area=potential_photosynthetic_area, Cdensity=1)
    # Calculate C available to the plant
    intercepted_PAR = f_intercepted_PAR(initial_leaf_area=leaf_area0, k=k, e_amax=e_amax)
    C_available = f_C_available(PAR_intercepted=intercepted_PAR, RUE=RUE, Cstock=3, time_step=time_step)






#N_seg_available=N_available(segment_length=1,segment_radius=1, Nsoil=10, N_remobilized=2, time_step=1)
#print(N_seg_available)

