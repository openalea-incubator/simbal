from openalea.mtg import *
from openalea.mtg.turtle import *

from numpy import *
import numpy as np
from pandas import *

import mtg_plant as mtg_plant

#######################################

def update_MTG(g, outputs):
    """
    Update the MTG 
    :param dict outputs: {'param1': { vid1: , vid2, ...}, 'param2': { vid1: , vid2, ...}}
    """
     # add the missing property
    for param in outputs.keys():
        if param not in g.properties():
            g.add_property(param)
        # update the MTG
        g.property(param).update(outputs[param])

###########################################

def vegetative_rapeseed(
    n_nodes: int = 10,
    node_length: float= 0.04,
    ):
    """ Build a plant with n_nodes internodes
    """

    above = g = MTG()
    vid = g.add_component(g.root, label='Plant')
    vid = g.add_component(vid, label='Axis')

    pid = g.add_component(vid, label='Phytomer')
    vid = g.add_component(pid, label='Internode')
    g.node(vid).NodeLength = node_length
    #g.node(vid).short = True
    leaf_id = g.add_child(vid, edge_type='+', label='Leaf')
    
    for i in range(n_nodes-1):
        pid = g.add_child(pid, edge_type='<', label='Phytomer')
        vid, pid = g.add_child_and_complex(vid, complex=pid, edge_type='<', label='Internode')
        g.node(vid).NodeLength = n_nodes
        #g.node(vid).short = True

        leaf_id = g.add_child(vid, edge_type='+', label='Leaf')

    return mtg_plant.Plant(name='rapeseed', above = above)

def vegetative_rapeseed_DJ(
    DJ: int = 950,
    node_length: float= 0.02,
    phyllochrone = 0.015

    ):
    """ Build an mtg of a single plant. The number of internodes is based on its age (DD) and phyllochrone 
    """
    n_nodes=round(DJ*phyllochrone)
    above = g = MTG()
    vid = g.add_component(g.root, label='Plant')
    vid = g.add_component(vid, label='Axis')

    pid = g.add_component(vid, label='Phytomer')
    vid = g.add_component(pid, label='Internode')
    g.node(vid).NodeLength = node_length
    #g.node(vid).short = True
    leaf_id = g.add_child(vid, edge_type='+', label='Leaf')
    
    for i in range(n_nodes-1):
        pid = g.add_child(pid, edge_type='<', label='Phytomer')
        vid, pid = g.add_child_and_complex(vid, complex=pid, edge_type='<', label='Internode')
        g.node(vid).NodeLength = node_length
        #g.node(vid).short = True

        leaf_id = g.add_child(vid, edge_type='+', label='Leaf')

    return mtg_plant.Plant(name='rapeseed', above = above)


##Mtg builder with mutliple plants
# functions from Adel used in vegetative_rapeseed_field()
def adel_regular(nb_plants, nb_rank, dx, dy, nx=None):
    
    if nx is None:
        nx = int(nb_plants/nb_rank)
    ny = nb_rank
    domain = ((0,0),(nx*dx, ny*dy))
    return [(i*dx+dx/2., j*dy+dy/2., 0.) for j in range(ny) for i in range(nx)], domain

def agronomicplot(nplants=None, length=None, width=None, sowing_density=None, 
                  plant_density=None, inter_row=None, noise=0, convunit=100, center_scene=True):
    """ Returns the number of plants, the positions, the domain (scene units), the domain area (square meter) and the conversion coefficient for meter to scene unit (1/convunit) of a micro-plot specified with agronomical variables
    nplants is the number of plants of the micro-plot (mandatory if length or width is None)
    length (m) is plot dimension along row direction (mandatory if nplants is None)
    width (m) is plot dimension perpendicular to row direction (mandatory if nplants is None)
    sowing density is the density of seeds sawn (mandatory)
    plant_density is the density of plants that are present (after loss due to bad emergence, early death...) (mandatory)
    inter_row (m) is for the  distance between rows (mandatory)
    noise (%), indicates the precision of the sowing for the inter plant spacing
    convunit is the conversion factor from meter to scene unit
    center_scene allows to center the position around origin. If False, the scene is in the x+,y+ sector, the origin being at the lower left corner of the domain
    
    Rows are parallel to x-axis
    Length and Width are adjusted to produce a canopy centered in its domain and compliant with infinitisation
    """
    from operator import itemgetter
    import random
    
    
    inter_plant = 1. / inter_row / sowing_density
    
    if nplants is None: # length and width are NOT None
        nrow = max(1, int(float(width) / inter_row))
        plant_per_row = max(1, int(float(length) / inter_plant))
        nplants = nrow * plant_per_row
        
    else: # nplants is NOT None ; length and width are None 
        import math
        nrow = max(1, int(math.sqrt(nplants)))
        
    dx = inter_plant * convunit
    dy = inter_row * convunit
    positions, domain = adel_regular(nplants, nrow, dx, dy)
    n_emerged = int(nplants * plant_density / sowing_density)
    positions = random.sample(positions, n_emerged)
    # sorting by ranks
    positions = sorted(positions, key= itemgetter(1,0))
    domain_area = abs(domain[1][0] - domain[0][0]) / convunit * abs(domain[1][1] - domain[0][1]) / convunit
    if center_scene:
        xc = float(domain[1][0] + domain[0][0]) / 2
        yc = float(domain[1][1] + domain[0][1]) / 2
        positions = [(x - xc, y - yc, z) for x, y, z in positions]
        domain = ((domain[0][0] - xc, domain[0][1] - yc), (domain[1][0] - xc, domain[1][1] - yc))
    
    return n_emerged, positions, domain, domain_area, 1. / convunit

def vegetative_rapeseed_field(
    DJ: int = 950,
    node_length: float= 0.02,
    phyllochrone = 0.015,
    coord = [(0,0,0)]
    
    ):
    """ Build a field with nrows(coord) plants that have n_nodes internodes

    """
    n_nodes=round(DJ*phyllochrone)
    #print('plantes à',n_nodes,'entrenoeuds')
    above = g = MTG()

    for pl in coord:
       
        #print(pl)
        vid = g.add_component(g.root, label='Plant',edge_type='/',position=pl) #vid = vertex_id
        vid = g.add_component(vid, label='Axis')

        pid = g.add_component(vid, label='Phytomer') #pid= phytomer_id
        vid = g.add_component(pid, label='Internode') 
        g.node(vid).NodeLength = node_length
        #g.node(vid).short = True
        leaf_id = g.add_child(vid, edge_type='+', label='Leaf')
        
        for i in range(n_nodes-1):
            pid = g.add_child(pid, edge_type='<', label='Phytomer')
            vid, pid = g.add_child_and_complex(vid, complex=pid, edge_type='<', label='Internode')
            g.node(vid).NodeLength = node_length
            #g.node(vid).short = True

            leaf_id = g.add_child(vid, edge_type='+', label='Leaf')

    return mtg_plant.Plant(name='rapeseed', above = above)   

def bell_shaped_dist(total_area=1, nb_phy=15, rmax=.7, skewness=5):
    """ OLD, use generate_normalized_bell_curve_points() instead
    returns the leaf area of individual leaves along a bell shaped model, as a proportion of the total leaf area """

    r = np.linspace(1./nb_phy, 1, nb_phy)
    k = skewness
    relative_surface = np.exp(-k / rmax * ( 2 * (r - rmax)**2 + (r - rmax)**3))
    leaf_area = relative_surface / relative_surface.sum() * total_area
    return leaf_area.tolist()

def generate_normalized_bell_curve_points(mean, std_dev, amplitude, N):
    """ returns the leaf area of individual leaves along a bell shaped model, as a proportion of the total leaf area """
    x_values = np.linspace(mean - 3 * std_dev, mean + 3 * std_dev, N)
    y_values = amplitude * np.exp(-0.5 * ((x_values - mean) / std_dev)**2)
    
    # Normalize y values to ensure the sum is 1
    normalized_y_values = y_values / np.sum(y_values)
    
    points = list(normalized_y_values)
    return points

def phenotype(g,
    leaf_max: float=0.6587057,
    skew: float=5,  #
    sd: float= 0.246378, 
    amp: float=0.1188,
    phyllot : float = 137.5,
    ins_angle : float = 60,
    leaf_angle : float = 80, #à mesurer expérimentallement 
    total_surface: float=300, 
    coeff_width: float =0.831562, 
    coeff_petiole: float =0.51017, #limb length to petiole length , 0.31017 length/whole leag length
    stem_d: float=0.1,
    coeff_petiole_d : float=0.5
    ): 
    """
    applies phenotypical attributes to an MTG based on parameters. 
        leaf_max : [float between 0 and 1] Expectation of the bell curve => Position of the biggest leaf along the stem.
        sd: [float] Standard deviation of the bell curve => repartition of the leaf surface
        amp: [float] amplitude of the bell curve 
        phyllot: [float, degrees] phyllotaxy angle 
        ins_angle: [float, degrees] insertion angle of the petiole
        leaf_angle: [float, degrees] angle of the flat of the leaf
        total_surface: [int, cm2] surface to be shared following the bell proportion
        coeff_width: [float] Leaf Width / Leaf Length ratio (limb only)
        coeff_petiole: [float] Petiole Length / Leaf Length 
        stem_d: [float, cm] stem diameter
        coeff_petiole_d : [float] petiole diameter / stem diameter ration 
    """
    sortie=dict()
    plants=[k for k, v in g.properties()['label'].items()  if v=='Plant']
    leaves=[k for k, v in g.properties()['label'].items()  if v=='Leaf']
    internodes=[k for k, v in g.properties()['label'].items()  if v=='Internode']
    
    #surface=bell_shaped_dist(1,int(len(leaves)/len(plants)),leaf_max,skew) Old functio
    surface=generate_normalized_bell_curve_points(mean=leaf_max,std_dev=sd,amplitude=amp,N=int(len(leaves)/len(plants)))
    sortie['SurfaceRepartition']=dict(zip(leaves,surface*len(leaves)))
    sortie['LeafSurface']=dict(zip(leaves,[element * total_surface for element in surface]*len(leaves)))


    sortie['LeafLength']= dict(zip(leaves,
             [2*(element/(coeff_width*pi))**(0.5) for element in sortie['LeafSurface'].values()]*len(leaves)))

    sortie['LeafWidth']=dict(zip(leaves,
                     [element *coeff_width for element in sortie['LeafLength'].values()]*len(leaves)))
    
    sortie['SurfTheo']=dict(zip(leaves,[pi*(element/2)*(coeff_width*element/2) for element in sortie['LeafLength'].values()]*len(leaves)))

    sortie['PetioleLength']=dict(zip(leaves,[element * coeff_petiole for element in sortie['LeafLength'].values()]*len(leaves)))
    
    sortie['InsertionAngle']=dict(zip(leaves,[ins_angle]*len(leaves)))
    sortie['LeafAngle']=dict(zip(leaves,[leaf_angle]*len(leaves)))
    sortie['Phyllotaxy']=dict(zip(leaves,[phyllot]*len(leaves)))

    sortie['StemDiam']=dict(zip(internodes,[stem_d]*len(internodes)))
    sortie['PetioleDiam']=dict(zip(leaves,[stem_d*coeff_petiole_d]*len(leaves)))

    update_MTG(g,sortie)

def metabolism(g,
    fall : bool = False,
    growth : bool = True,
    photosynthesis = 0.1,
    par = 100,
    respiration = 0.1,
    potential_photosynthesis= 0.1,
    potential_area= 0.1,
    nitrogen_source : bool = True,
    nitrogen_sink : bool = True,
    nitrogen_content = 10,
    nitrogen_demand = 0.1,
    nitrogen_storage= 0.1,
    nitrogen_import= 0.1,
    nitrogen_export= 0.1,
    rue = 0.1,
    weight_import= 0.1,
    weight_export= 0.1,
    carbon_content= 0.1,
    carbon_demand = 0.1,
    carbon_deficit = 0.1,
    carbon_import= 0.1,
    carbon_export= 0.1,
    carbon_spend= 0.1,
    carbon_growth = 0.1,
    carbon_storage = 0.1,
    carbon_nitrogen_ratio= 0.1,
    nitrogen_fall_threshold_empirical= 0.1,
    PAR_fall_threshold_empirical= 0.1,
    carbon_fall_threshold_mecanism= 0.1,
    nitrogen_fall_threshold_mecanism = 0.1): 
    """
    """
    sortie=dict()

    leaves=[k for k, v in g.properties()['label'].items()  if v=='Leaf']

    sortie['Growth']=dict(zip(leaves,[growth]*len(leaves)))

    sortie['PAR']=dict(zip(leaves,[par]*len(leaves)))
    sortie['Photosynthesis']=dict(zip(leaves,[photosynthesis]*len(leaves)))
    sortie['Respiration']=dict(zip(leaves,[respiration]*len(leaves)))

    sortie['PotentialPhotosynthesis']=dict(zip(leaves,[potential_photosynthesis]*len(leaves)))
    sortie['PotentialArea']=dict(zip(leaves, [potential_area]*len(leaves)))

    sortie['NitrogenSource']=dict(zip(leaves,[nitrogen_source]*len(leaves)))
    sortie['NitrogenSink']=dict(zip(leaves,[nitrogen_sink]*len(leaves)))
    
    sortie['NitrogenContent']=dict(zip(leaves,[nitrogen_content]*len(leaves)))
    sortie['NitrogenDemand']=dict(zip(leaves,[nitrogen_demand]*len(leaves)))
    sortie['NitrogenImport']=dict(zip(leaves,[nitrogen_import]*len(leaves)))
    sortie['NitrogenExport']=dict(zip(leaves,[nitrogen_export]*len(leaves)))
    sortie['NitrogenStorage']=dict(zip(leaves,[nitrogen_storage]*len(leaves)))
    sortie['RUE']=dict(zip(leaves,[rue]*len(leaves)))


    sortie['PonderationExport']=dict(zip(leaves,[weight_export]*len(leaves)))
    sortie['PonderationImport']=dict(zip(leaves,[weight_import]*len(leaves)))
    
    sortie['CarbonContent']=dict(zip(leaves,[carbon_content]*len(leaves)))
    sortie['CarbonDemand']=dict(zip(leaves,[carbon_demand]*len(leaves)))
    sortie['CarbonDeficit']=dict(zip(leaves,[carbon_deficit]*len(leaves)))
    sortie['CarbonImport']=dict(zip(leaves,[carbon_import]*len(leaves)))
    sortie['CarbonExport']=dict(zip(leaves,[carbon_export]*len(leaves)))
    sortie['CarbonSpend']=dict(zip(leaves,[carbon_spend]*len(leaves)))
    sortie['CarbonGrowth']=dict(zip(leaves,[carbon_growth]*len(leaves)))
    sortie['CarbonStorage']=dict(zip(leaves,[carbon_storage]*len(leaves)))
    sortie['CarbonNitrogenRatio']=dict(zip(leaves,[carbon_nitrogen_ratio]*len(leaves)))

    sortie['Fall']= dict(zip(leaves,[fall]*len(leaves)))
    sortie['NirogenFallThresholdEmpirical']=dict(zip(leaves,[nitrogen_fall_threshold_empirical]*len(leaves)))
    sortie['NirogenFallThresholdMecanism']=dict(zip(leaves,[nitrogen_fall_threshold_mecanism]*len(leaves)))
    sortie['CarbonFallThresholdMecanism']=dict(zip(leaves,[carbon_fall_threshold_mecanism]*len(leaves)))
    sortie['PARFallThresholdEmpirical']=dict(zip(leaves,[PAR_fall_threshold_empirical]*len(leaves)))

    update_MTG(g,sortie)

