from openalea.mtg import *


class Plant:
    """ Define a plant as an above and below ground MTG.

    Some usefull functions can be defined as proxies.

    Authors
    -------
    Arman GRUMEL and Laurène PERTHAME
        Colzette model

    """
    def __init__(self, name='Rapessed', above=None, below=None):
        self.name = name
        self.above = self.init_above_mtg() if above is None else above
        self.below = self.init_below_mtg() if below is None else below

    @staticmethod
    def init_above_mtg():
        """ Simple function to create a simple aerial MTG

        Authors
        -------
        Arman GRUMEL and Laurène PERTHAME
            Colzette model

        """
        g = MTG()
        
        # Add a plant (scale is 1)
        vid = g.add_component(g.root, label='Plant')
        vid = g.add_component(vid, label='Axis')
        vid = g.add_component(vid, label='Phytomer')
        
        # Last scale
        vid = g.add_component(vid, label='Internode')
        lid = g.add_child(vid, edge_type='+', label='Leaf')
        bid = g.add_component(vid, edge_type='<', label='Bud')

        fat_mtg(g)

        return g

    @staticmethod
    def init_below_mtg():
        """ Simple function to create a simple belowground MTG or RootSystem Architecture

        Authors
        -------
        Arman GRUMEL and Laurène PERTHAME
            Colzette model

        """
        g = MTG()
        
        # Add a plant (scale is 1)
        vid = g.add_component(g.root, label='RSA')
        vid = g.add_component(vid, label='Root')
        
        # Last scale
        vid = g.add_component(vid, label='Segment')
    
        fat_mtg(g)

        return g

def vegetative(
        n_nodes,
        node_length,
        name = 'rapeseed'):
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
        g.node(vid).NodeLength = node_length
        #g.node(vid).short = True

        leaf_id = g.add_child(vid, edge_type='+', label='Leaf')

    return Plant(name=name, above = above)

def vegetative_DJ(
    DJ: int = 950,
    node_length: float= 0.02,
    phyllochrone = 0.015,
    name = 'Rapeseed'

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

    return Plant(name, above = above)


##Mtg builder with mutliple plants
# functions from Adel used in vegetative_rapeseed_field()
def adel_regular(nb_plants, nb_rank, dx, dy, nx=None):
    
    if nx is None:
        nx = int(nb_plants/nb_rank)
    ny = nb_rank
    domain = ((0,0),(nx*dx, ny*dy))
    return [(i*dx+dx/2., j*dy+dy/2., 0.) for j in range(ny) for i in range(nx)], domain

def agronomicplot(nplants, 
                  length, 
                  width, 
                  sowing_density, 
                  plant_density, 
                  inter_row,
                  nrow, 
                  noise=0, 
                  convunit=100, 
                  center_scene=True
                  ):
    """ Returns the number of plants, the positions, the domain (scene units), the domain area (square meter) and the conversion coefficient for meter to scene unit (1/convunit) of a micro-plot specified with agronomical variables
    nplants is the number of plants of the micro-plot (mandatory if length or width is None)
    length (m) is plot dimension along row direction (mandatory if nplants is None)
    width (m) is plot dimension perpendicular to row direction (mandatory if nplants is None)
    sowing density is the density of seeds sawn (mandatory) #surfacic density ?
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

    sowing_density = plant_density

    inter_plant = 1. / inter_row / sowing_density
    
    if nplants is None: # length and width are NOT None
        nrow = max(1, int(float(width) / inter_row))
        plant_per_row = max(1, int(float(length) / inter_plant))
        nplants = nrow * plant_per_row
        
    else: # nplants is NOT None ; length and width are None 
        #import math
        #nrow = max(1, int(math.sqrt(nplants)))
        nrow = nrow
        
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

def vegetative_field(
    DJ,
    node_length,
    phyllochrone,
    coord,
    name):
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

    return Plant(name=name, above = above)   

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