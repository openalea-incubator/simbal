from . import mtg_rapeseed
from . import photosynthesis as ph
from . import aptos as apt


class Main :
    def _init_(self) :
        self.mtg = mtg_rapeseed.Rapseed()
        above = self.mtg.above() 
        self.g = above
    
    def run(self, n=10) :
        g = self.g
        leaves = [ph.Leaf(vid) for vid in g if g.label(vid)=='Leaf']

        for i in range(n) :
            for leaf in leaves:
                PAR_total,N_pool, N_potential_shoot_offer = leaf.N_potential_shoot_offer(g,...)
                plant_N_by_lightunit = ph.leaf.leaf_N_potential_offer
                g = ph.leaf.leaf_area
                g = ph.leaf.leaf_photosynthesis
                g = ph.leaf.potential_leaf_area
                g = ph.leaf.potential_leaf_light_interception
                g = ph.leaf.leaf_potential_photosynthesis
                g = ph.leaf.FSPM_Allirand_law
                g = ph.leaf.C_potential_leaf_demand
                g = ph.leaf.C_potential_leaf_offer
            
                g = apt.plant.leaf_status
                g = apt.plant.C_pool_budget

                g = apt.leaf.leaf_senescence
                g = apt.leaf.leaf_necrotic_area

                g = apt.leaf.fall

                return g, PAR_total,N_pool, N_potential_shoot_offer, PAR_total, plant_N_by_lightunit


if __name__ == "__main__" :
    """Run the unit test suite 

    """
    import doctest
    doctest.testmod()

