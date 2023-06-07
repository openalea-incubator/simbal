from . import mtg_rapeseed
from . import photosynthesis as ph
from . import aptos as cnb


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
                PAR_total,N_pool, N_potential_shoot_offer = leaf.define_N_potential_shoot_offer(g,...)
                plant_N_by_lightunit = ph.leaf.define_leaf_N_potential_offer
                g = ph.leaf.define_leaf_area
                g = ph.leaf.define_leaf_photosynthesis
                g = ph.leaf.define_potential_leaf_area
                g = ph.leaf.define_potential_leaf_light_interception
                g = ph.leaf.define_leaf_potential_photosynthesis
                g = ph.leaf.define_FSPM_Allirand_law
                g = ph.leaf.define_C_potential_leaf_demand
                g = ph.leaf.define_C_potential_leaf_offer
            
                g = cnb.plant.define_leaf_status
                g = cnb.plant.define_C_pool_budget

                g = cnb.leaf.define_leaf_senescence
                g = cnb.leaf.define_leaf_necrotic_area

                g = cnb.leaf.define_fall

                return g, PAR_total,N_pool, N_potential_shoot_offer, 
                    PAR_total, plant_N_by_lightunit


if __name__ == "__main__" :
    """Run the unit test suite 

    """
    import doctest
    doctest.testmod()

