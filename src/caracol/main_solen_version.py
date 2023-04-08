from caracol import mtg_rapeseed
from caracol import photosynthesis as ph
from caracol import CN_plant_balance as cnp
from caracol import CN_leaf_balance as cnl
from caracol import fall as fll


class main:
    def _init_(self):
        mtg = mtg_rapeseed
        above = mtg.above() 

        g = above
        
        n = 10
    
    def main(_init_):
        for i in range(n=_init_.n, g =_init_.g ):
            PAR_total,N_pool, N_potential_shoot_offer = ph.leaf.define_N_potential_shoot_offer
            plant_N_by_lightunit = ph.leaf.define_leaf_N_potential_offer
            g = ph.leaf.define_leaf_area
            g = ph.leaf.define_leaf_photosynthesis
            g = ph.leaf.define_potential_leaf_area
            g = ph.leaf.define_potential_leaf_light_interception
            g = ph.leaf.define_leaf_potential_photosynthesis
            g = ph.leaf.define_FSPM_Allirand_law
            g = ph.leaf.define_C_potential_leaf_demand
            g = ph.leaf.define_C_potential_leaf_offer
            
            g = cnp.plant.define_leaf_status
            g = cnp.plant.define_C_pool_budget

            g = cnl.leaf.define_leaf_senescence
            g = cnl.leaf.define_leaf_necrotic_area

            g = fll.leaf.define_fall

            return g, PAR_total,N_pool, N_potential_shoot_offer, PAR_total, plant_N_by_lightunit

           


