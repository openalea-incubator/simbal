from caracol import mtg_senescence
from caracol import fall as snc

scenario = mtg_senescence.PlantScenario
g = above = scenario.init_above_mtg()

n = 10
N_potential_plant_offer = 1
N_potential_root_demand = 1
RUE = 1
e_amax = 1
k = 1



for i in range (n):

    PAR_total,N_pool, N_potential_shoot_offer = N_potential_shoot_offer(N_potential_plant_offer, N_potential_root_demand)
    g = snc.leaf_area(g)
    g = snc.photosynthetic_leaf_area(g)
    g = snc.potential_leaf_PAR(g, k, e_amax)
    g = snc.leaf_photosynthesis(g, RUE)
    g, plant_N_by_Lightunit = snc.leaf_N_potential_offer(g, PAR_total, N_potential_shoot_offer)
    g = snc.potential_leaf_surface(g)
    g = snc.potential_leaf_light_interception(g, k, e_amax)
    g = snc.leaf_potential_photosynthesis(g, RUE)
    g = snc.FSPM_Allirand_law(g)
    g = snc.C_potential_leaf_demand(g)
    g = snc.C_potential_leaf_offer(g)

    g = snc.leaf_status(g)
    g, pool_C, PAR_total, C_plant_hunger = snc.C_pool_budget(g)
    g, plant_area, plant_N, pool_N, N_plant_demand, N_plant_deficit, N_plant_hunger, N_plant_remobilized = snc.leaf_senescence(g, N_potential_shoot_offer)
    g = snc.leaf_necrotic_area(g)

    g, plant_area, plant_N, pool_N, plant_N_by_lightunit, N_plant_demand, N_plant_deficit, N_plant_hunger, N_plant_remobilized, senescence_initiation, leaf_fall = snc.turnover_plant_budget(g)






