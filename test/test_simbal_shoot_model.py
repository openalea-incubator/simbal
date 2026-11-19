from openalea.simbal import simbal_shoot_model

def test_simple():
    shoot_g = simbal_shoot_model.create_shoot_MTG(number_of_phytomers=2,
                                      initial_internode_height=0.0001,
                                      initial_internode_radius=0.0003,
                                      initial_leaf_length=0.003,
                                      initial_labile_C_concentration=0,
                                      initial_labile_N_concentration=0
                                      )

    d = shoot_g.property('potential_radius')
    assert next(iter(d.values())) == 0.0003