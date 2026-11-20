import sys

import openalea.simbal.simbal_parameters as param
from openalea.simbal.simbal_orchestrator import Wholeplant, computing_shoot_root_exchange_rates

EPS = sys.float_info.epsilon

def test_computing_shoot_root_exchange_rates():

    plant = Wholeplant()
    param.C_diffusion_coefficient = 2.5e-6
    plant.shoot_labile_C_pool = 0.1
    plant.shoot_labile_C_deficit = 0.05
    plant.root_labile_C_pool = 0.5
    plant.root_labile_C_deficit = 0.01
    plant.shoot_struct_mass = 0.01
    plant.exchange_section = 0.002
    plant.root_struct_mass = 0.05
    computing_shoot_root_exchange_rates(plant, computing_C_exchange=True)
    assert (plant.net_C_transfer_from_shoot_to_root_rate - (-2.4e-8)) <= EPS

    param.N_diffusion_coefficient = 2.5e-6
    plant.shoot_labile_N_pool = 0.1
    plant.shoot_labile_N_deficit = 0.05
    plant.root_labile_N_pool = 0.5
    plant.root_labile_N_deficit = 0.01
    plant.shoot_struct_mass = 0.01
    plant.exchange_section = 0.002
    plant.root_struct_mass = 0.05
    computing_shoot_root_exchange_rates(plant, computing_N_exchange=True)
    assert (plant.net_N_transfer_from_root_to_shoot_rate - 2.4e-8) <= EPS
