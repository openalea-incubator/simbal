from simbal import CN_balance as cnb
import pytest

#Test function
def test_functions():
    assert cnb.leaf.define_leaf_senescence(0, 0) == 0
    assert cnb.leaf.define_leaf_necrotic_area(0, 0) == 0
    assert cnb.leaf.define_fall(0, 0) == 0
    assert cnb.plant.define_leaf_status(0, 0) == 0
    assert cnb.plant.define_C_pool_budget(0, 0) == 0


#Test function parametrize
@pytest.mark.parametrize(
    "x,a,b,expected", [(0, 0, 0, 0), (0, 0, 0, 0)]
    )
def test_functions(x, a, b, expected):
    assert cnb.leaf.define_leaf_senescence(x, a, b) == expected
    assert cnb.leaf.define_leaf_necrotic_area(x, a, b) == expected
    assert cnb.leaf.define_fall(x, a, b) == expected
    assert cnb.plant.define_leaf_status(x, a, b) == expected
    assert cnb.plant.define_C_pool_budget(x, a, b) == expected

