from simbal import plant_mtg as plt

def test1():
    rapeseed = plt.Plant(name='Rapeseed')
    above = rapeseed.above
    below = rapeseed.below

    assert len(above) == 7
    assert len(below) == 4

    