from caracol import mtg_rapeseed
from caracol import photosynthesis as ph
from caracol import CN_plant_balance as cnp
from caracol import CN_leaf_balance as cnl
from caracol import fall as fll


class main:
    def _init_(self):
        mtg = mtg_rapeseed
        above = mtg.above()
        below = mtg.below() 

        g = above
        
        n = 10
    
    def main(_init_, g):
        for i in range(n=_init_.n):
            ph.leaf
            cnp.plant
            cnl.leaf
            fll.leaf

            return g

           


