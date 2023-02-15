from caracol import mtg_senescence
from caracol import senescence as snc

scenario = mtg_senescence.PlantScenario
g = above = scenario.init_above_mtg()

max_scale = g.max_scale()
# Initial conditions
# define here the parameters 
def init(g):
    # PAR
    PAR1, PAR2 = 200, 450
    # surface
    s1, s2 = 5., 2.
    # C
    C1, C2 = 10, 10
    # N
    N1, N2 = 10, 10

    leaves = [lid 
        for lid in g.vertices(scale=max_scale) 
            if g.class_name(lid)=='Leaf'
            ]
    assert len(leaves)==2

    l1, l2 = leaves
    l1, l2 = g.node(l1), g.node(l2)

    l1.PAR = PAR1
    l2.PAR = PAR2

    l1.surface = s1    
    l2.surface = s2

    l1.C = C1    
    l2.C = C2

    l1.N = N1    
    l2.N = N2

    l1.growth = False
    l2.growth = True

    return g

g = init(g)

gN = g.property('N')
gC = g.property('C')
for i in range(10):
    
    leaves = [g.node(lid) 
        for lid in g.vertices(scale=max_scale) 
            if g.class_name(lid)=='Leaf'
            ]

    # growing leaves
    mature = [l for l in leaves if not l.growth]
    growing = [l for l in leaves if l.growth]

    pool_C = sum(l.C for l in leaves)
    pool_N = sum(l.N for l in leaves)
    PAR_total = sum(l.PAR for l in leaves)

    for l in leaves:
        l.weight = l.PAR / PAR_total
        l.C = pool_C * l.weight 
        l.N = pool_N * l.weight 

    print('Step ', i, gN, gC)

    

