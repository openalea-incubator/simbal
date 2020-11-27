def aerial_plant(N0):
    """
    < : successor
    + : ramification
    N0: the initial quantity of N available for the seedling
    """
    g = mtg.MTG()
    vid = g.add_component(g.root, label='P')
    vid = g.add_component(vid, label='A')
    vid = g.add_component(vid, label='U')
    bud = g.add_child(vid, label = 'U', N = N0)

    stem = g.add_component(vid, label='S', N = N0) # organ
    leaf = g.add_child(stem, label = 'L', N=N0)

    return g