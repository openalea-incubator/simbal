

class plant:
        def _init_(self):
            """plant scale turnover module initial parameters
          
            """
            self.C
            self.PAR
            self.leaf_fall
            self.C_potential_offer
            self.C_potential_leaf_demand
            self.C_producted
            self.C_storage

            N_potential_shoot_offer = None
            
            return self

        def define_leaf_status(g, self):
            """Evaluation of the leaf state to compute his value (weigth) in the pool 
    
    Parameters
    ----------
    g: MTG
        the above part

    Returns
    -------
    g : MTG
        the above part

    State Variables (on MTG)
    ------------------------
    PAR : float
        leaf light intercepted
    C_potential_leaf_offer : float
        leaf potential C content
    C_potential_leaf_demand : float
        leaf potential C demand
    dC : float
        C deficit in the leaf
    growth : bool
        growth state
    sink : bool
        organe state
    source : bool
        organe state
    storage : bool
        C storage value
    
    Hypothesis
    ----------
    H7 Source/sink status according to leaves c stress

    H8 A leaf cannot be a source and a sink organ at the same time step

    References
    ----------
    Khudairi (1970)
    Wingler (2004)
    """
            max_scale = g.max_scale()

            leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) if g.class_name(lid)=='Leaf']

            for l in leaves :
                l.dC = self.C_potential_leaf_offer - self.C_potential_leaf_demand

                if l.dC < 0:
                    l.storage = False
                    l.growth = False
                    l.sink = True
                    l.source = False
    
                else:
                    l.storage = True
                    l.growth = True
                    l.sink = False
                    l.source = True

            return g

 
        def define_C_pool_budget(g, self):
            """C distribution between leaves acording to their status and value (weight)
    
    Parameters
    ----------
    g: MTG
        the above part

    Returns
    -------
    g : MTG
        the above part
    pool_C : float
        C labile available in the plant
    PAR_total :
        total light intercepted by the plant
    C_plant_hunger :
        C plant stress state

    State Variables (on MTG)
    ------------------------
    sink : bool
        organe state
    leaf_fall : bool
        viability state of the leaf
    C_storage : float
        ????
    C_producted : float
        leaf producted biomass
    PAR : float
        leaf light intercepted
    weight : float
        value of the leaf participation to the pool accordint to his PAR
    C : float
        leaf carbon content
    C_exported : float
        C exported from the source leaf
    C_imported : float
        C imported to the sink leaf
    C_deficit : float
        leaf C stress
    growth : bool
        growth state
    C_potential_demand : float
        carbon potential leaf demand
    
    Hypothesis
    ----------
    H9 leaf participation to the C pool according to their PAR (higher PAR = greater weight)
    """
            max_scale = g.max_scale()

            leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) if g.class_name(lid)=='Leaf']
            alive = [l for l in leaves if not self.leaf_fall]
            sink = [l for l in alive if l.sink]

            pool_C = sum(self.C_storage + self.C_producted for l in alive)
            PAR_total = sum(self.PAR for l in alive)

            for l in alive :
                l.weight = self.PAR / PAR_total
                l.C_exported = pool_C * l.weight
                l.C = self.C - l.C_exported

                if pool_C < 0 :
                    C_plant_hunger = True
                 
                    for l in sink :
                        l.weight = self.PAR / PAR_total
                        l.C_imported = pool_C * l.weight
                    
                        l.dC = l.leaf_status()
                        l.C_deficit = l.dC + l.C_imported

                else :
                    C_plant_hunger = False
                
                    for l in sink :
                        l.growth = True
                        l.C = self.C_potential_leaf_demand
                    
    
            return g, pool_C, PAR_total, C_plant_hunger



        