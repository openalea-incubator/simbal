class leaf:
    def _init_(self):
        """Leaf scale turnover module initial parameters
          
        """
        self.N
        self.N_potential_offer
        self.C
        self.PAR
        self.photosynthetic_area
        self.senescence
        self.leaf_fall
        self.C_potential_offer
        self.C_potential_leaf_demand
        self.C_producted
        self.C_storage
        self.area_by_N_unit

        N_potential_shoot_offer = None
            
        return self


    def define_leaf_senescence(_init_, g, self,
                               N_potential_shoot_offer) :
        """Leaf senescence mecanism and leaf falling
    
    Parameters
    ----------
    g: MTG
        the above part
    N_potential_plant_offer : float
        plant potential nitrogen content

    Returns
    -------
    g : MTG
        the above part
    plant_area : float
        plant total area 
    PAR_total :
        total light intercepted by the plant
    plant_N : float
        plant nitrogen content
    pool_N : float
        senescence leaves nitrogen pool
    plant_N_by_lightunit : float
        potential nitrogen quantity by light unit
    N_plant_demand : float
        plant nitrogen demand
    N_plant_deficit : float
        plant nitrogen deficit
    N_plant_hunger :
        N plant stress state
    N_plant_remobilized : float
        plant remobilized nitrogen

    State Variables (on MTG)
    ------------------------
    growth : bool
        growth state
    leaf_fall : bool
        viability state of the leaf
    senescence : bool
        leaf development status
    photosynthetic_area : float
        photosynthetic leaf surface
    N_potential_offer : float
        potential leaf nitrogen content
    N_deficit : float
        leaf nitrogen deficit
    weight : float
        value of the leaf participation to the pool accordint to his PAR
    N_exported : float
        leaf exported nitrogen from the source/senescent leaves to the N pool
    N : float
        leaf nitrogen content
    C : float
        leaf carbon content
    N_senescnence_threshold : float
        leaf senescence initiation nitrogen threshold according to the leaf state (MTG)
    C_senescnence_threshold : float
        leaf senescence initiation carbon threshold according to the leaf state (MTG)
    PAR : float
        leaf light intercepted
    
    Hypothesis
    ----------
    H10 Senescence initiation of leaf choose according to their PAR (lowest PAR = sacrificial leaf) // Samurai law

    H11 Senescent leaf participation to the N pool according to their PAR (higher PAR = greater weight)

    References
    ----------
    Skulachev (2002)
    Thomas (2013)
    """
        max_scale = g.max_scale()
    
        leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) if g.class_name(lid)=='Leaf']
        alive = [l for l in leaves if not self.leaf_fall]
        source = [l for l in alive if self.growth]
        senescence = [l for l in alive if self.senescence]

        plant_area = sum(self.photosynthetic_area for l in alive)

        for l in alive :

            plant_N = sum(self.N for l in alive)
            pool_N = sum(self.N for l in source if self.senescence)
            PAR_total = N_potential_shoot_offer
        
            l.N_demand = self.N_potential_offer
            l.N_deficit = self.N - self.N_demand


            N_plant_demand = sum(l.N_deficit for l in alive)
            N_plant_deficit = pool_N - N_plant_demand

            if N_plant_deficit > 0 :
                N_plant_hunger = False
                N_plant_remobilized = N_plant_demand

                for l in senescence :
                    l.weight = self.PAR / PAR_total
                    l.N_exported = pool_N * l.weight
                    l.N = self.N - l.N_exported
                    l.C = self.C
                   
                for l in alive:
                    l.N = self.N_potential_offer
                    l.C = self.C

            else :
                N_plant_hunger = True

                alive.sort(key = lambda l:l.PAR)

                for l in alive[0]:
                    l.senescence = True
                    l.N_senescence_threshold = self.N
                    l.C_senescence_threshold = self.C

                    for l in senescence :
                        l.weight = self.PAR / PAR_total
                        l.N_exported = pool_N * l.weight
                        l.N = l.N_senescence_threshold - l.N_exported
                        l.C = l.C_senescence_threshold

                    N_plant_remobilized = sum(l.N_exported for l in senescence)

                    senescence.sort(key = lambda l:l.PAR)

                    for l in senescence[0] : 
                        l.leaf_fall = True
                        l.N_leaf_fall = l.N - l.N_exported
                        l.C_leaf_fall = l.C
                        l.N = l.N_leaf_fall
                        l.C = l.C_leaf_fall

                    survivors = [l for l in alive if not l.leaf_fall]

                    for l in survivors:
                        l.weight = l.PAR / PAR_total
                        l.N_imported = pool_N * l.weight 
                        l.N_deficit = self.N + l.N_imported
                

            return g, plant_area, plant_N, pool_N, N_plant_demand, N_plant_deficit, N_plant_hunger, N_plant_remobilized
        
        
    def define_leaf_necrotic_area(g, self):
        """Leaf senescence necrotic surface after N remobilization (yellowing area)
    
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
    
    leaf_fall : bool
        viability state of the leaf
    senescence : bool
        leaf development status
    area_by_Nunit : float
        surface nitrogen use efficiency ?
    N_exported : float
        leaf exported nitrogen from the source/senescent leaves to the N pool
    necrotic_area : float
        leaf necrotic surface
        """
   
        max_scale = g.max_scale()
    
        leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) if g.class_name(lid)=='Leaf']
        alive = [l for l in leaves if not l.leaf_fall]
        senescence = [l for l in alive if l.senescence]

        for l in senescence:
            l.necrotic_area = l.N_exported * self.area_by_N_unit

        return g



