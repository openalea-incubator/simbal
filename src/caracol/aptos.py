


class leaf :
    def _init_(self) :
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
        leaf senescence initiation nitrogen threshold according to the leaf 
        state (MTG)
    C_senescnence_threshold : float
        leaf senescence initiation carbon threshold according to the leaf 
        state (MTG)
    PAR : float
        leaf light intercepted
    
    Hypothesis
    ----------
    H10 Senescence initiation of leaf choose according to their PAR 
        (lowest PAR = sacrificial leaf) // Samurai law

    H11 Senescent leaf participation to the N pool according to their PAR 
        (higher PAR = greater weight)

    References
    ----------
    Skulachev (2002)
    Thomas (2013)

    Examples
    --------
    >>> define_leaf_senescence()

    """
        max_scale = g.max_scale()
        self = l
    
        leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) 
                  if g.class_name(lid) =='Leaf']
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
                   
                for l in alive :
                    l.N = self.N_potential_offer
                    l.C = self.C

            else :
                N_plant_hunger = True

                alive.sort(key=lambda l : l.PAR)

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

                    senescence.sort(key=lambda l : l.PAR)

                    for l in senescence[0] : 
                        l.leaf_fall = True
                        l.N_leaf_fall = l.N - l.N_exported
                        l.C_leaf_fall = l.C
                        l.N = l.N_leaf_fall
                        l.C = l.C_leaf_fall

                    survivors = [l for l in alive if not l.leaf_fall]

                    for l in survivors :
                        l.weight = l.PAR / PAR_total
                        l.N_imported = pool_N * l.weight 
                        l.N_deficit = self.N + l.N_imported
                

            return g, plant_area, plant_N, pool_N, N_plant_demand, 
        N_plant_deficit, N_plant_hunger, N_plant_remobilized
         
    def define_leaf_necrotic_area(g, self):
        """Leaf senescence necrotic surface after N remobilization 
        (yellowing area)
        
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
            leaf exported nitrogen from the source/senescent leaves 
            to the N pool
        necrotic_area : float
            leaf necrotic surface

        Examples
        --------
        >>> define_leaf_necrotic_area()
        
        """
   
        max_scale = g.max_scale()
        self = l
    
        leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) 
                  if g.class_name(lid) =='Leaf']
        alive = [l for l in leaves if not l.leaf_fall]
        senescence = [l for l in alive if l.senescence]

        for l in senescence :
            l.necrotic_area = l.N_exported * self.area_by_N_unit

        return g
    
    def define_fall(_init_, g, self):
        """Leaf fall mecanism 
        
        Parameters
        ----------
        g : MTG
            the above part

        Returns
        -------
        g : MTG
            the above part
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
        senescence_initiation : bool
            scenescence initiation indicator
        leaf_fall : bool
            leaf falling indicator

        State Variables (on MTG)
        ------------------------
        N : float
            leaf nitrogen content
        C : float
            leaf carbon content
        PAR : float
            leaf light intercepted
        C_potential : float
            potential leaf carbon demand
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
        leaf_fall : bool
            viability state of the leaf
        weight : float
            value of the leaf participation to the pool accordint to his PAR
        C_exported : float
            C exported from the source leaf
        C_imported : float
            C imported to the sink leaf
        C_deficit : float
            leaf C stress
        senescence : bool
            leaf development status
        area : float
            leaf area (surface)
        N_potential : float
            potential leaf nitrogen content
        N_deficit : float
            leaf nitrogen deficit
        N_exported : float
            leaf exported nitrogen from the source/senescent leaves 
            to the N pool
        N_senescnence_threshold : float
            leaf senescence initiation nitrogen threshold according 
            to the leaf state (MTG)
        C_senescnence_threshold : float
            leaf senescence initiation carbon threshold according 
            to the leaf state (MTG)
        N_leaf_fall : float
            fall leaf nitrogen content
        C_leaf_fall : float
            fall leaf carbon content
        N_imported : float
            leaf imported nitrogen from the N pool (senescent leaves)

        Hypothesis
        ----------
        H4 Nitrogen plant distribution accordind to PAR (light)

        H5 Nutrient contents and PAR are homogeneous in the leaf area

        H7 Source/sink status according to leaves c stress

        H8 A leaf cannot be a source and a sink organ at the same time step

        H10 Senescence initiation of leaf choose according to their PAR 
            (lowest PAR = sacrificial leaf) // Samurai law

        H11 Senescent leaf participate to the N pool according to their PAR 
            (higher PAR = greater weight)

        H12 Carbon plant stress definit nitrogen remobilization plant 
            necessities (photosynthetic system turnover to optimise the biomass 
            production according to the environmental factors)

        H13 Fall of unnecessary senescent leaves 
            (no N remobilization at t time step)

        References
        ----------
        Ono (1996),
        Takashima (2004),
        Hikosaka (2005),
        Yasumura (2007),
        Allirand (2009),

        Khudairi (1970),
        Aerts and Caluwe (1994),
        Paul and Driscoll (1997),
        Wingler (2004),
        Asim (2022),

        Skulachev (2002),
        Thomas (2013),

        Thornley (1998),
        Kull (2002),
        Bertheloot(2008)

        Examples
        --------
        >>> define_fall()
        
        """
            
        max_scale = g.max_scale()
        self = l
        
        leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) 
                  if g.class_name(lid) =='Leaf']
        alive = [l for l in leaves if not self.leaf_fall]
        leaf_fall = [l for l in alive if self.senescence]

        for l in leaf_fall :
            l.leaf_fall = True

            survivors = [l for l in alive if not l.leaf_fall]

            for l in survivors :
                l.weight = l.PAR / _init_.PAR_total
                l.N_imported = _init_.pool_N * l.weight 
                l.N_deficit = self.N + l.N_imported
                l.N = l.N_deficit
                        
            alive = survivors   

        return g, self


class plant:
    def _init_(self) :
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

    def define_leaf_status(g, self) :
        """Evaluation of the leaf state to compute his value (weigth) 
        in the pool 
    
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

    Examples
    --------
    >>> define_leaf_status()

    """
        max_scale = g.max_scale()

        leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) 
                  if g.class_name(lid) =='Leaf']

        for l in leaves :
            l.dC = self.C_potential_leaf_offer - self.C_potential_leaf_demand

            if l.dC < 0 :
                l.storage = False
                l.growth = False
                l.sink = True
                l.source = False
    
            else :
                l.storage = True
                l.growth = True
                l.sink = False
                l.source = True

        return g

 
    def define_C_pool_budget(g, self) :
        """C distribution between leaves acording to their status 
        and value (weight)
        
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
        H9 leaf participation to the C pool according to their PAR 
            (higher PAR = greater weight)

        Examples
        --------
        >>> define_C_pool_budget()

        """
        max_scale = g.max_scale()
        self = l

        leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) 
                  if g.class_name(lid) == 'Leaf']
        alive = [l for l in leaves if not self.leaf_fall]
        sink = [l for l in alive if l.sink]

        pool_C = sum(self.C_storage + self.C_producted for self in alive)
        PAR_total = sum(self.PAR for self in alive)

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

# FACTORY

class CNBalance(leaf, plant) :
    def CN_balance(version='default') :
        """Leaf and plant nutrient balance interface

        """
        if version == 'default': 
            return plant(), leaf()
            
        elif version == 'no':
            return leaf()
            
        else:
            return plant()
        
if __name__ == "__main__" :
    """Run the unit test suite 

    """
    import doctest
    doctest.testmod()
