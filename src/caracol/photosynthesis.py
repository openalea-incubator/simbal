from math import (exp, pi)

class leaf:
    def _init_(self):
        """Leaf photosynthesis module initial parameters
          
        """
        self.N
        self.C
        self.leaf_fall
        self.senescence
        self.necrotic_area
        self.PAR
        self.length
        self.radius

        N_potential_plant_offer = None
        N_potential_root_demand = None

        k = None
        e_amax = None
        RUE = None

        return self


    def define_N_potential_shoot_offer(g, self,
                                N_potential_plant_offer: float, 
                                N_potential_root_demand: float):
        """Calculation of the potential shoot demand (gN) in function of N available to the plant (gN) and potential root N demand (gN)

        Parameters
        ----------
        N_potential_shoot_offer: float
            N offer to the plant in gN
        N_potential_root_demand: float
            N potential root N demand in gN
    
        Returns
        -------
        N_potential_plant_offer: float
            N shoot potential offer in gN
    
        Hypothesis
        ----------
        H1 Roots take precedence on the nitrogen uptake

        """
        max_scale = g.max_scale()

        leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) if g.class_name(lid)=='Leaf']
        alive = [l for l in leaves if not self.leaf_fall and self.senescence]
        senescence = [l for l in alive if self.senescence]

        N_shoot = N_potential_plant_offer - N_potential_root_demand
        PAR_total = sum(self.PAR for l in leaves)

        for l in alive:
            l.weight = self.PAR / PAR_total
            l.N_available = self.N * l.weight
    
        for l in senescence:
            l.N_available = self.N
        
        N_pool = sum(l.N_available for l in alive) 
        N_potential_shoot_offer = N_shoot + N_pool
    
        return PAR_total,N_pool, N_potential_shoot_offer


    def define_leaf_area(g, self):
        """Calculation of the mtg leaf surface 

        Parameters
        ----------
        g : MTG
            the above part
    
        Returns
        -------
        g : MTG
            the above part

        State Variables (on MTG)
        ------------------------
        length : float
        radius : float
        area : float
            leaf surface acording to length and radius
        """
        max_scale = g.max_scale()

        leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) if g.class_name(lid)=='Leaf']
        alive = [l for l in leaves if not self.leaf_fall]

        for l in alive :
            l.area = 2 * pi * self.length  * self.radius
    
        return g
    

    def define_photosynthetic_leaf_area(g, self):
        """Calculation of the photosynthetic leaf area (m2) 

        Parameters
        ----------
        g : MTG
            the above part
    
        Returns
        -------
        g : MTG
            the above part

        State Variables (on MTG)
        ------------------------
        area : float
            leaf surface acording to length and radius
        necrotic_area : float
            non photosynthetic leaf area
        photosynthetic_area : float
            green leaf area
    
        Hypothesis
        ----------
        H2 leaf necrotic surface doesn't have photosynthetic activities anymore
        """
        max_scale = g.max_scale()

        leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) if g.class_name(lid)=='Leaf']
        alive = [l for l in leaves if not self.leaf_fall]

        for l in alive:
            l.photosynthetic_area = l.area - self.necrotic_area

        return g


    def define_leaf_photosynthesis(g, self,
                            RUE: float):
        """Transform radiation intercepted into biomass by the leaf

   Parameters
    ----------
    g : MTG
        the above part
    RUE : float
        ?????   
    
    Returns
    -------
    g : MTG
        the above part

    State Variables (on MTG)
    ------------------------
    PAR : float
        leaf intercepted light
    C_producted : float
        leaf producted biomass by photosynthesis 
    
    Hypothesis
    ----------
    H3 The Monteith law is valid at the leaf scale (with MTG)
    """
        max_scale = g.max_scale()

        leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) if g.class_name(lid)=='Leaf']
        alive = [l for l in leaves if not self.leaf_fall]

        for l in alive:
            l.C_producted = RUE * self.PAR

        return g


    def define_leaf_N_potential_offer(g, self,
                               PAR_total,
                               N_potential_shoot_offer):
        """Calculation of leaves N potential demand (in gN) in function of the N by light unit necessary

    Parameters
    ----------
    g : MTG
        the above part
    N_potential_plant_offer : float
        N shoot potential offer in gN

    Returns
    -------
    g : MTG
        the above part
    PAR_total : float
        total leaves PAR
    plant_N_by_lightunit : float
        nitrogen light use efficiency ?

    State Variables (on MTG)
    ------------------------
    PAR : float
        leaf intercepted light
    N_potential_offer : float
        nitrogen potential leaf offer
    
    Hypothesis
    ----------
    H4 Nitrogen plant distribution accordind to PAR (light)

    H5 Nutrient contents and PAR are homogeneous in the leaf area

    References
    ----------
    Takashima (2004)
    Hikosaka (2005)
    Yasumura (2007)
    """
        max_scale = g.max_scale()

        leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) if g.class_name(lid)=='Leaf']
        alive = [l for l in leaves if not self.leaf_fall]

        plant_N_by_lightunit = N_potential_shoot_offer / PAR_total

        for l in alive:
            l.N_potential_offer = plant_N_by_lightunit * self.PAR
   

        return g, PAR_total, plant_N_by_lightunit


    def define_potential_leaf_area(g, self):
        """Calculation of potential formed leaf area  (cm2) in function of potential N leaf offer (gN) and the coefficient of the relationship of leaf area - leaf N quantity (cm2.g-1)

    Parameters
    ----------
    g : MTG
        the above part

    Returns
    -------
    g : MTG
        the above part

    State Variables (on MTG)
    ------------------------
    N : float
        nitrogen leaf content
    area : float
        leaf surface
    area_by_Nunit : float
        surface nitrogen use efficiency ?
    N_potential_offer : float
        nitrogen potential leaf offer
    potential_area : float
        potential leaf surface
    
    Hypothesis
    ----------
    H5 Nutrient contents and PAR are homogeneous in the leaf area
    """
        max_scale = g.max_scale()

        leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) if g.class_name(lid)=='Leaf']
        alive = [l for l in leaves if not self.leaf_fall]

        for l in alive:
            l.area_by_Nunit = l.photosynthetic_area / self.N
            l.potential_area = l.area_by_Nunit * l.N_potential_offer
        
        return g

    def define_potential_leaf_light_interception(g, self,
                                          k: int,
                                          e_amax: float):
        """Intercepted light by the leaf

   Parameters
    ----------
    g : MTG
        the above part
    k: int
        extinction coefficient
    Light : int
        global incident light

    Returns
    -------
    g : MTG
        the above part
    
    State Variables (on MTG)
    ------------------------
    potential_area : float
        potential leaf surface
    PAR_potential : float
        potential leaf intercepted light
    
    Hypothesis
    ----------
    H3 The Beer-Lambert law is valid at the leaf scale (with MTG) 
    
    References
    ----------
    Bertheloot et al. 2008
    """
        max_scale = g.max_scale()

        leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) if g.class_name(lid)=='Leaf']
        alive = [l for l in leaves if not self.leaf_fall]

        for l in alive:
            l.PAR_potential = k * e_amax * (1 - exp(k * l.potential_area))
  
        return  g


    def define_leaf_potential_photosynthesis(g, self,
                                      RUE: float):
        """Transform radiation intercepted into biomass by the leaf

    Parameters
    ----------
    g : MTG
        the above part
    RUE : float
        radiation use efficiency ?

    Returns
    -------
    g : MTG
        the above part
    
    State Variables (on MTG)
    ------------------------
    PAR_potential : float
        potential leaf intercepted light
    C_potential_producted : float
        potential leaf producted biomass
    
    Hypothesis
    ----------
    H3 The Beer-Lambert law is valid at the leaf scale (with MTG)
    """
        max_scale = g.max_scale()

        leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) if g.class_name(lid)=='Leaf']
        alive = [l for l in leaves if not self.leaf_fall]

        for l in alive:
            l.C_potential_producted = RUE * l.PAR_potential

        return g


    def define_FSPM_Allirand_law(g, self):
        """Allirand law ; relationship between the light intercepted by the leaf and the nitrogen content

    Parameters
    ----------
    g : MTG
        the above part

    Returns
    -------
    g : MTG
        the above part

    State Variables (on MTG)
    ------------------------
    PAR_potential : float
        potential leaf intercepted light
    N_potential_offer : float
        nitrogen potential leaf offer
    FSPM_Allirand : float
        relationship between the light intercepted by the leaf and the nitrogen content
    
    Hypothesis
    ----------
    H6 The Allirand law is valid at the leaf scale (with MTG)

    References
    ----------
    Allirand (2009)
    """
        max_scale = g.max_scale()

        leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) if g.class_name(lid)=='Leaf']
        alive = [l for l in leaves if not self.leaf_fall]

        for l in alive:
            l.FSPM_Allirand = l.PAR_potential * l.N_potential

        return g


    def define_C_potential_leaf_demand(g, self):
        """Calculation of C potential leaf demand (in gC) in function of the potential photosynthetic area and leaf N density
    
    Parameters
    ----------
    g : MTG
        the above part

    Returns
    -------
    g : MTG
        the above part
    
    State Variables (on MTG)
    ------------------------
    area : float
        leaf surface
    potential_are : flooat
        potential leaf surface
    C_by_areaunit : float
        carbon surface use efficiency ?
    C_potential_demand : float
        leaf potential carbon demand
    
    Hypothesis
    ----------
    H5 Nutrient contents and PAR are homogeneous in the leaf area
    """
        max_scale = g.max_scale()

        leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) if g.class_name(lid)=='Leaf']
        alive = [l for l in leaves if not self.leaf_fall]

        for l in alive:
            l.C_by_areaunit = self.C / l.area
            l.C_potential_demand = l.C_by_areaunit * l.potential_area
        
        return g


    def define_C_potential_leaf_offer(g, self):
        """Calculation of C available to the leaf (in gC), in function of C uptake (in gC) 
    
    Parameters
    ----------
    g : MTG
        the above part

    Returns
    -------
    g : MTG
        the above part
    
    State Variables (on MTG)
    ------------------------
    C_producted : float
        leaf producted biomass
    C_potential_offer : float
        leaf potential carbon offer
    """
        max_scale = g.max_scale()

        leaves = [g.node(lid) for lid in g.vertices(scale=max_scale) if g.class_name(lid)=='Leaf']
        alive = [l for l in leaves if not self.leaf_fall]

        for l in alive:
            l.C_potential_offer = self.C + l.C_producted - l.C_respiration         

        return g

