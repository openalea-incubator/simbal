from CN_plant_balance import plant
  
    
class ICNbalance(plant):
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

        def define(_init_, self):
            """Plant scale turnover mecanism 
        
            Parameters
            ----------
            g: MTG
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
                leaf exported nitrogen from the source/senescent leaves to the N pool
            N_senescnence_threshold : float
                leaf senescence initiation nitrogen threshold according to the leaf state (MTG)
            C_senescnence_threshold : float
                leaf senescence initiation carbon threshold according to the leaf state (MTG)
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

            H10 Senescence initiation of leaf choose according to their PAR (lowest PAR = sacrificial leaf) // Samurai law

            H11 Senescent leaf participate to the N pool according to their PAR (higher PAR = greater weight)

            H12 Carbon plant stress definit nitrogen remobilization plant necessities (photosynthetic system turnover to optimise the biomass production according to the environmental factors)

            H13 Fall of unnecessary senescent leaves (no N remobilization at t time step)

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
            Bertheloot(2008),
        
            """
            self.define_leaf_status()
            self.define_C_pool_budget()
            
            return self


class NoCNbalance(plant):
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

        def define(_init_, self):
            """Plant scale turnover mecanism 
        
            Parameters
            ----------
            g: MTG
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
                leaf exported nitrogen from the source/senescent leaves to the N pool
            N_senescnence_threshold : float
                leaf senescence initiation nitrogen threshold according to the leaf state (MTG)
            C_senescnence_threshold : float
                leaf senescence initiation carbon threshold according to the leaf state (MTG)
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

            H10 Senescence initiation of leaf choose according to their PAR (lowest PAR = sacrificial leaf) // Samurai law

            H11 Senescent leaf participate to the N pool according to their PAR (higher PAR = greater weight)

            H12 Carbon plant stress definit nitrogen remobilization plant necessities (photosynthetic system turnover to optimise the biomass production according to the environmental factors)

            H13 Fall of unnecessary senescent leaves (no N remobilization at t time step)

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
            Bertheloot(2008),
        
            """
            self.define_leaf_status = None
            self.define_C_pool_budget = None

            return self



#CN BALANCE PLUGIN
class CNBalance(ICNbalance, NoCNbalance):
        def CN_plant_balance(version = 'default'):
            """Plant nutrients balance plugin
            
            """
            if version == 'default': 
                return ICNbalance()
            
            elif version == 'no':
                return NoCNbalance()
            
            else:
                return NoCNbalance()
