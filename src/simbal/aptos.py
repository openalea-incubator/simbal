from openalea.mtg import *

from numpy import *
import numpy as np
from pandas import *
import pandas as pd

import mtg_rapeseed as mtg_rapeseed

###################################

def update_MTG(g, outputs):
    """Update the MTG after run methods.

    Parameters
    ----------
    g: MTG 
        MTG where properties are recorded

    outputs: dict
        New properties values after run methods
    
    Returns
    -------
    g: MTG 
        MTG where properties are recorded

    Examples
    --------
    >>> outputs = {'param1': { vid1: , vid2, ...}, 'param2': { vid1: , vid2, ...}}

    Authors
    -------
    Armand GRUMEL and Laurène PERTHAME
        Colzette model

    """
    print("Update the MTG...")

    # add the missing property
    for param in outputs.keys():
        if param not in g.properties():
            g.add_property(param)

        # update the MTG
        g.property(param).update(outputs[param])

########################################################################################
class Photosynthesis:
    """Photosynthesis methods.

    Authors
    -------
    Solen FARRA

    """
    def __init__(self, g, photosynthesis = None) :

        self.g = g
        self.photosynthesis = self.leaf_photosynthesis(g) if photosynthesis is None else photosynthesis

        return
    
    # This function calculates the carbon produced by a leaf:
    @classmethod
    def leaf_photosynthesis(self, g, alpha = 0.1) :
        """Calculation of the actual quantity of carbon produced of each organ.
    
        Parameters
        ---------
        g: MTG 
            MTG where properties are recorded
        
        alpha: MTG proporties
            quantical yield (slope coefficient) from the light response-curve of photosynthesis

        PAR: MTG proporties
            Incident light on the organ (W/m2)
        
        Respiration: MTG proporties
            Carbon quantity spend for respiration activity (g)
        
        Returns
        ------
        g: MTG 
            MTG where properties are recorded

        Authors
        -------
        Solen FARRA

        References
        ----------
        Drouet (1998)

        """
        print("Run the actual photosynthesis...")

        sortie=dict()

        alive = [k for k, v in g.properties()['Fall'].items() if v == False]
        par = [v for k, v in g.properties()['PAR'].items()]
        respiration = [v for k, v in g.properties()['Respiration'].items()]

        sortie['Photosynthesis']= dict(zip(alive,[(alpha*a)-b for a, b in zip(par, respiration)]))

        return update_MTG(g,sortie)

class PotentialGrowth:
    """Growth methods. 

    Authors
    -------
    Solen FARRA

    """
    def __init__(self, g, potential_area = None, potential_photosynthesis = None, area = None) :
        
        self.g = g
        self.potential_area = self.leaf_potential_area(g) if potential_area is None else potential_area
        self.potential_photosynthesis = self.leaf_potential_photosynthesis(g) if potential_photosynthesis is None else potential_photosynthesis

        return 
    
    # This function calculates the potential growth of the leaf surface according to RUE:
    @classmethod
    def leaf_potential_area(self, g) :
        """Estimation of the potential leaf area based on a radiation use efficiency of each organ.
    
        Parameters
        ---------
        g: MTG 
            MTG where properties are recorded
    
        RUE: MTG proporties
            Radiation Use Efficiency coefficient

        PAR : MTG proporties
            Incident light on the organ (W/m2)
        
        Returns
        ------
        g: MTG 
            MTG where properties are recorded

        Authors
        -------
        Solen FARRA

        """
        print("Run a potential area...")

        sortie=dict()

        alive = [k for k, v in g.properties()['Fall'].items()  if v == False]
        rue = [v for k, v in g.properties()['RUE'].items()]
        par = [v for k, v in g.properties()['PAR'].items()]

        sortie['PotentialArea']= dict(zip(alive,[a*b for a, b in zip(rue, par)]))

        return update_MTG(g,sortie)
    
    # This function calculates the potential net photosynthesis according to the potential leaf area:
    @classmethod
    def leaf_potential_photosynthesis(self, g,
                                      epsilon_a : float = 0.1, 
                                      epsilon_b : float = 0.2 ):
        """Calculation of a potential (optimum) quantity of carbon produced of each organ.
    
        Parameters
        ---------
        g: MTG 
            MTG where properties are recorded
        
        alpha: MTG proporties
            quantical yield (slope coefficient) from the light response-curve of photosynthesis

        PAR: MTG proporties
            Incident light on the organ (W/m2)
        
        epsilon_b: float
            Biological efficiency coefficient (g/MJ)

        epsilon_a: float
            Absorption efficiency coefficient

        Returns
        ------
        g: MTG 
            MTG where properties are recorded

        Authors
        -------
        Solen FARRA

        References
        ----------
        Monteith (1981)
            Penman-Monteith equation

        """
        print("Run a potential photosynthesis...")

        sortie=dict()

        alive = [k for k, v in g.properties()['Fall'].items()  if v == False]
        par = [v for k, v in g.properties()['PAR'].items()]

        sortie['PotentialPhotosynthesis']= dict(zip(alive,[((a*epsilon_a)*epsilon_b) for a in par]))

        return update_MTG(g,sortie)
    
class Competition:
    """Nutrient competition status methods. 

    Authors
    -------
    Solen FARRA

    """
    def __init__(self, g, 
                 carbon_deficit = None, source_N = None, sink_N = None) :
    
        self.g = g
        self.carbon_deficit = self.leaf_carbon_deficit(g) if carbon_deficit is None else carbon_deficit
        self.source_N = self.leaf_nitrogen_competition_status(g) if source_N is None else source_N
        self.sink_N = self.leaf_nitrogen_competition_status(g) if sink_N is None else sink_N
    
        return 
    
    @classmethod
    def leaf_carbon_deficit(self, g):
        """Calculation of the carbon deficit according to the actual photosynthesis and the carbon stress for each organ.
    
        Parameters
        ---------
        g: MTG 
            MTG where properties are recorded
    
        carbon demand: MTG properties
            Actual quantity of carbon produced by photosynthesis (g)
        
        potential_area: MTG properties
            Carbon quantity spend for nitrogen redistribution among organs (g)

        area: MTG properties
            Actual leaf area 
        
        Returns
        ------
        g: MTG 
            MTG where properties are recorded

        Authors
        -------
        Solen FARRA

        """
        print("Run the actual carbon deficit...")

        sortie=dict()

        alive = [k for k, v in g.properties()['Fall'].items()  if v == False]
        photosynthesis = [v for k, v in g.properties()['Photosynthesis'].items()]
        carbon_demand = [v for k, v in g.properties()['CarbonDemand'].items()]

        sortie['CarbonDeficit']= dict(zip(alive,[a-b for a, b in zip(photosynthesis, carbon_demand)]))

        return update_MTG(g,sortie)
    
    @classmethod
    def leaf_nitrogen_competition_status(self, g):
        """Evaluation of the nitrogen source-sink relationship between organs.
    
        Parameters
        ---------
        g: MTG 
            MTG where properties are recorded
    
        carbon demand: MTG properties
            Actual quantity of carbon produced by photosynthesis (g)
        
        potential_area: MTG properties
            Carbon quantity spend for nitrogen redistribution among organs (g)

        area: MTG properties
            Actual leaf area 
        
        Returns
        ------
        g: MTG 
            MTG where properties are recorded

        Authors
        -------
        Solen FARRA

        """
        print("Run the actual nitrogen competition status of each organ...") 

        sortie=dict()

        alive = [k for k, v in g.properties()['Fall'].items()  if v == False]
        carbon_deficit =[v for k, v in g.properties()['CarbonDeficit'].items()]
        nitrogen_source =[v for k, v in g.properties()['NitrogenSource'].items()]
        nitrogen_sink =[v for k, v in g.properties()['NitrogenSink'].items()]

        for v in carbon_deficit :
            if v < 0 : 
                sortie['NitrogenSource']= dict(zip(alive,[a == True for a in nitrogen_source]))
                sortie['NitrogenSink']= dict(zip(alive,[a == False for a in nitrogen_sink]))
            
            else :
                sortie['NitrogenSource']= dict(zip(alive,[a == False for a in nitrogen_source]))
                sortie['NitrogenSink']= dict(zip(alive,[a == True for a in nitrogen_sink]))  

        return update_MTG(g,sortie)

class Nitrogen:
    """Nitrogen exchange methods. 

    Authors
    -------
    Solen FARRA

    """
    def __init__(self, g,
                 nitrogen_demand = None,
                 weight_export = None,
                 weight_import = None,
                 nitrogen_import = None,
                 nitrogen_export = None,
                 nitrogen_storage = None,
                 RUE = None,
                 carbon_expend = None,
                 nitrogen_content = None,
                 carbon_nitrogen_ratio = None) :
    
        self.g = g

        self.nitrogen_demand = self.leaf_nitrogen_demand(g) if nitrogen_demand is None else nitrogen_demand
        self.weight_export = self.N_exchange_within_shoots(g) if weight_export is None else weight_export
        self.weight_import = self.N_exchange_within_shoots(g) if weight_import is None else weight_import
        self.nitrogen_import = self.N_exchange_within_shoots(g) if nitrogen_import is None else nitrogen_import
        self.nitrogen_export = self.N_exchange_within_shoots(g) if nitrogen_export is None else nitrogen_export
        self.nitrogen_storage = self.N_exchange_within_shoots(g) if nitrogen_storage is None else nitrogen_storage
        self.RUE = self.leaf_RUE(g) if RUE is None else RUE
        self.carbon_expend = self.leaf_senescence_carbon_cost(g) if carbon_expend is None else carbon_expend
        self.nitrogen_content = self.leaf_nitrogen_content(g) if nitrogen_content is None else nitrogen_content
        self.carbon_nitrogen_ratio = self.leaf_carbon_nitrogen_ratio(g) if carbon_nitrogen_ratio is None else carbon_nitrogen_ratio
    
        return 

    # This function calculates the leaf nitrogen demand according to the potential leaf area:
    @classmethod
    def leaf_nitrogen_demand(self, g):
        """Calculation of the nitrogen demand of each organ.
    
        Parameters
        ---------
        g: MTG 
            MTG where properties are recorded
        
        potential_area: MTG proporties
            Potential organ area according to the RUE (m2)

        N_content : MTG proporties
            Actual nitrogen content of the oran (g)
        
        Returns
        ------
        g: MTG 
            MTG where properties are recorded

        Authors
        -------
        Solen FARRA

        """
        print("Run the actual nitrogen demand...")

        sortie=dict()

        alive=[k for k, v in g.properties()['Fall'].items()  if v == False]
        potential_area = [v for k, v in g.properties()['PotentialArea'].items()]
        nitrogen_content = [v for k, v in g.properties()['NitrogenContent'].items()]

        sortie['NitrogenDemand']= dict(zip(alive,[a-b for a, b in zip(potential_area, nitrogen_content)]))

        return update_MTG(g,sortie)

    # This function calculates the nitrogen remobilised quantity according to ponderation,
    # the redistribution of N among leaves according to ponderation,
    # and the final leaf nitrogen content:
    @classmethod
    def N_exchange_within_shoots(self, g,
                                 PAR: float = 100.1) :
        """Leaf senescence mecanism 
        
        Parameters
        ----------
        g: MTG
            the above part of the plant

        PAR : float
            Global incident light on each organs (W/m2)
        
        PAR : MTG proporties
            Incident light on the organ (W/m2)
            
        N_content : MTG proporties
            Nitrogen content of each organs (g)

        nitrogen_demand : MTG proporties
            Nitrogen content of each organs (g)

        Returns
        -------
        g : MTG
            the above part of the plant

        References
        ----------
        Lemaire et al. (2008)

        Examples
        --------
        >>> leaf_senescence() = 

        """
        print("Run the actual nitrogen re-allocation...")

        sortie=dict()

        senescent =[k for k, v in g.properties()['NitrogenSource'].items()  if v == True]
        par =[v for k, v in g.properties()['PAR'].items()]
        sink =[k for k, v in g.properties()['NitrogenSink'].items() if v=='True']
        nitrogen_content =[v for k, v in g.properties()['NitrogenContent'].items()]
        nitrogen_demand =[v for k, v in g.properties()['NitrogenDemand'].items()]
        weight_import =[v for k, v in g.properties()['PonderationImport'].items()]

        for  k, v in g.properties()['Fall'].items() :
            if v == False :
                for k, v in g.properties()['NitrogenSource'].items() :
                    if v == True :
                        sortie['PonderationExport']= dict(zip(senescent,[a/PAR for a in par]))

                        weight_export =[v for k, v in g.properties()['PonderationExport'].items()]
                        sortie['NitrogenExport']= dict(zip(senescent,[a*b for a, b in zip(nitrogen_content,weight_export)]))

                        pool = [v for k, v in g.properties()['NitrogenExport'].items()]
                        pool_N = sum(pool)

                    else :
                    
                        sortie['PonderationImport']= dict(zip(sink,[a/PAR for a in par]))

                        weight_export =[v for k, v in g.properties()['PonderationExport'].items()]
                        sortie['NitrogenExport']= dict(zip(senescent,[a*b for a, b in zip(nitrogen_content,weight_export)]))

                        for k, v in g.properties()['NitrogenDemand'].items() :
                            if v < 0 :

                                sortie['NitrogenImport']= dict(zip(sink,[a*pool_N for a in weight_import]))
                            
                            else :

                                sortie['NitrogenStorage']= dict(zip(sink,[a for a in nitrogen_demand]))
                                sortie['NitrogenContent']= dict(zip(sink,[a for a in nitrogen_content]))

        return update_MTG(g,sortie)
            

    # This function calculates the final radiation use efficiency leaf coefficient:
    @classmethod
    def leaf_RUE(self, g):
        """Calculation of the radiation use efficiency for each organs.
    
        Parameters
        ---------
        g: MTG 
            MTG where properties are recorded
    
        PAR : MTG proporties
            Incident light on the organ (W/m2)
            
        N_content : MTG proporties
            Nitrogen content of each organs (g) 
        
        Returns
        ------
        g: MTG 
            MTG where properties are recorded

        Authors
        -------
        Solen FARRA

        References
        ----------
        Luquet et al. (2006)
            EcoMeristem model

        """
        print("Run the actual RUE...")

        sortie=dict()

        alive=[k for k, v in g.properties()['Fall'].items()  if v == False]
        par = [v for k, v in g.properties()['PAR'].items()]
        nitrogen_content = [v for k, v in g.properties()['NitrogenContent'].items()]

        sortie['RUE']= dict(zip(alive,[a/b for a, b in zip(nitrogen_content, par)]))

        return update_MTG(g,sortie)
    

    # This function calculates the carbon cost of the redistribution of N among leaves:
    @classmethod
    def leaf_senescence_carbon_cost(self, g,
                                     kr : float = 0.08):
        """Calculation of the carbon cost of nitrogen exchange among organs.
    
        Parameters
        ---------
        g: MTG 
            MTG where properties are recorded
            
        N_import : MTG proporties
            Nitrogen quantity import of each organs (g) 
        
        N_content : MTG proporties
            Actual nitrogen content of each organs (g)
        
        kr : float
            Carbon cost coefficient of nitrogen resorption for each organs (kg/C/m2)
        
        Returns
        ------
        g: MTG 
            MTG where properties are recorded

        Authors
        -------
        Solen FARRA

        References
        ----------
        Fisher et al. (2010)
        Brzostek et al. (2014)
            FUN 2.0 model

        """
        print("Run the carbon cost of nitrogen re-allocation...")

        sortie = dict()

        alive = [k for k, v in g.properties()['Fall'].items()  if v == False]
        nitrogen_import = [v for k, v in g.properties()['NitrogenImport'].items()]
        nitrogen_content = [v for k, v in g.properties()['NitrogenContent'].items()]

        sortie['CarbonSpend']= dict(zip(alive,[(kr/a)*b for a, b in zip(nitrogen_content, nitrogen_import)]))

        return update_MTG(g,sortie)
    
    @classmethod
    def leaf_nitrogen_content(self, g):
        """Calculation of the final nitrogen content after exchange among organs.
    
        Parameters
        ---------
        g: MTG 
            MTG where properties are recorded
            
        N_import : MTG proporties
            Nitrogen quantity import of each organs (g) 
        
        N_export : MTG proporties
            Nitrogen quantity export of each organs (g) 
        
        Returns
        ------
        g: MTG 
            MTG where properties are recorded

        Authors
        -------
        Solen FARRA

        References
        ----------
        Fisher et al. (2010)
        Brzostek et al. (2014)
            FUN 2.0 model

        """
        print("Run the nitrogen content...")

        sortie = dict()

        alive = [k for k, v in g.properties()['Fall'].items()  if v == False]
        senescent =[k for k, v in g.properties()['NitrogenSource'].items()  if v=='True']
        nitrogen_import = [v for k, v in g.properties()['NitrogenImport'].items()]
        nitrogen_export = [v for k, v in g.properties()['NitrogenExport'].items()]
        nitrogen_content = [v for k, v in g.properties()['NitrogenContent'].items()]

        for k, v in g.properties()['Fall'].items() :
            if v == False :

                sortie['NitrogenContent']= dict(zip(alive,[a-b for a, b in zip(nitrogen_content, nitrogen_import)]))

                for k, v in g.properties()['NitrogenSource'].items() :
                    if v=='True' :

                        sortie['NitrogenContent']= dict(zip(senescent,[a-b for a, b in zip(nitrogen_content, nitrogen_export)]))

        return update_MTG(g,sortie)

    @classmethod
    def leaf_carbon_nitrogen_ratio(self, g):
        """Evaluation of the C:N ratio of each organ for growth regulation and shoot:root balance.
    
        Parameters
        ---------
        g: MTG 
            MTG where properties are recorded
            
        N_content : MTG proporties
           Final nitrogen content of each organs (g)

        carbon_growth : MTG properties
            Carbon quantity left to growth after nitrogen exchange among organs.
        
        Returns
        ------
        g: MTG 
            MTG where properties are recorded

        Authors
        -------
        Solen FARRA

        References
        ----------
        Fisher et al. (2010)
        Brzostek et al. (2014)
            FUN 2.0 model

        """
        print("Run the C:N ratio...")

        sortie = dict()

        alive = [k for k, v in g.properties()['Fall'].items()  if v == False]
        carbon_growth = [v for k, v in g.properties()['CarbonGrowth'].items()]
        nitrogen_import = [v for k, v in g.properties()['NitrogenImport'].items()]

        sortie['CarbonNitrogenRatio']= dict(zip(alive,[a/b for a, b in zip(carbon_growth, nitrogen_import)]))

        return update_MTG(g,sortie)
    
class Carbon:
    """Carbon balance methods. 

    Authors
    -------
    Solen FARRA

    """
    def __init__(self, g, carbon_growth = None, carbon_demand = None, carbon_export = None, 
               carbon_import = None, growth = None) :

        self.g = g
        self.carbon_growth = self.leaf_carbon_growth(g) if carbon_growth is None else carbon_growth
        self.carbon_demand = self.leaf_carbon_demand(g) if carbon_demand is None else carbon_demand
        self.carbon_export = self.C_exchange_within_shoots(g) if carbon_export is None else carbon_export
        self.carbon_import = self.C_exchange_within_shoots(g) if carbon_import is None else carbon_import
        self.growth = self.C_exchange_within_shoots(g) if growth is None else growth
            
        return 
    
    # This function calculates the carbon finally available for growth after nutrient redistribution:
    @classmethod
    def leaf_carbon_growth(self, g) :
        """Calculation of a potential (optimum) quantity of carbon produced of each organ.
    
        Parameters
        ---------
        g: MTG 
            MTG where properties are recorded
    
        photosynthesis: MTG properties
            Actual quantity of carbon produced by photosynthesis (g)
        
        carbon_spend: MTG properties
            Carbon quantity spend for nitrogen redistribution among organs (g)
        
        Returns
        ------
        g: MTG 
            MTG where properties are recorded

        Authors
        -------
        Solen FARRA

        """
        print("Run the carbon available for growth...")

        sortie = dict()

        alive = [k for k, v in g.properties()['Fall'].items()  if v == False]
        photosynthesis = [v for k, v in g.properties()['Photosynthesis'].items()]
        carbon_spend = [v for k, v in g.properties()['CarbonSpend'].items()]

        sortie['CarbonGrowth']= dict(zip(alive,[a-b for a, b in zip(photosynthesis, carbon_spend)]))

        return update_MTG(g,sortie)
    
    # This function calculates the carbon demand according to the potential carbon produced associated with the potential leaf area:
    @classmethod
    def leaf_carbon_demand(self, g) :
        """Calculation of a potential (optimum) quantity of carbon produced of each organ.
    
        Parameters
        ---------
        g: MTG 
            MTG where properties are recorded
        
        potential_photosynthesis: MTG properties
            Potential (optimum) quantity of carbon produced by photosynthesis according to a potential leaf area (g)
        
        carbon_growth: MTG properties
            Carbon quantity left for growth of each organs (g)
        
        Returns
        ------
        g: MTG 
            MTG where properties are recorded

        Authors
        -------
        Solen FARRA

        """
        print("Run the carbon demand...")

        sortie = dict()

        alive = [k for k, v in g.properties()['Fall'].items()  if v == False]
        potential_photosynthesis = [v for k, v in g.properties()['PotentialPhotosynthesis'].items()]
        carbon_growth = [v for k, v in g.properties()['CarbonGrowth'].items()]

        sortie['CarbonDemand']= dict(zip(alive,[a-b for a, b in zip(potential_photosynthesis, carbon_growth)]))

        return update_MTG(g,sortie)
    
    @classmethod
    def C_exchange_within_shoots(self, g) :
        """Carbon storage mecanism and carbon stress for growth state for each organs.
        
        Parameters
        ----------
        g: MTG
            the above part of the plant

        C_storage: MTG properties
            Carbon quantity storage after carbon balance
        
        weight_export: MTG properties
            Ponderation coefficient of each organs for nutrient pool export

        weight_import: MTG properties
            Ponderation coefficient of each organs for nutrient pool import

        carbon_growth : MTG proporties
            Carbon quantity left to growth after nitrogen exchange among leaves (g)

        Returns
        -------
        g : MTG
            the above part of the plant

        References
        ----------
        Lemaire et al. (2008)

        """
        print("Run the re-allocation of storage carbon among organs...")

        sortie=dict()

        senescent =[k for k, v in g.properties()['NitrogenSource'].items()  if v == True]
        growth =[v for k, v in g.properties()['Growth'].items()]
        sink =[k for k, v in g.properties()['NitrogenSink'].items() if v == True]
        carbon_growth =[k for k, v in g.properties()['CarbonGrowth'].items() if v ==True]
        carbon_import =[k for k, v in g.properties()['CarbonImport'].items() if v ==  True]
        carbon_storage =[k for k, v in g.properties()['CarbonStorage'].items() if v == True]
        weight_import =[v for k, v in g.properties()['PonderationImport'].items()]
        weight_export =[v for k, v in g.properties()['PonderationExport'].items()]

        for  k, v in g.properties()['Fall'].items() :
            if v == False :

                sortie['CarbonExport']= dict(zip(senescent,[a*b for a, b in zip(carbon_storage, weight_export)]))

                pool =[v for k, v in g.properties()['CarbonExport'].items()]
                pool_C = sum(pool)

                for  k, v in g.properties()['NitrogenSink'].items() :
                    if v == True :
                        for k, v in g.properties()['CarbonGrowth'].items() :
                            if v < 0 :
                    
                                sortie['CarbonImport']= dict(zip(sink,[pool_C/a for a in weight_import]))
                                sortie['CarbonGrowth']= dict(zip(sink,[a+b for a, b in zip(carbon_growth, carbon_import)]))

                            for  k, v in g.properties()['CarbonGrowth'].items() :
                                if v > 0 : 

                                    sortie['Growth']= dict(zip(sink,[a == True for a in growth]))
                            
                            else :

                                sortie['Growth']= dict(zip(sink,[a == False for a in growth]))

                    else :
                        
                        sortie['Growth']= dict(zip(sink,[a == True for a in growth]))

class Growth:
    """Growth methods. 

    Authors
    -------
    Solen FARRA

    """
    def __init__(self, g, area = None) :

        self.g = g
        self.area = self.leaf_area(g) if area is None else area
    
        return 

    @classmethod
    def leaf_area(self, g) :
        """leaf area extend according to final carbon content.
        
        Parameters
        ----------
        g: MTG
            the above part of the plant

        carbon_demand : MTG proporties
            Nitrogen content of each organs (g)

        potential_area : MTG proporties
            Potential area of each organs according to a radiation use efficiency (m2)

        carbon_growth : MTG proporties
            Carbon quantity left to growth after nitrogen exchange among leaves (g)

        carbon_content : MTG proporties 
            Carbon content of each organs (g)

        Returns
        -------
        g : MTG
            the above part of the plant

        """
        print("Run the leaf area evolution...")

        sortie=dict()

        growing =[k for k, v in g.properties()['Growth'].items()  if v == True]
        potential_area =[v for k, v in g.properties()['PotentialArea'].items()]
        area =[v for k, v in g.properties()['LeafSurface'].items()]
        carbon_content =[v for k, v in g.properties()['CarbonContent'].items()]
        carbon_growth =[v for k, v in g.properties()['CarbonGrowth'].items()]

        for k, v in g.properties()['CarbonDemand'].items() :
            if v < 0 : 
                sortie['LeafSurface']= dict(zip(growing,[a for a in potential_area]))
            
            else :
                sortie['LeafSurface']= dict(zip(growing,[(a*b)/c for a, b, c in zip(carbon_growth, area, carbon_content)]))  

        return update_MTG(g,sortie)
    
class Fall:
    """Leaves falling methods. 

    Authors
    -------
    Solen FARRA

    """
    def __init__(self, g, fall = None, nitrogen_fall_threshold = None, nitrogen_fall_threshold_empirical = False,
               PAR_fall_threshold_empirical = False, carbon_fall_threshold = None) :

        self.g = g
        self.fall = self.leaf_falling(g) if fall is None else fall
        self.nitrogen_fall_threshold_empirical = self.leaf_falling_empirical_threshold(g) if nitrogen_fall_threshold_empirical is None else nitrogen_fall_threshold
        self.nitrogen_fall_threshold = self.leaf_falling_mecanism_threshold(g) if nitrogen_fall_threshold is None else nitrogen_fall_threshold_empirical
        self.PAR_fall_threshold_empirical = self.leaf_falling_empirical_threshold(g) if PAR_fall_threshold_empirical is None else carbon_fall_threshold
        self.carbon_fall_threshold = self.leaf_falling_mecanism_threshold(g) if carbon_fall_threshold is None else PAR_fall_threshold_empirical
            
        return 

    # This function considers whether a leaf should fall or not:
    @classmethod
    def leaf_falling(self, g):
        """Leaf fall mecanism.
    
        Parameters
        ---------
        g: MTG 
            MTG where properties are recorded
    
        carbon_growth : MTG proporties
            Carbon left for growth after nitrogen exchange accoring to photosynthesis of each organs (g) 
        
        Returns
        ------
        g: MTG 
            MTG where properties are recorded

        Authors
        -------
        Solen FARRA

        References
        ----------
        Field (1983)

        """
        print("Run organs fall...")

        sortie = dict()

        alive = [k for k, v in g.properties()['Fall'].items()  if v == False]
        senescent = [k for k, v in g.properties()['NitrogenSource'].items()  if v ==True]
        fall = [v for k, v in g.properties()['Fall'].items()]
        
        for k, v in g.properties()['Fall'].items() :
            if v ==False :
                for k, v in g.properties()['CarbonGrowth'].items() :
                    if v < 0 :
                        sortie['Fall']= dict(zip(senescent,[a == True for a in fall]))
                    
                    else :
                        sortie['Fall']= dict(zip(senescent,[a == False for a in fall]))

        return update_MTG(g,sortie)
    
    # This function considers the final carbon and nitrogen mecanistical threshold for leaf fall:
    @classmethod
    def leaf_falling_mecanism_threshold(self, g):
        """Leaf fall nutrient threshold estimation.
    
        Parameters
        ---------
        g: MTG 
            MTG where properties are recorded
    
        N_content : MTG proporties
            Actual nitrogen content of each organs (g)
        
        C_content : MTG propoerties
            Actual carbon content of each organs (g)

        Returns
        ------
        g: MTG 
            MTG where properties are recorded

        Authors
        -------
        Solen FARRA

        References
        ----------
        Field (1983)

        """
        print("Run fall thresholds (mecanistical hypothesis)...")

        sortie = dict()

        fall = [v for k, v in g.properties()['Fall'].items() if v == True]
        carbon_content = [v for k, v in g.properties()['CarbonContent'].items()]
        nitrogen_content = [v for k, v in g.properties()['NitrogenContent'].items()]

        for k, v in g.properties()['Fall'].items() :
            if v == True:
                
                sortie['CarbonFallThresholdMecanism']= dict(zip(fall,[a for a in carbon_content]))
                sortie['NitrogenFallThresholdMecanism']= dict(zip(fall,[a for a in nitrogen_content]))

        return update_MTG(g,sortie)
    
    # This function considers the final nitrogen and light empirical threshold for leaf fall:
    @classmethod
    def leaf_falling_empirical_threshold(self, g):
        """Leaf fall empirical (light) threshold estimation.
    
        Parameters
        ---------
        g: MTG 
            MTG where properties are recorded
    
        N_content : MTG proporties
            Actual nitrogen content of each organs (g)

        PAR : MTG proporties
            Incident light on the organ (W/m2)
        
        Returns
        ------
        g: MTG 
            MTG where properties are recorded

        Authors
        -------
        Solen FARRA

        References
        ----------
        Thomas (2013)

        """
        print("Run fall thresholds (empirical hypothesis)...")

        sortie = dict()

        fall = [v for k, v in g.properties()['Fall'].items() if v == True]
        par = [v for k, v in g.properties()['PAR'].items()]
        nitrogen_content = [v for k, v in g.properties()['NitrogenContent'].items()]

        for k, v in g.properties()['Fall'].items() :
            if v == True :
                
                sortie['PARFallThresholdEmpirical']= dict(zip(fall,[a for a in par]))
                sortie['NitrogenFallThresholdEmpirical']= dict(zip(fall,[a for a in nitrogen_content]))

        return update_MTG(g,sortie)
    