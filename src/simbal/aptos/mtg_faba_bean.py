from openalea.mtg import *
from openalea.mtg.turtle import *

from numpy import *
import numpy as np
from pandas import *

import mtg_main as mtg

#######################################

def update_MTG(g, outputs):
    """
    Update the MTG 
    :param dict outputs: {'param1': { vid1: , vid2, ...}, 'param2': { vid1: , vid2, ...}}
    """
     # add the missing property
    for param in outputs.keys():
        if param not in g.properties():
            g.add_property(param)
        # update the MTG
        g.property(param).update(outputs[param])

###########################################

def faba_bean_phenotype(g,
                        node_length: float = 1,
                        leaf_max: float=0.6587057,
                        skew: float=5,  #
                        sd: float= 0.246378, 
                        amp: float=0.1188,
                        phyllot : float = 150,
                        ins_angle : float = 50,
                        leaf_angle : float = 80, #à mesurer expérimentallement 
                        total_surface: float=100, 
                        coeff_width: float =0.631562, 
                        coeff_petiole: float =0.71017, #limb length to petiole length , 0.31017 length/whole leag length
                        stem_d: float=0.1,
                        coeff_petiole_d : float=0.5):

    sortie = dict()

    keys = [k for k, v in g.properties()['Species'].items() if v =='Faba bean']
      
    l = [k for k, v in g.properties()['label'].items()  if v =='Leaf']
    p = [k for k, v in g.properties()['label'].items()  if v=='Plant']
    i = [k for k, v in g.properties()['label'].items()  if v=='Internode']

    leaves = [k for k in keys if k in l]
    plants = [k for k in keys if k in p]
    internodes  = [k for k in keys if k in i]

    surface=mtg.generate_normalized_bell_curve_points(mean=leaf_max,std_dev=sd,amplitude=amp,N=int(len(leaves)/len(plants)))
    sortie['SurfaceRepartition']=dict(zip(leaves,surface*len(leaves)))
    sortie['LeafSurface']=dict(zip(leaves,[element * total_surface for element in surface]*len(leaves)))

    sortie['LeafLength']= dict(zip(leaves,
            [2*(element/(coeff_width*pi))**(0.5) for element in sortie['LeafSurface'].values()]*len(leaves)))

    sortie['LeafWidth']=dict(zip(leaves,
                    [element *coeff_width for element in sortie['LeafLength'].values()]*len(leaves)))
    
    sortie['SurfTheo']=dict(zip(leaves,[pi*(element/2)*(coeff_width*element/2) for element in sortie['LeafLength'].values()]*len(leaves)))

    sortie['PetioleLength']=dict(zip(leaves,[1]*len(leaves)))
    
    sortie['InsertionAngle']=dict(zip(leaves,[ins_angle]*len(leaves)))
    sortie['LeafAngle']=dict(zip(leaves,[leaf_angle]*len(leaves)))
    sortie['Phyllotaxy']=dict(zip(leaves,[phyllot]*len(leaves)))

    sortie['StemDiam']=dict(zip(internodes,[stem_d]*len(internodes)))
    sortie['PetioleDiam']=dict(zip(leaves,[stem_d*coeff_petiole_d]*len(leaves)))

    sortie['NodeLength']=dict(zip(keys,[node_length]*len(keys)))

    return update_MTG(g,sortie)


    
   

def faba_bean_metabolism(g,
    fall : bool = False,
    growth : bool = True,
    photosynthesis = 0.2,
    par = 100,
    respiration = 0.2,
    potential_photosynthesis= 0.2,
    potential_area= 0.2,
    nitrogen_source : bool = True,
    nitrogen_sink : bool = True,
    nitrogen_content = 20,
    nitrogen_demand = 0.2,
    nitrogen_storage= 0.2,
    nitrogen_import= 0.2,
    nitrogen_export= 0.2,
    rue = 0.2,
    weight_import= 0.2,
    weight_export= 0.2,
    carbon_content= 0.2,
    carbon_demand = 0.2,
    carbon_deficit = 0.2,
    carbon_import= 0.2,
    carbon_export= 0.2,
    carbon_spend= 0.2,
    carbon_growth = 0.2,
    carbon_storage = 0.2,
    carbon_nitrogen_ratio= 0.2,
    nitrogen_fall_threshold_empirical= 0.2,
    PAR_fall_threshold_empirical= 0.2,
    carbon_fall_threshold_mecanism= 0.2,
    nitrogen_fall_threshold_mecanism = 0.2): 
    """
    """
    sortie=dict()

    for k, v in g.properties()['Species'].items() :
        if v == 'Faba bean': 
    
            leaves=[k for k, v in g.properties()['label'].items()  if v=='Leaf']

            sortie['Growth']=dict(zip(leaves,[growth]*len(leaves)))

            sortie['PAR']=dict(zip(leaves,[par]*len(leaves)))
            sortie['Photosynthesis']=dict(zip(leaves,[photosynthesis]*len(leaves)))
            sortie['Respiration']=dict(zip(leaves,[respiration]*len(leaves)))

            sortie['PotentialPhotosynthesis']=dict(zip(leaves,[potential_photosynthesis]*len(leaves)))
            sortie['PotentialArea']=dict(zip(leaves, [potential_area]*len(leaves)))

            sortie['NitrogenSource']=dict(zip(leaves,[nitrogen_source]*len(leaves)))
            sortie['NitrogenSink']=dict(zip(leaves,[nitrogen_sink]*len(leaves)))
            
            sortie['NitrogenContent']=dict(zip(leaves,[nitrogen_content]*len(leaves)))
            sortie['NitrogenDemand']=dict(zip(leaves,[nitrogen_demand]*len(leaves)))
            sortie['NitrogenImport']=dict(zip(leaves,[nitrogen_import]*len(leaves)))
            sortie['NitrogenExport']=dict(zip(leaves,[nitrogen_export]*len(leaves)))
            sortie['NitrogenStorage']=dict(zip(leaves,[nitrogen_storage]*len(leaves)))
            sortie['RUE']=dict(zip(leaves,[rue]*len(leaves)))


            sortie['PonderationExport']=dict(zip(leaves,[weight_export]*len(leaves)))
            sortie['PonderationImport']=dict(zip(leaves,[weight_import]*len(leaves)))
            
            sortie['CarbonContent']=dict(zip(leaves,[carbon_content]*len(leaves)))
            sortie['CarbonDemand']=dict(zip(leaves,[carbon_demand]*len(leaves)))
            sortie['CarbonDeficit']=dict(zip(leaves,[carbon_deficit]*len(leaves)))
            sortie['CarbonImport']=dict(zip(leaves,[carbon_import]*len(leaves)))
            sortie['CarbonExport']=dict(zip(leaves,[carbon_export]*len(leaves)))
            sortie['CarbonSpend']=dict(zip(leaves,[carbon_spend]*len(leaves)))
            sortie['CarbonGrowth']=dict(zip(leaves,[carbon_growth]*len(leaves)))
            sortie['CarbonStorage']=dict(zip(leaves,[carbon_storage]*len(leaves)))
            sortie['CarbonNitrogenRatio']=dict(zip(leaves,[carbon_nitrogen_ratio]*len(leaves)))

            sortie['Fall']= dict(zip(leaves,[fall]*len(leaves)))
            sortie['NirogenFallThresholdEmpirical']=dict(zip(leaves,[nitrogen_fall_threshold_empirical]*len(leaves)))
            sortie['NirogenFallThresholdMecanism']=dict(zip(leaves,[nitrogen_fall_threshold_mecanism]*len(leaves)))
            sortie['CarbonFallThresholdMecanism']=dict(zip(leaves,[carbon_fall_threshold_mecanism]*len(leaves)))
            sortie['PARFallThresholdEmpirical']=dict(zip(leaves,[PAR_fall_threshold_empirical]*len(leaves)))

            return update_MTG(g,sortie)
        
       

