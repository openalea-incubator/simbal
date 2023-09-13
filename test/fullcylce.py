"""
Full simulation cycle script
"""

import sys
from openalea.mtg import *
from openalea.mtg.turtle import *
from openalea.plantgl import all as pgl
from openalea.plantgl.all import *
from openalea.plantgl.all import Vector3, Color3, Scaled, Viewer
import time
#from numpy import *

import pandas 

import openalea.plantgl.all as pgl
from alinea.caribu.CaribuScene import CaribuScene

from colzette_lib import *
from semis import * 

tps1 = time.process_time()

#grâce à l'hypercube on a une matrice avec 
# 1 param par colonne 
# 1 combinaison de paramètres par ligne

#On itere sur toutes les lignes en excutant la séquence suivante à chaque fois
# 1. crée le mtg correspondant aux paramètres
# 2. faire tourner Caribu dessus 
# 3. retourner : Ei -> energie incidente sur chaque triangle
#                Eabs -> energie absorbée sur chaque triangle
#                area -> surface de chaque triangle
#--> raw contient les valeurs pour chaque triangle. agg les valeurs agglomérées par objet (ex: ts les triangles d'une feuille)

# Le dataframe 'Results' collecte l'ensemble des variables intermédiaires utiles et des variables de sorties

PlantAge=950 #plant age in DD 
numofseeds=25 #number of plants to simulate
pattern='Homogeneous' # pattern can be 'Rows' or 'Homogeneous'   /!\ if numofseeds=1, 'Rows' is applied by default

outputfilename='Results_'+str(PlantAge)+'DD_'+str(numofseeds)+pattern[0]


data = pandas.read_csv("/Users/armangrumel/Documents/CEC ECOSYS/COLZA/docs_stats/ValeursParam.csv",sep=";")

### paramètres visuels : 
# couleur des elements 
colors = dict()
colors['Internode'] = (100, 100, 80)
colors['Petiole'] = (100, 100, 80)
colors['Leaf'] = (9, 82, 40)

# forme des feuilles

#leafshape = NurbsPatch([[(-0.0241355, -0.0549754, -0.00804807, 1), (-0.00877065, -0.00933309, 0.00100194, 1), (0.0221632, 0.0405477, -0.00898557, 1), (-0.00223729, -0.00247516, 0.000385988, 1)], 
#                       [(0.314084, -0.472921, 0.392759, 1), (0.000957513, 0.000400347, 0.33376, 1), (-0.00092555, 0.00517653, 0.332921, 1), (0.408856, 0.483731, 0.413612, 1)], 
#                       [(0.329337, -0.501345, 0.730665, 1), (0.000326354, 0.0164005, 0.665257, 1), (-0.000118909, 0.019707, 0.66817, 1), (0.414001, 0.558038, 0.74373, 1)], 
#                       [(0.00223079, 0.0223592, 0.999352, 1), (0.00112035, 0.028812, 1.0005, 1), (-0.000767699, 0.0327184, 0.999658, 1), (-0.0011293, 0.0361733, 1.00279, 1)]])

leafshape = unileaflet()

##


##Cycle 
Results=pandas.DataFrame()
"""
Results['SurfProp']=list() #Surface proportion repartition
Results['LeafS']=list()
Results['LeafL']=list()
Results['LeafW']=list()
Results['RealLeafS']=list()
Results['Ei']=list()
Results['Eabs']=list()
Results['TheoLeafS']=list()
Results['RatiosWL']=list()
"""
progress=0 

print(len(data),"combinations to run. Let's get started !")

for index,row in data.iterrows() :

    Temporary_Res=dict()

    progress+=1
    print('Progress : ' , progress, '/', len(data))
    #1 create an MTG for a plant at a certain age (DD), for a given phyllochron (n_leaves/DD)
    
    if pattern=='Rows' or numofseeds==1:
        parcelle = agronomicplot(nplants=numofseeds,length=1, width=1,sowing_density=35,inter_row=0.45,plant_density=35)
        #parcelle[0] : nb de plantes, 
        #parcelle[1] : coordonnées, 
        #parcelle[2] : limites du domaine (-x,-y,x,y),
        #parcelle[3] : surface en m2 du domaine, 
        #parcelle[4] : coefficient de conversion de mètre à l'unitée utilisée dans la scène
            
            
        g = vegetative_rapeseed_field(PlantAge,phyllochrone=0.015,coord=parcelle[1]).above 

    elif pattern=='Homogeneous': 
        parcelle = generate_centered_diagonal_mesh_points(numofseeds)
        g = vegetative_rapeseed_field(PlantAge,phyllochrone=0.015,coord=multiply(parcelle[1],100)).above 
        

    #2 attribute architectural properties to the MTG 
    phenotype(g,
            total_surface=300,    
            leaf_max = row["LargLeafPos"],
            sd = row["LS_spread"],
            amp = row['Bell_amp'],
            ins_angle = row["Ins_angle"],
            coeff_width = row["coefLenWid"],
            #coeff_petiole= row["coefLenPeti"], 
            )

    Temporary_Res = Temporary_Res| g.properties()
   
    #TO DO (optional) : generate 3D archi and save picture 



    #3 CARIBU
    can = scene3d(g)
    zenith=[(100,(0,0,-1))]
    #opts = list(map(str, [data_path('par.opt'), data_path('nir.opt')]))

    cs = CaribuScene(scene=can,light=zenith,pattern=parcelle[2])
    raw,agg=cs.run(infinite=True,simplify=True,direct=False)

    Temporary_Res = Temporary_Res | agg
    Temporary_Res=pandas.DataFrame.from_dict(Temporary_Res,orient='columns')

    Results=pandas.concat([Results,Temporary_Res])
    filtered_results = Results[Results['label']=='Leaf']

   



tps2 = time.process_time()
print('Execution time: ', tps2 - tps1)
 
#Results.to_csv('Colzetteresults2.csv')
filtered_results.to_csv(outputfilename+'.csv')
