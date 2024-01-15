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
from alinea.caribu.data_samples import data_path

from colzette_lib import *

tps1 = time.process_time()

#grâce à l'hypercube on a une matrice avec 
# 1 param par colonne 
# 1 combinaison de paramètres par ligne

#On itere sur toutes les lignes en excutant la séquence suivante à chaque fois
# 1. crée le mtg correspondant aux paramètres
# 2. faire tourner Caribu dessus 
# 3. retourner : par absorbé total (pour tout le couvert),
#               (par absorbé par individu ?), 
#                par absorbé par métamere 

# garder ça dans un grand tableau d'acceuil des résultats 




data = pandas.read_csv("/Users/armangrumel/Documents/CEC ECOSYS/COLZA/docs_stats/NonCorrelated100.csv",sep=";")

### paramètres visuels : 
# couleur des elements 
colors = dict()
colors['Internode'] = (100, 100, 80)
colors['Petiole'] = (100, 100, 80)
colors['Leaf'] = (9, 82, 40)

# forme des feuilles
leafshape = NurbsPatch([[(-0.0241355, -0.0549754, -0.00804807, 1), (-0.00877065, -0.00933309, 0.00100194, 1), (0.0221632, 0.0405477, -0.00898557, 1), (-0.00223729, -0.00247516, 0.000385988, 1)], 
                        [(0.314084, -0.472921, 0.392759, 1), (0.000957513, 0.000400347, 0.33376, 1), (-0.00092555, 0.00517653, 0.332921, 1), (0.408856, 0.483731, 0.413612, 1)], 
                        [(0.329337, -0.501345, 0.730665, 1), (0.000326354, 0.0164005, 0.665257, 1), (-0.000118909, 0.019707, 0.66817, 1), (0.414001, 0.558038, 0.74373, 1)], 
                        [(0.00223079, 0.0223592, 0.999352, 1), (0.00112035, 0.028812, 1.0005, 1), (-0.000767699, 0.0327184, 0.999658, 1), (-0.0011293, 0.0361733, 1.00279, 1)]])



##Cycle 

LeafS=list()
sumLight=list()

for index,row in data.iterrows() :
    #1 create an MTG for a plant at a certain age (DD), for a given phyllochron (n_leaves/DD)
    g=vegetative_rapeseed_DJ(950,phyllochrone=0.015).above 
    #2 attribute architectural properties to the MTG 
    phenotype(g,
            leaf_max = row["LargLeafPos"],
            skew = row["LS_spread"],
            ins_angle = row["Ins_angle"],
            coeff_leaf = row["coefSurfLen"],
            coeff_width = row["coefLenWid"],
            coeff_petiole= row["coefLenPeti"],
            )
    
    #optional : generate 3D archi and save picture 

    #3 CARIBU

    can = scene3d(g)
    sky = str(data_path('Turtle16soc.light'))
    #opts = list(map(str, [data_path('par.opt'), data_path('nir.opt')]))
    pattern = str(data_path('filter.8'))
    
    cs = CaribuScene(scene=can,pattern=pattern)
    raw,agg=cs.run(simplify=True)

    #cs.plot(raw['Ei'],display=True)
    LeafS.append(__builtins__.sum(g.properties()['LeafWidth'].values())*__builtins__.sum(g.properties()['LeafLength'].values()))
    #print(LS)
    


    sumLight.append(__builtins__.sum(agg['Ei'].values()))
    #print(sumLight)
    
from matplotlib import pyplot

print(LeafS)
print(sumLight)
pyplot.plot(LeafS,sumLight,'bo') 
pyplot.show()

 
tps2 = time.process_time()
print('Execution time: ', tps2 - tps1)
