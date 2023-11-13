from openalea.mtg import *
import openalea.plantgl.all as pgl
from openalea.plantgl.all import *
from openalea.plantgl.all import Viewer

from numpy import *
from pandas import *
import pandas as pd
import time

import mtg_rapeseed 
import mtg_visualization 
import aptos 
import files_management as files

import alinea.caribu 
from alinea.caribu.CaribuScene import CaribuScene


class PlotMTG:
    """Run MTG methods.

    Authors
    -------
    Solen FARRA

    """
    print("Choose a method to create the MTG...")
    
    def __init__(self, csv_filename, input_dataframe = None) :
        
        self.csv_filename = self.create_MTG_from_csv_file() if csv_filename is None else input_dataframe
        self.input_dataframe = self.create_MTG_from_scratch() if input_dataframe is None else csv_filename

        return 

    # OPTION 1: We create a MTG from its known properties (NOT COMPLETELY FINISHED!):
    @staticmethod
    def create_MTG_from_csv_file():
        """Creation of the MTG from a initial csv file with properties datas.
        
        Returns
        ------
        g: MTG 
            MTG where properties are recorded

        Authors
        -------
        Frederic REES

        """
        print("Create the MTG from csv file...")

        csv_filename='MTG_00003.csv'

        # We first read the csv file where all the properties of each vertex has been previously recorded:
        dataframe = pd.read_csv(csv_filename, sep=',', header=0)
    
        # We initialize an empty MTG:
        g = MTG()

        # We define the list of vertices and properties to be added to the MTG:
        list_of_vid = dataframe['node_index']
        list_of_properties = list(dataframe.columns.values)
        list_of_properties.remove('node_index')

        # We cover each new vertex to be added to the MTG:
        for vid in list_of_vid:
            # We get the previous element on which to add the new element:
            n = g.node(vid-1)
            # We now define the new element as the child of the previous element:
            n = n.add_child()
            # For this new element, we cover all the properties defined in the csv file:
            for property in list_of_properties:
                # For the specific property, we get the single value corresponding to the current vertex ID:
                selected_series = dataframe[property].loc[dataframe['node_index'] == vid].iloc[0]
                # NB: we need to add 'iloc[0]' as the previous expression otherwise returns a series, not a single value.
        
        return g

    # OPTION 2: We create before a generic MTG:
    @staticmethod
    def create_plant_MTG_from_scratch(n_nodes, node_length):
        """Creation of the MTG from zero with a generic MTG pipeline and properties functions.
        
        Returns
        ------
        g: MTG 
            MTG where properties are recorded

        Authors
        -------
        Frederic REES

        """
        print("Create the MTG from scratch...")
        # We read the environemental/initial conditions:
        #input_dataframe = pd.read_excel(r'C:\Users\FARRA\AppData\Local\GitHubDesktop\app-3.2.3\src\caracol\src\caracol\APTOS_input_1.xlsx', header=0, sheet_name="scenario")

        #We plot the MTG:
        g = mtg_rapeseed.vegetative_rapeseed(n_nodes, node_length).above
        mtg_rapeseed.phenotype(g)
        mtg_rapeseed.metabolism(g)

        files.recording_MTG_properties(g, file_name = 'MTG_0000.csv')

        return g
    
    def create_field_MTG_from_scratch(numofseeds,
                                      nplants,
                                      length,
                                      width,
                                      sowing_density,
                                      inter_row,
                                      nrow,
                                      plant_density,
                                      DJ,
                                      node_length,
                                      phyllochrone
                                      ):
        """Creation of the MTG from zero with a generic MTG pipeline and properties functions.
        
        Returns
        ------
        g: MTG 
            MTG where properties are recorded

        Authors
        -------
        Frederic REES

        """
        print("Create the field MTG from scratch...")
        # We read the environemental/initial conditions:
        #input_dataframe = pd.read_excel(r'C:\Users\FARRA\AppData\Local\GitHubDesktop\app-3.2.3\src\caracol\src\caracol\APTOS_input_1.xlsx', header=0, sheet_name="scenario")

        #We plot the MTG:
        parcelle = mtg_rapeseed.agronomicplot(numofseeds,
                                      nplants,
                                      length,
                                      width,
                                      sowing_density,
                                      inter_row,
                                      nrow,
                                      plant_density)
             
        g = mtg_rapeseed.vegetative_rapeseed_field( DJ,
                                           node_length,
                                           phyllochrone,
                                           coord=parcelle[1]).above
        mtg_rapeseed.phenotype(g)
        mtg_rapeseed.metabolism(g)

        files.recording_MTG_properties(g, file_name = 'MTG_0000.csv')

        return g, parcelle
    
    

    @staticmethod
    def create_MTG_from_Colzette(PlantAge,
                                  phylloc,
                                  numofseeds,
                                  pattern,
                                  inter_row,
                                  plant_density) :
        """Creation of the MTG base on Colzette model methodology
        
        Authors
        -------
        Arman GRUMEL & Laurène PERTHAME
            Colzette model
        
        """

        #PlantAge plant age in DD 
        #phylloc plant phyllochrone Leaf/DD
        #numofseeds number of plants to simulate
        #pattern  pattern can be 'Rows' or 'Homogeneous'   /!\ if numofseeds=1, 'Rows' is applied by default

        outputfilename ='Results_'+str(round(PlantAge*phylloc))+'L_'+str(numofseeds)+pattern[0]

        data = pd.read_csv("ValeursParam.csv",sep=";")

        ##Cycle 
        Results = pd.DataFrame()
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
        progress = 0 

        print(len(data),"combinations to run. Let's get started !")

        for index,row in data.iterrows() :

            Temporary_Res = dict()

            progress += 1
            print('Progress : ' , progress, '/', len(data))
            
            #1 create an MTG for a plant at a certain age (DD), for a given phyllochron (n_leaves/DD)
            if pattern=='Rows' or numofseeds==1:
                parcelle = mtg_rapeseed.agronomicplot(nplants = numofseeds,
                                                    length=1,
                                                    width=1,
                                                    sowing_density = plant_density,
                                                    inter_row = inter_row,
                                                    plant_density = plant_density)
                #parcelle[0] : nb de plantes, 
                #parcelle[1] : coordonnées, 
                #parcelle[2] : limites du domaine (-x,-y,x,y),
                #parcelle[3] : surface en m2 du domaine, 
                #parcelle[4] : coefficient de conversion de mètre à l'unitée utilisée dans la scène
                    
                    
                g = mtg_rapeseed.vegetative_rapeseed_field(PlantAge,
                                                        phyllochrone = phylloc,
                                                        coord=parcelle[1]).above 

            elif pattern == 'Homogeneous': 
                parcelle = mtg_visualization.generate_centered_diagonal_mesh_points(numofseeds)
                g = mtg_rapeseed.vegetative_rapeseed_field(PlantAge,
                                                           phyllochrone = phylloc,
                                                           coord = multiply(parcelle[1],100)).above 
                

            #2 attribute architectural properties to the MTG 
            mtg_rapeseed.phenotype(g,
                                   skew = 5,
                                   phyllot = 137.5,
                                   leaf_angle = 80,
                                   coeff_petiole = 0.51017, #limb length to petiole length , 0.31017 length/whole leag length
                                   stem_d = 0.1,
                                   coeff_petiole_d = 0.5,
                                   total_surface=300,
                                   leaf_max = row["LargLeafPos"],
                                   sd = row["LS_spread"],
                                   amp = row['Bell_amp'],
                                   ins_angle = row["Ins_angle"],
                                   coeff_width = row["coefLenWid"]
                                   )
            
            #2bis attribute metabolic properties to the MTG 
            mtg_rapeseed.metabolism(g,
                                    fall = False,
                                    growth = True,
                                    photosynthesis = 0.1,
                                    par = 100,
                                    respiration = 0.1,
                                    potential_photosynthesis= 0.1,
                                    potential_area= 0.1,
                                    nitrogen_source = True,
                                    nitrogen_sink = True,
                                    nitrogen_content = 10,
                                    nitrogen_demand = 0.1,
                                    nitrogen_storage= 0.1,
                                    nitrogen_import= 0.1,
                                    nitrogen_export= 0.1,
                                    rue = 0.1,
                                    weight_import= 0.1,
                                    weight_export= 0.1,
                                    carbon_content= 0.1,
                                    carbon_demand = 0.1,
                                    carbon_deficit = 0.1,
                                    carbon_import= 0.1,
                                    carbon_export= 0.1,
                                    carbon_spend= 0.1,
                                    carbon_growth = 0.1,
                                    carbon_storage = 0.1,
                                    carbon_nitrogen_ratio= 0.1,
                                    nitrogen_fall_threshold_empirical= 0.1,
                                    PAR_fall_threshold_empirical= 0.1,
                                    carbon_fall_threshold_mecanism= 0.1,
                                    nitrogen_fall_threshold_mecanism = 0.1
                                    )


            Temporary_Res = Temporary_Res | g.properties()
        
            #TO DO (optional) : generate 3D archi and save picture 

            #3 CARIBU
            can = mtg_visualization.scene3d(g)
            zenith=[(100,(0,0,-1))]
            #opts = list(map(str, [data_path('par.opt'), data_path('nir.opt')]))

            cs = CaribuScene(scene=can,light=zenith,pattern=parcelle[2])
            raw,agg=cs.run(infinite=True,simplify=True,direct=False)

            Temporary_Res = Temporary_Res | agg

            Temporary_Res=pd.DataFrame.from_dict(Temporary_Res,orient='columns')

            Results=pd.concat([Results,Temporary_Res])
            filtered_results = Results[Results['label']=='Leaf']
            
            Results.to_csv('Colzetteresults2.csv')
            filtered_results.to_csv(outputfilename+'.csv')

            return g

class Simulation:
    """Run a simulation according a configure scenario.

    Authors
    -------
    Solen FARRA
    Frederic REES

    """
    print("Choose and init a simulation scenario...")
    # We define an initial and final step for covering a certain period:
    def __init__(self,g, default = None) :
    
        self.g
        self.default = self.main() if default is None else default
       
    @staticmethod
    def main(g, initial_step, final_step):
        """Run a default simulation.

        Parameters
        ----------
        g: MTG 
            MTG where properties are recorded

        initial_step : int
            first time step of the simulation
        
        finl_step : int
            end time step of the simulation
        
        Returns
        -------
        g: MTG 
            MTG where properties are recorded

        Authors
        -------
        Solen FARRA
        Frederic REES

        """
        print("")
        print("Start simulation...")

        tps1 = time.process_time()

        global parcelle

        # For each time step covering the simulation period:
        for step in range(initial_step, final_step+1):
            print("")
            print("For step", step,":")

            # We successively call the different class and methods:
            aptos.Photosynthesis(g, parcelle = parcelle,
                                  photosynthesis = None, 
                                  caribu = None)
            
            #aptos.PotentialGrowth(g,
                            #potential_photosynthesis = None,
                            #potential_area = None, 
                            #area = None)
            
            #aptos.Competition(g, 
                        #carbon_deficit = None,
                        #source_N = None, 
                        #sink_N = None)
            
            #aptos.Nitrogen(g,
                     #nitrogen_demand = None,
                     #weight_import = None,
                     #weight_export = None,
                     #nitrogen_import =None,
                     #nitrogen_export = None,
                     #nitrogen_storage = None,
                     #carbon_expend = None,
                     #nitrogen_content = None,
                     #RUE = None,
                     #carbon_nitrogen_ratio = None)
            
            #aptos.Carbon(g, 
                   #carbon_growth = None,
                   #carbon_demand = None,
                   #carbon_import = None,
                   #carbon_export = None)
            
            #aptos.Growth(g, area = None)
            
            #aptos.Fall(g, 
                 #fall = None,
                 #nitrogen_fall_threshold = None,
                 #nitrogen_fall_threshold_empirical = False,
                 #PAR_fall_threshold_empirical = False,
                 #carbon_fall_threshold = None)
            
            #cs = mtg_visualization.scene3d(g, select_visitor = mtg_visualization.ColzetteVisitor)
            # We view the plant MTG with PlantGL:
            Viewer.display(mtg_visualization.scene3d(g, select_visitor = mtg_visualization.ColzetteVisitor))
            #cs = mtg_visualization.scene3d(g, select_visitor = mtg_visualization.ColzetteVisitor)
            #cs.plot(g.properties()['Ei'],display=True)

            # We recording a snapshot a the MTG:
            plot_name = 'plot_%.5d.png' % step
            pgl.Viewer.saveSnapshot(plot_name)

            # We record the main properties of the plant:
            properties_filename = 'MTG_%.5d.csv' % step
            files.recording_MTG_properties(g, properties_filename)

            return g

        print("")
        print("Simulation is done!")

        tps2 = time.process_time()
        print('Execution time: ', tps2 - tps1)

        return g

        
#RUN
g, parcelle = PlotMTG.create_field_MTG_from_scratch(numofseeds = 18,
                                                    nplants = 18,
                                                    sowing_density = 35,
                                                    inter_row = 0.2,
                                                    nrow = 3,#raw by box (PHOCUS)
                                                    plant_density = 35,
                                                    DJ = 950,
                                                    node_length = 0.02,
                                                    phyllochrone = 0.015,
                                                    length = None,
                                                    width = None)#if nplants is not None, keep always length and width = None

#g = PlotMTG.create_MTG_from_Colzette(PlantAge = 950,
                                  #phylloc = 0.015,
                                  #numofseeds = 9,
                                  #pattern = 'Rows',
                                  #inter_row = 0.2,
                                  #plant_density = 35)

#g = PlotMTG.create_MTG_from_scratch(n_nodes = 10, 
                                    # node_length = 0.02)

Simulation.main(g, initial_step = 1, final_step = 3)

