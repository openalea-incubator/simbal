from openalea.mtg import *
import openalea.plantgl.all as pgl

from numpy import *
from pandas import *
import pandas as pd

import mtg_rapeseed as mtg_rapeseed
import mtg_visualization as mtg_visualization
import aptos as aptos
import files_management as files


class ScratchMTG:
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
    def create_MTG_from_scratch():
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
        g = mtg_rapeseed.vegetative_rapeseed().above
        mtg_rapeseed.phenotype(g)
        mtg_rapeseed.metabolism(g)

        # We initialize the content of N in each leaf/internode:
        #for vid in g.vertices_iter(scale=1):
            #n = g.node(vid)
            # n.N_content = input_dataframe['Initial_N_content'].loc[vid-1]
            #n.N_content =input_dataframe['Initial_N_content'].loc[input_dataframe['element_ID'] == vid]

        # We initialize the PAR in each leaf/internode:
        #for vid in g.vertices_iter(scale=1):
            #n = g.node(vid)
            #n.PAR = input_dataframe['PAR'].loc[input_dataframe['element_ID'] == vid]

        # We initialize a new property of the MTG g that corresponds to the status "Fallen / Not fallen":
        #for vid in g.vertices_iter(scale=1):
            #n = g.node(vid)
            #n.fallen = False

        # We record the properties in its initial state:
        files.recording_MTG_properties(g, file_name = 'MTG_0000.csv')

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
    def __init__(self, default = None) :
    
        self.default = self.main() if default is None else default
       
    @staticmethod
    def main(initial_step, final_step):
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
        print("Start simulation...")

        # For each time step covering the simulation period:
        for step in range(initial_step, final_step+1):
            print("")
            print("For step", step,":")

            # We successively call the different class and methods:
            aptos.Photosynthesis(g, photosynthesis = None)
            
            aptos.PotentialGrowth(g,
                            potential_photosynthesis = None,
                            potential_area = None, 
                            area = None)
            
            aptos.Competition(g, 
                        carbon_deficit = None,
                        source_N = None, 
                        sink_N = None)
            
            aptos.Nitrogen(g,
                     nitrogen_demand = None,
                     weight_import = None,
                     weight_export = None,
                     nitrogen_import =None,
                     nitrogen_export = None,
                     nitrogen_storage = None,
                     carbon_expend = None,
                     nitrogen_content = None,
                     RUE = None,
                     carbon_nitrogen_ratio = None)
            
            aptos.Carbon(g, 
                   carbon_growth = None,
                   carbon_demand = None,
                   carbon_import = None,
                   carbon_export = None)
            
            aptos.Growth(g, area = None)
            
            aptos.Fall(g, 
                 fall = None,
                 nitrogen_fall_threshold = None,
                 nitrogen_fall_threshold_empirical = False,
                 PAR_fall_threshold_empirical = False,
                 carbon_fall_threshold = None)
            
            # We view the plant MTG with PlantGL:
            pgl.Viewer.display(mtg_visualization.scene3d(g))

            # We recording a snapshot a the MTG:
            plot_name = 'plot_%.5d.png' % step
            pgl.Viewer.saveSnapshot(plot_name)

            # We record the main properties of the plant:
            properties_filename = 'MTG_%.5d.csv' % step
            files.recording_MTG_properties(g, properties_filename)

        print("")
        print("Simulation is done!")

        
#RUN
g = ScratchMTG.create_MTG_from_scratch()
Simulation.main(initial_step = 1, final_step = 3)
