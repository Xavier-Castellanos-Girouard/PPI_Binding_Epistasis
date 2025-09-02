# Xavier Castellanos-Girouard

# Date First Created: June 17th 2024
# Date Last Modified: June 17th 2024


#### Import Libraries ####

import numpy as np
import pandas as pd
import networkx as nx
import random
import gc
import pickle
from functools import reduce

#### Define Functions ####

# This function takes the main path for the whole projet (string) and the path
# specific to the file (strings). The paths are combined to import GI module
# table (Pandas DataFrame), which contains module IDs and the ORFs within each
# module. A Formatted pandas Dataframe is returned.

def import_GI_modules(main_path, file_path, ORF_LUT_dict):
    
    # Import GI Module Table
    GI_Modules = pd.read_excel((main_path + file_path))
    
    # Assign proper names to columns
    GI_Modules.columns = list(GI_Modules.iloc[0].values)
    
    # Remove row containing column names. Reset index
    GI_Modules = GI_Modules.drop(index = 0, axis = 0, inplace = False).reset_index(drop = True)
    
    # Remove duplicated ORFs
    GI_Modules = GI_Modules[~GI_Modules.duplicated(['Systematic array ORF  name', 'Pathway/complex level (PCC > 0.4)'])]
    
    GI_Modules = GI_Modules[~GI_Modules.duplicated(['Systematic array ORF  name'])].reset_index(drop = True)
    
    ## Convert yeast ORFs to Integer IDs using Lookup Table
    GI_Modules['Systematic array ORF  name'] = GI_Modules['Systematic array ORF  name'].apply(lambda x: ORF_LUT_dict[x])
    
    return(GI_Modules)


def sort_GI_ORFs(main_path, file_path, GI_ORF_arr, ORF_LUT_dict):
    
    # Import Module Dataframe
    GI_Modules = import_GI_modules(main_path, file_path, ORF_LUT_dict = ORF_LUT_dict)
    
    # Get modules from module dataframe. Make list unique.
    GI_Modules_arr = np.unique(GI_Modules['Pathway/complex level (PCC > 0.4)'])
    
    # Make empty dictionary with keys as cluster number
    GI_module_dict = {key: [] for key in GI_Modules_arr}
    
    # Make 'None' cluster to store connector ORFs
    GI_module_dict['None'] = []
    
    # Make an array of Module ORFs
    GI_Module_ORFs = np.unique(GI_Modules['Systematic array ORF  name'])
    
    # Sort ORFs into dictionnary of GI modules
    for ORF in GI_ORF_arr:
    
        # If the ORF is a module:
        if ORF in GI_Module_ORFs:
            
            # Find cluster to which the ORF belongs:
            cluster = GI_Modules[GI_Modules['Systematic array ORF  name'] == ORF]['Pathway/complex level (PCC > 0.4)'].values
            
            # Make sure ORF belongs to only one cluster
            assert(len(cluster) == 1)
            cluster = cluster[0] # Get value from array containing 1 element
            
            # If the ORF has not yet been sorted into this cluster:
            if ORF not in GI_module_dict[cluster]:
                GI_module_dict[cluster].append(ORF) # Add ORF to cluster in dict
            
            else:
                continue
        
        # If ORF is not in a cluster (not in a module):
        elif ORF not in GI_Module_ORFs:
            GI_module_dict['None'].append(ORF) # Add to non-module list in dict
    
    return(GI_module_dict)


def import_GI_network_dataframe(main_path, 
        file_path, 
        ORF_LUT_dict, 
        GI_sign, 
        thresh_pval, 
        thresh_score):
    
    # Combine main path and file path to make complete path
    complete_path = main_path + file_path
    
    # Import Network Dataframe
    GI_Epsi_network_DF = pd.read_csv(complete_path)
    
    # Filter interactions by significance threshold
    GI_Epsi_network_DF = GI_Epsi_network_DF[GI_Epsi_network_DF['pval'] < 0.05]
    
    ## Filter interactions by sign and score threshold
    
    # For Negative Networks
    if GI_sign == "Neg":
        GI_Epsi_network_DF = GI_Epsi_network_DF[GI_Epsi_network_DF['scores'] < thresh_score]
    
    # For Positive Networks
    elif GI_sign == "Pos":
        GI_Epsi_network_DF = GI_Epsi_network_DF[GI_Epsi_network_DF['scores'] > thresh_score]
    
    # For bidirectional Networks
    elif GI_sign == "Both":
        GI_Epsi_network_DF = GI_Epsi_network_DF[
            (GI_Epsi_network_DF['scores'] > thresh_score[0]) # Apply positive threshold
            | (GI_Epsi_network_DF['scores'] < thresh_score[1])] # Apply negative threshold
    
    # Convert scores to distances
    GI_Epsi_network_DF['scores'] = 1/abs(GI_Epsi_network_DF['scores'])
    
    ## Convert Yeast ORFs to Integer IDs using Lookup Table
    GI_Epsi_network_DF['ORF_query'] = GI_Epsi_network_DF['ORF_query'].apply(lambda x: ORF_LUT_dict[x])
    GI_Epsi_network_DF['ORF_array'] = GI_Epsi_network_DF['ORF_array'].apply(lambda x: ORF_LUT_dict[x])
    
    return(GI_Epsi_network_DF)


def link_ORF_2_GI_Module(GI_module_dict):
    
    ORF_list = np.array(reduce(lambda x, y: x+y, list(GI_module_dict.values())))
    
    GI_ORF_2_module_dict = dict.fromkeys(ORF_list)
    
    for module in list(GI_module_dict.keys()):
        
        for ORF in list(GI_ORF_2_module_dict.keys()):
            if ORF in GI_module_dict[module]:
                GI_ORF_2_module_dict[ORF] = module
    
    return(GI_ORF_2_module_dict)
    

def randomize_GI_networks(main_path, n, ORF_LUT_dict):
    
    ## Import Negative GI Network
    print("Importing Negative GI Network...")
    GI_negEpsi_network_DF = import_GI_network_dataframe(main_path, 
        file_path = "data/costanzo_2016_longer_withoutReps.csv",
        ORF_LUT_dict = ORF_LUT_dict,
        GI_sign = "Neg",
        thresh_pval = 0.05,
        thresh_score = -0.16)
    
    ## Create an array of all unique ORFs
    GI_ORF_arr = np.unique(GI_negEpsi_network_DF['ORF_query'].tolist()+GI_negEpsi_network_DF['ORF_array'].tolist())
    
    ## Sort GI ORFs into Modules
    print("Sorting GI ORFs into Modules...")
    GI_module_dict = sort_GI_ORFs(
        main_path,
        file_path = "data/Costanzo2016_DataFileS6.xlsx",
        GI_ORF_arr = GI_ORF_arr,
        ORF_LUT_dict = ORF_LUT_dict)
    
    ## For each interaction. Determine what modules the ORFs are in.
    
    # Make a dictionary that maps ORFs to their modules
    GI_ORF_2_module_dict = link_ORF_2_GI_Module(GI_module_dict)
    
    # Add module information to query ORFs
    GI_negEpsi_network_DF['ModuleID_query'] = GI_negEpsi_network_DF['ORF_query'].apply(lambda x: GI_ORF_2_module_dict[x])
    
    # Add module information to array ORFs
    GI_negEpsi_network_DF['ModuleID_array'] = GI_negEpsi_network_DF['ORF_array'].apply(lambda x: GI_ORF_2_module_dict[x])
    
    
    ## Get dataframe of intra module interactions
    intramodule_GI_negEpsi_network_DF = GI_negEpsi_network_DF[(GI_negEpsi_network_DF['ModuleID_query'] != 'None') & (GI_negEpsi_network_DF['ModuleID_array'] != 'None')].copy()
    
    # Make graph from dataframe
    intramodule_GI_negEpsi_network_G = nx.from_pandas_edgelist(intramodule_GI_negEpsi_network_DF,
        source = 'ORF_array', 
        target = 'ORF_query')
    
    ## Get dataframe of all other interactions
    Connector_GI_negEpsi_network_DF = GI_negEpsi_network_DF[~((GI_negEpsi_network_DF['ModuleID_query'] != 'None') & (GI_negEpsi_network_DF['ModuleID_array'] != 'None'))].copy()
    
    # Clear memory unused dataframes
    del intramodule_GI_negEpsi_network_DF
    del GI_negEpsi_network_DF
    
    gc.collect()
    
    
    ## Create degree preserved randomized networks, where intra-module interactions are stable.
    for i in range(n):
        randomized_graph = intramodule_GI_negEpsi_network_G.copy()
        randomized_connector_DF = Connector_GI_negEpsi_network_DF.copy()
        
        ## Randomize edges 
        randomized_connector_DF['ORF_array'] = randomized_connector_DF['ORF_array'].sample(frac=1, replace = False).tolist()
        
        Connector_GI_negEpsi_network_G = nx.from_pandas_edgelist(randomized_connector_DF,
            source = 'ORF_query',
            target = 'ORF_array')
            
        randomized_graph.add_edges_from([e for e in Connector_GI_negEpsi_network_G.edges])
        
        with open((main_path+f"data/Controls/DegPres_Randomized_GI_Networks/random_GI_network_{i}_G.pickle"),'wb') as handle:
            pickle.dump(randomized_graph, handle)
    
    

### Main ####

# Note: We want to create randomized networks. However, because we are looking
#       only at connectors, interactions within modules have to stay identical.

def main():
    
    ## Main Path
    main_path = "/home/xaviercg/scratch/"
    
    ## Import LookUp Table
    file_path = "data/encode_LUT_dict.pickle"
    with open((main_path+file_path), "rb") as handle:
        ORF_LUT_dict = pickle.load(handle)
    
    ## Generate Random GI networks
    randomize_GI_networks(main_path = main_path, 
        n = 10000, 
        ORF_LUT_dict = ORF_LUT_dict)
    
    
    
if __name__ == "__main__":
    main()
