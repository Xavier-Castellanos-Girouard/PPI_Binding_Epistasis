# Xavier Castellanos-Girouard

# Modified From a Script Date First Created: Oct 22 2023
# Date Last Modified: June 6th 2024


#### Import Libraries ####

import numpy as np
import pandas as pd
import networkx as nx
import pickle
from functools import reduce
from itertools import combinations
import gc



#### Define Functions ####


# Function takes pandas series where elements are strings describing lists.
# Converts string into actual list. This list is returned.
def convert_ListString_to_List(pandas_series):
    
    # Remove brackets, apostrophes and white space. Then split string by comma.
    pandas_series = pandas_series.str.replace("[", "").str.replace("]", "").str.replace("'", "").str.replace("'", "").str.replace(" ", "").str.split(",")
    
    return(pandas_series)

def boolean_indexing(v, fillval, dtype):
    try:
        lens = np.array([len(item) for item in v])
        mask = lens[:,None] > np.arange(lens.max())
        out = np.full(mask.shape,fillval, dtype = dtype)
        out[mask] = np.concatenate(v)
        
    except ValueError:
        out = np.empty(shape=(0,0), dtype = dtype)
    
    return out

# This function takes the main path for the whole projet (string) and the path
# specific to the file (strings). The paths are combined to import GI/PPI 
# overlap table (Pandas DataFrame). Dataframe is formatted to convert strings 
# to lists. A Formatted Pandas Dataframe is returned.

def import_overlap_table(main_path, file_path):
    
    # Combine main path and file path to make complete path
    complete_path = main_path + file_path
    
    # Import GI/PPI overlap table
    overlap_table = pd.read_csv(complete_path, index_col = 0)
    
    # Convert to proper datatypes: String type to list type.
    GI_PPI_optimal_module_overlap_DF['GI_Module_ORFs'] = convert_ListString_to_List(GI_PPI_optimal_module_overlap_DF['GI_Module_ORFs'])
    
    GI_PPI_optimal_module_overlap_DF['PPI_Module_ORFs'] = convert_ListString_to_List(GI_PPI_optimal_module_overlap_DF['PPI_Module_ORFs'])
    
    return(overlap_table)


def import_GI_network(main_path, 
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
    
    # Construct Graph
    GI_Epsi_network_G = nx.from_pandas_edgelist(
        df = GI_Epsi_network_DF,
        source = "ORF_query",
        target = "ORF_array",
        edge_attr = ['scores'])
    
    del GI_Epsi_network_DF
    gc.collect()
    
    return(GI_Epsi_network_G)


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
    

def trim_shortest_path(shortest_path, GI_module_dict, module_GI_arr):
    trimmed_shortest_path = shortest_path.copy()
    
    for node in trimmed_shortest_path.copy(): # For every ORF in Shortest Path
        if node in module_GI_arr: # If source ORF is in a module
            trimmed_shortest_path.remove(node) # Remove ORF from shortest path
    
    return(trimmed_shortest_path)

def get_module_pair_connectors(source_module_ID, 
    target_module_ID,
    source_module_ORFs,
    target_module_ORFs,
    GI_module_dict,
    GI_module_arr,
    GI_APSP_dict,
    GI_SP_dict):
    
    for source_ORF in source_module_ORFs: # For all nodes in source module
        for target_ORF in target_module_ORFs: # For all nodes in target module
            
            # Get Shortest Path between source and target
            shortest_path = GI_APSP_dict[source_ORF][target_ORF]
            
            # Trim Shortest Path
            shortest_path_trimmed = trim_shortest_path(shortest_path,
                GI_module_dict,
                GI_module_arr)
            
            if len(shortest_path_trimmed) > 0:
                    
                # Store Shortest Path
                GI_SP_dict[source_module_ID][target_module_ID]['Shortest_Path'].append(np.array(shortest_path, dtype = np.int16))
                    
                # Store Trimmed Shortest Path
                GI_SP_dict[source_module_ID][target_module_ID]['Trimmed_Shortest_Path'].append(np.array(shortest_path_trimmed, dtype = np.int16))
                    
                # Store Source and Target ORFs
                GI_SP_dict[source_module_ID][target_module_ID]['Source_Node'].append(source_ORF)
                GI_SP_dict[source_module_ID][target_module_ID]['Target_Node'].append(target_ORF)
            


def find_GI_SP_connectors(GI_APSP_dict, GI_module_dict, GI_module_arr):
    
    ## Make a dictionary to conserve best shortest paths:
    # Note: Some nested keys are redundant. They will be removed later.
    GI_clusters_unique_ls = list(GI_module_dict.keys())[:-1]
    GI_pairwise_module_shortest_paths_dict = {
        key1: {
            key2: {'Source_Node': [], 
                'Target_Node': [], 
                'Shortest_Path': [], 
                'Trimmed_Shortest_Path':[]
                } for key2 in GI_clusters_unique_ls
            } for key1 in GI_clusters_unique_ls
        }
    
    # Make zip object of GI Modules and Genes in GI modules. 
    # Note: [:-1] excludes the 'None' item which contains connectors 
    GI_Module_ORFs_zip = zip(list(GI_module_dict.keys())[:-1], list(GI_module_dict.values())[:-1])
    GI_Module_ORFs_list = list(GI_Module_ORFs_zip)
    GI_Module_ORFs_pairwise_comb = list(combinations(GI_Module_ORFs_list , 2))
    
    print(len(GI_Module_ORFs_pairwise_comb))
    
    ## Find shortest path for every pairwise combination of nodes from two distict modules
    for pair in GI_Module_ORFs_pairwise_comb:
        source_module_ID = pair[0][0] # Retrieve Source Module ID
        target_module_ID = pair[1][0] # Retrieve Target Module ID
        
        source_module_ORFs = pair[0][1] # Retrieve Source Module ORFs
        target_module_ORFs = pair[1][1] # Retrieve Target Module ORFs
        
        # Retreive all shortest path connectors between module pair
        get_module_pair_connectors(source_module_ID,
            target_module_ID,
            source_module_ORFs,
            target_module_ORFs,
            GI_module_dict,
            GI_module_arr,
            GI_APSP_dict,
            GI_pairwise_module_shortest_paths_dict)
    
    return(GI_pairwise_module_shortest_paths_dict)


def clean_GI_shortest_path_connectors(SP_connector_dict):
    ## Convert all to array and clean
    for i in SP_connector_dict.keys():
        for j in SP_connector_dict.keys():
            
            if i == j:
                del SP_connector_dict[i][j]
                continue
            
            # If array is empty
            if len(SP_connector_dict[i][j]['Source_Node']) == 0:
                del SP_connector_dict[i][j]
                continue
            
            SP_connector_dict[i][j]['Source_Node'] = np.array(SP_connector_dict[i][j]['Source_Node'], dtype = np.int16)
            SP_connector_dict[i][j]['Target_Node'] = np.array(SP_connector_dict[i][j]['Target_Node'], dtype = np.int16)
            
            SP_connector_dict[i][j]['Shortest_Path'] = boolean_indexing(SP_connector_dict[i][j]['Shortest_Path'], 
                fillval = -1,
                dtype = np.int16)
        
            SP_connector_dict[i][j]['Trimmed_Shortest_Path'] = boolean_indexing(SP_connector_dict[i][j]['Trimmed_Shortest_Path'],
                fillval = -1,
                dtype = np.int16)
        


### Main ####
def main():
    
    ## Main Path used for the whole project.
    main_path = "/home/xavier/Desktop/Cell_interactome_stoichiometry/Yeast_Epistasis_Stoichiometries/Third_Run/"
    
    ## Import LookUp Table
    file_path = "Python_Connector_Overlap/results/encode_LUT_dict.pickle"
    with open((main_path+file_path), "rb") as handle:
        ORF_LUT_dict = pickle.load(handle)
    
    ## Import Negative GI Network
    print("Importing Negative GI Network...")
    GI_negEpsi_network_G = import_GI_network(main_path, 
        file_path = "Python_GI_Network_Formatting/data/Costanzo_GI_processed/costanzo_2016_longer_withoutReps.csv",
        ORF_LUT_dict = ORF_LUT_dict,
        GI_sign = "Neg",
        thresh_pval = 0.05,
        thresh_score = -0.16)
    
    ## Create an array of all unique ORFs
    GI_ORF_arr = np.unique(GI_negEpsi_network_G.nodes)
    
    
    ## Sort GI ORFs into Modules
    print("Sorting GI ORFs into Modules...")
    GI_module_dict = sort_GI_ORFs(
        main_path,
        file_path = "Python_Module_Overlap/data/Costanzo2016_DataFileS6.xlsx",
        GI_ORF_arr = GI_ORF_arr,
        ORF_LUT_dict = ORF_LUT_dict)
    
    # Get list of nodes that are in a module
    GI_module_arr = np.array(reduce(lambda x, y: x+y, list(GI_module_dict.values())[:-1]))
    
    
    ## Perform All Pairs Shortest Paths Algorithm on Negative GI Graph
    print("Performing All Pairs Shortest Path Search...")
    GI_negEpsi_APSP_dict = dict(nx.all_pairs_shortest_path(GI_negEpsi_network_G))
    
    ## Clear memory
    del GI_negEpsi_network_G
    gc.collect()
    
    # Get Shortest Path connectors between modules.
    # Note: Every unique node combination between two modules will generate
    #       a shortest path.
    print("Retreiving Shortest Paths Connectors...")
    negGI_pairwise_module_shortest_paths_dict = find_GI_SP_connectors(
        GI_negEpsi_APSP_dict,
        GI_module_dict,
        GI_module_arr)
    
    
    ## Clean shortest path dict.
    # Note: Lists are converted to arrays.
    # Note: Empty dictionary entries are removed.
    print("Cleaning Shortest Paths...")
    clean_GI_shortest_path_connectors(negGI_pairwise_module_shortest_paths_dict)
    
    print("Exporting...")
    ## Export dictionary as pickle file
    #file_path = "Python_Connector_Overlap/results/GI_pairwise_module_shortest_paths_dict.pickle"
    #with open(main_path + file_path, "wb") as handle:
    #    pickle.dump(negGI_pairwise_module_shortest_paths_dict, handle)
    
    print("Done...")

if __name__ == "__main__":
    main()
