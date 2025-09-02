# Xavier Castellanos-Girouard
# Date First Created: June 22nd 2024
# Date Last Modified: June 22nd 2024


#### Import Libraries ####

import numpy as np
import pandas as pd
import pickle
from functools import reduce
import sys

#### Define Functions ####

def make_module_conversion_table(main_path, file_path):
    
    # Import GI-PPI module overlap table
    GI_PPI_Module_Overlap_DF = pd.read_csv(main_path+file_path, index_col = 0)
    
    GI_PPI_Module_Overlap_DF = GI_PPI_Module_Overlap_DF.replace({np.nan: "None"})
    
    ## Create a dictionary that maps GI module IDs to PPI module IDs
    GI_Module_ls = [int(ID) for ID in GI_PPI_Module_Overlap_DF['Cluster_ID'].tolist() if ID != 'None']
    
    conversion_dict = dict.fromkeys(GI_Module_ls)
    
    for key in list(conversion_dict.keys()):
        index = GI_PPI_Module_Overlap_DF['Cluster_ID'].tolist().index(int(key))
        conversion_dict[key] = GI_PPI_Module_Overlap_DF['CPX_ID'].tolist()[index]
    
    return(conversion_dict)



def find_shortestpath_overlap(PPI_shortest_path_dict, GI_shortest_path_dict, conversion_dict):
    
    # Construct a dictionary to store the shortest paths overlaps
    shortest_path_overlap_dict = {
        key1: {
            key2: [] for key2 in GI_shortest_path_dict[key1].keys()
            } for key1 in GI_shortest_path_dict.keys()
        }
    
    for key1 in GI_shortest_path_dict.keys():
        for key2 in GI_shortest_path_dict[key1].keys():
            print(key1, key2)
            
            # Not all GI modules will be in the conversion dict (e.g. no match with PPI module)
            # This if statement takes care of this exception
            if ((key1 not in list(conversion_dict.keys())) | (key2 not in list(conversion_dict.keys()))):
                continue
            
            # If the GI modules do not map to a PPI module, skip.
            if (conversion_dict[key1] == 'None') | (conversion_dict[key2] == 'None'):
                continue
            
            # Retrieve shortest paths from both networks
            GI_paths = GI_shortest_path_dict[key1][key2]['Trimmed_Shortest_Path']
            PPI_paths = PPI_shortest_path_dict[conversion_dict[key1]][conversion_dict[key2]]['Trimmed_Shortest_Path']
            
            # For every GI and PPI path combination
            for GI_path in GI_paths:
                for PPI_path in PPI_paths:
                
                    # Find overlap in both shortest paths
                    GI_PPI_Overlap = list(set(GI_path).intersection(PPI_path))
                    
                    # Remove all instances of -1 (filler for staggered 2d np.arrays)
                    GI_PPI_Overlap = np.array([x for x in GI_PPI_Overlap.copy() if x != -1], dtype = np.int16)
                    # If there are common nodes (overlap) conserve this in overlap dict
                    if len(GI_PPI_Overlap) >= 1:
                        shortest_path_overlap_dict[key1][key2].append(GI_PPI_Overlap)
                    
                    #if len(GI_PPI_Overlap) > 1:
                    #    print(GI_PPI_Overlap)
    
    return(shortest_path_overlap_dict)


def clean_shortestpath_overlap(shortest_path_overlap_dict, GI_shortest_path_dict):
    
    shortest_path_overlap_cleaned_dict = {
        key1: {
            key2: [] for key2 in GI_shortest_path_dict[key1].keys()
            } for key1 in GI_shortest_path_dict.keys()
        }
    
    for key1 in shortest_path_overlap_dict.keys():
        
        # Skip if there are no entries for this GI clusters (case of 112th cluster).
        if len(shortest_path_overlap_dict[key1]) == 0:
            del shortest_path_overlap_cleaned_dict[key1]
            continue
        
        # If this key contains no shortest paths at all, remove all its nested entries.
        # Note: length > 0 is necessary to avoid error from 112th cluster which is completely empty.
        if (len(reduce(lambda x, y: x+y, shortest_path_overlap_dict[key1].values())) == 0):
            del shortest_path_overlap_cleaned_dict[key1]
            continue
        
        for key2 in shortest_path_overlap_dict[key1].keys():
            if len(shortest_path_overlap_dict[key1][key2]) == 0:
                del shortest_path_overlap_cleaned_dict[key1][key2]
            
            #elif len(shortest_path_overlap_dict[key1][key2]) == 1:
            #    shortest_path_overlap_cleaned_dict[key1][key2] = shortest_path_overlap_dict[key1][key2].copy()
            # Note: If entry is == 1. Then array stays the same because no duplicates are possible.
            
            elif len(shortest_path_overlap_dict[key1][key2]) >= 1:
                shortest_path_overlap_cleaned_dict[key1][key2] = np.unique(shortest_path_overlap_dict[key1][key2].copy())
    
    return(shortest_path_overlap_cleaned_dict)
    

def convert_overlap_dict_to_df(shortest_path_overlap_clean_dict, conversion_dict):
    
    ## Make a list of lists containing shortest path overlap data
    shortest_path_data_ls = []
    for key1 in shortest_path_overlap_clean_dict.keys():
        for key2 in shortest_path_overlap_clean_dict[key1].keys():
            for path in shortest_path_overlap_clean_dict[key1][key2]:
                shortest_path_data_ls.append([key1, key2, path.copy()])
    
    
    ## Convert list of lists to dataframe
    Shortest_Path_Overlap_DF = pd.DataFrame(shortest_path_data_ls, 
        columns = ['Source_GI_Cluster', 'Target_GI_Cluster', 'Shortest_Path_Overlap'])
    
    ## Add Protein Complex IDs
    Shortest_Path_Overlap_DF['Source_CPX_ID'] = Shortest_Path_Overlap_DF['Source_GI_Cluster'].apply(lambda x: conversion_dict[x])
    Shortest_Path_Overlap_DF['Target_CPX_ID'] = Shortest_Path_Overlap_DF['Target_GI_Cluster'].apply(lambda x: conversion_dict[x])
    
    
    return(Shortest_Path_Overlap_DF)
    

#### Main Loop ####

def main():
    main_path = "/home/xaviercg/scratch/"
    
    iter_num = sys.argv[1]
    
    
    with open((main_path + f"data/Controls/Random_PPI_Shortest_Path_Connectors/PPI_pairwise_module_shortest_paths_{iter_num}_dict.pickle"), "rb") as handle:
        PPI_shortest_path_dict = pickle.load(handle)
    
    with open((main_path + f"data/Controls/Random_GI_Shortest_Path_Connectors/GI_pairwise_module_shortest_paths_{iter_num}_dict.pickle"), "rb") as handle:
        GI_shortest_path_dict = pickle.load(handle)
    
    
    conversion_dict = make_module_conversion_table(main_path,
        file_path = "data/GI_PPI_optimal_module_overlap_clean.csv")
    
    print(conversion_dict)
    
    print("Identifying Shortest Path Overlap...")
    shortest_path_overlap_dict = find_shortestpath_overlap(
        PPI_shortest_path_dict = PPI_shortest_path_dict,
        GI_shortest_path_dict = GI_shortest_path_dict,
        conversion_dict = conversion_dict)
    
    print("Clean Overlap Dictionary...")
    shortest_path_overlap_clean_dict = clean_shortestpath_overlap(
        shortest_path_overlap_dict = shortest_path_overlap_dict,
        GI_shortest_path_dict = GI_shortest_path_dict)
    
    print("Convert Overlap Dictionary to Dictionary...")
    Shortest_Path_Overlap_DF = convert_overlap_dict_to_df(
        shortest_path_overlap_clean_dict = shortest_path_overlap_clean_dict,
        conversion_dict = conversion_dict)
    
    with open((main_path+f"data/Controls/Random_Shortest_Path_Connector_Overlap/GI_PPI_pairwise_connector_overlap_{iter_num}_DF.pickle"),"wb") as handle:
        pickle.dump(Shortest_Path_Overlap_DF, handle)
    
if __name__ == "__main__":
    main()
