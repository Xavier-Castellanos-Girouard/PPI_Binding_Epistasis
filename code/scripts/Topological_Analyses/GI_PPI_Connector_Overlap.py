# Xavier Castellanos-Girouard

# Date Last Modified: June 22th 2024
# Date First Created: June 22th 2024

#### Import libraries ####

import numpy as np
import pandas as pd
import pickle
from tqdm import tqdm
from functools import reduce


#### Import Connectors and Module Overlap data ####

with open("../results/IntermediateFiles/PPI_pairwise_module_shortest_paths_dict.pickle", "rb") as handle:
    PPI_shortest_path_dict = pickle.load(handle)
    
with open("../results/IntermediateFiles/GI_pairwise_module_shortest_paths_dict.pickle", "rb") as handle:
    GI_shortest_path_dict = pickle.load(handle)
    

## Import GI-PPI module overlap table
GI_PPI_Module_Overlap_DF = pd.read_csv("../results/Topological_Analyses/GI_PPI_optimal_module_overlap_clean.csv", index_col = 0)

## Create a dictionary that maps GI module IDs to PPI module IDs
GI_Module_ls = [int(ID) for ID in GI_PPI_Module_Overlap_DF['Cluster_ID'].tolist() if ID != 'None']

conversion_dict = dict.fromkeys(GI_Module_ls)

for key in list(conversion_dict.keys()):
    index = GI_PPI_Module_Overlap_DF['Cluster_ID'].tolist().index(key)
    conversion_dict[key] = GI_PPI_Module_Overlap_DF['CPX_ID'].tolist()[index]


## Import LUT
with open("../results/IntermediateFiles/decode_LUT_dict.pickle", "rb") as handle:
    LUT = pickle.load(handle)



#### Check for overlap in GI and PPI connectors ####

shortest_path_overlap = {key1: {key2: [] for key2 in GI_shortest_path_dict[key1].keys()} for key1 in GI_shortest_path_dict.keys()}

for key1 in tqdm(GI_shortest_path_dict.keys()):
    for key2 in GI_shortest_path_dict[key1].keys():
        
        if ((key1 not in list(conversion_dict.keys())) | (key2 not in list(conversion_dict.keys()))):
            continue
        
        if (conversion_dict[key1] == 'None') | (conversion_dict[key2] == 'None'):
            continue
        
        GI_paths = GI_shortest_path_dict[key1][key2]['Trimmed_Shortest_Path']
        PPI_paths = PPI_shortest_path_dict[conversion_dict[key1]][conversion_dict[key2]]['Trimmed_Shortest_Path']
        
        for GI_path in GI_paths:
            for PPI_path in PPI_paths:
                GI_PPI_Overlap = list(set(GI_path).intersection(PPI_path))
                GI_PPI_Overlap = np.array([x for x in GI_PPI_Overlap.copy() if x != -1], dtype = np.int16)
                if len(GI_PPI_Overlap) >= 1:
                    shortest_path_overlap[key1][key2].append(GI_PPI_Overlap)
                
                if len(GI_PPI_Overlap) > 1:
                    print(GI_PPI_Overlap)

## Clean shortest path dictionary.
# 1. If there are no shortest paths between two modules, remove dict entry
# 2. If there are identical instances of shortest paths between two modules. Reduce to one instance.

shortest_path_overlap_cleaned = {key1: {key2: [] for key2 in GI_shortest_path_dict[key1].keys()} for key1 in GI_shortest_path_dict.keys()}

for key1 in shortest_path_overlap.keys():
    #print(key1)
    
    # Skip if there are no entries for this GI clusters (case of 112th cluster).
    if len(shortest_path_overlap[key1]) == 0:
        del shortest_path_overlap_cleaned[key1]
        continue
    
    # If this key contains no shortest paths at all, remove all its nested entries.
    # Note: length > 0 is necessary to avoid error from 112th cluster which is completely empty.
    if (len(reduce(lambda x, y: x+y, shortest_path_overlap[key1].values())) == 0):
        del shortest_path_overlap_cleaned[key1]
        continue
    
    for key2 in shortest_path_overlap[key1].keys():
        if len(shortest_path_overlap[key1][key2]) == 0:
            del shortest_path_overlap_cleaned[key1][key2]
        
        #elif len(shortest_path_overlap[key1][key2]) == 1:
        #    shortest_path_overlap_cleaned[key1][key2] = shortest_path_overlap[key1][key2].copy()
        # Note: If entry is == 1. Then array stays the same because no duplicates are possible.
        
        elif len(shortest_path_overlap[key1][key2]) >= 1:
            shortest_path_overlap_cleaned[key1][key2] = np.unique(shortest_path_overlap[key1][key2].copy())

all_connectors = []
for key1 in shortest_path_overlap_cleaned.keys():
    for key2 in shortest_path_overlap_cleaned[key1].keys():
        if len(shortest_path_overlap_cleaned[key1][key2]) > 0:
            all_connectors.extend(shortest_path_overlap_cleaned[key1][key2])

all_unique_connectors = list(np.unique(all_connectors.copy()))


#### Convert Clean OVerlap to Dataframe ####

## Make a list of lists containing shortest path overlap data
shortest_path_data_ls = []
for key1 in shortest_path_overlap_cleaned.keys():
    for key2 in shortest_path_overlap_cleaned[key1].keys():
        for path in shortest_path_overlap_cleaned[key1][key2]:
            shortest_path_data_ls.append([key1, key2, path.copy()])

     
## Convert list of lists to dataframe
Shortest_Path_Overlap_DF = pd.DataFrame(shortest_path_data_ls, columns = ['Source_GI_Cluster', 'Target_GI_Cluster', 'Shortest_Path_Overlap'])

## Add Protein Complex IDs
Shortest_Path_Overlap_DF['Source_CPX_ID'] = Shortest_Path_Overlap_DF['Source_GI_Cluster'].apply(lambda x: conversion_dict[x])
Shortest_Path_Overlap_DF['Target_CPX_ID'] = Shortest_Path_Overlap_DF['Target_GI_Cluster'].apply(lambda x: conversion_dict[x])

# Convert ORF integer IDs to ORF systematic IDs (characters)
Shortest_Path_Overlap_DF['Shortest_Path_Overlap'] = Shortest_Path_Overlap_DF['Shortest_Path_Overlap'].apply(lambda x: LUT[x])


## Export
Shortest_Path_Overlap_DF.to_csv("../results/Topological_Analyses/Shortest_Path_Connector_Overlap.csv")
