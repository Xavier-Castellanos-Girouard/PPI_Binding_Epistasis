# Xavier Castellanos-Girouard

# Date First Created: June 11 2024
# Date Last Modified: June 11 2024

#### Import Libraries ####
import pandas as pd
import numpy as np
import networkx as nx
import pickle
from tqdm import tqdm
from functools import reduce
from itertools import combinations
import gc
import time


def Convert_ListString_to_List(pandas_series):
    pandas_series = pandas_series.str.replace("[", "", regex=False).str.replace("]", "", regex=False).str.replace("'", "", regex=False).str.replace("'", "", regex=False).str.replace(" ", "", regex=False).str.split(",")
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


def import_module_overlap_table(main_path, file_path):
    overlap_table = pd.read_csv((main_path + file_path), index_col = 0)
    
    # Convert to proper datatypes
    overlap_table['GI_Module_ORFs'] = Convert_ListString_to_List(overlap_table['GI_Module_ORFs'])
    overlap_table['PPI_Module_ORFs'] = Convert_ListString_to_List(overlap_table['PPI_Module_ORFs'])
    
    return(overlap_table)



def import_PPI_network(main_path, file_path):
    
    modular_PPI_network_DF = pd.read_csv((main_path + file_path), index_col = 0)
    
    modular_PPI_network_DF = modular_PPI_network_DF.replace({np.nan: "None"})
    
    # Convert PPI Complex strings to list
    modular_PPI_network_DF['source_Complex'] = modular_PPI_network_DF['source_Complex'].str.split("|")
    modular_PPI_network_DF['target_Complex'] = modular_PPI_network_DF['target_Complex'].str.split("|")
    
    # Reset index
    modular_PPI_network_DF = modular_PPI_network_DF.reset_index(drop = True)
    
    return(modular_PPI_network_DF)
    

def sort_PPI_loop(PPI_DF, module_dict, locus_col, complex_col):
    
    # Pair ORF to list of complexes it is part of:
    Interactor_assign = zip(PPI_DF[locus_col], PPI_DF[complex_col])
    
    for ORF, CPX_ls in Interactor_assign: # For every ORF and list of complexes
        for CPX in CPX_ls: # For every complex in list of complexes
            if (ORF not in module_dict[CPX]): # If the Complexes exists (i.e. ORF is in a complex) and ORF is not in dict
                module_dict[CPX].append(ORF) # Assign ORF to complex in dictionnary


def sort_PPI_ORFs(Modular_PPI_Network_DF):
    ## Get list of PPI modules
    
    # Initiate list:
    PPI_module_list = []
    
    # Flatten list from interactorA complexes, and add to list
    PPI_module_list.extend([CPX for sublist in Modular_PPI_Network_DF["source_Complex"] for CPX in sublist])
    
    # Flatten list from interactorB complexes, and add to list
    PPI_module_list.extend([CPX for sublist in Modular_PPI_Network_DF["target_Complex"] for CPX in sublist])
    
    # Get unique list
    PPI_module_list = list(np.unique(PPI_module_list))
    
    # Make dictionary keys from list items
    PPI_module_dict = {key: [] for key in PPI_module_list}
    
    # * Note: 'None' value is purpousfully kept in dict, it serves as a list of connectors
    
    ## Sort ORFs from source column
    sort_PPI_loop(Modular_PPI_Network_DF, PPI_module_dict, "source_locus", "source_Complex")
    
    ## Sort ORFs from interactorB column
    sort_PPI_loop(Modular_PPI_Network_DF, PPI_module_dict, "target_locus", "target_Complex")
    
    ## Convert lists to numpy arrays
    for key in PPI_module_dict.keys():
        PPI_module_dict[key] = np.array(PPI_module_dict[key], dtype = np.int16)
    
    return(PPI_module_dict)


def trim_shortest_path(shortest_path, PPI_module_dict, module_PPI_arr):
    
    for node in shortest_path.copy(): # For every ORF in Shortest Path
        if node in module_PPI_arr: # If source ORF is in a module
            shortest_path.remove(node) # Remove ORF from shortest path
    
    return(shortest_path)


def get_PPI_module_pair_connectors(shortest_path_dict, 
    source_module_ORFs, 
    target_module_ORFs,
    source_module_ID,
    target_module_ID,
    PPI_module_dict,
    PPI_module_arr,
    Modular_PPI_Network_G):
    
    for source_ORF in source_module_ORFs: # For all nodes in source module
        for target_ORF in target_module_ORFs: # For all nodes in target module
            
            
            try:
                # Find shortest path between source and target node. 
                # Returns generator object
                shortest_paths_gen = nx.all_shortest_paths(Modular_PPI_Network_G, 
                    source_ORF, 
                    target_ORF)
                    
                for PPI_Shortest_path_ls in shortest_paths_gen: # Iterate through generator object
                    
                    shortest_path_trimmed = trim_shortest_path(PPI_Shortest_path_ls,
                        PPI_module_dict,
                        PPI_module_arr)
                    
                    # Only keep non-empty shortest paths
                    if len(shortest_path_trimmed) > 0:
                        # Store shortest path information. 
                        # Save Shortest Paths as array to save memory space.
                        shortest_path_dict[source_module_ID][target_module_ID]['Source_Node'].append(source_ORF)
                        shortest_path_dict[source_module_ID][target_module_ID]['Target_Node'].append(target_ORF)
                        shortest_path_dict[source_module_ID][target_module_ID]['Trimmed_Shortest_Path'].append(shortest_path_trimmed)
                
            except:
                pass
                
    # Optimize Shortest Paths for Memory
    shortest_path_dict[source_module_ID][target_module_ID]['Source_Node'] = np.array(shortest_path_dict[source_module_ID][target_module_ID]['Source_Node'], dtype = np.int16)
    shortest_path_dict[source_module_ID][target_module_ID]['Target_Node'] = np.array(shortest_path_dict[source_module_ID][target_module_ID]['Target_Node'], dtype = np.int16)
    shortest_path_dict[source_module_ID][target_module_ID]
                
                ## Make an 2D array from list. Fill empty spaces with "-1".
    shortest_path_dict[source_module_ID][target_module_ID]['Trimmed_Shortest_Path'] = boolean_indexing(shortest_path_dict[source_module_ID][target_module_ID]['Trimmed_Shortest_Path'],
        fillval = np.int16(-1),
        dtype = np.int16)


def find_PPI_SP_connectors(module_overlap_DF, 
    PPI_module_dict, 
    PPI_module_arr, 
    Modular_PPI_Network_G):
    
    ## Get pairwise combination of PPI Modules and their ORFs
    # Note: Here only complexes that overlap with GI modules are looked into.
    #       There is no need to explore the rest.
    #       [:-1] is used to exclude "None" module (containing connectors)
    overlap_complexes_arr = np.array([CPX for CPX in module_overlap_DF['CPX_ID'].tolist() if CPX != "None"])
    #print(PPI_module_dict[overlap_complexes_arr[0]], overlap_complexes_arr[0])
    overlap_complex_ORFs_ls = [PPI_module_dict[CPX] for CPX in overlap_complexes_arr]
    PPI_Module_ORFs_zip = zip(overlap_complexes_arr, overlap_complex_ORFs_ls)
    
    ## Get all combinations of protein modules
    PPI_Module_ORFs_list = list(PPI_Module_ORFs_zip)
    PPI_Module_ORFs_pairwise_comb = list(combinations(PPI_Module_ORFs_list, 2))
    
    
    # Make a dictionary to conserve shortest paths:
    PPI_modules_unique_arr = np.unique(overlap_complexes_arr)

    # NOTE: SHORTEST PATHS SHOULD BE APPENDED TO LIST IN KEY2 FOR PPI NETWORKS. This is because there are multiple shortest paths
    PPI_pairwise_module_shortest_paths_dict = {
        key1: {
            key2: {'Source_Node':[], 
                'Target_Node':[], 
                'Trimmed_Shortest_Path':[]
                } for key2 in PPI_modules_unique_arr
            } for key1 in PPI_modules_unique_arr
        }
    
    
    for pair in tqdm(PPI_Module_ORFs_pairwise_comb): # For every pair of PPI modules
        source_module_ID = pair[0][0] # Retrieve Source Module ID
        target_module_ID = pair[1][0] # Retrieve Target Module ID
        
        
        source_module_ORFs = pair[0][1] # Retrieve Source Module ORFs
        target_module_ORFs = pair[1][1] # Retrieve Target Module ORFs
        
        #print(source_module_ID, " : ", target_module_ID)
        #print(source_module_ORFs, " : ", target_module_ORFs)
        
        get_PPI_module_pair_connectors(PPI_pairwise_module_shortest_paths_dict,
            source_module_ORFs,
            target_module_ORFs,
            source_module_ID,
            target_module_ID,
            PPI_module_dict,
            PPI_module_arr,
            Modular_PPI_Network_G)
        
    return(PPI_pairwise_module_shortest_paths_dict)




### Main ####
def main():
    
    ## Main Path used for the whole project.
    main_path = "/home/xavier/Desktop/Cell_interactome_stoichiometry/Yeast_Epistasis_Stoichiometries/Third_Run/"
    
    ## Import Lookup Table
    file_path = "Python_Connector_Overlap/results/encode_LUT_dict.pickle"
    with open((main_path+file_path), "rb") as handle:
        ORF_LUT_dict = pickle.load(handle)
    
    ## Import Module Overlap Table
    GI_PPI_optimal_module_overlap_DF = import_module_overlap_table(main_path,
        file_path = "Python_Module_Overlap/results/GI_PPI_optimal_module_overlap_clean.csv")
    
    ## Import PPI network dataframe
    modular_PPI_network_DF = import_PPI_network(main_path, 
        file_path = "Python_Module_Overlap/results/Yeast_PPI_Network_CPX.csv")
    
    ## Convert yeast ORFs to Integer IDs using Lookup Table
    modular_PPI_network_DF['source_locus'] = modular_PPI_network_DF['source_locus'].apply(lambda x: ORF_LUT_dict[x])
    modular_PPI_network_DF['target_locus'] = modular_PPI_network_DF['target_locus'].apply(lambda x: ORF_LUT_dict[x])
    
    ## Sort ORFs into modules
    PPI_module_dict = sort_PPI_ORFs(modular_PPI_network_DF)
    
    ## Import PPI Network
    print("Constructing PPI Network...")
    PPI_network_G = nx.from_pandas_edgelist(modular_PPI_network_DF,
        source = "source_locus",
        target = "target_locus")
    
    ## Free Memory
    del modular_PPI_network_DF
    gc.collect()
    
    ## Get list of all ORFs in PPI network
    #all_PPI_ORFs = np.array(list(PPI_network_G.nodes), dtype = np.int16)
    
    ## Get list of all module nodes
    PPI_module_arr = np.unique([x for key in PPI_module_dict.keys() for x in PPI_module_dict[key] if key != "None"])
    
    ## Get Shortest paths
    PPI_pairwise_module_shortest_paths_dict = find_PPI_SP_connectors(
        GI_PPI_optimal_module_overlap_DF,
        PPI_module_dict, 
        PPI_module_arr,
        PPI_network_G)
    
    
    print("Exporting...")
    ## Export dictionary as pickle file
    file_path = "Python_Connector_Overlap/results/PPI_pairwise_module_shortest_paths_dict.pickle"
    with open(main_path + file_path, "wb") as handle:
        pickle.dump(PPI_pairwise_module_shortest_paths_dict, handle)
    
    print("Done...")


if __name__ == "__main__":
    main()
