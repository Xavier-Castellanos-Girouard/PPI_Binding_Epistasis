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
    
    ## Sort ORFs from target column
    sort_PPI_loop(Modular_PPI_Network_DF, PPI_module_dict, "target_locus", "target_Complex")
    
    ## Convert lists to numpy arrays
    for key in PPI_module_dict.keys():
        PPI_module_dict[key] = np.array(PPI_module_dict[key], dtype = np.int16)
    
    return(PPI_module_dict)


    
def randomize_PPI_networks(main_path, file_path, n, ORF_LUT_dict):
    
    modular_PPI_network_DF = import_PPI_network(main_path, file_path)
    
    ## Convert yeast ORFs to Integer IDs using Lookup Table
    modular_PPI_network_DF['source_locus'] = modular_PPI_network_DF['source_locus'].apply(lambda x: ORF_LUT_dict[x])
    modular_PPI_network_DF['target_locus'] = modular_PPI_network_DF['target_locus'].apply(lambda x: ORF_LUT_dict[x])
    
    ## Sort ORFs into modules
    PPI_module_dict = sort_PPI_ORFs(modular_PPI_network_DF)
    
    ## Get dataframe of intra module interactions
    intramodule_PPI_network_DF = modular_PPI_network_DF[modular_PPI_network_DF['Interaction_Type']=='IntraModule'].copy()
    
    # Make graph from dataframe
    intramodule_PPI_network_G = nx.from_pandas_edgelist(intramodule_PPI_network_DF,
        source = 'source_locus', 
        target = 'target_locus')
    
    ## Get dataframe of all other interactions
    Connector_PPI_network_DF = modular_PPI_network_DF[(modular_PPI_network_DF['Interaction_Type']!='IntraModule')].copy()
    
    # Clear memory unused dataframes
    del intramodule_PPI_network_DF
    del modular_PPI_network_DF
    
    gc.collect()
    
    
    ## Create degree preserved randomized networks, where intra-module interactions are stable.
    for i in range(n):
        randomized_graph = intramodule_PPI_network_G.copy()
        randomized_connector_DF = Connector_PPI_network_DF.copy()
        
        ## Randomize edges 
        randomized_connector_DF['target_locus'] = randomized_connector_DF['target_locus'].sample(frac=1, replace = False).tolist()
        
        Connector_PPI_network_G = nx.from_pandas_edgelist(randomized_connector_DF,
            source = 'source_locus',
            target = 'target_locus')
            
        randomized_graph.add_edges_from([e for e in Connector_PPI_network_G.edges])
        
        with open((main_path+f"data/Controls/DegPres_Randomized_PPI_Networks/random_PPI_network_{i}_G.pickle"),'wb') as handle:
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
    
    ## Generate Random PPI networks
    randomize_PPI_networks(main_path = main_path, 
        file_path = "data/Yeast_PPI_Network_CPX.csv",
        n = 10000, 
        ORF_LUT_dict = ORF_LUT_dict)
    
    
    
if __name__ == "__main__":
    main()
