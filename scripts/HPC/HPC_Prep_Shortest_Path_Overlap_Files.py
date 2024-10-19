# Xavier Castellanos-Girouard

# Date First Created: June 19 2024
# Date Last Modified: June 19 2024


#### Libraries ####
import pandas as pd
import numpy as np
import networkx as nx
import pickle
import gc


#### Define Functions ####

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



def Convert_ListString_to_List(pandas_series):
    pandas_series = pandas_series.str.replace("[", "", regex=False).str.replace("]", "", regex=False).str.replace("'", "", regex=False).str.replace("'", "", regex=False).str.replace(" ", "", regex=False).str.split(",")
    return(pandas_series)




def sort_PPI_loop(PPI_DF, module_dict, locus_col, complex_col):
    
    # Pair ORF to list of complexes it is part of:
    Interactor_assign = zip(PPI_DF[locus_col], PPI_DF[complex_col])
    
    for ORF, CPX_ls in Interactor_assign: # For every ORF and list of complexes
        for CPX in CPX_ls: # For every complex in list of complexes
            if (ORF not in module_dict[CPX]): # If the Complexes exists (i.e. ORF is in a complex) and ORF is not in dict
                module_dict[CPX].append(ORF) # Assign ORF to complex in dictionnary



def import_PPI_network(main_path, file_path):
    
    modular_PPI_network_DF = pd.read_csv((main_path + file_path), index_col = 0)
    
    modular_PPI_network_DF = modular_PPI_network_DF.replace({np.nan: "None"})
    
    # Convert PPI Complex strings to list
    modular_PPI_network_DF['source_Complex'] = modular_PPI_network_DF['source_Complex'].str.split("|")
    modular_PPI_network_DF['target_Complex'] = modular_PPI_network_DF['target_Complex'].str.split("|")
    
    # Reset index
    modular_PPI_network_DF = modular_PPI_network_DF.reset_index(drop = True)
    
    return(modular_PPI_network_DF)



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


def import_module_overlap_table(main_path, file_path):
    overlap_table = pd.read_csv((main_path + file_path), index_col = 0)
    
    overlap_table = overlap_table.replace({np.nan: "None"})
    
    # Convert to proper datatypes
    overlap_table['GI_Module_ORFs'] = Convert_ListString_to_List(overlap_table['GI_Module_ORFs'])
    overlap_table['PPI_Module_ORFs'] = Convert_ListString_to_List(overlap_table['PPI_Module_ORFs'])
    
    return(overlap_table)



#### Main ####


def main():
    ## Main Path used for the whole project.
    main_path = "/home/xaviercg/scratch/"
    
    ## Import LookUp Table
    file_path = "data/encode_LUT_dict.pickle"
    with open((main_path+file_path), "rb") as handle:
        ORF_LUT_dict = pickle.load(handle)
    
    
    modular_PPI_network_DF = import_PPI_network(main_path, file_path = "data/Yeast_PPI_Network_CPX.csv")
    
    ## Convert yeast ORFs to Integer IDs using Lookup Table
    modular_PPI_network_DF['source_locus'] = modular_PPI_network_DF['source_locus'].apply(lambda x: ORF_LUT_dict[x])
    modular_PPI_network_DF['target_locus'] = modular_PPI_network_DF['target_locus'].apply(lambda x: ORF_LUT_dict[x])
    
    ## Sort ORFs into modules
    PPI_module_dict = sort_PPI_ORFs(modular_PPI_network_DF)
    
    ## Import Module Overlap Table
    GI_PPI_optimal_module_overlap_DF = import_module_overlap_table(main_path,
        file_path = "data/GI_PPI_optimal_module_overlap_clean.csv")
    
    
    ## Import Negative GI Network
    print("Importing Negative GI Network...")
    GI_negEpsi_network_G = import_GI_network(main_path, 
        file_path = "data/costanzo_2016_longer_withoutReps.csv",
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
        file_path = "data/Costanzo2016_DataFileS6.xlsx",
        GI_ORF_arr = GI_ORF_arr,
        ORF_LUT_dict = ORF_LUT_dict)
    
    print("Exporting...")
    
    # Export as pickle
    with open(main_path+"data/GI_negEpsi_network_G.pickle", "wb") as handle:
        pickle.dump(GI_negEpsi_network_G, handle)
        
    with open(main_path+"data/GI_module_dict.pickle", "wb") as handle:
        pickle.dump(GI_module_dict, handle)
    
    with open(main_path+"data/PPI_module_dict.pickle", "wb") as handle:
        pickle.dump(PPI_module_dict, handle)
    
    with open(main_path+"data/GI_PPI_optimal_module_overlap_clean_DF.pickle", "wb") as handle:
        pickle.dump(GI_PPI_optimal_module_overlap_DF, handle)
    
    print("Done")

if __name__ == "__main__":
    main()


