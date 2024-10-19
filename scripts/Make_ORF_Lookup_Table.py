# Xavier Castellanos-Girouard

# Date First Created: June 11 2024
# Date Last Modified: June 11 2024

#### Import Libraries ####
import pandas as pd
import numpy as np
import pickle
  
### Main ####
def main():
    
    main_path = "/home/xavier/Desktop/Cell_interactome_stoichiometry/Yeast_Epistasis_Stoichiometries/Third_Run/"
    
    
    # Import GI Module Table
    file_path = "Python_Module_Overlap/data/Costanzo2016_DataFileS6.xlsx"
    GI_modules_DF = pd.read_excel((main_path + file_path))
    
    # Assign proper names to columns
    GI_modules_DF.columns = list(GI_modules_DF.iloc[0].values)

    # Remove row containing column names. Reset index
    GI_modules_DF = GI_modules_DF.drop(index = 0, axis = 0, inplace = False).reset_index(drop = True)
    
    
    # Import GI network table
    file_path = "Python_GI_Network_Formatting/data/Costanzo_GI_processed/costanzo_2016_longer_withoutReps.csv"
    EpsiGI_network_DF = pd.read_csv((main_path + file_path))
    
    
    # Import PPI network table
    file_path = "Python_Module_Overlap/results/Yeast_PPI_Network_CPX.csv"
    modular_PPI_network_DF = pd.read_csv((main_path+file_path), index_col = 0)
    
    # Make an array of all unique ORFs from different sources
    GI_ORF_list = []
    GI_ORF_list.extend(np.unique(EpsiGI_network_DF['ORF_query']))
    GI_ORF_list.extend(np.unique(EpsiGI_network_DF['ORF_array']))
    GI_ORF_list.extend(GI_modules_DF['Systematic array ORF  name'])
    GI_ORF_list.extend(modular_PPI_network_DF['source_locus'])
    GI_ORF_list.extend(modular_PPI_network_DF['target_locus'])
    GI_ORF_arr = np.array(np.unique([GI_ORF_list]), dtype = "<U9")
    
    ## Sort Alphabetically
    GI_ORF_sorted_arr = np.sort(GI_ORF_arr)
    
    LUT_indexes = np.arange(len(GI_ORF_sorted_arr), dtype = np.int16)
    
    # ORF 2 index
    encode_LUT_dict = {GI_ORF_sorted_arr[key]: key for key in LUT_indexes}
    
    # index 2 ORF
    decode_LUT_dict = {key: GI_ORF_sorted_arr[key] for key in LUT_indexes}
    
    file_path = "Python_Connector_Overlap/results/encode_LUT_dict.pickle"
    with open((main_path+file_path), "wb") as handle:
        pickle.dump(encode_LUT_dict, handle)
    
    file_path = "Python_Connector_Overlap/results/decode_LUT_dict.pickle"
    with open((main_path+file_path), "wb") as handle:
        pickle.dump(decode_LUT_dict, handle)
    
if __name__ == "__main__":
    main()
