# Xavier Castellanos-Girouard
# Date First Created: June 4 2024
# Date Last Modified: June 4 2024


#### Import modules ####

import numpy as np
import pandas as pd
import pickle


#### Import and Concatenate DataFrames ####

# Import ID to ORF lookup table
with open(f"/home/xaviercg/scratch/data/decode_LUT_dict.pickle", "rb") as handle:
    decode_LUT_dict = pickle.load(handle) # Open pickle file containing DF

DataFrame_list = [] # Initiate list of dataframes

for i in range(0, 999, 1): # For all 1000 controls
    with open(f"/home/xaviercg/scratch/data/Controls/Random_Shortest_Path_Connector_Overlap/GI_PPI_pairwise_connector_overlap_{i}_DF.pickle", "rb") as handle:
        Current_DataFrame = pickle.load(handle) # Open pickle file containing DF
    
    Current_DataFrame['Shortest_Path_Overlap'] = Current_DataFrame['Shortest_Path_Overlap'].apply(lambda x: decode_LUT_dict[x])
    DataFrame_list.append(Current_DataFrame) # Add Dataframe to dataframe list

Master_DF = pd.concat(DataFrame_list) # Concatenate all dataframes

Master_DF.to_csv("/home/xaviercg/scratch/data/Controls/GI_PPI_ConnectorOverlap_Controls.csv")

