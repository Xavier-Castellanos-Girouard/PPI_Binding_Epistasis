# Xavier Castellanos-Girouard
# Date First Created: June 4 2024
# Date Last Modified: June 4 2024


#### Import modules ####

import numpy as np
import pandas as pd
import random
import pickle


#### Import and Concatenate DataFrames ####

DataFrame_list = [] # Initiate list of dataframes

for i in range(0, 10000, 1): # For all 10000 controls
    with open(f"/home/xaviercg/scratch/data/Controls/Random_GI_PPI_Module_Match/GI_PPI_ModuleMatch_{i}_DF.pickle", "rb") as file:
        Current_DataFrame = pickle.load(file) # Open pickle file containing DF
    Current_DataFrame['iter_ID'] = str(i)
    DataFrame_list.append(Current_DataFrame) # Add Dataframe to dataframe list

Master_DF = pd.concat(DataFrame_list) # Cancatenate all dataframes

Master_DF.to_csv("/home/xaviercg/scratch/data/Controls/GI_PPI_ModuleOverlap_Controls.csv")

