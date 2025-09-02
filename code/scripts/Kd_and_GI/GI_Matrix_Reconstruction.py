# Xavier Castellanos-Girouard

# Reconstruction of GI matrix using PPIs

# Date First Created: March 17 2024
# Date Last Modified: May 28 2024


#### Import Libraries ####

import numpy as np
import pandas as pd
import networkx as nx
#import seaborn as sns
import scipy.optimize as optimization
import scipy as sp
#import matplotlib.pyplot as plt
#from matplotlib.colors import LogNorm, Normalize
from itertools import combinations
from itertools import repeat


#### Import Data ####

Kd_GI_Network_DF = pd.read_csv("../results/Kd_and_GI/Yeast_Kd_GI.csv", index_col=0)


#### Find the Optimal Number of Quantiles for Distance vs. GI Fit ####

def exp_func(dist, a, r): # Exponential fit
    exponential = a*np.exp(r*dist)
    return exponential

def lin_func(dist, a): # Linear fit
    lin = a*dist
    return lin

# Define function that will serve as basis for optimization
def OptiQuantile(Binding_GI_DF, GI_col, Dist_col, denom_exp, nquant):
    
    # Calculate an Association score
    Binding_GI_DF['Association'] = 1/Binding_GI_DF[Dist_col]**denom_exp
    
    # Sort by Distance
    Binding_GI_DF = Binding_GI_DF.sort_values(by = ['Association'], inplace = False)
    
    # Assign Deciles
    Binding_GI_DF['AssociationQuantile']= pd.qcut(Binding_GI_DF['Association'], q = nquant, labels = False)
    Binding_GI_DF['AssociationQuantile'] = Binding_GI_DF['AssociationQuantile'] + 1 # Quantiles will now start at 1
    
    
    # Get mean distance and mean scores per decile
    Binding_GI_Quantile_DF = Binding_GI_DF.groupby(['AssociationQuantile']).agg({'Association': 'mean', GI_col: 'mean'})
    
    ## Fit curve relating distance to GI (exponential fit)
    popt, pcov = optimization.curve_fit(exp_func, Binding_GI_Quantile_DF['Association'], abs(Binding_GI_Quantile_DF[GI_col]), p0=[0.1, 0.1])
    
    # Calculate Root Mean Square Errors (RMSE)
    expoFit_RMSE = np.mean((abs(Binding_GI_Quantile_DF[GI_col]) - exp_func(Binding_GI_Quantile_DF['Association'], *popt)**2))
    
    
    ## Fit curve relating distance to GI (exponential fit)
    popt, pcov = optimization.curve_fit(lin_func, Binding_GI_Quantile_DF['Association'], abs(Binding_GI_Quantile_DF[GI_col]), p0=[0.1])
    
    # Calculate Root Mean Square Errors (RMSE)
    linFit_RMSE = np.mean((abs(Binding_GI_Quantile_DF[GI_col]) - lin_func(Binding_GI_Quantile_DF['Association'], *popt)**2))
    
    return([expoFit_RMSE, linFit_RMSE])

## Set distances for shortest paths
# Logarithmize the distances. Stronger the interaction, lower the distance
Kd_GI_Network_DF['Distance'] = 1/(-np.log10(Kd_GI_Network_DF['Kd']))

## Only keep negative genetic interactions
Kd_NegGI_DF = Kd_GI_Network_DF[Kd_GI_Network_DF['scores'] < 0].copy()

# Association is the 
#Kd_NegGI_DF['Association'] = 1/Kd_NegGI_DF['Distance']

## Compute fits for different numbers of quantiles
Opti_RMSE_ls = list(map(OptiQuantile, repeat(Kd_GI_Network_DF[Kd_GI_Network_DF['scores'] < 0].copy()), repeat('scores'), repeat('Distance'), repeat(1), list(range(5, 1001))))

## Make list of values for Exponential fit and Linear fit
exp_RMSE_ls = [x[0] for x in Opti_RMSE_ls]
lin_RMSE_ls = [x[1] for x in Opti_RMSE_ls]

## View Distribution of RMSE
#sns.scatterplot(x = list(range(5, 505)), y = lin_RMSE_ls[0:500])
## Note: Graph suggests a 25 quantiles is optimal

## Export Fitting data

GI_Kd_distance_fit_DF = pd.DataFrame({'exponentialFit_RMSE': exp_RMSE_ls, 'linearFit_RMSE': lin_RMSE_ls, 'quantiles': list(range(5, 1001))})
GI_Kd_distance_fit_DF.to_csv("../results/Kd_and_GI/GI_Kd_distance_fit.csv")


#### Retrieve Optimal Equation ####

## NOTE: Exponential Fit Seems to work well. A local minima is observed around q = 25. Use this value.


## Get equation relating distance to negative epistasis
# Sort by Distance
Kd_NegGI_DF['Association'] = 1/Kd_NegGI_DF['Distance']
Kd_NegGI_DF = Kd_NegGI_DF.sort_values(by = ['Association'], inplace = False)

# Assign Deciles
Kd_NegGI_DF['AssociationQuantile']= pd.qcut(Kd_NegGI_DF['Association'], q = 25, labels = False)
Kd_NegGI_DF['AssociationQuantile'] = Kd_NegGI_DF['AssociationQuantile'] + 1 # Deciles will now start at 1

# Get mean distance and mean scores per decile
Kd_NegGI_Quantile_DF = Kd_NegGI_DF.groupby(['AssociationQuantile']).agg({'Association': 'mean', 'scores': 'mean'})

# Export Binned Table
Kd_NegGI_Quantile_DF.to_csv("../results/Kd_and_GI/Quantile_Dist_GI_Table.csv")

popt, pcov = optimization.curve_fit(exp_func, Kd_NegGI_Quantile_DF['Association'], abs(Kd_NegGI_Quantile_DF['scores']), p0=[0.1, 0.1])
#plt.plot(Kd_NegGI_Quantile_DF['Association'], abs(Kd_NegGI_Quantile_DF['scores']), 'o')
#plt.plot(Kd_NegGI_Quantile_DF['Association'], exp_func(Kd_NegGI_Quantile_DF['Association'], *popt), '-')
#plt.show()
#print(popt[0], popt[1])

#print(np.mean((abs(Kd_NegGI_Quantile_DF['scores']) - exp_func(Kd_NegGI_Quantile_DF['Association'], *popt)**2)))


#### Construct Network and Extract Shortest Paths ####

## Construct Weighted PPI Network
Kd_GI_Network_G = nx.from_pandas_edgelist(Kd_GI_Network_DF, "source", "target", ["Distance"])

## Get all pair shortest paths lengths
Kd_GI_APSP_len_it = nx.shortest_path_length(Kd_GI_Network_G, source=None, target=None, weight='Distance', method='dijkstra')

# Extract adjacencies from adjacency list generator function
# Note: adjacency[0] will be a source node, adjacency[1] will be a dictionary with target nodes as keys and weight as entry
Kd_GI_APSP_len_dict = {adjacency[0]: adjacency[1] for adjacency in Kd_GI_APSP_len_it}

## Make 2D array representing all pairwise shortest path distances
node_list = np.array([key for key in Kd_GI_APSP_len_dict.keys()])

# Note: Data in the array will be ordered according to the node list
#           i.e., the ith element of node list will be the ith and jth in the array (where i==j)

## Initiate 2D array
Kd_GI_APSP_len_arr = np.empty(shape=(len(node_list), len(node_list)))
Kd_GI_APSP_len_arr.fill(np.nan) # Fill with distances of 1 (highest; weakest interaction)


## Fill Array with shortest path distances
for i, ORF in enumerate(node_list):
    adjacency_dict = Kd_GI_APSP_len_dict[ORF]
    
    for key in adjacency_dict.keys():
        j = np.where(node_list == key)[0][0]
        Kd_GI_APSP_len_arr[i,j] = adjacency_dict[key]
        Kd_GI_APSP_len_arr[j,i] = adjacency_dict[key]
        
#### Make Reconstructed GI Matrix Using Shortest Paths Distances ####

## Convert to similarity matrix (1/d; 'mod distance') and apply fitted equation
Kd_GI_Sim_arr = 1/Kd_GI_APSP_len_arr
Kd_GI_inferred_arr = popt[0]*np.exp(popt[1]*Kd_GI_Sim_arr)
Kd_GI_inferred_arr[Kd_GI_inferred_arr == np.inf] = np.nan


## Conserve Reconstructed Matrix as Dataframe
Kd_GI_inferred_arr_DF = pd.DataFrame(Kd_GI_inferred_arr)
Kd_GI_inferred_arr_DF.columns = node_list
Kd_GI_inferred_arr_DF.index = node_list
Kd_GI_inferred_arr_DF.index.name = "Gene"

## Export
Kd_GI_inferred_arr_DF.to_csv("../results/IntermediateFiles/Kd_inferred_GI.tsv", sep = "\t")
