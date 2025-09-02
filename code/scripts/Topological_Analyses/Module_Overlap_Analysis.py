# Xavier Castellanos-Girouard

# Date First Created: Oct 22 2023
# Date Last Modified: Jun 26 2024

#### Import Libraries ####

import numpy as np
import pandas as pd
#import seaborn as sns
import matplotlib.pyplot as plt
import networkx as nx
import pickle
from itertools import repeat
from functools import reduce
from itertools import combinations
from itertools import product
from itertools import chain
import random
from tqdm import tqdm
import scipy
import sys
import statistics
import time


#### Define Functions ####

# Calculate Jaccard index for two sets of elements
def Jaccard(a, b):
    set_a = set(a)
    set_b = set(b)
    intersect_ab = set_a.intersection(set_b)
    union_ab = set_a.union(set_b)
    
    try:
        Jaccard_index = len(intersect_ab)/len(union_ab)
    
    except:
        Jaccard_index = np.nan
    
    return(Jaccard_index)



#### Import and Format Data ####

### PPI Network Dataframes
# Import PPI network
Modular_PPI_Network_DF = pd.read_csv("../results/Topological_Analyses/Yeast_PPI_Network_CPX.csv", index_col = 0)

# Convert PPI Complex strings to list
Modular_PPI_Network_DF['source_Complex'] = Modular_PPI_Network_DF['source_Complex'].str.split("|")
Modular_PPI_Network_DF['target_Complex'] = Modular_PPI_Network_DF['target_Complex'].str.split("|")

# Reset index
Modular_PPI_Network_DF = Modular_PPI_Network_DF.reset_index(drop = True)

## GI Network Dataframes
GI_Epsi_network_DF = pd.read_csv("../results/Stoichiometry_and_GI/costanzo_2016_longer_withoutReps.csv")


### Import other network attributes
GI_Modules = pd.read_excel("../data/Modules/Costanzo2016_DataFileS6.xlsx")

# Assign proper names to columns
GI_Modules.columns = list(GI_Modules.iloc[0].values)

# Remove row containing column names. Reset index
GI_Modules = GI_Modules.drop(index = 0, axis = 0, inplace = False).reset_index(drop = True)

# Remove duplicates
GI_Modules = GI_Modules[~GI_Modules.duplicated(['Systematic array ORF  name', 'Pathway/complex level (PCC > 0.4)'])]
GI_Modules = GI_Modules[~GI_Modules.duplicated(['Systematic array ORF  name'])].reset_index(drop = True)



#### Sort PPI ORFs into modules ####


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

PPI_module_dict

# * Note: 'None' value is purpousfully kept in dict, it serves as a list of connectors


## Sort ORFs from interactorA column

# Pair ORF to list of complexes it is part of:
InteractorA_assign = zip(Modular_PPI_Network_DF["source_locus"], Modular_PPI_Network_DF["source_Complex"])

for ORF, CPX_ls in InteractorA_assign: # For every ORF and list of complexes
    for CPX in CPX_ls: # For every complex in list of complexes
        if (ORF not in PPI_module_dict[CPX]): # If the Complexes exists (i.e. ORF is in a complex)
            PPI_module_dict[CPX].append(ORF) # Assign ORF to complex in dictionnary

## Sort ORFs from interactorB column

# Pair ORF to list of complexes it is part of:
InteractorB_assign = zip(Modular_PPI_Network_DF["target_locus"], Modular_PPI_Network_DF["target_Complex"])

for ORF, CPX_ls in InteractorB_assign: # For every ORF and list of complexes
    for CPX in CPX_ls: # For every complex in list of complexes
        if (ORF not in PPI_module_dict[CPX]): # If the Complexes exists (i.e. ORF is in a complex) and ORF is not in dict
            PPI_module_dict[CPX].append(ORF) # Assign ORF to complex in dictionnary
            


#### Sort GI ORFs into clusters/modules ####

# Add ORF from correlation network
GI_ORF_list = []
GI_ORF_list.extend(GI_Epsi_network_DF['ORF_query'])
GI_ORF_list.extend(GI_Epsi_network_DF['ORF_array'])

# Make list unique
GI_ORF_list = list(np.unique(GI_ORF_list))


# Get list for only EpsiGI
all_EpsiGI_ORFs = []
all_EpsiGI_ORFs.extend(GI_Epsi_network_DF['ORF_query'])
all_EpsiGI_ORFs.extend(GI_Epsi_network_DF['ORF_array'])

all_EpsiGI_ORFs = list(np.unique(all_EpsiGI_ORFs))


## Get list of GI clusters/modules

# Get modules from module dataframe. Make list unique
GI_Modules_list = list(np.unique(GI_Modules['Pathway/complex level (PCC > 0.4)']))

# Make empty dictionary with keys as cluster number
GI_module_dict = {key: [] for key in GI_Modules_list}

# Make 'None' cluster to store connector ORFs
GI_module_dict['None'] = []


# Sort ORFs into dictionnary of GI modules
for ORF in GI_ORF_list:
    if ORF in list(GI_Modules['Systematic array ORF  name']): # If the ORF is a module
        
        # Find cluster to which the ORF belongs:
        cluster = GI_Modules[GI_Modules['Systematic array ORF  name'] == ORF]['Pathway/complex level (PCC > 0.4)'].values
        
        # Make sure ORF belongs to only one cluster
        assert(len(cluster) == 1)
        cluster = cluster[0] # Get value from array containing 1 item
        
        
        # If the ORF has not yet been sorted into this cluster:
        if ORF not in GI_module_dict[cluster]:
            GI_module_dict[cluster].append(ORF) # Add ORF to cluster in dict
        
        else:
            continue
        
    # If ORF is not in a cluster (not in a module):
    elif ORF not in list(GI_Modules['Systematic array ORF  name']):
        GI_module_dict['None'].append(ORF) # Add to non-module list in dict



#### Map GI modules to PPI modules ####

### Find Overlap between individual GI modules and PPI modules
GI_PPI_module_Overlap_matrix = pd.DataFrame(columns = list(PPI_module_dict.keys()), index = list(GI_module_dict.keys()))

# For every cell in matrix, calculate jaccard index between GI module (index) and PPI module (column)
for i in list(GI_module_dict.keys()):
    for j in list(PPI_module_dict.keys()):
        GI_PPI_module_Overlap_matrix.loc[i,j] = Jaccard(PPI_module_dict[j], GI_module_dict[i])

# Convert all values to numeric
for col_name in GI_PPI_module_Overlap_matrix.columns:
    GI_PPI_module_Overlap_matrix[col_name] = pd.to_numeric(GI_PPI_module_Overlap_matrix[col_name])


### Find best PPI complex for each GI cluster

# Get list of complexes, in the same order as columns of matrix
PPI_modules_matrix_list = GI_PPI_module_Overlap_matrix.columns.tolist()

# For every GI cluster, get the column index for which overlap is maximal
GI_PPI_optimal_overlap_DF = pd.DataFrame(GI_PPI_module_Overlap_matrix.apply(lambda x: x.tolist().index(max(x)), axis = 1), columns = ["CPX_index"])

# Get complex ID from index
GI_PPI_optimal_overlap_DF['CPX_ID'] = [PPI_modules_matrix_list[index] for index in GI_PPI_optimal_overlap_DF['CPX_index']]

# Make column for Jaccard index
GI_PPI_optimal_overlap_DF['Jaccard_index'] = GI_PPI_optimal_overlap_DF.apply(lambda x: GI_PPI_module_Overlap_matrix.loc[x.name, x['CPX_ID']], axis = 1)

# Make list of common ORF between GI and PPI modules for each pair
GI_PPI_optimal_overlap_DF['Common_Module_ORFs'] = [set(PPI_module_dict[CPX_ID]).intersection(set(GI_module_dict[cluster])) for CPX_ID, cluster in zip(GI_PPI_optimal_overlap_DF['CPX_ID'].tolist(), list(GI_PPI_optimal_overlap_DF.index))]



#### Add PPI and GI Modules ORFs ####

GI_PPI_optimal_overlap_DF['Cluster_ID'] = GI_PPI_optimal_overlap_DF.index
GI_PPI_optimal_overlap_DF = GI_PPI_optimal_overlap_DF.reset_index(drop = True)

# Add list of ORFs for GI modules
GI_PPI_optimal_overlap_DF['GI_Module_ORFs'] = [GI_module_dict[key_] for key_ in GI_PPI_optimal_overlap_DF['Cluster_ID']]

# Add list of ORFs for PPI modules
GI_PPI_optimal_overlap_DF['PPI_Module_ORFs'] = [PPI_module_dict[key_] for key_ in GI_PPI_optimal_overlap_DF['CPX_ID']]

## Export
GI_PPI_optimal_overlap_DF.to_csv('../results/Topological_Analyses/GI_PPI_optimal_module_overlap.csv')


GI_PPI_optimal_overlap_DF = GI_PPI_optimal_overlap_DF[GI_PPI_optimal_overlap_DF['Common_Module_ORFs'].apply(lambda x: len(x) > 0)]
GI_PPI_optimal_overlap_DF = GI_PPI_optimal_overlap_DF[GI_PPI_optimal_overlap_DF['CPX_ID'] != 'None']
GI_PPI_optimal_overlap_DF

## Export
GI_PPI_optimal_overlap_DF.to_csv('../results/Topological_Analyses/GI_PPI_optimal_module_overlap_clean.csv')



#### Randomization to get p-values for individual module pairings ####

Cluster_ID_iter_ls = []
Jaccard_index_iter_ls = []
pool_of_PPI_Module_ORFs_stable = [x for list_ in PPI_module_dict.values() for x in list_]


for i in tqdm(range(0,1000)):
    
    pool_of_PPI_Module_ORFs_random = pool_of_PPI_Module_ORFs_stable.copy()
    
    random.shuffle(pool_of_PPI_Module_ORFs_random)
    
    ## Map GI modules to PPI modules 
    ### Find Overlap between individual GI modules and PPI modules
    GI_PPI_module_Overlap_matrix_iter = pd.DataFrame(columns = list(PPI_module_dict.keys()), index = list(GI_module_dict.keys()))
    
    
    k = 0
    # For every cell in matrix, calculate jaccard index between GI module (index) and PPI module (column)
    for j in list(PPI_module_dict.keys()):
        #print(j)
        size = len(PPI_module_dict[j])
        random_PPI_module = pool_of_PPI_Module_ORFs_random[k:k+size]
        k = k+size
        for i in list(GI_module_dict.keys()):
            #print(i)
            GI_PPI_module_Overlap_matrix_iter.loc[i,j] = Jaccard(random_PPI_module, GI_module_dict[i])
            #print(k, size, random_PPI_module, Jaccard(random_PPI_module, GI_module_dict[i]))
        #time.sleep(60)
    
    # Convert all values to numeric
    for col_name in GI_PPI_module_Overlap_matrix_iter.columns:
        GI_PPI_module_Overlap_matrix_iter[col_name] = pd.to_numeric(GI_PPI_module_Overlap_matrix_iter[col_name])
        
    ### Find best PPI complex for each GI cluster
    
    # Get list of complexes, in the same order as columns of matrix
    PPI_modules_matrix_list_iter = GI_PPI_module_Overlap_matrix_iter.columns.tolist()
    
    # For every GI cluster, get the column index for which overlap is maximal
    GI_PPI_optimal_overlap_DF_iter = pd.DataFrame(GI_PPI_module_Overlap_matrix_iter.apply(lambda x: x.tolist().index(max(x)), axis = 1), columns = ["CPX_index"])
    
    # Get complex ID from index
    GI_PPI_optimal_overlap_DF_iter['CPX_ID'] = [PPI_modules_matrix_list_iter[index] for index in GI_PPI_optimal_overlap_DF_iter['CPX_index']]
    
    # Make column for Jaccard index
    GI_PPI_optimal_overlap_DF_iter['Jaccard_index'] = GI_PPI_optimal_overlap_DF_iter.apply(lambda x: GI_PPI_module_Overlap_matrix_iter.loc[x.name, x['CPX_ID']], axis = 1)
    
    # Make list of common ORF between GI and PPI modules for each pair
    #GI_PPI_optimal_overlap_DF_iter['Common_Module_ORFs'] = [set(PPI_module_dict[CPX_ID]).intersection(set(GI_module_dict[cluster])) for CPX_ID, cluster in zip(GI_PPI_optimal_overlap_DF_iter['CPX_ID'].tolist(), list(GI_PPI_optimal_overlap_DF_iter.index))]
    
    Cluster_ID_iter_ls.extend(GI_PPI_optimal_overlap_DF_iter.index.copy().tolist())
    Jaccard_index_iter_ls.extend(GI_PPI_optimal_overlap_DF_iter['Jaccard_index'].copy().tolist())
    

randomized_module_overlap_DF = pd.DataFrame({"Cluster_ID":Cluster_ID_iter_ls, "Jaccard_index":Jaccard_index_iter_ls})

## Export
randomized_module_overlap_DF.to_csv("../results/Topological_Analyses/GI_PPI_randomized_ModuleOverlap.csv")
