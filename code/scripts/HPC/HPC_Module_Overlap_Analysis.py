# Xavier Castellanos-Girouard
# Adapted from a Script First Created in: Oct 22 2023
# Date Last Modified: June 3 2024


#### Import modules ####

import numpy as np
import pandas as pd
from itertools import product
from itertools import chain
import random
import sys
import statistics
import pickle

random.seed(sys.argv[1])

ctrl_number = sys.argv[1]

#### Functions ####

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
Modular_PPI_Network_DF = pd.read_csv("/home/xaviercg/scratch/data/Yeast_PPI_Network_CPX.csv", index_col = 0)

Modular_PPI_Network_DF = Modular_PPI_Network_DF.replace({np.nan: "None"})


# Convert PPI Complex strings to list
Modular_PPI_Network_DF['source_Complex'] = Modular_PPI_Network_DF['source_Complex'].str.split("|")
Modular_PPI_Network_DF['target_Complex'] = Modular_PPI_Network_DF['target_Complex'].str.split("|")


# Reset index
Modular_PPI_Network_DF = Modular_PPI_Network_DF.reset_index(drop = True)

### GI Network Dataframes
GI_Epsi_network_DF = pd.read_csv("/home/xaviercg/scratch/data/costanzo_2016_longer_withoutReps.csv")


## Import other network attributes
GI_Modules = pd.read_excel("/home/xaviercg/scratch/data/Costanzo2016_DataFileS6.xlsx")

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



## Make a dictionnary storing GI cluster sizes
GI_cluster_sizes = {}
for key in list(GI_module_dict.keys()):
    size = len(GI_module_dict[key])
    GI_cluster_sizes[key] = size


## Make useful lists
all_GI_ORFs = list(chain(*list(GI_module_dict.values()))) # All ORFs in GI networks
all_GI_cluster_IDs = list(GI_module_dict.keys()) # List of all GI modules (cluster ID)
all_PPI_cluster_IDs = [x for x in list(PPI_module_dict.keys()) if x != "None"] # List of all PPI modules (CPX ID)


## Get list of all combinations of GI to PPI module pairings (will be used for calculating overlap)
module_pair_combinations = [list(i) for i in product(all_GI_cluster_IDs, all_PPI_cluster_IDs)]




#### Make Dataframe of Randomized GI Clusters Optimally Mapped to PPI Complexes ####

random.shuffle(all_GI_ORFs) # Randomize ORF Order
    
random_GI_cluster = dict.fromkeys(all_GI_cluster_IDs) # Create a dictionary that will hold random GI clusters
    
i = 0 # i is the index from list of ORFs
for clusterID in all_GI_cluster_IDs:
    cluster_size = GI_cluster_sizes[clusterID] # Size of the current randomized cluster
    random_GI_ORFs = all_GI_ORFs[i:(i+cluster_size)] # Select slice of ORF list equivalent in size to 'cluster_size'
    random_GI_cluster[clusterID] = random_GI_ORFs # Assign randomly selected ORFs to random cluster
    i += cluster_size # Increment i
    
    
### Find Overlap between individual GI modules and PPI modules
GI_PPI_module_Overlap_matrix_rand = pd.DataFrame(columns = all_PPI_cluster_IDs, index = all_GI_cluster_IDs) # Only counts the number of overlapping genes
GI_PPI_module_Jaccard_matrix_rand = pd.DataFrame(columns = all_PPI_cluster_IDs, index = all_GI_cluster_IDs) # Collects Jaccard indexes
    
    
# For every cell in matrix, calculate jaccard index between GI module (index) and PPI module (column)
for pair in module_pair_combinations:
    GI_module = pair[0]
    PPI_module = pair[1]
    GI_PPI_module_Jaccard_matrix_rand.loc[pair[0],pair[1]] = Jaccard(PPI_module_dict[pair[1]], random_GI_cluster[pair[0]])
    #GI_PPI_module_Overlap_matrix_rand.loc[pair[0],pair[1]] = len(set(PPI_module_dict[pair[1]]).intersection(random_GI_cluster[pair[0]]))
    #GI_PPI_module_Overlap_matrix_rand.loc[pair[1],pair[0]] = len(set(PPI_module_dict[pair[1]]).intersection(random_GI_cluster[pair[0]]))


#with open(f'/home/xaviercg/scratch/data/Controls/Random_GI_PPI_Module_Overlap/GI_PPI_ModuleOverlap_{ctrl_number}_DF.pickle', 'wb') as handle:
#        pickle.dump(GI_PPI_module_Overlap_matrix_rand, handle)

## Find best PPI complex for each GI cluster

# For every GI cluster, get the column index for which overlap is maximal
GI_PPI_optimal_Jaccard_DF_rand = pd.DataFrame(GI_PPI_module_Jaccard_matrix_rand.apply(lambda x: x.tolist().index(max(x)), axis = 1), columns = ["CPX_index"])

# Get complex ID from index
GI_PPI_optimal_Jaccard_DF_rand['CPX_ID'] = [all_PPI_cluster_IDs[index] for index in GI_PPI_optimal_Jaccard_DF_rand['CPX_index']]

# Make column for Jaccard index
GI_PPI_optimal_Jaccard_DF_rand['Jaccard_index'] = GI_PPI_optimal_Jaccard_DF_rand.apply(lambda x: GI_PPI_module_Jaccard_matrix_rand.loc[x.name, x['CPX_ID']], axis = 1)

# Make list of common ORF between GI and PPI modules for each pair
GI_PPI_optimal_Jaccard_DF_rand['Common_Module_ORFs'] = [set(PPI_module_dict[CPX_ID]).intersection(set(GI_module_dict[cluster])) for CPX_ID, cluster in zip(GI_PPI_optimal_Jaccard_DF_rand['CPX_ID'].tolist(), list(GI_PPI_optimal_Jaccard_DF_rand.index))]

#GI_PPI_optimal_Jaccard_DF_rand = GI_PPI_optimal_Jaccard_DF_rand[GI_PPI_optimal_Jaccard_DF_rand['Common_Module_ORFs'].apply(lambda x: len(x) > 0)]

#GI_PPI_optimal_Jaccard_DF_rand = GI_PPI_optimal_Jaccard_DF_rand[GI_PPI_optimal_Jaccard_DF_rand['CPX_ID'] != 'None']

with open(f'/home/xaviercg/scratch/data/Controls/Random_GI_PPI_Module_Match/GI_PPI_ModuleMatch_{ctrl_number}_DF.pickle', 'wb') as handle:
        pickle.dump(GI_PPI_optimal_Jaccard_DF_rand, handle)
