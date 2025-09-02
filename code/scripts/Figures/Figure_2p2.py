# Xavier Castellanos-Girouard

# Visualizing Reconstructed of GI matrix

# Date First Created: May 21 2024
# Date Last Modified: May 21 2024


## Import Libraries

import numpy as np
import pandas as pd
import networkx as nx
import seaborn as sns
import scipy.optimize as optimization
import scipy as sp
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from itertools import combinations
from itertools import repeat


#### Import and Format Data ####

## Import GI data
Costanzo_longer_DF = pd.read_csv("../results/Stoichiometry_and_GI/costanzo_2016_longer_withoutReps.csv", index_col = 0)

# Combined AB and BA to construct Matrix
Costanzo_longer_DF1 = Costanzo_longer_DF[['ORF_query', 'ORF_array', 'scores']].copy()
Costanzo_longer_DF2 = Costanzo_longer_DF[['ORF_query', 'ORF_array', 'scores']].copy()
Costanzo_longer_DF2 = Costanzo_longer_DF2.rename(columns={"ORF_query": "ORF_array", "ORF_array": "ORF_query"})

Costanzo_longer_DF = pd.concat([Costanzo_longer_DF1, Costanzo_longer_DF2], axis = 0, ignore_index = True)

Costanzo_longer_DF = Costanzo_longer_DF.reset_index(drop = True)

# Pivot Wider
Costanzo_longer_DF = Costanzo_longer_DF.pivot(index = "ORF_query", columns = "ORF_array", values = "scores")
#Costanzo_longer_DF

# Clear memory
Costanzo_longer_DF1 = None
Costanzo_longer_DF2 = None

Kd_GI_inferred_clustered_DF = pd.read_csv("../results/Kd_and_GI/Clustered_Kd_inferred_GI.csv", index_col = 0)

Kd_GI_inferred_clustered_DF = Kd_GI_inferred_clustered_DF.fillna(0, inplace = False)

Kd_GI_inferred_clustered_DF.columns = [x.replace(".", "-") for x in Kd_GI_inferred_clustered_DF.columns]

#print(Kd_GI_inferred_clustered_DF.head(5))
#print(Costanzo_longer_DF.head(5))


#### Common ORFs to filter ####

PPI_ORF_ls = list(np.unique(Kd_GI_inferred_clustered_DF.index))
GI_ORF_ls = list(np.unique(Costanzo_longer_DF.index))
common_ORFs = [x for x in PPI_ORF_ls if x in GI_ORF_ls]

## Filter GI dataset (only keep those that appear in Kd network)
indices_subset = Costanzo_longer_DF.index[Costanzo_longer_DF.index.isin(common_ORFs)]
columns_subset = Costanzo_longer_DF.columns[Costanzo_longer_DF.columns.isin(common_ORFs)]
Costanzo_longer_DF = Costanzo_longer_DF.loc[indices_subset, columns_subset]


## Filter PPI dataset
indices_subset = Kd_GI_inferred_clustered_DF.index[Kd_GI_inferred_clustered_DF.index.isin(common_ORFs)]
columns_subset = Kd_GI_inferred_clustered_DF.columns[Kd_GI_inferred_clustered_DF.columns.isin(common_ORFs)]
Kd_GI_inferred_clustered_DF = Kd_GI_inferred_clustered_DF.loc[indices_subset, columns_subset]



#### Viz of Kd PPI Network ####


## Full matrix
fig, ax = plt.subplots(figsize=(360, 360))
g = sns.heatmap(Kd_GI_inferred_clustered_DF, cmap='viridis', norm=LogNorm(), cbar=False)
g.set_facecolor('xkcd:black')
plt.savefig("../results/Figures/Main_Figures/GI_reconstructed_full.png")
plt.close()


## Selected Matrix
#index_start = np.where(Kd_GI_inferred_clustered_DF.index.values == 'YNL299W')[0][0]
index_start = np.where(Kd_GI_inferred_clustered_DF.index.values == 'YOL151W')[0][0]
#index_end = np.where(Kd_GI_inferred_clustered_DF.index.values == 'YOR039W')[0][0]
index_end = np.where(Kd_GI_inferred_clustered_DF.index.values == 'YLR424W')[0][0]
fig, ax = plt.subplots(figsize=(15, 15))
g = sns.heatmap(Kd_GI_inferred_clustered_DF.iloc[index_start:index_end, index_start:index_end], cmap='viridis', norm=LogNorm(), cbar=False, yticklabels=False,xticklabels=False)
g.set_facecolor('xkcd:black')
plt.savefig("../results/Figures/Main_Figures/GI_reconstructed_partial.png")
#plt.show()



#### Compare to GI network ####

## Initiate GI array (will be filled)
GI_array = np.empty(Kd_GI_inferred_clustered_DF.shape)
GI_array.fill(np.nan)

## Get a list of ORFs. This defined the order of genes in the rows and columns of 2d array
ORF_arr = np.array(Kd_GI_inferred_clustered_DF.columns)

for i in Costanzo_longer_DF.index:
    for j in Costanzo_longer_DF.columns:
        if np.isnan(Costanzo_longer_DF.loc[i,j]):
            continue
        
        if (i in ORF_arr) and (j in ORF_arr):
            i_pos, = np.where(ORF_arr == i)
            j_pos, = np.where(ORF_arr == j)
            GI_array[i_pos[0], j_pos[0]] = Costanzo_longer_DF.loc[i,j]
            GI_array[j_pos[0], i_pos[0]] = Costanzo_longer_DF.loc[i,j]
            
GI_DF = pd.DataFrame(GI_array)
GI_DF.index = ORF_arr
GI_DF.columns = ORF_arr


fig, ax = plt.subplots(figsize=(15, 15))
g = sns.heatmap(GI_DF.iloc[index_start:index_end, index_start:index_end], center=0, vmin = -0.16, vmax = 0.16, cmap=sns.diverging_palette(220, 80, l=70, s = 200,center="dark", as_cmap=True), cbar=False, yticklabels=False,xticklabels=False)
g.set_facecolor('xkcd:black')
plt.savefig("../results/Figures/Main_Figures/GI_experimental_partial.png")
#plt.show()


#### Export Legends/Colorbars ####

#colobar_l = sns.diverging_palette(220, 80, l=70, s = 200,center="dark", as_cmap=True)
