# Jupyter Notebook script generously shared by Anastasia Baryshnikova


# Script makes a minimal Genetic Interaction set.

# Date First Created: April 5 2023
# Date Last Modified: April 16 2024

## Import libraries
import pandas as pd
import numpy as np
import re
from collections import defaultdict

## Set directories
data_folder = './data/GI_data/'
output_folder = './results/'


source_folder = data_folder
datasets = ['NxN','ExE','ExN_NxE']



#### Load all 3 datasets ####

source_folder = data_folder
datasets = ['NxN','ExE','ExN_NxE']

costanzo_2016 = {}
for dataset in datasets:
    print(dataset)
    data_file = source_folder + 'SGA_' + dataset + '.txt'
    costanzo_2016[dataset] = pd.read_table(data_file, delimiter='\t', header=0)
    
    
## First, make matrices for each of the datasets
costanzo_2016_matrix_scores = {}
costanzo_2016_matrix_pvals = {}

for dataset in datasets:
    print(dataset)
    # "Longer" Dataframe with "Array Strain" along x and "Query Strain" along y, GI scores as elements:
    costanzo_2016_matrix_scores[dataset] = costanzo_2016[dataset].pivot(index='Query Strain ID', 
                                                                          columns='Array Strain ID', 
                                                                          values='Genetic interaction score (ε)')
    
    # "Longer" Dataframe with "Array Strain" along x and "Query Strain" along y, p-values as elements:
    costanzo_2016_matrix_pvals[dataset] = costanzo_2016[dataset].pivot(index='Query Strain ID',
                                                                         columns='Array Strain ID',
                                                                         values ='P-value')
    
    # Print dimensions of dataframe:
    print('%d x %d' % (costanzo_2016_matrix_scores[dataset].shape[0],costanzo_2016_matrix_scores[dataset].shape[1]))


## Merge the 3 datasets into 1 (& keep the score with the lowest p-value)

queries = []
arrays = []

for dataset in datasets:
    queries = np.append(queries, costanzo_2016_matrix_scores[dataset].index.values) # Append row names (Query Strain)
    arrays = np.append(arrays, costanzo_2016_matrix_scores[dataset].columns.values) # Append col names (Array Strain)


queries = np.unique(queries)
arrays = np.unique(arrays)



costanzo_2016_scores = pd.DataFrame(data=np.full([len(queries), len(arrays)], np.nan), index=queries, columns=arrays)
costanzo_2016_pvals = pd.DataFrame(data=np.full([len(queries), len(arrays)], np.nan), index=queries, columns=arrays)


stacked_scores = []
stacked_pvals = []
for idx, dataset in enumerate(datasets):
    #print(dataset)
    [t0,_] = costanzo_2016_matrix_scores[dataset].align(costanzo_2016_scores, join='right')
    stacked_scores.append(t0)
    
    print(t0)
    
    [t0,_] = costanzo_2016_matrix_pvals[dataset].align(costanzo_2016_pvals, join='right')
    stacked_pvals.append(t0)
    


stacked_scores = np.dstack((stacked_scores[0],stacked_scores[1],stacked_scores[2]))
stacked_pvals = np.dstack((stacked_pvals[0],stacked_pvals[1],stacked_pvals[2]))


stacked_pvals = np.ma.masked_array(stacked_pvals, mask=np.isnan(stacked_pvals))


stacked_pvals_ix_min = np.argmin(stacked_pvals, axis=2)



m,n = stacked_pvals.shape[:2]
c,r = np.mgrid[:m,:n]
pvals_min = stacked_pvals[(c,r,stacked_pvals_ix_min)].filled(fill_value=np.nan)
scores_min = stacked_scores[(c,r,stacked_pvals_ix_min)]


costanzo_2016_scores = pd.DataFrame(data = scores_min, index=queries, columns=arrays)
costanzo_2016_pvals = pd.DataFrame(data = pvals_min, index=queries, columns=arrays)


## Now, make a unique query & array details dataframes
q = {}
a = {}
for dataset in datasets:
    
    print(dataset)
    
    q[dataset] = pd.DataFrame({'strain_id': costanzo_2016[dataset]['Query Strain ID'], 
                               'allele': costanzo_2016[dataset]['Query allele name']})
    a[dataset] = pd.DataFrame({'strain_id': costanzo_2016[dataset]['Array Strain ID'], 
                               'allele': costanzo_2016[dataset]['Array allele name']})

costanzo_2016_queries = pd.concat(q, axis=0, join='outer', ignore_index=True)
costanzo_2016_arrays = pd.concat(a, axis=0, join='outer', ignore_index=True)



costanzo_2016_queries.drop_duplicates(inplace=True)
costanzo_2016_arrays.drop_duplicates(inplace=True)



# Get the ORF info
pattern = re.compile('^Y[A-P][RL][0-9]{3}[CW](-[A-H])*')
costanzo_2016_queries['orf'] = [pattern.match(q).group(0) for q in costanzo_2016_queries['strain_id']]
costanzo_2016_arrays['orf'] = [pattern.match(a).group(0) for a in costanzo_2016_arrays['strain_id']]




costanzo_2016_queries.index = costanzo_2016_queries['strain_id']
costanzo_2016_arrays.index = costanzo_2016_arrays['strain_id']




# Drop the strains with suppressors
pattern_strain = re.compile('^Y[A-P][RL][0-9]{3}[CW](-[A-H])*_S')
pattern_allele = re.compile('.*supp')

queries_notsup = [q for q in costanzo_2016_queries['strain_id'] if not pattern_strain.match(q)]
arrays_notsup = [a for a in costanzo_2016_arrays['allele'] if not pattern_allele.match(a)]

costanzo_2016_queries = costanzo_2016_queries.loc[queries_notsup,:]

costanzo_2016_arrays = costanzo_2016_arrays.loc[costanzo_2016_arrays['allele'].isin(arrays_notsup),:]


## Check the uniqueness of the allele names
costanzo_2016_queries.loc[costanzo_2016_queries['allele'].duplicated(),:]


# Renaming the only duplicated allele in the query dataset
costanzo_2016_queries.loc[costanzo_2016_queries['strain_id'] == 'YHR058C_tsq1019','allele'] = 'med6_ts_ab1'


costanzo_2016_arrays.loc[costanzo_2016_arrays['allele'].duplicated(),:]



# Make sure that the same allele name corresponds to the same ORF across queries and arrays
tmp = pd.merge(costanzo_2016_queries, costanzo_2016_arrays, how='outer', on='allele', suffixes=['_q','_a'])


tmp['orf_match'] = ((tmp['orf_q'] == tmp['orf_a']) | (pd.isnull(tmp['orf_q'])) | (pd.isnull(tmp['orf_a'])))


np.sum(~tmp['orf_match'])


## Remove the data for all queries/arrays that were filtered out
costanzo_2016_scores = costanzo_2016_scores.reindex(index=costanzo_2016_queries.index, columns=costanzo_2016_arrays.index)
costanzo_2016_pvals = costanzo_2016_pvals.reindex(index=costanzo_2016_queries.index, columns=costanzo_2016_arrays.index)



## Constructing pairwise interaction format dataframes

# Reset column level name and index level name to be non-identical
costanzo_2016_pvals.columns.name = "strain_id_array"
costanzo_2016_pvals.index.name = "strain_id_query"

# Do the same for scores
costanzo_2016_scores.columns.name = "strain_id_array"
costanzo_2016_scores.index.name = "strain_id_query"


# Make longer, easier to scan for duplicated this way.
costanzo_2016_pvals_longer = costanzo_2016_pvals.melt(ignore_index=False).reset_index()
costanzo_2016_scores_longer = costanzo_2016_scores.melt(ignore_index=False).reset_index()


# Combine scores and pvalues into one dataframe
costanzo_2016_longer = costanzo_2016_pvals_longer
costanzo_2016_longer = costanzo_2016_longer.rename(columns={'value': 'pval'})

costanzo_2016_longer['scores'] = costanzo_2016_scores_longer['value']

costanzo_2016_longer = costanzo_2016_longer.dropna(axis = 0).reset_index(drop = True)


# Clear variables for memory issues
costanzo_2016_scores = None
costanzo_2016_pvals = None
costanzo_2016_pvals_longer = None
costanzo_2016_scores_longer = None
costanzo_2016_pvals = None
costanzo_2016_scores = None
costanzo_2016 = None


# Remove strain IDs
costanzo_2016_longer['ORF_query'] = costanzo_2016_longer['strain_id_query'].str.replace("_.*", "")
costanzo_2016_longer['strain_query'] = costanzo_2016_longer['strain_id_query'].str.replace(".*_", "")

costanzo_2016_longer['ORF_array'] = costanzo_2016_longer['strain_id_array'].str.replace("_.*", "")
costanzo_2016_longer['strain_array'] = costanzo_2016_longer['strain_id_array'].str.replace(".*_", "")

## Full strain ID no longer needed
costanzo_2016_longer = costanzo_2016_longer.drop('strain_id_query', axis = 1)
costanzo_2016_longer = costanzo_2016_longer.drop('strain_id_array', axis = 1)


## Remove duplicate ORF pairs by keeping the one with the lowest p-value

# Sort by p-values
costanzo_2016_longer = costanzo_2016_longer.sort_values(['pval'], ascending=True, axis = 0).reset_index(drop = True)

# Combine
costanzo_2016_longer_asSet = list(costanzo_2016_longer[['ORF_query', 'ORF_array']].apply(frozenset, axis=1))
costanzo_2016_longer['ORFs_asSet'] = costanzo_2016_longer_asSet

costanzo_2016_longer = costanzo_2016_longer.drop_duplicates(subset=['ORFs_asSet'], keep = 'first')

costanzo_2016_longer = costanzo_2016_longer.drop('ORFs_asSet', axis = 1)
outputfile = output_folder + "costanzo_2016_longer_withoutReps.csv"
costanzo_2016_longer.to_csv(outputfile, index=False)

