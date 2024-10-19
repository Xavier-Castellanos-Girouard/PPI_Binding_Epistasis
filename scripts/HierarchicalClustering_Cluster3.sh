#!/bin/bash

# Xavier Castellanos-Girouard

# Use Cluster3 to perform Hierarchical Clustering

# Date First Created: March 15 2024
# Date Last Modified: May 22 2024

# Execute in ./scripts directory

## Use cluster3 with pearson correlation and average linkage hierarchical clustering

cluster3 -f ../results/Kd_inferred_GI.tsv -g 2 -e 2 -m a
