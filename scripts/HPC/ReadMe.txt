#### Notes on HPC scripts ####

# All HPC scripts are run on Python 3.9.6

## 1. Creating Randomized Networks

# There are two Run scripts which should be executed with the 'sbatch' command: 
Create_GI_Connector_DegPres_RandNet_Run.sh
Create_PPI_Connector_DegPres_RandNet_Run.sh

# These will run several instances of the following python scripts:
HPC_Create_GI_Connector_DegPres_RandNet.py
HPC_Create_PPI_Connector_DegPres_RandNet.py

# The python scripts will generate degree-preserved randomized networks.

## 2. Performing the Module Overlap Analysis

# The following Launch script should be executed using the 'bash' command:
GI_PPI_TotalOverlapSig_Launch.sh

# Which will launch several instances of the following Run script
GI_PPI_TotalOverlapSig_Run.sh

# Which will run several instances of the following python script:
HPC_Module_Overlap_Analysis.py

# The resulting csv files can then be concatenated into a single csv file using the following script:
HPC_Module_Overlap_Analysis_Concat.py

## 3. Finding the Shortest Path Connectors between Modules

# The following python script needs to be run to prepare the files for Connector analysis:
HPC_Prep_Shortest_Path_Overlap_Files.py

# Next, the following Launch scripts should be executed using the 'bash' command:
Find_GI_Shortest_Path_Connectors_Launch.sh
Find_PPI_Shortest_Path_Connectors_Launch.sh

# Which will launch several instances of the following Run scripts:
Find_GI_ShortestPath_Connectors_Run.sh
Find_PPI_ShortestPath_Connectors_Run.sh

# Which will run several instances of the following python scripts:
HPC_Find_GI_ShortestPath_Connectors.py
HPC_Find_PPI_ShortestPath_Connectors.py

## 4. Performing Shortest Path Overlap Analysis

# The following Launch scripts should be executed using the 'bash' command:
Identify_Connector_Overlap_Launch.sh

# Which will Launch several instances of the following Run script:
Identify_Connector_Overlap_Run.sh

# Which will run several instances of the following python script:
HPC_Identify_Connector_Overlap.py

# The resulting csv files can then be concatenated into a single csv file using the following script:
HPC_Connector_Overlap_Analysis_Concat.py
