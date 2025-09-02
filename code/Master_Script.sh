#### Master Script ####

## Executing this bash script will run all the analysis scripts in the correct order.
## Output will generate all necessary files to recreate figures and results in the manuscript.
## HPC scripts are not included in this pipeline.

## Note that directories are hardcoded in each bash, python and R scripts.
## Note that scripts assume all necessary packages and libraries are installed.


## Command unique to codeocean. Installs Cluster3 with program files in local directory
#echo "##### Installing Cluster3 #####"
#bash ./CmdLine_Program_Install/Install_Cluster3.sh


## Create directory to store intermediate files. Solely to be used as input for other scripts.
## These output files are not used for figures or intended to be interpreted.
mkdir ../results/IntermediateFiles

## Make directories for results
mkdir ../results/Stoichiometry_and_GI
mkdir ../results/Kd_and_GI
mkdir ../results/Paralog_Analyses
mkdir ../results/Topological_Analyses

#### PPI Stoichiometries vs GI Analyses ####

## 1. Processing of GI data
echo "##### Processing GI data #####"
python3.9 scripts/Stoichiometry_and_GI/Make_Costanzo_2016_dataset_minimal.py

## 2. Calculating yeast Interaction Stoichiometries from MS data
echo "##### Calculating Yeast Interaction Stoichiometries #####"
Rscript scripts/Stoichiometry_and_GI/Calc_Interaction_Stoichiometry.R

## 3. Merging PPI Stoichiometries and GIs

# Yeast
echo "##### Merging Yeast Stoichiometries and GIs #####"
Rscript scripts/Stoichiometry_and_GI/Yeast_Stoichiometries_vs_GI.R

# Human
echo "##### Merging Human Stoichiometries and GIs #####"
Rscript scripts/Stoichiometry_and_GI/Human_Stoichiometries_vs_GI.R

#### PPI Kd vs GI Analyses ####

## 4. Collecting Kd values from literature
echo "##### Extract Kd values and PDB IDs #####"
Rscript scripts/Kd_and_GI/PDB_bind_Extract.R

echo "##### Convert PDB IDs to Uniprot IDs #####"
# Note: This cannot be run on codeocean because of internet restriction. Pregenerated file is used instead.
python3.9 scripts/Kd_and_GI/PDB2Uniprot.py

echo "##### Convert PDB IDs to systematic yeast IDs #####"
Rscript scripts/Kd_and_GI/Kd_PDB2Sys.R

## 5. Calculating Kd values using PCA values
echo "##### Estimate Kd values from PCA data #####"
Rscript scripts/Kd_and_GI/Kd_PCA_Estimation.R

## 6. Estimating and benchmarking Kd values from Mass Spectrometry
echo "##### Estimate Kd values from MS data and Benchmark #####"
Rscript scripts/Kd_and_GI/Calc_Kd_HumanYeast.R

## 7. Analysis of Relation between Kd and GIs

echo "##### Analysis of Relation between Kd and GIs ####"
Rscript scripts/Kd_and_GI/Kd_vs_Epistasis.R

## 8. Reconstruction of GI matrix using Kd-weighted PPI network data

# Reconstruct GI matrix from kd values
echo "##### Reconstruct GI matrix from Kd values #####"
python3.9 scripts/Kd_and_GI/GI_Matrix_Reconstruction.py

# Perform Hierarchical clustering on reconstructed Matrix
echo "##### Hierarchical clustering of reconstructed GI matrix #####"
bash scripts/Kd_and_GI/HierarchicalClustering_Cluster3.sh

# Convert cluster files to csv
echo "##### Convert cluster files to csv format #####"
Rscript scripts/Kd_and_GI/CDT2DF.R


## 9. Paralog Analysis
mkdir ../results/IntermediateFiles/Paralog_proteinSeq_ClustalW_input

mkdir ../results/IntermediateFiles/reverseAlign_nucl_FastaFiles

mkdir ../results/IntermediateFiles/reverseAlign_output_Files

mkdir ../results/IntermediateFiles/ClustalW_Output

# Paralog annotation with duplication type, Kd, stoichiometries and sequence
echo "##### Paralog annotation #####"
Rscript "scripts/Paralog_Analyses/Paralog_Annot.R"

# Protein alignment with Clustal
echo "##### Protein Alignment with Clustal v2.1 #####"

cd ../results/IntermediateFiles/
bash ../../code/scripts/Paralog_Analyses/Execute_ProtAlign_ClustalW.sh
cd -

# Perform analysis related to paralog divergence
echo "##### Paralog divergence analysis #####"
Rscript "scripts/Paralog_Analyses/Paralog_Divergence_Analysis.R"


## 10. Module Overlap Analysis
echo "##### Creating Modular PPI Network #####"
Rscript scripts/Topological_Analyses/R_Create_Modular_PPI_Network.R

echo "##### Find overlap in GI and PPI modules #####"
python3.9 scripts/Topological_Analyses/Module_Overlap_Analysis.py

## 11. Connector Overlap Analysis
echo "##### Creating Lookup Tables for Yeast Systematic IDs #####"
python3.9 scripts/Topological_Analyses/Make_ORF_Lookup_Table.py

echo "##### Finding Connectors in PPI network #####"
python3.9 scripts/Topological_Analyses/Find_PPI_ShortestPath_Connectors.py

echo "##### Finding Connectors in GI network #####"
python3.9 scripts/Topological_Analyses/Find_GI_ShortestPath_Connectors.py

echo "##### Finding Overlap between GI and PPI connectors #####"
python3.9 scripts/Topological_Analyses/GI_PPI_Connector_Overlap.py

## 12. Create Figures

mkdir ../results/Figures
mkdir ../results/Figures/Main_Figures
mkdir ../results/Figures/Extended_Data_Figures

echo "##### Creating Main Figures #####"

Rscript "scripts/Figures/Figure_1.R"
Rscript "scripts/Figures/Figure_2p1.R"
python3.9 "scripts/Figures/Figure_2p2.py"
Rscript "scripts/Figures/Figure_3.R"
Rscript "scripts/Figures/Figure_4.R"

echo "##### Creating Extended Data Figures #####"

Rscript "scripts/Figures/Ext_Data_Figure_2.R"
Rscript "scripts/Figures/Ext_Data_Figure_3.R"
Rscript "scripts/Figures/Ext_Data_Figure_4.R"
Rscript "scripts/Figures/Ext_Data_Figure_5.R"
Rscript "scripts/Figures/Ext_Data_Figure_6.R"
Rscript "scripts/Figures/Ext_Data_Figure_7.R"
Rscript "scripts/Figures/Ext_Data_Figure_8.R"
Rscript "scripts/Figures/Ext_Data_Figure_9.R"
Rscript "scripts/Figures/Ext_Data_Figure_10.R"
Rscript "scripts/Figures/Ext_Data_Figure_11.R"
Rscript "scripts/Figures/Ext_Data_Figure_12.R"
