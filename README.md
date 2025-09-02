# PPI_Binding_Epistasis
This Repository contains all scripts and data related to the manuscript on protein-protein interaction binding and epistasis <br>
<br>
All R scripts were run on version 4.4 ; BiocManager is used for package management <br>
All python scripts were run on version 3.9 <br>
Scripts for analyses performed on the High Performance Computing (HPC) clusters are specifically written for the Slurm Workload Manager. <br>
<br>
See data/ReadMe for information on data that is too heavy to store on GitHub <br>
Below, *Italic text* Indicates a dependency i.e., an input file created using another script <br>
Executing the analysis in the order described below should not yield any dependency errors <br>
<br>
A Master Script (Master_Script.sh) and dockerfile is provided to run the entire pipeline from start to finish <br>
<br>
## 1. Stoichiometry and GI
### Yeast Genetic Interaction data processing script
**Methods section "Collecting yeast and human genetic and protein-protein interaction data"** <br>
This script takes 3 files containing genetic interaction scores and p-values from Costanzo et al. (Science, 2016) and merges them into a single dataframe. The dataset contains no duplicated interactions: for every unique pair of genes, the interaction with the lowest p-value is retained.  <br>
 <br>
&ensp; Script Name: Make_Constanzo_2016_dataset_minimal.py <br>
&ensp; Input files: SGA_ExE.txt ; SGA_ExN_NxE.txt ; SGA_NxN.txt <br>
&ensp; Output files: costanzo_2016_longer_withoutReps.csv <br>
<br>
### Yeast Stoichiometries data processing script
**Methods section "Collecting yeast and human genetic and protein-protein interaction data"** <br>
This script takes the preprocessed mass-spectrometry datatable from Michaelis et al. (Nature, 2023) to calculate interaction stoichiometries (and abundance stoichiometries). <br>
<br>
&ensp; Script Name: Calc_Interaction_Stoichiometry.R <br>
&ensp; Input files: Abundances_copies_per_cell.xlsx ; preprocessedProteinDataForCorrelation.txt <br>
&ensp; Output files: Yeast_Interactome_Stoich.csv <br>
<br>
### Analysis of Relation between Stoichiometries and GIs
**Methods section "Analysis of genetic interactions and protein-protein interaction stoichiometries"** <br>
These scripts integrate interaction stoichiometries and genetic interaction datasets. This is done separately for yeast and human.<br>
<br>
&ensp; Script Name: Yeast_Stoichiometries_vs_GI.R <br>
&ensp; Input files: *costanzo_2016_longer_withoutReps.csv* ; *Yeast_Interactome_Stoich.csv* <br>
&ensp; Output files: Yeast_Stoich_GI.csv ; Yeast_Stoich_GI_wRegions.csv <br>
<br>
&ensp; Script Name: Human_Stoichiometries_vs_GI.R <br>
&ensp; Input files: GI_Horlbeck.xlsx ; opencell-protein-abundance.csv ; opencell-protein-interactions.csv ; Stoichiometries_Hein.xlsx <br>
&ensp; Output files: Human_Interactome_Stoich.csv ; Human_Stoich_GI.csv ; Human_Stoich_GI_wRegions.csv <br>
<br>
## 2. Kd and GI
### Collecting Kd values from literature
**Methods section "Validating estimated dissociation constants using values from literature"**
These three scripts do the following: extract PPI information from the PDBbind database, files, convert PDB IDs to Uniprot IDs, converts the Uniprot IDs to systematic gene names and creates clean tables. <br>
<br>
&ensp; Script Name: PDB_bind_Extract.R <br>
&ensp; Input files: INDEX_general_PP.2020.txt ; PDBbind_vs2020/PDB_files/* <br>
&ensp; Output files: Kd_PDB_list.csv <br>
<br>
&ensp; Script Name: PDB2Uniprot.py <br>
&ensp; Input files: *Kd_PDB_list.csv* <br>
&ensp; Output files: PDB2Uniprot.csv <br>
<br>
&ensp; Script Name: Kd_PDB2Sys.R <br>
&ensp; Input files: *PDB2Uniprot.csv* ; *Kd_PDB_list.csv* ; uniprotkb_database_HGNC_2024_02_28.tsv <br>
&ensp; Output files: Yeast_Kd_literature.csv ; Human_Kd_literature.csv <br>
### Calculating Kd values using PCA values
**Methods section "Comparison of mass spectrometry Kd estimates in yeast using PCA Kd estimates"** <br>
This script estimates Kd values from a systematic Protein-fragment complementation screen (Tarassov et al., Science, 2008). <br>
<br>
Kd_PCA_Estimation.R <br>
&ensp; Input files: Abundances_copies_per_cell.xlsx ; Abundace_Levy_2014.xls ; Tarassov_2008_PPI.xls ; network_1.tsv <br>
&ensp; Output files: Estimated_Kd2_PCA_intensitiesNorm_LinFit_MeanCopyNumberUsed.csv <br>
### Estimating and benchmarking Kd values from Mass Spectrometry 
**Methods section "Estimating protein-protein interaction dissociation constants" and "Validating estimated dissociation constants using values from literature"**
These scripts estimate Kd values for PPIs in Human and Yeast, and compares them to values from the litterature. <br>
<br>
&ensp; Script Name: Calc_Kd_HumanYeast.R <br>
&ensp; Input files: *Human_Interactome_Stoich.csv* ; Hein_TableS3.xlsx ; opencell-protein-abundance.csv ; *Human_Kd_literature.csv* ; *Yeast_Interactome_Stoich.csv* ; *Yeast_Kd_literature.csv* ; *Estimated_Kd2_PCA_intensitiesNorm_LinFit_MeanCopyNumberUsed.csv* <br>
&ensp; Output files: Human_Kd_Est_vs_Lit.csv ; Human_Estimated_Kd.csv ; Yeast_Estimated_Kd.csv ; Yeast_Kd_Est_vs_PCA_Kd.csv <br>
### Analysis of Relation between Kd and GIs
**Methods section "Reconstructing a GI network matrix from a PPI binding network"** <br>
This script integrates protein-protein interaction Kd values and genetic interactions. The yeast and human analyses are independent but contained in the same script. <br>
<br>
&ensp; Calc_Kd_HumanYeast.R <br>
&ensp; Input files: *costanzo_2016_longer_withoutReps.csv* ; *Yeast_Estimated_Kd.csv* ; GI_Horlbeck.xlsx ; *Human_Estimated_Kd* <br>
&ensp; Output files: Yeast_Kd_GI.csv ; Yeast_Mean_GI_KdDeciles.csv ; Human_Kd_GI.csv ; Human_Mean_GI_KdQuintiles.csv <br>
### Predicting GI matrix from PPI binding data
**Methods section "Reconstructing a GI network matrix from a PPI binding network"** <br>
This script uses an approach similar to that of Rives & Galitsky (PNAS, 2003) on the yeast PPI network weighted by dissocation constants (Kd) to predict a genetic interaction network. Cluster3 is used for Hierarchical clustering, and the output is processed by CDT2DF.R <br>
<br>
&ensp; Script Name: GI_Matrix_Reconstruction.py <br>
&ensp; Input files: *Yeast_Kd_GI* <br>
&ensp; Output files: Kd_inferred_GI.tsv <br>
<br>
&ensp; Script Name: HierarchicalClustering_Cluster3.sh <br>
&ensp; Input files: *Kd_inferred_GI.tsv* <br>
&ensp; Output files: Kd_inferred_GI.cdt ; Kd_inferred_GI.gtr ; Kd_inferred_GI.atr <br>
<br>
&ensp; Script Name: CDT2DF.R <br>
&ensp; Input files: *Kd_inferred_GI.cdt* ; *Kd_inferred_GI.gtr* ; *Kd_inferred_GI.atr* <br>
&ensp; Output files: Clustered_Kd_inferred_GI.csv <br>
<br>
## 3. Paralog Analyses
### Analysis of dissociation constants in yeast paralogs
**Methods section "Analysis of dissociation constants in yeast paralogs"** <br>
These scripts compute information on paralog genes for each PPI and integrates these with interaction stoichiometry, Kd, and GI scores are incorported. ClustalW is used to align paralog sequences; these are analysed in Paralog_Divergence_Analysis.R. <br>
<br>
&ensp; Script Name: Paralog_Annot.R <br>
&ensp; Input files: Kellis_S8_DupGenes/*.txt ; msb201082-sup-0002.xls ; *Yeast_Kd_GI.csv* ; *Yeast_Stoich_GI.csv* <br>
&ensp; Output files: Paralog_DuplicationType.csv ; Stoich_GI_ParalogType.csv ; Paralog_proteinSeq_ClustalW_input/*.fasta
<br>
&ensp; Script Name: Execute_ProtAlign_ClustalW.sh <br>
&ensp; Input files:  Paralog_proteinSeq_ClustalW_input/*.fasta <br>
&ensp; Output files: ClustalW_Output/*.aln<br>
<br>
&ensp; Script Name: Paralog_Divergence_Analysis.R <br>
&ensp; Input files:  *Paralog_DuplicationType.csv* ; *Kd_GI_ParalogType.csv* ; orf_coding_all.fasta <br>
&ensp; Output files: all_paralog_kaks.csv ; Kd_GI_Divergence.csv<br>
<br>
## 4. PPI and GI topology Analyses
### Integrate Module information (ComplexPortal) to PPI network
**Methods section "Mapping GI Modules to PPI Modules"** <br>
This script annotates PPIs with information on protein complexes (ComplexPortal). <br>
<br>
&ensp; Script Name: R_Create_Modular_PPI_Network.R <br>
&ensp; Input files: ComplexPortal_559292.tsv ; The_Yeast_Interactome_Edges.csv <br>
&ensp; Output files: Yeast_PPI_Network_CPX.csv <br>
<br>
### Map PPI to GI modules
**Methods section "Mapping GI Modules to PPI Modules"** <br>
This script maps GI modules to PPI modules and determines the statistical significance  of this overlap. <br>
<br>
&ensp; Script Name: Module_Overlap_Analysis.py <br>
&ensp; Input files: *Yeast_PPI_Network_CPX.csv* ; *costanzo_2016_longer_withoutReps.csv* ; Costanzo2016_DataFileS6.xlsx <br>
&ensp; Output files: GI_PPI_randomized_ModuleOverlap.csv <br>
<br>
### Find Shortest Path connectors in GI and PPI networks
**Methods section "Identifying Nodes Connecting GI and PPI Modules"** <br>
These scripts find shortest path connectors in GI and PPI networks. <br>
<br>
&ensp; Script Name: Find_GI_ShortestPath_Connectors.py <br>
&ensp; Input files: *costanzo_2016_longer_withoutReps.csv* ; *Costanzo2016_DataFileS6.xlsx* <br>
&ensp; Output files: GI_pairwise_module_shortest_paths_dict.pickle <br>
<br>
&ensp; Script Name: Find_PPI_ShortestPath_Connectors.py <br>
&ensp; Input files: *GI_PPI_optimal_module_overlap_clean.csv* ; *Yeast_PPI_Network_CPX.csv* <br>
&ensp; Output files: PPI_pairwise_module_shortest_paths_dict.pickle <br>
<br>
### Find common Shortest Path connectors that in GI and PPI networks
**Methods section "Identifying Nodes Connecting GI and PPI Modules"** <br>
This script finds common Shortest Path connectors that in GI and PPI networks. <br>
<br>
&ensp; Script Name: Find_GI_ShortestPath_Connectors.py <br>
&ensp; Input files: *GI_PPI_optimal_module_overlap_clean.csv* ; *GI_pairwise_module_shortest_paths_dict.pickle* ; *PPI_pairwise_module_shortest_paths_dict.pickle* <br>
&ensp; Output files: Shortest_Path_Connector_Overlap.csv <br>
<br>
