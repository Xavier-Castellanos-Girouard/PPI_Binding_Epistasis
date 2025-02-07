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
## Yeast Genetic Interaction data processing script
**Methods section "Collecting yeast and human genetic and protein-protein interaction data"** <br>
Make_Constanzo_2016_dataset_minimal.py <br>
&ensp; Input files: SGA_ExE.txt ; SGA_ExN_NxE.txt ; SGA_NxN.txt <br>
&ensp; Output files: costanzo_2016_longer_withoutReps.csv <br>
<br>
## Yeast Stoichiometries data processing script
**Methods section "Collecting yeast and human genetic and protein-protein interaction data"** <br>
Calc_Interaction_Stoichiometry.R <br>
&ensp; Input files: Abundances_copies_per_cell.xlsx ; preprocessedProteinDataForCorrelation.txt <br>
&ensp; Output files: Yeast_Interactome_Stoich.csv <br>
<br>
## Analysis of Relation between Stoichiometries and GIs
**Methods section "Analysis of genetic interactions and protein-protein interaction stoichiometries"** <br>
Yeast_Stoichiometries_vs_GI.R <br>
&ensp; Input files: *costanzo_2016_longer_withoutReps.csv* ; *Yeast_Interactome_Stoich.csv* <br>
&ensp; Output files: Yeast_Stoich_GI.csv ; Yeast_Stoich_GI_wRegions.csv <br>
<br>
Human_Stoichiometries_vs_GI.R <br>
&ensp; Input files: GI_Horlbeck.xlsx ; opencell-protein-abundance.csv ; opencell-protein-interactions.csv ; Stoichiometries_Hein.xlsx <br>
&ensp; Output files: Human_Interactome_Stoich.csv ; Human_Stoich_GI.csv ; Human_Stoich_GI_wRegions.csv <br>
## Collecting Kd values from literature
**Methods section "Validating estimated dissociation constants using values from literature"**
<br>
PDB_bind_Extract.R <br>
&ensp; Input files: INDEX_general_PP.2020.txt ; PDBbind_vs2020/PDB_files/* <br>
&ensp; Output files: Kd_PDB_list.csv <br>
PDB2Uniprot.py <br>
&ensp; Input files: *Kd_PDB_list.csv* <br>
&ensp; Output files: PDB2Uniprot.csv <br>
Kd_PDB2Sys.R <br>
&ensp; Input files: *PDB2Uniprot.csv* ; *Kd_PDB_list.csv* ; uniprotkb_database_HGNC_2024_02_28.tsv <br>
&ensp; Output files: Yeast_Kd_literature.csv ; Human_Kd_literature.csv <br>
