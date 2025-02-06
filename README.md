# PPI_Binding_Epistasis
All scripts and data related to the manuscript on protein-protein interaction binding and epistasis <br>
All R scripts were run on version 4.4 <br>
All python scripts were run on version 3.8 <br>
See data/ReadMe for information on data that is too heavy to store on GitHub <br>
*Italic text* Indicates a dependency i.e., an input file created using another script
<br>
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

## Analysis of Relation between Stoichiometries and GIs
**Methods section "Analysis of genetic interactions and protein-protein interaction stoichiometries"** <br>
Yeast_Stoichiometries_vs_GI.R <br>
&ensp; Input files: costanzo_2016_longer_withoutReps.csv ; *Yeast_Interactome_Stoich.csv* <br>
&ensp; Output files:  <br>
<br>
Human_Stoichiometries_vs_GI.R <br>
&ensp; Input files:  <br>
&ensp; Output files:  <br>


