## Xavier Castellanos-Girouard

## Calculating Interaction Stoichiometry using Dataset Given by Andre Michaelis

## Date First Created: Jan 12 2024
## Date Last Modified: April 18 2024

#### Import Libraries ####

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(seqinr)
library(data.table)
library(readxl)

#### Import Data ####

## Import cellular copy numbers from Ho et al. (2018, Cell systems)
Yeast_protein_abundances <- 
  readxl::read_xlsx("/home/xavier/Desktop/Cell_interactome_stoichiometry/Yeast_Epistasis_Stoichiometries/Third_Run/R_Calc_Stoichiometries/data/Abundances_copies_per_cell.xlsx")

colnames(Yeast_protein_abundances) <- Yeast_protein_abundances[2,]
Yeast_protein_abundances <- Yeast_protein_abundances[3:nrow(Yeast_protein_abundances),]

Yeast_protein_abundances$`Median molecules per cell` <- as.numeric(Yeast_protein_abundances$`Median molecules per cell`)


## Import processed data given 
processed_Data <- read.table("/home/xavier/Desktop/Cell_interactome_stoichiometry/Yeast_Epistasis_Stoichiometries/Third_Run/R_Calc_Stoichiometries/data/preprocessedProteinDataForCorrelation.txt", sep = "\t", header = TRUE)


#### Substract Background data ####

processed_Data_normalized <- processed_Data
processed_Data_normalized[,1:8294] <- 2^processed_Data[,1:8294]

## background substraction:
processed_Data_normalized[,1:8294] <- processed_Data_normalized[,1:8294]-apply(processed_Data_normalized[,1:8294], MARGIN = 1, FUN = median)


#### Divide by number of theoretical peptides ####

# Get list of all detected proteins in mass spec
detected_proteins <- unique(unlist(strsplit(processed_Data_normalized$Majority.protein.IDs, split = ";")))

# Import Yeast full uniprot database containing protein sequences. 
Uniprot_reference_sequences <- 
  read.fasta(file = "/home/xavier/Desktop/Cell_interactome_stoichiometry/Yeast_Epistasis_Stoichiometries/Third_Run/R_Calc_Stoichiometries/data/uniprotkb_proteome_UP000002311_2023_11_30.fasta",
             seqtype = c("AA"))

# Get uniprot ID from fasta sequence titles
Uniprot_reference_sequences_uniprot_ID <- 
  str_match(names(Uniprot_reference_sequences), 
            pattern = "\\|.*\\|")

# Remove vertical bar from strings
Uniprot_reference_sequences_uniprot_ID <- str_remove_all(Uniprot_reference_sequences_uniprot_ID, pattern = "\\|")

names(Uniprot_reference_sequences) <- Uniprot_reference_sequences_uniprot_ID

# Filter uniprot database, only keep entries that were detected in experiment.
Uniprot_reference_sequences <- Uniprot_reference_sequences[Uniprot_reference_sequences_uniprot_ID %in% detected_proteins]

# From protein sequence, make continuous string of characters instead of list of characters
Uniprot_reference_sequences <- 
  lapply(Uniprot_reference_sequences,
         FUN = paste0,
         collapse = "")

# Perform LysC digestion on sequences
Uniprot_reference_sequences_dig <-
  lapply(Uniprot_reference_sequences,
         str_split,
         pattern = "(?<=[K])")

# Keep digests of length 7-30
for (i in 1:length(Uniprot_reference_sequences_dig)){
  curr_digests <- Uniprot_reference_sequences_dig[[i]][[1]]
  curr_digests <- curr_digests[(str_length(curr_digests) >= 7) & (str_length(curr_digests) <= 30)]
  Uniprot_reference_sequences_dig[[i]][[1]] <- curr_digests
}

theoretical_peptide_count <- lapply(Uniprot_reference_sequences_dig,
                                    FUN = function(x){return(length(x[[1]]))})


# Assign number of theoretical peptides to table
processed_Data_normalized$num_theoretical_peptide <- NA

# For loop to do the assignment
for (i in 1:nrow(processed_Data_normalized)){ # for every Protein ID detected in Mass spec
  uniprot <- processed_Data_normalized$Majority.protein.IDs[i] # Retrieve protein ID
  uniprot <- unlist(str_split(uniprot, pattern = ";")) # If multiple protein IDs, make vector of IDs
  #print(uniprot)
  
  
  if (length(uniprot)==1){ # IF only one ID, get number of theoretical peptides
    #print("1")
    if (!is.null(theoretical_peptide_count[[uniprot]])){
      processed_Data_normalized$num_theoretical_peptide[i] = theoretical_peptide_count[[uniprot]]
    }
  } else if (length(uniprot)>1){ # If more than one ID, get max number of theoretical peptides ** maybe use sum here?
    #print("2")
    max_ = 0
    for (j in uniprot){
      #print(i)
      
      if (is.null(theoretical_peptide_count[[j]])){next} 
      
      else if (theoretical_peptide_count[[j]]>max_){
        max_ = theoretical_peptide_count[[j]]
      }
    }
    processed_Data_normalized$num_theoretical_peptide[i] = max_
  }
}

#### Calculate Stoichiometries ####


## Get a vector of Bait Intensities that match with column names (length should be 8294)
# NA is returned if bait is not detected.
Bait_intensities <- c()
col_ORFs <- c()
for (col_name in colnames(processed_Data_normalized[,1:8294])){
  #print(col_name)
  split_col <- str_split(col_name, pattern = "_") # Split col name
  col_ORF <- split_col[[1]][1] # Bait ORF ID
  #col_rep <- split_col[[1]][2] # Replicate number
  
  # Replace dashed with dots because of character restriction in column names
  col_ORF <- str_replace(col_ORF, pattern = "-", replacement = "\\.")
  
  col_ORFs <- c(col_ORFs, col_ORF)
  
  
  # Search for Bait ORF in detected proteins: 
  index <- which(processed_Data_normalized$ORF==col_ORF)
  
  # Get Bait Intensity
  bait_intensity <- processed_Data_normalized[index, col_name]
  
  # If no intensity is found (bait not detected in pulldown) store NA in vector
  if (length(bait_intensity) == 0){
    Bait_intensities <- c(Bait_intensities, NA)
  } else if (length(bait_intensity) == 1){ # If intensity is found, store value.
    Bait_intensities <- c(Bait_intensities, bait_intensity)
    }
}


# Calculate Stoichiometries
Stoichiometries_DF <- processed_Data_normalized
Stoichiometries_DF[,1:8294] <- mapply(processed_Data_normalized[,1:8294], FUN = `/`, Bait_intensities)

#### Map stoichiometries to significant interactions ####

## Import statistically significant interactions
Yeast_Interactome_Sig <- read.csv("/home/xavier/Desktop/Cell_interactome_stoichiometry/Yeast_Epistasis_Stoichiometries/Third_Run/R_Calc_Stoichiometries/data/The_Yeast_Interactome_Edges.csv", sep = ";")


## Make functions for the mapping

# Map stoichiometries (source prot to bait prot; target prot to prey prot)
mapStoichiometries <- function(index, IntDF, StoichDF, rep){
  bait <- IntDF$source[index] # Get bait ORF
  prey <- IntDF$target[index] # Get source ORF
  
  # Replace dashed with dots because of character restriction in column names
  bait <- stringr::str_replace(bait, pattern = "-", replacement = "\\.")
  
  prey <- unlist(str_split(prey, pattern = ";")) # Some preys correspond to multiple genes
  for (j in prey){
    if(length(which(StoichDF$ORF == j))>=1){ # If Prey is detected in the assay
      
      prey_index <- which(StoichDF$ORF == j) # Find Index of the prey proteingroup
      Stoich <- StoichDF[prey_index, paste0(bait, "_", rep)] # Get the Interaction Stoichiometry for the bait::prey pair
      
      if(!is.null(Stoich)){return(unlist(unname(Stoich)))} # If there was no Interaction Stoichiometry value, return NA
      
    }
  }
  return(NA)
}


## Apply Functions
Yeast_Interactome_Sig$Stoichiometry1 <- 
  sapply(1:nrow(Yeast_Interactome_Sig),
        FUN = mapStoichiometries,
        IntDF = Yeast_Interactome_Sig,
        StoichDF = Stoichiometries_DF,
        rep = 1)

Yeast_Interactome_Sig$Stoichiometry2 <- 
  sapply(1:nrow(Yeast_Interactome_Sig),
         FUN = mapStoichiometries,
         IntDF = Yeast_Interactome_Sig,
         StoichDF = Stoichiometries_DF,
         rep = 2)

## Make new column for combined stoichometries
Yeast_Interactome_Sig$IntStoichiometry <- NA

## If Stoichiometry 1 is negative, take stoichiometry 2
Yeast_Interactome_Sig$IntStoichiometry[(Yeast_Interactome_Sig$Stoichiometry1 < 0) & (Yeast_Interactome_Sig$Stoichiometry2 > 0) & !is.na(Yeast_Interactome_Sig$Stoichiometry2)] <- 
  Yeast_Interactome_Sig$Stoichiometry2[(Yeast_Interactome_Sig$Stoichiometry1 < 0) & (Yeast_Interactome_Sig$Stoichiometry2 > 0) & !is.na(Yeast_Interactome_Sig$Stoichiometry2)]

## If stoichiometry 2 is negative, take stoichiometry 1
Yeast_Interactome_Sig$IntStoichiometry[(Yeast_Interactome_Sig$Stoichiometry2 < 0) & (Yeast_Interactome_Sig$Stoichiometry1 > 0) & !is.na(Yeast_Interactome_Sig$Stoichiometry2)] <- 
  Yeast_Interactome_Sig$Stoichiometry1[(Yeast_Interactome_Sig$Stoichiometry2 < 0) & (Yeast_Interactome_Sig$Stoichiometry1 > 0) & !is.na(Yeast_Interactome_Sig$Stoichiometry2)]

## If both stoichiometries are positive, take the mean
Yeast_Interactome_Sig$IntStoichiometry[(Yeast_Interactome_Sig$Stoichiometry1 > 0) & (Yeast_Interactome_Sig$Stoichiometry2 > 0) & !is.na(Yeast_Interactome_Sig$Stoichiometry2)] <- 
  (Yeast_Interactome_Sig$Stoichiometry1[(Yeast_Interactome_Sig$Stoichiometry1 > 0) & (Yeast_Interactome_Sig$Stoichiometry2 > 0) & !is.na(Yeast_Interactome_Sig$Stoichiometry2)] + 
     Yeast_Interactome_Sig$Stoichiometry2[(Yeast_Interactome_Sig$Stoichiometry1 > 0) & (Yeast_Interactome_Sig$Stoichiometry2 > 0) & !is.na(Yeast_Interactome_Sig$Stoichiometry2)])/2

# Remove NA values
Yeast_Interactome_Stoich <- Yeast_Interactome_Sig[!is.na(Yeast_Interactome_Sig$IntStoichiometry),]

# Select only useful columns
Yeast_Interactome_Stoich <- Yeast_Interactome_Stoich %>% dplyr::select(source, target,IntStoichiometry)

#### Calculate Abundance Stoichiometry ####

# Add copy numbers for source
Yeast_Interactome_Stoich <-
  merge(x = Yeast_Interactome_Stoich,
        y = Yeast_protein_abundances[, c(1,5)],
        by.x = "source",
        by.y = "Systematic Name",
        all.x = TRUE,
        all.y = FALSE)

colnames(Yeast_Interactome_Stoich)[4] <- "source_copy_number"

# Add copy numbers for target
Yeast_Interactome_Stoich <-
  merge(x = Yeast_Interactome_Stoich,
        y = Yeast_protein_abundances[, c(1,5)],
        by.x = "target",
        by.y = "Systematic Name",
        all.x = TRUE,
        all.y = FALSE)

colnames(Yeast_Interactome_Stoich)[5] <- "target_copy_number"

## Calculate Abundance Stoichiometry
Yeast_Interactome_Stoich$Abundance_Stoichiometry <- Yeast_Interactome_Stoich$target_copy_number/Yeast_Interactome_Stoich$source_copy_number


#### Export data ####

## Reorganize and rename columns to make them neat

# Rename Interaction Stoichiometry column
colnames(Yeast_Interactome_Stoich)[3] <- "Interaction_Stoichiometry"

# Order columns
Yeast_Interactome_Stoich <-
  Yeast_Interactome_Stoich %>%
  dplyr::select(source, target, Abundance_Stoichiometry, Interaction_Stoichiometry)

# Export
write.csv(Yeast_Interactome_Stoich, "/home/xavier/Desktop/Cell_interactome_stoichiometry/Yeast_Epistasis_Stoichiometries/Third_Run/R_Calc_Stoichiometries/results/Yeast_Interactome_Stoich.csv")

