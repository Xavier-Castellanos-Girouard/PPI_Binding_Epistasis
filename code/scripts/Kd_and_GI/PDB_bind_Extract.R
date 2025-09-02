# Xavier Castellanos-Girouard
# 
# Getting a set of PPI dissociation constants from PDBbind database
#
# Date First Created: June 12 2023
# Date last updated: April 23 2024

#### Import libraries ####
library(stringr)
library(dplyr)
library(tidyr)

#### Import data #### 

# PDBbind
index_data_DF <- 
  read.csv("../data/PDBbind_vs2020/INDEX_general_PP.2020.txt",
           sep = "\t")

#### Functions ####

## Extract Binding information from PDBbind index file. 
## Input is PDBbind index file (dataframe).
## Output is dataframe with PDB ID, type of binding measurement, binding value, and units
Extract_PDBbind_index <-
  function(index_DF){
    ## Extract PDB ID
    PDB_ID <- 
      stringr::str_extract(index_DF[6:nrow(index_DF),1],
                           pattern = "^....")
    
    ## Extract all quantitative data related to binding
    Binding_Info <- 
      stringr::str_extract(index_DF[6:nrow(index_DF),1],
                           pattern = "(IC50|Kd|Ki)(=|~|<|>)[^M]*")
    
    ## Get the Binding value type (IC50 or Kd or Ki)
    Binding_value_type <-
      stringr::str_extract(Binding_Info,
                           pattern = "(IC50|Kd|Ki)")
    
    ## Get the binding value itself
    Binding_value <-
      stringr::str_match(Binding_Info,
                         pattern = "(IC50|Kd|Ki)(=|~|<|>)\\s*(.*)\\s*(n|p|u|f|m)$")
    
    Binding_value <- as.numeric(as.vector(Binding_value[,4]))
    
    ## Get binding unit
    Binding_unit <- 
      stringr::str_match(Binding_Info,
                         pattern = "[a-z]$")
    
    Binding_unit <- paste0(Binding_unit, "M")
    
    ## Put all into a common dataframe 
    BindingData_DF <- 
      data.frame(ID = PDB_ID,
                 Type = Binding_value_type,
                 Value = Binding_value,
                 Unit = Binding_unit)
    
    return(BindingData_DF)
  }

# The following function converts molar concentrations to uM. Input should be 
# two vectors of equal size; one containing numeric values of concentrations 
# (vals_vec; numeric vector type), the other containing units of concentration 
# (unit_vec; character vector type). Output should be vector of numeric values.
convert_conc_to_M <- 
  function(vals_vec, unit_vec){
    vals_femto_bool <- (unit_vec == "fM")
    vals_pico_bool <- (unit_vec == "pM")
    vals_nano_bool <- (unit_vec == "nM")
    vals_micro_bool <- (unit_vec == "uM")
    vals_milli_bool <- (unit_vec == "mM")
    
    vals_vec[vals_femto_bool] <- vals_vec[vals_femto_bool]*10^-15
    vals_vec[vals_pico_bool] <- vals_vec[vals_pico_bool]*10^-12
    vals_vec[vals_nano_bool] <- vals_vec[vals_nano_bool]*10^-9
    vals_vec[vals_micro_bool] <- vals_vec[vals_micro_bool]*10^-6
    vals_vec[vals_milli_bool] <- vals_vec[vals_milli_bool]*10^-3
    
    return(vals_vec)
  }

## Extracting name of molecule, species or other from vector of strings
PDB_str_extract <-
  function(str_vector, row_pattern, str_regex, isSynonym){
    
    # Extract name (molecule, species etc.) from vector
    name_str <-
      str_vector[grepl(pattern = row_pattern, x = str_vector)] %>% # Find row with str type (e.g. MOLECULE:)
      stringr::str_squish() %>% # Replace whitespace with space, remove outer spaces
      stringr::str_match(pattern = str_regex) %>% # Extract name
      stringr::str_remove(pattern = "^:(\\s?)") %>% # Remove unwanted colon
      stringr::str_remove(pattern = "(\\s?);$") # Remove unwanted semi-colon
    
    
    if (isSynonym){ # If extracted string is a synonym
      name_str <-
        name_str %>%
        stringr::str_split(pattern = ",\\s") %>% # Split synonyms into list
        unlist() %>% # Unlist to remove list index
        as.character() # Convert to character vector
    }
    
    return(name_str)
  }


## Extracting names of interacting proteins from PDB file
# Input is PDB identifier and directory to PDBbind PDB files
getPDB_PPI <- 
  function(pdb_ID, PDB_dir){
    
    # Set path to PDB
    pdb_path <- paste0(PDB_dir, pdb_ID, ".ent.pdb")
    
    ## Import PBD file (line = dataframe row), cols not distinguished.
    pdb_file <- read.delim(pdb_path, sep = "\n")
    
    
    ### Only keep compound and source sections
    
    # Get vector containing information on compound information
    pdb_COMPND <- 
      as.character(pdb_file[grepl(pattern = "^COMPND",
                                  x = as.character(pdb_file[, 1])),1])
    
    # Get boolean vector for position of source information
    pdb_SOURCE <- 
      as.character(pdb_file[grepl(pattern = "^SOURCE", 
                                  x = as.character(pdb_file[, 1])), 1])
    
    # If there are more than 3 molecules (proteins), discard the file
    not_2Mol <- (2 != sum(grepl(pattern = "MOL_ID:", x = pdb_COMPND)))
    if (not_2Mol){return(NA)}
    
    
    ### Get name of proteins (and synonyms)
    
    # Get indices indicating start of molecule information
    mol_indices <- which(grepl(pattern = "MOL_ID:", x = pdb_COMPND) == TRUE)
    
    # Extract info for both molecules separately
    mol1 <- pdb_COMPND[mol_indices[1]:(mol_indices[2]-1)]
    mol2 <- pdb_COMPND[mol_indices[2]:length(pdb_COMPND)]
    
    ## Get molecule name for molecule 1 and molecule 2
    mol1_name <- PDB_str_extract(mol1, "MOLECULE:", ":\\s[A-Za-z0-9][A-Za-z0-9\\s\\-\\(\\)]*;?", FALSE)
    mol2_name <- PDB_str_extract(mol2, "MOLECULE:", ":\\s[A-Za-z0-9][A-Za-z0-9\\s\\-\\(\\)]*;?", FALSE)
    
    ## Get synonyms for molecule 1 and molecule 2
    mol1_synonyms <- PDB_str_extract(mol1, "SYNONYM", ":\\s[A-Za-z0-9][A-Za-z0-9\\s\\-,\\(\\)]*;?", TRUE)
    mol2_synonyms <- PDB_str_extract(mol2, "SYNONYM", ":\\s[A-Za-z0-9][A-Za-z0-9\\s\\-,\\(\\)]*;?", TRUE)
    
    
    ### Get ID and scientific name of species
    
    # Get indices indicating start of each species information
    spec_indices <- which(grepl(pattern = "MOL_ID:", x = pdb_SOURCE) == TRUE)
    
    # Extract info for both molecules separately
    spec1 <- pdb_SOURCE[spec_indices[1]:(spec_indices[2]-1)]
    spec2 <- pdb_SOURCE[spec_indices[2]:length(pdb_SOURCE)]
    
    ## Get scientific name of species for molecule 1 and molecule 2
    spec1_name <- PDB_str_extract(spec1, "[^:]\\sORGANISM_SCIENTIFIC", ":\\s[A-Za-z0-9][A-Za-z0-9\\s\\-\\.\\(\\)]*;?", FALSE)
    spec2_name <- PDB_str_extract(spec2, "[^:]\\sORGANISM_SCIENTIFIC", ":\\s[A-Za-z0-9][A-Za-z0-9\\s\\-\\.\\(\\)]*;?", FALSE)
    
    ## Get Taxonomy ID for molecule 1 and molecule 2
    spec1_taxID <- PDB_str_extract(spec1, "ORGANISM_TAXID", ":\\s[A-Za-z0-9][A-Za-z0-9\\s\\-\\.\\(\\)]*;?", FALSE)
    spec2_taxID <- PDB_str_extract(spec2, "ORGANISM_TAXID", ":\\s[A-Za-z0-9][A-Za-z0-9\\s\\-\\.\\(\\)]*;?", FALSE)
    
    # Catch exceptions
    if (length(mol1_name) == 0){mol1_name = NA}
    if (length(mol2_name) == 0){mol2_name = NA}
    if (length(mol1_synonyms) == 0){mol1_synonyms = NA}
    if (length(mol2_synonyms) == 0){mol2_synonyms = NA}
    if (length(spec1_name) == 0){spec1_name = NA}
    if (length(spec2_name) == 0){spec2_name = NA}
    if (length(spec1_taxID) == 0){spec1_taxID = NA}
    if (length(spec2_taxID) == 0){spec2_taxID = NA}
    
    # Construct dataframe from extracted information
    pdb_dataframe <-
      data.frame(PDB_ID = pdb_ID, 
                 protein1_name = mol1_name,
                 protein2_name = mol2_name,
                 protein1_synonyms = paste(mol1_synonyms, collapse = ","),
                 protein2_synonyms = paste(mol2_synonyms, collapse = ","),
                 species1_name = spec1_name,
                 species2_name = spec2_name,
                 species1_TaxID = spec1_taxID,
                 species2_TaxID = spec2_taxID)
    
    return(pdb_dataframe)
  }

#### Main ####

## Extract relevant data from PDBbind index file
Binding_DF <- Extract_PDBbind_index(index_data_DF)
Binding_DF <- Binding_DF[(Binding_DF$Type != "IC50") & (Binding_DF$Type != "Ki" ), ]

# Convert Kd values to uM
Binding_DF$Value <- convert_conc_to_M(Binding_DF$Value, Binding_DF$Unit)

# Modify units to uM to reflect conversion
Binding_DF$Unit <- "M"

## Extract PPI name and species information using PDB file names
pdb_PPI_DF_list <- 
  sapply(unique(Binding_DF$ID),
         FUN = getPDB_PPI,
         PDB_dir = "../data/PDBbind_vs2020/PDB_files/")

# Remove entries where information is not extractable (NA)
pdb_PPI_DF_list <- pdb_PPI_DF_list[!is.na(pdb_PPI_DF_list)]

# Combine PDB dataframes into one master dataframe
pdb_PPI_DF <- do.call(rbind, pdb_PPI_DF_list)

### Merge kd values to pdb files
pdb_PPI_Kd_DF <-
  merge(pdb_PPI_DF,
        Binding_DF,
        by.x = "PDB_ID",
        by.y = "ID",
        all.x = FALSE,
        all.y = FALSE)

### Only keep instances where proteins are either both from yeast or both from human

# Find instances where interactor 1 is from yeast
spec1_isYeast <- grepl(pattern = "CEREVISIAE", x = pdb_PPI_Kd_DF$species1_name)
# Find instances where interactor 2 is from yeast
spec2_isYeast <- grepl(pattern = "CEREVISIAE", x = pdb_PPI_Kd_DF$species2_name)

# Filter for human PPIs
pdb_PPI_Kd_DF <-
  pdb_PPI_Kd_DF %>%
  dplyr::filter(((species1_TaxID=="9606") & (species2_TaxID=="9606")) | (spec1_isYeast & spec2_isYeast))


### Quality control by retaining values within a range that is reasonable for most methods
pdb_PPI_Kd_DF <- 
  pdb_PPI_Kd_DF %>%
  dplyr::filter((Value >= 10^-9) & (Value <= 10^-3))

### Export for webscraping
write.csv(pdb_PPI_Kd_DF, "../results/IntermediateFiles/Kd_PDB_list.csv")
