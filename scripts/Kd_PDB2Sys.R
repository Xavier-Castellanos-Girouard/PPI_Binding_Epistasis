# Xavier Castellanos-Girouard
# 
# Match PDB to systematic name via Uniprot IDs
#
# Date First Created: April 23 2024
# Date last updated: April 29 2024

#### Import libraries ####

library(dplyr)
library(tidyr)
library(org.Sc.sgd.db)

#### Import data ####

dir <- "/home/xavier/Desktop/Cell_interactome_stoichiometry/Yeast_Epistasis_Stoichiometries/Third_Run/"

# Import PDB to Uniprot conversion table for PDBbind Kd dataset
PDB_2_Uniprot_DF <- read.csv(paste0(dir, "R_Calc_Kd_Benchmark/data/IntermediateFiles/PDB2Uniprot.csv"), row.names = 1)
row.names(PDB_2_Uniprot_DF) <- NULL # Reset row names

# Import PDB bind PPI file with extracted information
pdb_kd_PPI_DF <- read.csv(paste0(dir, "R_Calc_Kd_Benchmark/data/IntermediateFiles/Kd_PDB_list.csv"), row.names = 1)
row.names(pdb_kd_PPI_DF) <- NULL # Reset row names

# Import Uniprot to HGNC gene symbol conversion table
human_Uniprot2HGNC_DF <- read.csv(paste0(dir, "R_Calc_Kd_Benchmark/data/uniprotkb_database_HGNC_2024_02_28.tsv"), sep = "\t")

# Format Uniprot to HGNC conversion table
human_Uniprot2HGNC_DF <- 
  human_Uniprot2HGNC_DF %>%
  tidyr::separate_longer_delim(cols = "Gene.Names",
                               delim = " ")

#### Yeast ####

### Filter for yeast interactions

# Find instances where interactor 1 is from yeast
spec1_isYeast <- grepl(pattern = "CEREVISIAE", x = pdb_kd_PPI_DF$species1_name)
# Find instances where interactor 2 is from yeast
spec2_isYeast <- grepl(pattern = "CEREVISIAE", x = pdb_kd_PPI_DF$species2_name)

# Only keep interactions where both species were yeast
pdb_kd_PPI_yeast_DF <- pdb_kd_PPI_DF[spec1_isYeast & spec2_isYeast,]

## Merge Uniprot names to PDB
pdb_kd_PPI_yeast_DF <-
  merge(x = pdb_kd_PPI_yeast_DF,
        y = PDB_2_Uniprot_DF,
        by.x = "PDB_ID",
        by.y = "PDB_id",
        all.x = FALSE,
        all.y = FALSE)


### Convert Uniprot IDs to Systematic ORFs
bimap <- org.Sc.sgdUNIPROT # Get Uniprot to sys ORF bimap

# Get the Systematic ORF IDs that are mapped to a Uniprot ID
mapped_genes <- mappedkeys(bimap)

# Convert map to list. ORFs as entry names, Uniprot IDs as entry values
list_map <- as.list(bimap[mapped_genes])

# Some ORFs have multiple Uniprot IDs, collapse them into a single string
list_map <- sapply(list_map, FUN = paste0, collapse = ",")

# Construct a conversion dataframe
Uniprot_2_ORF_DF <- 
  data.frame(ORF_ID = names(list_map),
             Uniprot_ID = unlist(unname(list_map)))

# Separate ORF that have more than one Uniprot ID into new rows
Uniprot_2_ORF_DF <-
  Uniprot_2_ORF_DF %>%
  tidyr::separate_longer_delim(cols = "Uniprot_ID",
                               delim = ",")

# Add systematic ORF ids for interactor 1
pdb_kd_PPI_yeast_DF <-
  merge(x = pdb_kd_PPI_yeast_DF,
        y = Uniprot_2_ORF_DF,
        by.x = "protein1_UniprotID",
        by.y = "Uniprot_ID",
        all.x = FALSE,
        all.y = FALSE)
colnames(pdb_kd_PPI_yeast_DF)[15] <- "protein1_ORF"

# Add systematic ORF ids for interactor 2
pdb_kd_PPI_yeast_DF <-
  merge(x = pdb_kd_PPI_yeast_DF,
        y = Uniprot_2_ORF_DF,
        by.x = "protein2_UniprotID",
        by.y = "Uniprot_ID",
        all.x = FALSE,
        all.y = FALSE)
colnames(pdb_kd_PPI_yeast_DF)[16] <- "protein2_ORF"

## Rearrange columns and export
pdb_kd_PPI_yeast_DF <-
  pdb_kd_PPI_yeast_DF %>%
  dplyr::select(protein1_ORF,
                protein2_ORF,
                Type,
                Value,
                Unit,
                PDB_ID,
                protein1_UniprotID,
                protein2_UniprotID)

# Export
#write.csv(pdb_kd_PPI_yeast_DF, paste0(dir, "R_Calc_Kd_Benchmark/results/Yeast_Kd_literature.csv"))

#### Human ####

# Filter for human PPIs
human_pdb_kd_PPI_DF <-
  pdb_kd_PPI_DF %>%
  dplyr::filter((species1_TaxID=="9606") & (species2_TaxID=="9606"))

## Merge Kd DataFrame to Uniprot Dataframe
human_pdb_kd_PPI_DF <-
  merge(x = human_pdb_kd_PPI_DF,
        y = PDB_2_Uniprot_DF,
        by.x = "PDB_ID",
        by.y = "PDB_id",
        all.x = FALSE,
        all.y = FALSE)

## Map UniprotID to HGNC
human_pdb_kd_PPI_DF <-
  merge(x = human_pdb_kd_PPI_DF,
        y = human_Uniprot2HGNC_DF[,c(1,4)],
        by.x = "protein1_UniprotID",
        by.y = "Entry",
        all.x = FALSE,
        all.y = FALSE)

colnames(human_pdb_kd_PPI_DF)[15] <- "protein1_GeneName"

## Map UniprotID to HGNC
human_pdb_kd_PPI_DF <-
  merge(x = human_pdb_kd_PPI_DF,
        y = human_Uniprot2HGNC_DF[,c(1,4)],
        by.x = "protein2_UniprotID",
        by.y = "Entry",
        all.x = FALSE,
        all.y = FALSE)

colnames(human_pdb_kd_PPI_DF)[16] <- "protein2_GeneName"

### Export Dataframe

# Select only useful columns
human_pdb_kd_PPI_DF <-
  human_pdb_kd_PPI_DF %>%
  dplyr::select(protein1_GeneName,
                protein2_GeneName,
                Type,
                Value,
                Unit,
                PDB_ID,
                protein1_UniprotID,
                protein2_UniprotID)

# Export
#write.csv(human_pdb_kd_PPI_DF, paste0(dir, "R_Calc_Kd_Benchmark/results/Human_Kd_literature.csv"))

