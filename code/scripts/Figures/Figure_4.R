# Xavier Castellanos-Girouard

# Date First Created: July 8 2024
# Date Last Modified: July 9 2024

#### Import Libraries ####

library(dplyr)
library(tidyr)
library(org.Sc.sgd.db)

#### Import Data ####

Module_Overlap_DF <- read.csv('../results/Topological_Analyses/GI_PPI_optimal_module_overlap.csv', row.names = 1)
row.names(Module_Overlap_DF) <- NULL

Connector_Overlap_DF <- read.csv('../results/Topological_Analyses/Shortest_Path_Connector_Overlap.csv', row.names = 1)
row.names(Connector_Overlap_DF) <- NULL

#### Functions ####

# Convert string (from python list) to vector of ORFs
StringedList_to_vector <-
  function(string){
    new_string = stringr::str_replace(string = string, 
                                      pattern = "\\]", 
                                      replacement = "")
    
    new_string = stringr::str_replace(string = new_string, 
                                      pattern = "\\[", 
                                      replacement = "")
    
    new_string = stringr::str_replace_all(string = new_string,
                                          pattern = "'",
                                          replacement = "")
    
    new_string = stringr::str_replace_all(string = new_string,
                                          pattern = " ",
                                          replacement = "")
    
    new_string = stringr::str_split_1(string = new_string,
                                      pattern = ",")
    
    return(unname(unlist(new_string)))
  }

StringedSet_to_vector <-
  function(string){
    new_string = stringr::str_replace(string = string, 
                                      pattern = "\\}", 
                                      replacement = "")
    
    new_string = stringr::str_replace(string = new_string, 
                                      pattern = "\\{", 
                                      replacement = "")
    
    new_string = stringr::str_replace_all(string = new_string,
                                          pattern = "'",
                                          replacement = "")
    
    new_string = stringr::str_replace_all(string = new_string,
                                          pattern = " ",
                                          replacement = "")
    
    new_string = stringr::str_split_1(string = new_string,
                                      pattern = ",")
    
    return(unname(unlist(new_string)))
  }


ORF_Set_2_GeneName <- function(ORF_set){
  
  ORF_set <- unname(unlist(ORF_set))
  
  # If ORF set is a vector
  if (is.vector(ORF_set)){
    return(unname(unlist(mapped_genes_list[ORF_set])))
  }
  
  # Some ORF sets are empty
  else if (ORF_set == ""){
    return("")
  }
  
  # If ORF set is just one ORF (string type)
  else{
    return(mapped_genes_list[[ORF_set]])
  }
}

#### Format Dataframe Columns ####

## Modules
Module_Overlap_DF <-
  Module_Overlap_DF %>%
  dplyr::select(Cluster_ID, 
                CPX_ID, 
                Common_Module_ORFs, 
                GI_Module_ORFs, 
                PPI_Module_ORFs)

Module_Overlap_DF <-
  Module_Overlap_DF %>%
  dplyr::filter(Common_Module_ORFs != "set()") %>%
  dplyr::filter(CPX_ID != "None")

colnames(Module_Overlap_DF) <- 
  c("GI_Module_ID", "PPI_Module_ID",
    "Common_Module_ORFs", "GI_Module_ORFs", "PPI_Module_ORFs")

## Connectors
Connector_Overlap_DF <-
  Connector_Overlap_DF %>%
  dplyr::select(Source_GI_Cluster,
                Source_CPX_ID,
                Target_GI_Cluster,
                Target_CPX_ID,
                Shortest_Path_Overlap)

colnames(Connector_Overlap_DF) <- 
  c("Source_GI_Module_ID", "Source_PPI_Module_ID",
    "Target_GI_Module_ID", "Target_PPI_Module_ID", 
    "Shortest_Path_Overlap")


#### Get list for ORF 2 Name conversion ####

## Bimap interface:
Bimap_GeneName <- org.Sc.sgdGENENAME

## Get the gene names that are mapped to an ORF identifier
mapped_genes <- mappedkeys(Bimap_GeneName)

# Convert to a list
mapped_genes_list <- as.list(Bimap_GeneName[mapped_genes])

#### Convert Strings from python list to vector of ORFs ####

### For Module Overlaps
# Convert GI module ORFs
Module_Overlap_DF$GI_Module_ORFs <-
  sapply(Module_Overlap_DF$GI_Module_ORFs,
         FUN = StringedList_to_vector)

Module_Overlap_DF$PPI_Module_ORFs <-
  sapply(Module_Overlap_DF$PPI_Module_ORFs,
         FUN = StringedList_to_vector)

Module_Overlap_DF$Common_Module_ORFs <-
  sapply(Module_Overlap_DF$Common_Module_ORFs,
         FUN = StringedSet_to_vector)

#Module_Overlap_DFDifferent_ORFs <-
#  sapply(Module_Overlap_DF$Different_ORFs,
#         FUN = StringedList_to_vector)

### For Shortest Paths Overlaps
#Connector_Overlap_DF$GI_path <-
#  sapply(Connector_Overlap_DF$GI_path,
#         FUN = StringedList_to_vector)

#Connector_Overlap_DF$PPI_path <-
#  sapply(GConnector_Overlap_DF$PPI_path,
#         FUN = StringedList_to_vector)

#Connector_Overlap_DF$Path_Overlap <-
#  sapply(Connector_Overlap_DF$Path_Overlap,
#         FUN = StringedList_to_vector)


#### Convert ORFs to Gene Names ####

### For Module Overlaps
Module_Overlap_DF$GI_Module_ORFs <-
  sapply(Module_Overlap_DF$GI_Module_ORFs,
         ORF_Set_2_GeneName)

Module_Overlap_DF$PPI_Module_ORFs <-
  sapply(Module_Overlap_DF$PPI_Module_ORFs,
         ORF_Set_2_GeneName)

Module_Overlap_DF$Common_Module_ORFs <-
  sapply(Module_Overlap_DF$Common_Module_ORFs,
         ORF_Set_2_GeneName)

#Module_Overlap_DF$Different_ORFs <-
#  sapply(Module_Overlap_DF$Different_ORFs,
#         ORF_Set_2_GeneName)


### For shortest Paths Overlaps
#Connector_Overlap_DF$GI_path <-
#  sapply(Connector_Overlap_DF$GI_path,
#         ORF_Set_2_GeneName)

#Connector_Overlap_DF$PPI_path <-
#  sapply(Connector_Overlap_DF$PPI_path,
#         ORF_Set_2_GeneName)

Connector_Overlap_DF$Shortest_Path_Overlap_GeneName <-
  sapply(Connector_Overlap_DF$Shortest_Path_Overlap,
         ORF_Set_2_GeneName)

Connector_Overlap_DF$Shortest_Path_Overlap_GeneName <- as.character(Connector_Overlap_DF$Shortest_Path_Overlap_GeneName)


#### Convert Into Human Easily Readable Form ####

### For Module Overlaps

# Convert vectors into strings
Module_Overlap_DF$GI_Module_ORFs <-
  sapply(Module_Overlap_DF$GI_Module_ORFs,
         FUN = paste,
         collapse = ", ")

Module_Overlap_DF$PPI_Module_ORFs <-
  sapply(Module_Overlap_DF$PPI_Module_ORFs,
         FUN = paste,
         collapse = ", ")

Module_Overlap_DF$Common_Module_ORFs <-
  sapply(Module_Overlap_DF$Common_Module_ORFs,
         FUN = paste,
         collapse = ", ")

#Module_Overlap_DF$Different_ORFs <-
#  sapply(Module_Overlap_DF$Different_ORFs,
#         FUN = paste,
#         collapse = ", ")

### For Shortest Paths Overlaps
#Connector_Overlap_DF$GI_path <-
#  sapply(Connector_Overlap_DF$GI_path,
#         FUN = paste,
#         collapse = ", ")

#Connector_Overlap_DF$PPI_path <-
#  sapply(Connector_Overlap_DF$PPI_path,
#         FUN = paste,
#         collapse = ", ")


#### Export ####

write.csv(Module_Overlap_DF, file = "../results/Figures/Main_Figures/Module_Overlap_Readable.csv")
write.csv(Connector_Overlap_DF, file = "../results/Figures/Main_Figures/Connector_Overlap_Readable.csv")
