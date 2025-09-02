# Xavier Castellanos-Girouard

# Paralog type annotation with duplication type, Kd, stoichiometries and sequence

# Date First Created: March 14 2024 
# Date Last Modified: May 20 2024


#### Import Libraries ####

library(dplyr)
library(tidyr)
library(readxl)
library(seqinr)

#### Import data ####

## Import Kellis WGD paralogs
paralog_files <- list.files("../data/Paralogs/Kellis_S8_DupGenes") # Get names from pairwise alignment file names
paralog_files <- gsub(paralog_files, pattern = "\\.txt", replacement = "") # Remove file type string
paralog_DF <- data.frame(syst1 = paralog_files) # Make Dataframe containing paralog names

paralog_DF <- 
  paralog_DF %>% 
  tidyr::separate_wider_delim(cols = "syst1", delim = "_", names = c("syst1", "syst2"))
paralog_DF$DupType <- "WGD"

## Import VanderSluis et al. paralogs
SSD_paralog_DF <- read.csv("../data/Paralogs/msb201082-sup-0002.xls", sep = "\t") # From VanderSluis et al. Mol. Sys. Biol. 2010
colnames(SSD_paralog_DF)[c(1,2)] <- c("syst1", "syst2")
SSD_paralog_DF <- SSD_paralog_DF %>% dplyr::filter(WGD == FALSE)
SSD_paralog_DF$DupType <- "SSD"

MS_Kd_GI_DF <- read.csv("../results/Kd_and_GI/Yeast_Kd_GI.csv", row.names = 1)
MS_Stoich_GI_DF <-  read.csv("../results/Stoichiometry_and_GI/Yeast_Stoich_GI.csv", row.names = 1)

# Combine SSD and WGD tables
col_keeps <- c("syst1", "syst2", "DupType")
paralog_DF <- 
  rbind(paralog_DF[,colnames(paralog_DF) %in% col_keeps],
        SSD_paralog_DF[,colnames(SSD_paralog_DF) %in% col_keeps])

# Remove duplicates. Discrepencies in SSD vs WGD are resolved by the annotation
# given on SGD (as of March 2024). These all suggest WGD instead of SSD
paralog_DF <- paralog_DF[!duplicated(paralog_DF[,c(1,2)]),]

paralog_list <- unique(c(paralog_DF$syst1, paralog_DF$syst2))

write.csv(paralog_DF, "../results/IntermediateFiles/Paralog_DuplicationType.csv")

#### Annotation GI/Kd dataset with Duplication type (WGD or SSD) ####


### 1. Case where none of the interacting genes have a paralog

## Keep instances where none of the genes have paralogs
MS_Kd_GI_NoPara_DF <- 
  MS_Kd_GI_DF[!(MS_Kd_GI_DF$source %in% paralog_list) & 
                !(MS_Kd_GI_DF$target %in% paralog_list),]

MS_Kd_GI_NoPara_DF$Type <- "NoParalog"

# Add duplication type
MS_Kd_GI_NoPara_DF$Source_DupType <- "None"
MS_Kd_GI_NoPara_DF$Target_DupType <- "None"



### 2. Case where one (but not the other) gene has a paralog

## Keep instances where one or the other gene has a paralog, but not both.
MS_Kd_GI_OnePara_DF <- 
  MS_Kd_GI_DF[xor(MS_Kd_GI_DF$source %in% paralog_list,
                  MS_Kd_GI_DF$target %in% paralog_list),]

MS_Kd_GI_OnePara_DF$Type <- "OneParalog"

# Add duplicate type
Source_DupType_vec <- c()
Target_DupType_vec <- c()
for (i in 1:nrow(MS_Kd_GI_OnePara_DF)){
  ORF1 <- MS_Kd_GI_OnePara_DF$source[i] # Get source ORF at index i
  ORF2 <- MS_Kd_GI_OnePara_DF$target[i] # Get target ORF at index i
  
  #print(c(ORF1, ORF2))
  
  ## Get duplication types of source and target orfs using paralog table
  # For source
  match1_ORF1 <- paralog_DF$DupType[which(paralog_DF$syst1 == ORF1)]
  match2_ORF1 <- paralog_DF$DupType[which(paralog_DF$syst2 == ORF1)]
  # For target
  match1_ORF2 <- paralog_DF$DupType[which(paralog_DF$syst1 == ORF2)]
  match2_ORF2 <- paralog_DF$DupType[which(paralog_DF$syst2 == ORF2)]
  
  #print(c(match1, match2))
  
  if (length(match1_ORF1) > 0){ # If source has match, assign duptype to source and "none" type to target
    Source_DupType_vec <- c(Source_DupType_vec, match1_ORF1)
    Target_DupType_vec <- c(Target_DupType_vec, "None")
  } else if (length(match2_ORF1) > 0){
    Source_DupType_vec <- c(Source_DupType_vec, match2_ORF1)
    Target_DupType_vec <- c(Target_DupType_vec, "None")
  } else if (length(match1_ORF2) > 0){
    Source_DupType_vec <- c(Source_DupType_vec, "None")
    Target_DupType_vec <- c(Target_DupType_vec, match1_ORF2)
  } else if (length(match2_ORF2) > 0){
    Source_DupType_vec <- c(Source_DupType_vec, "None")
    Target_DupType_vec <- c(Target_DupType_vec, match2_ORF2)
  }
  
}

MS_Kd_GI_OnePara_DF$Source_DupType <- Source_DupType_vec
MS_Kd_GI_OnePara_DF$Target_DupType <- Target_DupType_vec



### 3. Case where both gene have paralogs

## Keep instances where both genes have paralogs (includes paralogous heteromers)
MS_Kd_GI_TwoPara_DF <- 
  MS_Kd_GI_DF[(MS_Kd_GI_DF$source %in% paralog_list) & 
                (MS_Kd_GI_DF$target %in% paralog_list),]

MS_Kd_GI_TwoPara_DF$Type <- "TwoParalog"

# Add duplicate type
Source_DupType_vec <- c()
Target_DupType_vec <- c()
for (i in 1:nrow(MS_Kd_GI_TwoPara_DF)){
  ORF1 <- MS_Kd_GI_TwoPara_DF$source[i] # Get source ORF at index i
  ORF2 <- MS_Kd_GI_TwoPara_DF$target[i] # Get target ORF at index i
  
  #print(c(ORF1, ORF2))
  
  ## Get duplication types of source and target orfs using paralog table
  # For source
  match1_ORF1 <- paralog_DF$DupType[which(paralog_DF$syst1 == ORF1)]
  match2_ORF1 <- paralog_DF$DupType[which(paralog_DF$syst2 == ORF1)]
  # For target
  match1_ORF2 <- paralog_DF$DupType[which(paralog_DF$syst1 == ORF2)]
  match2_ORF2 <- paralog_DF$DupType[which(paralog_DF$syst2 == ORF2)]
  
  #print(c(match1_ORF1, match2_ORF1, match1_ORF2, match2_ORF2))
  
  Source_DupType_vec <- c(Source_DupType_vec, c(match1_ORF1, match2_ORF1))
  Target_DupType_vec <- c(Target_DupType_vec, c(match1_ORF2, match2_ORF2))
  
}

MS_Kd_GI_TwoPara_DF$Source_DupType <- Source_DupType_vec
MS_Kd_GI_TwoPara_DF$Target_DupType <- Target_DupType_vec


### 4. Case where interacting genes are paralogs of each other (paralogous heteromers)


## Keep instances where pair of interactors are paralogous
MS_Kd_GI_DF_1 <- 
  merge(MS_Kd_GI_DF, 
        #paralog_DF[, c(2, 3, 7, 8)], # Diss Dataframe
        paralog_DF[, c(1,2,3)], # Kellis dataframe
        by.x = c("source", "target"),
        by.y = c("syst1", "syst2"),
        all.x = FALSE,
        all.y = FALSE)

MS_Kd_GI_DF_2 <- 
  merge(MS_Kd_GI_DF, 
        #paralog_DF[, c(2, 3, 7, 8)], # Diss Dataframe
        paralog_DF[, c(1,2,3)], # Kellis dataframe
        by.x = c("source", "target"),
        by.y = c("syst2", "syst1"),
        all.x = FALSE,
        all.y = FALSE)

MS_Kd_GI_IntPara_DF <-
  rbind(MS_Kd_GI_DF_1,
        MS_Kd_GI_DF_2)

# Add duptypes to individual ORFs
MS_Kd_GI_IntPara_DF$Source_DupType <- MS_Kd_GI_IntPara_DF$DupType
MS_Kd_GI_IntPara_DF$Target_DupType <- MS_Kd_GI_IntPara_DF$DupType
MS_Kd_GI_IntPara_DF <- MS_Kd_GI_IntPara_DF %>% dplyr::select(-c(DupType))


MS_Kd_GI_IntPara_DF$Type <- "IntParalog"


## Merge
MS_Kd_GI_AllPara_DF <- 
  rbind(MS_Kd_GI_NoPara_DF,
        MS_Kd_GI_OnePara_DF,
        MS_Kd_GI_IntPara_DF, # Kellis dataframe
        MS_Kd_GI_TwoPara_DF)

MS_Kd_GI_AllPara_DF <- 
  MS_Kd_GI_AllPara_DF[!duplicated(MS_Kd_GI_AllPara_DF[,c("source", "target")]),]


write.csv(MS_Kd_GI_AllPara_DF, "../results/Paralog_Analyses/Kd_GI_ParalogType.csv", row.names = FALSE)


#### Merge Stoichiometries and Paralogs ####


# Keep instances where one or the other gene has a paralog, but not both.
MS_Stoich_GI_OnePara_DF <- 
  MS_Stoich_GI_DF[xor(MS_Stoich_GI_DF$source %in% paralog_list,
                         MS_Stoich_GI_DF$target %in% paralog_list),]

MS_Stoich_GI_OnePara_DF$Type <- "OneParalog"

# Keep instances none of the genes have paralogs
MS_Stoich_GI_NoPara_DF <- 
  MS_Stoich_GI_DF[!(MS_Stoich_GI_DF$source %in% paralog_list) & 
                       !(MS_Stoich_GI_DF$target %in% paralog_list),]

MS_Stoich_GI_NoPara_DF$Type <- "NoParalog"

# Keep instances where both genes have paralogs (not necessarily paralogous between them)
MS_Stoich_GI_TwoPara_DF <- 
  MS_Stoich_GI_DF[(MS_Stoich_GI_DF$source %in% paralog_list) & 
                           (MS_Stoich_GI_DF$target %in% paralog_list),]

MS_Stoich_GI_TwoPara_DF$Type <- "TwoParalog"

# Keep instances where pair of interactors are paralogous
MS_Stoich_GI_DF_1 <- 
  merge(MS_Stoich_GI_DF, 
        #paralog_DF[, c(2, 3, 7, 8)], # Diss
        paralog_DF[, c(1, 2)], # Kellis
        by.x = c("source", "target"),
        by.y = c("syst1", "syst2"),
        all.x = FALSE,
        all.y = FALSE)

MS_Stoich_GI_DF_2 <- 
  merge(MS_Stoich_GI_DF, 
        #paralog_DF[, c(2, 3, 7, 8)], # Diss
        paralog_DF[, c(1,2)], # Kellis
        by.x = c("source", "target"),
        by.y = c("syst2", "syst1"),
        all.x = FALSE,
        all.y = FALSE)

MS_Stoich_GI_IntPara_DF <-
  rbind(MS_Stoich_GI_DF_1,
        MS_Stoich_GI_DF_2)

MS_Stoich_GI_IntPara_DF$Type <- "IntParalog"

## Merge
MS_Stoich_GI_AllPara_DF <- 
  rbind(MS_Stoich_GI_NoPara_DF,
        MS_Stoich_GI_OnePara_DF,
        MS_Stoich_GI_IntPara_DF,
        MS_Stoich_GI_TwoPara_DF)


MS_Stoich_GI_AllPara_DF <- 
  MS_Stoich_GI_AllPara_DF[!duplicated(MS_Stoich_GI_AllPara_DF[,c("source", "target")]),]

write.csv(MS_Stoich_GI_AllPara_DF, "../results/Paralog_Analyses/Stoich_GI_ParalogType.csv")



#### Compute sequence divergence ####

### Calculate Sequence Divergence
# Note:
#   1) Make fasta files with paralog amino acid sequences
#   2) Align paralog proteins using ClustalW
#   3) Make fasta files with paralog nucleotide sequences
#   4) Use Protein Alignment to make a nucleotide alignment
#   5) Calculate sequence divergence

# Import amino acid sequences of yeast ORFs
ORF_Prot_sequences <- 
  read.fasta(file = "../data/Paralogs/orf_trans_all.fasta",
             seqtype = c("AA"))

# From protein sequence, make continuous string of characters instead of list of characters
ORF_Prot_sequences <- 
  lapply(ORF_Prot_sequences,
         FUN = paste0,
         collapse = "")

## Proteins must be exported into fasta file then ClustalW is used to align pairs
for (i in 1:nrow(paralog_DF)){ # For every paralog pair
  paralog1 <- paralog_DF$syst1[i] # Get systematic ORF name of paralog 1
  paralog2 <- paralog_DF$syst2[i] # Get systematic ORF name of paralog 2
  
  paralog1_seq <- ORF_Prot_sequences[paralog1] # Get protein Seq of paralog 1
  paralog2_seq <- ORF_Prot_sequences[paralog2] # Get protein Seq of paralog 1
  
  #print(paste(paralog1, paralog2))
  
  # Write sequences in fasta file
  seqinr::write.fasta(sequences = list(paralog1_seq, paralog2_seq),
                      as.string = TRUE,
                      names = c(paralog1, paralog2),
                      file.out = paste0("../results/IntermediateFiles/Paralog_proteinSeq_ClustalW_input/", 
                                        paralog1,
                                        "_",
                                        paralog2,
                                        ".fasta"))
}

### ** At this point, ClustalW was used to perform pairwise sequence alignment
### ** Execute Execute_ProtAlign_ClustalW.sh
