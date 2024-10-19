# Xavier Castellanos-Girouard

# Date First Created: March 14 2024 
# Date Last Modified: May 20 2024


#### Import Libraries ####

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(readxl)
library(seqinr)

#### Import data ####

dir <- "/home/xavier/Desktop/Cell_interactome_stoichiometry/Yeast_Epistasis_Stoichiometries/Third_Run/"

## Import Kellis WGD paralogs
paralog_files <- list.files(paste0(dir, "R_Paralog_Analysis/data/Kellis_S8_DupGenes/nucleotide")) # Get names from pairwise alignment file names
paralog_files <- gsub(paralog_files, pattern = "\\.txt", replacement = "") # Remove file type string
paralog_DF <- data.frame(syst1 = paralog_files) # Make Dataframe containing paralog names

paralog_DF <- 
  paralog_DF %>% 
  tidyr::separate_wider_delim(cols = "syst1", delim = "_", names = c("syst1", "syst2"))
paralog_DF$DupType <- "WGD"

## Import VanderSluis et al. paralogs
SSD_paralog_DF <- read.csv(paste0(dir, "R_Paralog_Analysis/data/VanderSluis_2010_SuppTable1/msb201082-sup-0002.xls"), sep = "\t") # From VanderSluis et al. Mol. Sys. Biol. 2010
colnames(SSD_paralog_DF)[c(1,2)] <- c("syst1", "syst2")
SSD_paralog_DF <- SSD_paralog_DF %>% dplyr::filter(WGD == FALSE)
SSD_paralog_DF$DupType <- "SSD"

MS_Kd_GI_DF <- read.csv(paste0(dir, "R_Python_Kd_vs_GI/results/Yeast_Kd_GI.csv"), row.names = 1)
MS_Stoich_GI_DF <-  read.csv(paste0(dir, "R_Stoichiometries_vs_GI/results/Yeast_Stoich_GI.csv"), row.names = 1)

# Combine SSD and WGD tables
col_keeps <- c("syst1", "syst2", "DupType")
paralog_DF <- 
  rbind(paralog_DF[,colnames(paralog_DF) %in% col_keeps],
        SSD_paralog_DF[,colnames(SSD_paralog_DF) %in% col_keeps])

# Remove duplicates. Discrepencies in SSD vs WGD are resolved by the annotation
# given on SGD (as of March 2024). These all suggest WGD instead of SSD
paralog_DF <- paralog_DF[!duplicated(paralog_DF[,c(1,2)]),]

paralog_list <- unique(c(paralog_DF$syst1, paralog_DF$syst2))

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


#write.csv(MS_Kd_GI_AllPara_DF, paste0(dir, "R_Paralog_Analysis/results/Kd_GI_ParalogType.csv"))


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

#write.csv(MS_Stoich_GI_AllPara_DF, paste0(dir, "R_Paralog_Analysis/results/Stoich_GI_ParalogType.csv"))



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
  read.fasta(file = paste0(dir, "R_Paralog_Analysis/data/Yeast_ORF_seq/orf_trans_all.fasta"),
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
                      file.out = paste0(dir, "R_Paralog_Analysis/data/Paralog_proteinSeq_ClustalW_input/", 
                                        paralog1,
                                        "_",
                                        paralog2,
                                        ".fasta"))
}

### ** At this point, ClustalW was used to perform pairwise sequence alignment
### ** Execute Execute_ProtAlign_ClustalW.sh

### Read individual DNA sequences, and write fasta file containing both paralog sequences
ORF_Nucl_sequences <- 
  read.fasta(file = paste0(dir, "R_Paralog_Analysis/data/Yeast_ORF_seq/orf_coding_all.fasta"),
             seqtype = c("DNA"))

ORF_Nucl_sequences <- 
  lapply(ORF_Nucl_sequences,
         FUN = paste0,
         collapse = "")

for (i in 1:nrow(paralog_DF)){
  paralog1 <- paralog_DF$syst1[i]
  paralog2 <- paralog_DF$syst2[i]
  
  paralog1_seq <- ORF_Nucl_sequences[paralog1]
  paralog2_seq <- ORF_Nucl_sequences[paralog2]
  
  #print(paste(paralog1, paralog2))
  
  seqinr::write.fasta(sequences = list(paralog1_seq, paralog2_seq),
                      as.string = TRUE,
                      names = c(paralog1, paralog2),
                      file.out = paste0(dir, "R_Paralog_Analysis/data/reverseAlign_nucl_FastaFiles/", 
                                        paralog1,
                                        "_",
                                        paralog2,
                                        ".fasta"))
}

### Align nucleic acid sequences using amino acid sequences

bug <- c("YDL061C", "YLR388W", "YLR466C-B", "YLL066W-B")


for (i in 1:nrow(paralog_DF)){
  paralog1 <- paralog_DF$syst1[i]
  paralog2 <- paralog_DF$syst2[i]

  #print(paste(paralog1, paralog2))
  
  # Note: Some paralogs raise an error, skip these
  if ((paralog1 %in% bug) & (paralog2 %in% bug)){
    next
  }
  
  input_dir <- paste0(dir, "R_Paralog_Analysis/data/reverseAlign_nucl_FastaFiles/",
                      paralog1, 
                      "_",
                      paralog2,
                      ".fasta")
  
  prot_aln_dir <- paste0(dir, "R_Paralog_Analysis/data/ClustalW_Ouput/",
                         paralog1,
                         "_",
                         paralog2,
                         ".aln")
  
  output_dir <- paste0(dir, "R_Paralog_Analysis/data/reverseAlign_output_Files/",
                       paralog1,
                       "_",
                       paralog2,
                       ".fasta")
  
  seqinr::reverse.align(nucl.file = input_dir,
                        protaln.file = prot_aln_dir,
                        input.format = "clustal",
                        out.file = output_dir)
}



## Import alignments

Nucl_align_list <- list()

for (i in 1:nrow(paralog_DF)){
  paralog1 <- paralog_DF$syst1[i]
  paralog2 <- paralog_DF$syst2[i]
  
  #print(paste(paralog1, paralog2))
  
  # Note: Some paralogs cause a bug, skip these
  if ((paralog1 %in% bug) & (paralog2 %in% bug)){
    next
  }
  
  alignment <- 
    read.alignment(file = paste0(dir, "R_Paralog_Analysis/data/reverseAlign_output_Files/",
                                 paralog1,
                                 "_",
                                 paralog2,
                                 ".fasta"),
                   format = "fasta")
  
  Nucl_align_list[paste0(paralog1, "_", paralog2)] <- list(alignment)
}

### Calculate Ka and Ks (synonymous and non-synonymous substitutions)

kaks_list <- list()

for (paralogs in names(Nucl_align_list)){
  kaks_computation <- kaks(Nucl_align_list[[paralogs]])
  
  kaks_list[paralogs] <- list(kaks_computation)
}


## Make a DF from Ka data
kaks_paralog_names <- names(kaks_list)

paralog_kaks_DF <- do.call(rbind, kaks_list)

paralog_kaks_DF <- as.data.frame(paralog_kaks_DF)

rownames(paralog_kaks_DF) <- NULL

paralog_kaks_DF$paralogs <- kaks_paralog_names
paralog_kaks_DF$ka <- as.numeric(paralog_kaks_DF$ka)
paralog_kaks_DF$ks <- as.numeric(paralog_kaks_DF$ks)
paralog_kaks_DF$vka <- as.numeric(paralog_kaks_DF$vka)
paralog_kaks_DF$vks <- as.numeric(paralog_kaks_DF$vks)

#write.csv(paralog_kaks_DF, paste0(dir, "R_Paralog_Analysis/results/all_paralog_kaks.csv"))


#### Assign Sequence Divergences to Kd/GI dataset ####

## Assign sequence divergence to paralogs

MS_Kd_GI_AllPara_DF$source_Ka <- NA
MS_Kd_GI_AllPara_DF$target_Ka <- NA

paralog_list <- unlist(strsplit(paste0(names(kaks_list), collapse = "_"), split = "_"))

for (i in 1:nrow(MS_Kd_GI_AllPara_DF)){
  
  if (MS_Kd_GI_AllPara_DF$Type[i] == "NoParalog"){
    next
  }
  
  gene1 <- MS_Kd_GI_AllPara_DF$source[i]
  gene2 <- MS_Kd_GI_AllPara_DF$target[i]
  
  #print(paste(gene1, gene2))
  
  if ((gene1 %in% bug) | (gene2 %in% bug)){
    next
  }
  
  gene1_paralogs <- names(kaks_list)[which(grepl(gene1, names(kaks_list)) == TRUE)]
  gene2_paralogs <- names(kaks_list)[which(grepl(gene2, names(kaks_list)) == TRUE)]
  
  
  if((length(gene1_paralogs)>0) & (length(gene2_paralogs)>0)){ ## Paralogous heteromers; IntParalog
    if ((gene1_paralogs == gene2_paralogs)){
      MS_Kd_GI_AllPara_DF$source_Ka[i] <- kaks_list[[gene1_paralogs]]$ka
      MS_Kd_GI_AllPara_DF$target_Ka[i] <- kaks_list[[gene2_paralogs]]$ka
    }
  }
  
  if (MS_Kd_GI_AllPara_DF$Type[i] == "OneParalog"){
    if (gene1 %in% paralog_list){
      MS_Kd_GI_AllPara_DF$source_Ka[i] <- kaks_list[[gene1_paralogs]]$ka
      }
    if (gene2 %in% paralog_list){
      MS_Kd_GI_AllPara_DF$target_Ka[i] <- kaks_list[[gene2_paralogs]]$ka
    }
  }
  
  if (MS_Kd_GI_AllPara_DF$Type[i] == "TwoParalog"){
    MS_Kd_GI_AllPara_DF$source_Ka[i] <- kaks_list[[gene1_paralogs]]$ka
    MS_Kd_GI_AllPara_DF$target_Ka[i] <- kaks_list[[gene2_paralogs]]$ka
  }
}


### Export
#write.csv(MS_Kd_GI_AllPara_DF, paste0(dir, "R_Paralog_Analysis/results/Kd_GI_Divergence.csv"))

